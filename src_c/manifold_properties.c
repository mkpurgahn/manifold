// Copyright 2021 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Properties computation for the C11 Manifold port.

#include "manifold_impl.h"
#include <stdlib.h>
#include <string.h>

static int cmp_int_by_u32(const void *a_ptr, const void *b_ptr, void *ctx) {
  const uint32_t *morton = (const uint32_t *)ctx;
  int a = *(const int *)a_ptr;
  int b = *(const int *)b_ptr;
  if (morton[a] < morton[b]) return -1;
  if (morton[a] > morton[b]) return 1;
  return (a < b) ? -1 : (a > b) ? 1 : 0;
}

#ifdef __APPLE__
static int cmp_int_by_u32_apple(void *ctx, const void *a_ptr, const void *b_ptr) {
  return cmp_int_by_u32(a_ptr, b_ptr, ctx);
}
#endif

double manifold_impl_get_volume(const ManifoldImpl *impl) {
  double vol = 0.0;
  double comp = 0.0;  // Kahan compensation
  size_t numTri = manifold_impl_num_tri(impl);
  for (size_t tri = 0; tri < numTri; tri++) {
    ManifoldVec3 a = impl->vertPos.data[impl->halfedge.data[3 * tri].startVert];
    ManifoldVec3 b = impl->vertPos.data[impl->halfedge.data[3 * tri + 1].startVert];
    ManifoldVec3 c = impl->vertPos.data[impl->halfedge.data[3 * tri + 2].startVert];
    ManifoldVec3 crossP = vec3_cross(
        vec3_sub(b, a), vec3_sub(c, a));
    double v1 = vec3_dot(crossP, a) / 6.0;
    double t = vol + v1;
    comp += (vol - t) + v1;
    vol = t;
  }
  return vol + comp;
}

double manifold_impl_get_surface_area(const ManifoldImpl *impl) {
  double area = 0.0;
  double comp = 0.0;  // Kahan compensation
  size_t numTri = manifold_impl_num_tri(impl);
  for (size_t tri = 0; tri < numTri; tri++) {
    ManifoldVec3 a = impl->vertPos.data[impl->halfedge.data[3 * tri].startVert];
    ManifoldVec3 b = impl->vertPos.data[impl->halfedge.data[3 * tri + 1].startVert];
    ManifoldVec3 c = impl->vertPos.data[impl->halfedge.data[3 * tri + 2].startVert];
    ManifoldVec3 cross = vec3_cross(vec3_sub(b, a), vec3_sub(c, a));
    double a1 = vec3_length(cross) / 2.0;
    double t = area + a1;
    comp += (area - t) + a1;
    area = t;
  }
  return area + comp;
}

void manifold_impl_calculate_curvature(ManifoldImpl *impl, int gaussianIdx,
                                        int meanIdx) {
  if (manifold_impl_is_empty(impl)) return;
  if (gaussianIdx < 0 && meanIdx < 0) return;

  size_t numVert = manifold_impl_num_vert(impl);
  size_t numTri = manifold_impl_num_tri(impl);

  double *vertMeanCurvature = (double *)calloc(numVert, sizeof(double));
  double *vertGaussianCurvature = (double *)malloc(numVert * sizeof(double));
  double *vertArea = (double *)calloc(numVert, sizeof(double));
  double *degree = (double *)calloc(numVert, sizeof(double));

  for (size_t i = 0; i < numVert; i++) {
    vertGaussianCurvature[i] = MANIFOLD_TWO_PI;
  }

  for (size_t tri = 0; tri < numTri; tri++) {
    ManifoldVec3 edge[3];
    double edgeLength[3];
    for (int i = 0; i < 3; i++) {
      int sv = impl->halfedge.data[3 * tri + i].startVert;
      int ev = impl->halfedge.data[3 * tri + i].endVert;
      edge[i] = vec3_sub(impl->vertPos.data[ev], impl->vertPos.data[sv]);
      edgeLength[i] = vec3_length(edge[i]);
      if (edgeLength[i] > 0) edge[i] = vec3_scale(edge[i], 1.0 / edgeLength[i]);

      int ph = impl->halfedge.data[3 * tri + i].pairedHalfedge;
      if (ph >= 0 && ph < (int)impl->halfedge.len) {
        int neighborTri = ph / 3;
        if (neighborTri < (int)numTri && tri < numTri) {
          ManifoldVec3 n1 = impl->faceNormal.data[tri];
          ManifoldVec3 n2 = impl->faceNormal.data[neighborTri];
          double crossDotEdge = vec3_dot(vec3_cross(n1, n2), edge[i]);
          double sinVal = crossDotEdge;
          if (sinVal > 1.0) sinVal = 1.0;
          if (sinVal < -1.0) sinVal = -1.0;
          double dihedral = 0.25 * edgeLength[i] * asin(sinVal);
          vertMeanCurvature[sv] += dihedral;
          vertMeanCurvature[ev] += dihedral;
        }
      }
      degree[sv] += 1.0;
    }

    // Angles
    double phi[3];
    double d02 = vec3_dot(edge[2], edge[0]);
    double d01 = vec3_dot(edge[0], edge[1]);
    phi[0] = acos(fmax(-1.0, fmin(1.0, -d02)));
    // dot of -edge[0] and edge[1]
    phi[1] = acos(fmax(-1.0, fmin(1.0, -d01)));
    phi[2] = MANIFOLD_PI - phi[0] - phi[1];
    double area3 = edgeLength[0] * edgeLength[1] *
                   vec3_length(vec3_cross(edge[0], edge[1])) / 6.0;

    for (int i = 0; i < 3; i++) {
      int vert = impl->halfedge.data[3 * tri + i].startVert;
      vertGaussianCurvature[vert] -= phi[i];
      vertArea[vert] += area3;
    }
  }

  // Normalize
  for (size_t i = 0; i < numVert; i++) {
    if (vertArea[i] > 0) {
      double factor = degree[i] / (6.0 * vertArea[i]);
      vertMeanCurvature[i] *= factor;
      vertGaussianCurvature[i] *= factor;
    }
  }

  // Expand properties
  int oldNumProp = impl->numProp;
  int newNumProp = oldNumProp;
  if (gaussianIdx >= 0 && gaussianIdx + 1 > newNumProp)
    newNumProp = gaussianIdx + 1;
  if (meanIdx >= 0 && meanIdx + 1 > newNumProp)
    newNumProp = meanIdx + 1;

  size_t numPropVert = manifold_impl_num_prop_vert(impl);
  ManifoldVecDouble oldProperties = impl->properties;

  impl->properties = (ManifoldVecDouble)MANIFOLD_VEC_INIT;
  vec_double_resize(&impl->properties, (size_t)newNumProp * numPropVert);
  memset(impl->properties.data, 0, sizeof(double) * impl->properties.len);
  impl->numProp = newNumProp;

  // Track which propVerts we've already processed
  uint8_t *visited = (uint8_t *)calloc(numPropVert, sizeof(uint8_t));

  for (size_t tri = 0; tri < numTri; tri++) {
    for (int i = 0; i < 3; i++) {
      ManifoldHalfedge *he = &impl->halfedge.data[3 * tri + i];
      int vert = he->startVert;
      int propVert = he->propVert;
      if (propVert < 0 || (size_t)propVert >= numPropVert) continue;
      if (visited[propVert]) continue;
      visited[propVert] = 1;

      // Copy old properties
      for (int p = 0; p < oldNumProp; p++) {
        impl->properties.data[newNumProp * propVert + p] =
            oldProperties.data[oldNumProp * propVert + p];
      }

      if (gaussianIdx >= 0 && vert >= 0 && (size_t)vert < numVert) {
        impl->properties.data[newNumProp * propVert + gaussianIdx] =
            vertGaussianCurvature[vert];
      }
      if (meanIdx >= 0 && vert >= 0 && (size_t)vert < numVert) {
        impl->properties.data[newNumProp * propVert + meanIdx] =
            vertMeanCurvature[vert];
      }
    }
  }

  free(visited);
  vec_double_free(&oldProperties);
  free(vertMeanCurvature);
  free(vertGaussianCurvature);
  free(vertArea);
  free(degree);
}

void manifold_impl_set_properties(ManifoldImpl *impl, int numProp,
    void (*propFunc)(double *newProp, ManifoldVec3 pos, const double *oldProp, void *ctx),
    void *ctx) {
  int oldNumProp = impl->numProp;
  ManifoldVecDouble oldProperties = impl->properties;
  size_t numPropVert = manifold_impl_num_prop_vert(impl);

  if (numProp == 0) {
    impl->properties = (ManifoldVecDouble)MANIFOLD_VEC_INIT;
    impl->numProp = 0;
    vec_double_free(&oldProperties);
    return;
  }

  impl->properties = (ManifoldVecDouble)MANIFOLD_VEC_INIT;
  vec_double_resize(&impl->properties, (size_t)numProp * numPropVert);
  memset(impl->properties.data, 0, sizeof(double) * impl->properties.len);
  impl->numProp = numProp;

  if (propFunc == NULL) {
    vec_double_free(&oldProperties);
    return;
  }

  size_t numTri = manifold_impl_num_tri(impl);
  for (size_t tri = 0; tri < numTri; tri++) {
    for (int i = 0; i < 3; i++) {
      ManifoldHalfedge *he = &impl->halfedge.data[3 * tri + i];
      int vert = he->startVert;
      int propVert = he->propVert;
      if (propVert < 0 || (size_t)propVert >= numPropVert) continue;
      if (vert < 0 || (size_t)vert >= impl->vertPos.len) continue;

      const double *oldProp = (oldNumProp > 0 && oldProperties.data)
          ? &oldProperties.data[oldNumProp * propVert]
          : NULL;

      propFunc(&impl->properties.data[numProp * propVert],
               impl->vertPos.data[vert], oldProp, ctx);
    }
  }

  vec_double_free(&oldProperties);
}

// Decompose into connected components
int manifold_impl_decompose(const ManifoldImpl *impl, ManifoldImpl *components,
                             int maxComponents) {
  size_t numVert = manifold_impl_num_vert(impl);
  size_t numTri = manifold_impl_num_tri(impl);
  if (numVert == 0 || numTri == 0) return 0;

  // Union-Find
  int *parent = (int *)malloc(numVert * sizeof(int));
  int *rank_ = (int *)calloc(numVert, sizeof(int));
  for (size_t i = 0; i < numVert; i++) parent[i] = (int)i;

  // Find with path compression
  #define UF_FIND(x) uf_find(parent, x)
  // Inline iterative find
  int uf_find_val;
  (void)uf_find_val;

  for (size_t i = 0; i < impl->halfedge.len; i++) {
    ManifoldHalfedge he = impl->halfedge.data[i];
    if (he.startVert < 0 || he.endVert < 0) continue;
    if (he.startVert >= (int)numVert || he.endVert >= (int)numVert) continue;
    // Only process forward edges
    if (he.startVert >= he.endVert) continue;

    // Find roots
    int a = he.startVert, b = he.endVert;
    while (parent[a] != a) { parent[a] = parent[parent[a]]; a = parent[a]; }
    while (parent[b] != b) { parent[b] = parent[parent[b]]; b = parent[b]; }
    if (a != b) {
      if (rank_[a] < rank_[b]) { int t = a; a = b; b = t; }
      parent[b] = a;
      if (rank_[a] == rank_[b]) rank_[a]++;
    }
  }

  // Find all roots and assign component indices
  int *componentLabel = (int *)malloc(numVert * sizeof(int));
  int numComponents = 0;
  int *rootToComponent = (int *)malloc(numVert * sizeof(int));
  for (size_t i = 0; i < numVert; i++) rootToComponent[i] = -1;

  for (size_t i = 0; i < numVert; i++) {
    int r = (int)i;
    while (parent[r] != r) { parent[r] = parent[parent[r]]; r = parent[r]; }
    parent[i] = r;
    if (rootToComponent[r] == -1) {
      rootToComponent[r] = numComponents++;
    }
    componentLabel[i] = rootToComponent[r];
  }

  if (numComponents <= 1) {
    // Single component - just copy the whole thing
    if (maxComponents >= 1 && components) {
      manifold_impl_init(&components[0]);
      components[0].bBox = impl->bBox;
      components[0].epsilon = impl->epsilon;
      components[0].tolerance = impl->tolerance;
      components[0].numProp = impl->numProp;
      components[0].status = impl->status;
      components[0].vertPos = vec_vec3_copy(&impl->vertPos);
      components[0].halfedge = vec_halfedge_copy(&impl->halfedge);
      components[0].properties = vec_double_copy(&impl->properties);
      components[0].vertNormal = vec_vec3_copy(&impl->vertNormal);
      components[0].faceNormal = vec_vec3_copy(&impl->faceNormal);
      components[0].halfedgeTangent = vec_vec4_copy(&impl->halfedgeTangent);
      components[0].meshRelation.originalID = impl->meshRelation.originalID;
      components[0].meshRelation.meshIDtransform =
          vec_meshid_copy(&impl->meshRelation.meshIDtransform);
      components[0].meshRelation.triRef =
          vec_triref_copy(&impl->meshRelation.triRef);
    }
    free(parent);
    free(rank_);
    free(componentLabel);
    free(rootToComponent);
    return 1;
  }

  int outputCount = numComponents < maxComponents ? numComponents : maxComponents;

  for (int comp = 0; comp < outputCount; comp++) {
    // Build vertex remap: old index → new index
    int nVert = 0;
    int *vertOld2New = (int *)malloc(numVert * sizeof(int));
    for (size_t i = 0; i < numVert; i++) vertOld2New[i] = -1;
    for (size_t i = 0; i < numVert; i++) {
      if (componentLabel[i] == comp) {
        vertOld2New[i] = nVert++;
      }
    }

    // Gather face indices belonging to this component
    int *faceNew2Old = (int *)malloc(numTri * sizeof(int));
    int nFace = 0;
    for (size_t tri = 0; tri < numTri; tri++) {
      int sv = impl->halfedge.data[3 * tri].startVert;
      if (sv < 0 || (size_t)sv >= numVert) continue;
      if (componentLabel[sv] != comp) continue;
      faceNew2Old[nFace++] = (int)tri;
    }

    if (nFace == 0) {
      free(vertOld2New);
      free(faceNew2Old);
      manifold_impl_init(&components[comp]);
      continue;
    }

    manifold_impl_init(&components[comp]);
    components[comp].epsilon = impl->epsilon;
    components[comp].tolerance = impl->tolerance;

    // Gather vertex positions and normals
    vec_vec3_resize(&components[comp].vertPos, (size_t)nVert);
    if (impl->vertNormal.len >= numVert)
      vec_vec3_resize(&components[comp].vertNormal, (size_t)nVert);
    for (size_t i = 0; i < numVert; i++) {
      if (componentLabel[i] == comp) {
        components[comp].vertPos.data[vertOld2New[i]] = impl->vertPos.data[i];
        if (impl->vertNormal.len >= numVert)
          components[comp].vertNormal.data[vertOld2New[i]] = impl->vertNormal.data[i];
      }
    }

    // Copy properties wholesale (like C++ GatherFaces)
    if (impl->numProp > 0) {
      components[comp].numProp = impl->numProp;
      components[comp].properties = vec_double_copy(&impl->properties);
    }

    // Build faceOld2New for pairedHalfedge remapping
    int *faceOld2New = (int *)malloc(numTri * sizeof(int));
    for (size_t i = 0; i < numTri; i++) faceOld2New[i] = -1;
    for (int i = 0; i < nFace; i++) faceOld2New[faceNew2Old[i]] = i;

    // Gather halfedges and tangents
    int numHe = nFace * 3;
    vec_halfedge_resize(&components[comp].halfedge, (size_t)numHe);
    if (impl->halfedgeTangent.len > 0)
      vec_vec4_resize(&components[comp].halfedgeTangent, (size_t)numHe);
    for (int i = 0; i < nFace; i++) {
      int oldTri = faceNew2Old[i];
      for (int j = 0; j < 3; j++) {
        ManifoldHalfedge he = impl->halfedge.data[3 * oldTri + j];
        he.startVert = vertOld2New[he.startVert];
        he.endVert = vertOld2New[he.endVert];
        // Remap pairedHalfedge
        int oldPair = he.pairedHalfedge;
        if (oldPair >= 0) {
          int oldPairTri = oldPair / 3;
          int oldPairIdx = oldPair % 3;
          int newPairTri = faceOld2New[oldPairTri];
          he.pairedHalfedge = (newPairTri >= 0) ? (3 * newPairTri + oldPairIdx) : -1;
        }
        // propVert stays the same (indexes into global properties array)
        components[comp].halfedge.data[3 * i + j] = he;
        if (impl->halfedgeTangent.len > 0)
          components[comp].halfedgeTangent.data[3 * i + j] =
              impl->halfedgeTangent.data[3 * oldTri + j];
      }
    }

    // Gather triRef and meshIDtransform
    if (impl->meshRelation.triRef.len > 0) {
      vec_triref_resize(&components[comp].meshRelation.triRef, (size_t)nFace);
      for (int i = 0; i < nFace; i++) {
        components[comp].meshRelation.triRef.data[i] =
            impl->meshRelation.triRef.data[faceNew2Old[i]];
      }
    }
    for (size_t i = 0; i < impl->meshRelation.meshIDtransform.len; i++) {
      vec_meshid_push(&components[comp].meshRelation.meshIDtransform,
                      impl->meshRelation.meshIDtransform.data[i]);
    }

    // Gather face normals
    if (impl->faceNormal.len >= numTri) {
      vec_vec3_resize(&components[comp].faceNormal, (size_t)nFace);
      for (int i = 0; i < nFace; i++) {
        components[comp].faceNormal.data[i] =
            impl->faceNormal.data[faceNew2Old[i]];
      }
    }

    manifold_impl_calculate_bbox(&components[comp]);
    manifold_impl_set_epsilon(&components[comp], impl->epsilon, false);
    manifold_impl_sort_geometry(&components[comp]);

    free(vertOld2New);
    free(faceNew2Old);
    free(faceOld2New);
  }

  free(parent);
  free(rank_);
  free(componentLabel);
  free(rootToComponent);
  return outputCount;
}

// ===== Triangle Distance (from tri_dist.h) =====

static void edge_edge_dist(ManifoldVec3 *x, ManifoldVec3 *y,
                           ManifoldVec3 p, ManifoldVec3 a,
                           ManifoldVec3 q, ManifoldVec3 b) {
  ManifoldVec3 T = vec3_sub(q, p);
  double ADotA = vec3_dot(a, a);
  double BDotB = vec3_dot(b, b);
  double ADotB = vec3_dot(a, b);
  double ADotT = vec3_dot(a, T);
  double BDotT = vec3_dot(b, T);
  double Denom = ADotA * BDotB - ADotB * ADotB;
  double t = (Denom != 0.0)
      ? fmax(0.0, fmin((ADotT * BDotB - BDotT * ADotB) / Denom, 1.0))
      : 0.0;
  double u;
  if (BDotB != 0.0) {
    u = (t * ADotB - BDotT) / BDotB;
    if (u < 0.0) {
      u = 0.0;
      t = (ADotA != 0.0) ? fmax(0.0, fmin(ADotT / ADotA, 1.0)) : 0.0;
    } else if (u > 1.0) {
      u = 1.0;
      t = (ADotA != 0.0) ? fmax(0.0, fmin((ADotB + ADotT) / ADotA, 1.0)) : 0.0;
    }
  } else {
    u = 0.0;
    t = (ADotA != 0.0) ? fmax(0.0, fmin(ADotT / ADotA, 1.0)) : 0.0;
  }
  *x = vec3_add(p, vec3_scale(a, t));
  *y = vec3_add(q, vec3_scale(b, u));
}

double distance_tri_tri_squared(ManifoldVec3 p[3], ManifoldVec3 q[3]) {
  ManifoldVec3 Sv[3], Tv[3];
  Sv[0] = vec3_sub(p[1], p[0]);
  Sv[1] = vec3_sub(p[2], p[1]);
  Sv[2] = vec3_sub(p[0], p[2]);
  Tv[0] = vec3_sub(q[1], q[0]);
  Tv[1] = vec3_sub(q[2], q[1]);
  Tv[2] = vec3_sub(q[0], q[2]);

  bool shown_disjoint = false;
  double mindd = 1e300;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ManifoldVec3 cp, cq;
      edge_edge_dist(&cp, &cq, p[i], Sv[i], q[j], Tv[j]);
      ManifoldVec3 V = vec3_sub(cq, cp);
      double dd = vec3_dot(V, V);
      if (dd <= mindd) {
        mindd = dd;
        int id = i + 2;
        if (id >= 3) id -= 3;
        ManifoldVec3 Z = vec3_sub(p[id], cp);
        double a = vec3_dot(Z, V);
        id = j + 2;
        if (id >= 3) id -= 3;
        Z = vec3_sub(q[id], cq);
        double b = vec3_dot(Z, V);
        if (a <= 0.0 && b >= 0.0) return vec3_dot(V, V);
        if (a <= 0.0) a = 0.0;
        else if (b > 0.0) b = 0.0;
        if ((mindd - a + b) > 0.0) shown_disjoint = true;
      }
    }
  }

  ManifoldVec3 Sn = vec3_cross(Sv[0], Sv[1]);
  double Snl = vec3_dot(Sn, Sn);
  if (Snl > 1e-15) {
    double Tp[3];
    Tp[0] = vec3_dot(vec3_sub(p[0], q[0]), Sn);
    Tp[1] = vec3_dot(vec3_sub(p[0], q[1]), Sn);
    Tp[2] = vec3_dot(vec3_sub(p[0], q[2]), Sn);
    int index = -1;
    if (Tp[0] > 0.0 && Tp[1] > 0.0 && Tp[2] > 0.0) {
      index = (Tp[0] < Tp[1]) ? 0 : 1;
      if (Tp[2] < Tp[index]) index = 2;
    } else if (Tp[0] < 0.0 && Tp[1] < 0.0 && Tp[2] < 0.0) {
      index = (Tp[0] > Tp[1]) ? 0 : 1;
      if (Tp[2] > Tp[index]) index = 2;
    }
    if (index >= 0) {
      shown_disjoint = true;
      ManifoldVec3 qIdx = q[index];
      ManifoldVec3 V = vec3_sub(qIdx, p[0]);
      ManifoldVec3 Z = vec3_cross(Sn, Sv[0]);
      if (vec3_dot(V, Z) > 0.0) {
        V = vec3_sub(qIdx, p[1]);
        Z = vec3_cross(Sn, Sv[1]);
        if (vec3_dot(V, Z) > 0.0) {
          V = vec3_sub(qIdx, p[2]);
          Z = vec3_cross(Sn, Sv[2]);
          if (vec3_dot(V, Z) > 0.0) {
            ManifoldVec3 cp = vec3_add(qIdx, vec3_scale(Sn, Tp[index] / Snl));
            ManifoldVec3 diff = vec3_sub(cp, qIdx);
            return vec3_dot(diff, diff);
          }
        }
      }
    }
  }

  ManifoldVec3 Tn = vec3_cross(Tv[0], Tv[1]);
  double Tnl = vec3_dot(Tn, Tn);
  if (Tnl > 1e-15) {
    double Sp[3];
    Sp[0] = vec3_dot(vec3_sub(q[0], p[0]), Tn);
    Sp[1] = vec3_dot(vec3_sub(q[0], p[1]), Tn);
    Sp[2] = vec3_dot(vec3_sub(q[0], p[2]), Tn);
    int index = -1;
    if (Sp[0] > 0.0 && Sp[1] > 0.0 && Sp[2] > 0.0) {
      index = (Sp[0] < Sp[1]) ? 0 : 1;
      if (Sp[2] < Sp[index]) index = 2;
    } else if (Sp[0] < 0.0 && Sp[1] < 0.0 && Sp[2] < 0.0) {
      index = (Sp[0] > Sp[1]) ? 0 : 1;
      if (Sp[2] > Sp[index]) index = 2;
    }
    if (index >= 0) {
      shown_disjoint = true;
      ManifoldVec3 pIdx = p[index];
      ManifoldVec3 V = vec3_sub(pIdx, q[0]);
      ManifoldVec3 Z = vec3_cross(Tn, Tv[0]);
      if (vec3_dot(V, Z) > 0.0) {
        V = vec3_sub(pIdx, q[1]);
        Z = vec3_cross(Tn, Tv[1]);
        if (vec3_dot(V, Z) > 0.0) {
          V = vec3_sub(pIdx, q[2]);
          Z = vec3_cross(Tn, Tv[2]);
          if (vec3_dot(V, Z) > 0.0) {
            ManifoldVec3 cq = vec3_add(pIdx, vec3_scale(Tn, Sp[index] / Tnl));
            ManifoldVec3 diff = vec3_sub(pIdx, cq);
            return vec3_dot(diff, diff);
          }
        }
      }
    }
  }

  return shown_disjoint ? mindd : 0.0;
}

// ===== MinGap =====

typedef struct {
  const ManifoldImpl *self;
  const ManifoldImpl *other;
  double minDistSq;
} MinGapCtx;

MANIFOLD_UNUSED static void mingap_callback(int queryIdx, int leafIdx, void *ctx) {
  MinGapCtx *mg = (MinGapCtx *)ctx;
  ManifoldVec3 p[3], q[3];
  for (int j = 0; j < 3; j++) {
    p[j] = mg->self->vertPos.data[mg->self->halfedge.data[3 * queryIdx + j].startVert];
    q[j] = mg->other->vertPos.data[mg->other->halfedge.data[3 * leafIdx + j].startVert];
  }
  double d = distance_tri_tri_squared(p, q);
  if (d < mg->minDistSq) mg->minDistSq = d;
}

double manifold_impl_min_gap(const ManifoldImpl *self,
                              const ManifoldImpl *other,
                              double searchLength) {
  if (manifold_impl_is_empty(self) || manifold_impl_is_empty(other))
    return searchLength;

  size_t numTriOther = manifold_impl_num_tri(other);

  // Build BVH for other
  ManifoldVecBox faceBoxOther = {0};
  ManifoldVecU32 faceMortonOther = {0};
  manifold_impl_get_face_box_morton(other, &faceBoxOther, &faceMortonOther);

  // Expand boxes by searchLength
  for (size_t i = 0; i < faceBoxOther.len; i++) {
    ManifoldVec3 sl = manifold_vec3_splat(searchLength);
    faceBoxOther.data[i].min = vec3_sub(faceBoxOther.data[i].min, sl);
    faceBoxOther.data[i].max = vec3_add(faceBoxOther.data[i].max, sl);
  }

  // Sort morton codes and reorder boxes
  ManifoldVecInt sortedIdx = {0};
  vec_int_resize(&sortedIdx, numTriOther);
  for (size_t i = 0; i < numTriOther; i++) sortedIdx.data[i] = (int)i;
  // Sort by morton codes using qsort_r
  {
    const uint32_t *mortonCtx = faceMortonOther.data;
#ifdef __APPLE__
    qsort_r(sortedIdx.data, numTriOther, sizeof(int), (void *)mortonCtx,
            cmp_int_by_u32_apple);
#else
    qsort_r(sortedIdx.data, numTriOther, sizeof(int), cmp_int_by_u32,
            (void *)mortonCtx);
#endif
  }

  // Create sorted arrays
  ManifoldVecBox sortedBoxes = {0};
  ManifoldVecU32 sortedMorton = {0};
  vec_box_resize(&sortedBoxes, numTriOther);
  vec_u32_resize(&sortedMorton, numTriOther);
  for (size_t i = 0; i < numTriOther; i++) {
    sortedBoxes.data[i] = faceBoxOther.data[sortedIdx.data[i]];
    sortedMorton.data[i] = faceMortonOther.data[sortedIdx.data[i]];
  }

  // Build collider from other's faces
  ManifoldCollider collOther = manifold_collider_create(
      sortedBoxes.data, sortedMorton.data, numTriOther);

  // Query self's faces against other's BVH
  ManifoldVecBox faceBoxSelf = {0};
  ManifoldVecU32 faceMortonSelf = {0};
  manifold_impl_get_face_box_morton(self, &faceBoxSelf, &faceMortonSelf);

  MinGapCtx mgCtx;
  mgCtx.self = self;
  mgCtx.other = other;
  mgCtx.minDistSq = searchLength * searchLength;

  // We need to map back from sorted leaf indices to original tri indices
  // Use a wrapper callback
  typedef struct {
    MinGapCtx *base;
    const int *leafMap;
  } MappedMinGapCtx;

  MappedMinGapCtx mapped;
  mapped.base = &mgCtx;
  mapped.leafMap = sortedIdx.data;

  // Custom callback that maps leaf indices back
  // We can't easily use a closure, so inline the collider query
  size_t numTriSelf = manifold_impl_num_tri(self);
  for (size_t qi = 0; qi < numTriSelf; qi++) {
    ManifoldBox queryBox = faceBoxSelf.data[qi];
    if (collOther.internalChildren.len == 0) continue;
    int stack[64];
    int top = -1;
    int node = MANIFOLD_COLLIDER_ROOT;
    while (1) {
      int internal = collider_node2internal(node);
      int child1 = collOther.internalChildren.data[internal].first;
      int child2 = collOther.internalChildren.data[internal].second;
      bool overlap1 = manifold_box_overlaps(collOther.nodeBBox.data[child1], queryBox);
      bool overlap2 = manifold_box_overlaps(collOther.nodeBBox.data[child2], queryBox);
      bool traverse1 = false, traverse2 = false;
      if (overlap1 && collider_is_leaf(child1)) {
        int leafIdx = collider_node2leaf(child1);
        int origTri = sortedIdx.data[leafIdx];
        // Compute tri-tri distance
        ManifoldVec3 p[3], q[3];
        for (int j = 0; j < 3; j++) {
          p[j] = self->vertPos.data[self->halfedge.data[3 * qi + j].startVert];
          q[j] = other->vertPos.data[other->halfedge.data[3 * origTri + j].startVert];
        }
        double d = distance_tri_tri_squared(p, q);
        if (d < mgCtx.minDistSq) mgCtx.minDistSq = d;
      } else {
        traverse1 = overlap1 && collider_is_internal(child1);
      }
      if (overlap2 && collider_is_leaf(child2)) {
        int leafIdx = collider_node2leaf(child2);
        int origTri = sortedIdx.data[leafIdx];
        ManifoldVec3 p[3], q[3];
        for (int j = 0; j < 3; j++) {
          p[j] = self->vertPos.data[self->halfedge.data[3 * qi + j].startVert];
          q[j] = other->vertPos.data[other->halfedge.data[3 * origTri + j].startVert];
        }
        double d = distance_tri_tri_squared(p, q);
        if (d < mgCtx.minDistSq) mgCtx.minDistSq = d;
      } else {
        traverse2 = overlap2 && collider_is_internal(child2);
      }
      if (!traverse1 && !traverse2) {
        if (top < 0) break;
        node = stack[top--];
      } else {
        node = traverse1 ? child1 : child2;
        if (traverse1 && traverse2) stack[++top] = child2;
      }
    }
  }

  double result = sqrt(fmin(mgCtx.minDistSq, searchLength * searchLength));

  vec_box_free(&faceBoxOther);
  vec_u32_free(&faceMortonOther);
  vec_box_free(&sortedBoxes);
  vec_u32_free(&sortedMorton);
  vec_int_free(&sortedIdx);
  manifold_collider_free(&collOther);
  vec_box_free(&faceBoxSelf);
  vec_u32_free(&faceMortonSelf);

  return result;
}

// ===== CalculateNormals =====

void manifold_impl_calculate_normals(ManifoldImpl *impl, int normalIdx,
                                      double minSharpAngle) {
  manifold_impl_set_normals_smooth(impl, normalIdx, minSharpAngle);
}

bool manifold_impl_is_convex(const ManifoldImpl *impl) {
  // Convex shape must have genus of 0
  int chi = (int)manifold_impl_num_vert(impl) -
            (int)manifold_impl_num_edge(impl) +
            (int)manifold_impl_num_tri(impl);
  int genus = 1 - chi / 2;
  if (genus != 0) return false;

  // Compute scale-based tolerance for near-convex shapes
  ManifoldVec3 sz = vec3_sub(impl->bBox.max, impl->bBox.min);
  double extent = fmax(fmax(fabs(sz.x), fabs(sz.y)), fabs(sz.z));
  double tol = -extent * 1e-2;

  size_t nbEdges = impl->halfedge.len;
  for (size_t idx = 0; idx < nbEdges; idx++) {
    ManifoldHalfedge edge = impl->halfedge.data[idx];
    if (!manifold_halfedge_is_forward(&edge)) continue;

    ManifoldVec3 normal0 = impl->faceNormal.data[idx / 3];
    if (edge.pairedHalfedge < 0 ||
        (size_t)edge.pairedHalfedge >= impl->halfedge.len)
      return false;
    ManifoldVec3 normal1 = impl->faceNormal.data[edge.pairedHalfedge / 3];

    if (vec3_equal(normal0, normal1)) continue;

    ManifoldVec3 edgeVec = vec3_sub(impl->vertPos.data[edge.endVert],
                                     impl->vertPos.data[edge.startVert]);
    double dot = vec3_dot(edgeVec, vec3_cross(normal0, normal1));
    if (dot < tol) return false;
  }
  return true;
}
