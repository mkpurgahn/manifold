// Copyright 2021 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Properties computation for the C11 Manifold port.

#include "manifold_impl.h"
#include <stdlib.h>
#include <string.h>

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
    // Count verts and tris in this component
    int nVert = 0;
    int *vertOld2New = (int *)malloc(numVert * sizeof(int));
    for (size_t i = 0; i < numVert; i++) vertOld2New[i] = -1;
    for (size_t i = 0; i < numVert; i++) {
      if (componentLabel[i] == comp) {
        vertOld2New[i] = nVert++;
      }
    }

    ManifoldVecIVec3 triVerts = {0};
    for (size_t tri = 0; tri < numTri; tri++) {
      int sv = impl->halfedge.data[3 * tri].startVert;
      if (sv < 0 || (size_t)sv >= numVert) continue;
      if (componentLabel[sv] != comp) continue;

      ManifoldIVec3 tv;
      tv.x = vertOld2New[impl->halfedge.data[3 * tri].startVert];
      tv.y = vertOld2New[impl->halfedge.data[3 * tri + 1].startVert];
      tv.z = vertOld2New[impl->halfedge.data[3 * tri + 2].startVert];
      if (tv.x < 0 || tv.y < 0 || tv.z < 0) continue;
      vec_ivec3_push(&triVerts, tv);
    }

    manifold_impl_init(&components[comp]);
    components[comp].epsilon = impl->epsilon;
    components[comp].tolerance = impl->tolerance;

    // Gather verts
    for (size_t i = 0; i < numVert; i++) {
      if (componentLabel[i] == comp) {
        vec_vec3_push(&components[comp].vertPos, impl->vertPos.data[i]);
      }
    }

    ManifoldVecIVec3 emptyTriVert = {0};
    manifold_impl_create_halfedges(&components[comp], &triVerts, &emptyTriVert);
    manifold_impl_initialize_original(&components[comp]);
    manifold_impl_calculate_bbox(&components[comp]);
    manifold_impl_set_epsilon(&components[comp], impl->epsilon, false);
    manifold_impl_sort_geometry(&components[comp]);
    manifold_impl_set_normals_and_coplanar(&components[comp]);

    vec_ivec3_free(&triVerts);
    vec_ivec3_free(&emptyTriVert);
    free(vertOld2New);
  }

  free(parent);
  free(rank_);
  free(componentLabel);
  free(rootToComponent);
  return outputCount;
}
