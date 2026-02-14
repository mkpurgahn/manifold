// Copyright 2021 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Core implementation for the C11 Manifold port.

#include "manifold_impl.h"
#include <stdio.h>
#include <limits.h>

// Global mesh ID counter
uint32_t manifold_mesh_id_counter = 1;

uint32_t manifold_reserve_ids(uint32_t n) {
  uint32_t start = manifold_mesh_id_counter;
  manifold_mesh_id_counter += n;
  return start;
}

void manifold_impl_init(ManifoldImpl *impl) {
  impl->bBox = manifold_box_empty();
  impl->epsilon = -1.0;
  impl->tolerance = -1.0;
  impl->numProp = 0;
  impl->status = MANIFOLD_ERROR_NO_ERROR;
  impl->vertPos = (ManifoldVecVec3)MANIFOLD_VEC_INIT;
  impl->halfedge = (ManifoldVecHalfedge)MANIFOLD_VEC_INIT;
  impl->properties = (ManifoldVecDouble)MANIFOLD_VEC_INIT;
  impl->vertNormal = (ManifoldVecVec3)MANIFOLD_VEC_INIT;
  impl->faceNormal = (ManifoldVecVec3)MANIFOLD_VEC_INIT;
  impl->halfedgeTangent = (ManifoldVecVec4)MANIFOLD_VEC_INIT;
  impl->meshRelation.originalID = -1;
  impl->meshRelation.meshIDtransform = (ManifoldVecMeshIDEntry)MANIFOLD_VEC_INIT;
  impl->meshRelation.triRef = (ManifoldVecTriRef)MANIFOLD_VEC_INIT;
  impl->collider.nodeBBox = (ManifoldVecBox)MANIFOLD_VEC_INIT;
  impl->collider.nodeParent = (ManifoldVecInt)MANIFOLD_VEC_INIT;
  impl->collider.internalChildren = (ManifoldVecIntPair)MANIFOLD_VEC_INIT;
  impl->colliderBuilt = false;
}

void manifold_impl_free(ManifoldImpl *impl) {
  vec_vec3_free(&impl->vertPos);
  vec_halfedge_free(&impl->halfedge);
  vec_double_free(&impl->properties);
  vec_vec3_free(&impl->vertNormal);
  vec_vec3_free(&impl->faceNormal);
  vec_vec4_free(&impl->halfedgeTangent);
  vec_meshid_free(&impl->meshRelation.meshIDtransform);
  vec_triref_free(&impl->meshRelation.triRef);
  manifold_collider_free(&impl->collider);
}

void manifold_impl_make_empty(ManifoldImpl *impl, ManifoldError status) {
  impl->status = status;
  vec_vec3_clear(&impl->vertPos);
  vec_halfedge_clear(&impl->halfedge);
  vec_double_clear(&impl->properties);
  vec_vec3_clear(&impl->vertNormal);
  vec_vec3_clear(&impl->faceNormal);
  vec_vec4_clear(&impl->halfedgeTangent);
  vec_meshid_clear(&impl->meshRelation.meshIDtransform);
  vec_triref_clear(&impl->meshRelation.triRef);
}

void manifold_impl_calculate_bbox(ManifoldImpl *impl) {
  impl->bBox = manifold_box_empty();
  for (size_t i = 0; i < impl->vertPos.len; i++) {
    manifold_box_union_point(&impl->bBox, impl->vertPos.data[i]);
  }
}

void manifold_impl_set_epsilon(ManifoldImpl *impl, double minEpsilon,
                                bool useSingle) {
  if (minEpsilon < 0) {
    double precision = useSingle ? (double)1.1920929e-07f : MANIFOLD_PRECISION;
    impl->epsilon = fmax(precision * manifold_box_scale(impl->bBox), 0.0);
  } else {
    impl->epsilon = fmax(minEpsilon, MANIFOLD_PRECISION * manifold_box_scale(impl->bBox));
  }
  if (!isfinite(impl->epsilon)) impl->epsilon = -1.0;

  if (impl->tolerance < 0) {
    impl->tolerance = impl->epsilon;
  } else {
    impl->tolerance = fmax(impl->tolerance, impl->epsilon);
  }
}

bool manifold_impl_is_finite(const ManifoldImpl *impl) {
  for (size_t i = 0; i < impl->vertPos.len; i++) {
    if (!vec3_isfinite(impl->vertPos.data[i])) return false;
  }
  return true;
}

bool manifold_impl_is_manifold(const ManifoldImpl *impl) {
  size_t numTri = manifold_impl_num_tri(impl);
  if (numTri == 0) return true;

  // Check: every edge appears exactly twice with opposite orientations
  for (size_t i = 0; i < impl->halfedge.len; i++) {
    const ManifoldHalfedge *he = &impl->halfedge.data[i];
    if (he->pairedHalfedge < 0 || (size_t)he->pairedHalfedge >= impl->halfedge.len)
      return false;
    const ManifoldHalfedge *paired = &impl->halfedge.data[he->pairedHalfedge];
    if (he->startVert != paired->endVert || he->endVert != paired->startVert)
      return false;
  }
  return true;
}

bool manifold_impl_is_2manifold(const ManifoldImpl *impl) {
  // Check that each edge has exactly one pair
  if (!manifold_impl_is_manifold(impl)) return false;

  // Check that each vertex has a consistent fan of triangles
  size_t numVert = manifold_impl_num_vert(impl);
  for (size_t v = 0; v < numVert; v++) {
    // Find a halfedge starting at this vertex
    int startEdge = -1;
    for (size_t i = 0; i < impl->halfedge.len; i++) {
      if (impl->halfedge.data[i].startVert == (int)v) {
        startEdge = (int)i;
        break;
      }
    }
    if (startEdge < 0) continue; // unreferenced vertex

    // Walk around the vertex
    int current = startEdge;
    int count = 0;
    do {
      current = manifold_next_halfedge(impl->halfedge.data[current].pairedHalfedge);
      count++;
      if (count > (int)impl->halfedge.len) return false; // infinite loop
    } while (current != startEdge);
  }
  return true;
}

// Sort helpers
static int cmp_int_by_u32(const void *a_ptr, const void *b_ptr, void *ctx) {
  const uint32_t *morton = (const uint32_t *)ctx;
  int a = *(const int *)a_ptr;
  int b = *(const int *)b_ptr;
  if (morton[a] < morton[b]) return -1;
  if (morton[a] > morton[b]) return 1;
  return (a < b) ? -1 : (a > b) ? 1 : 0; // stable by index
}

// Platform-specific qsort_r
#ifdef __APPLE__
// macOS qsort_r has different signature
static int cmp_int_by_u32_apple(void *ctx, const void *a_ptr, const void *b_ptr) {
  return cmp_int_by_u32(a_ptr, b_ptr, ctx);
}
#endif

static void sort_int_by_u32(int *arr, size_t n, uint32_t *morton) {
#ifdef __APPLE__
  qsort_r(arr, n, sizeof(int), morton, cmp_int_by_u32_apple);
#else
  qsort_r(arr, n, sizeof(int), cmp_int_by_u32, morton);
#endif
}

void manifold_impl_sort_verts(ManifoldImpl *impl) {
  size_t numVert = manifold_impl_num_vert(impl);
  if (numVert == 0) return;

  ManifoldVecU32 vertMorton = vec_u32_create_n(numVert);
  for (size_t i = 0; i < numVert; i++) {
    ManifoldVec3 pos = impl->vertPos.data[i];
    if (isnan(pos.x)) {
      vertMorton.data[i] = 0xFFFFFFFFu;
    } else {
      vertMorton.data[i] = manifold_morton_code(pos, impl->bBox);
    }
  }

  ManifoldVecInt vertNew2Old = vec_int_create_n(numVert);
  vec_int_sequence(&vertNew2Old);

  sort_int_by_u32(vertNew2Old.data, numVert, vertMorton.data);

  manifold_impl_reindex_verts(impl, &vertNew2Old, numVert);

  // Remove NaN verts (sorted to end)
  size_t newNumVert = numVert;
  for (size_t i = 0; i < numVert; i++) {
    if (vertMorton.data[vertNew2Old.data[i]] == 0xFFFFFFFFu) {
      newNumVert = i;
      break;
    }
  }

  vec_int_resize(&vertNew2Old, newNumVert);
  vec_vec3_permute(&impl->vertPos, &vertNew2Old);

  if (impl->vertNormal.len == numVert) {
    vec_vec3_permute(&impl->vertNormal, &vertNew2Old);
  }

  vec_u32_free(&vertMorton);
  vec_int_free(&vertNew2Old);
}

void manifold_impl_reindex_verts(ManifoldImpl *impl,
                                  const ManifoldVecInt *vertNew2Old,
                                  size_t oldNumVert) {
  ManifoldVecInt vertOld2New = vec_int_create_n(oldNumVert);
  for (size_t i = 0; i < vertNew2Old->len; i++) {
    vertOld2New.data[vertNew2Old->data[i]] = (int)i;
  }

  bool hasProp = impl->numProp > 0;
  for (size_t i = 0; i < impl->halfedge.len; i++) {
    ManifoldHalfedge *edge = &impl->halfedge.data[i];
    if (edge->startVert < 0) continue;
    if ((size_t)edge->startVert >= oldNumVert ||
        (size_t)edge->endVert >= oldNumVert) continue;
    edge->startVert = vertOld2New.data[edge->startVert];
    edge->endVert = vertOld2New.data[edge->endVert];
    if (!hasProp) {
      edge->propVert = edge->startVert;
    }
  }

  vec_int_free(&vertOld2New);
}

void manifold_impl_compact_props(ManifoldImpl *impl) {
  if (impl->numProp == 0) return;

  int numProp = impl->numProp;
  size_t numVerts = impl->properties.len / numProp;
  ManifoldVecInt keep = vec_int_create_fill(numVerts, 0);

  for (size_t i = 0; i < impl->halfedge.len; i++) {
    int pv = impl->halfedge.data[i].propVert;
    if (pv >= 0 && (size_t)pv < numVerts) keep.data[pv] = 1;
  }

  ManifoldVecInt propOld2New = vec_int_create_fill(numVerts + 1, 0);
  for (size_t i = 0; i < numVerts; i++) {
    propOld2New.data[i + 1] = propOld2New.data[i] + keep.data[i];
  }

  ManifoldVecDouble oldProp = vec_double_copy(&impl->properties);
  int numVertsNew = propOld2New.data[numVerts];
  vec_double_resize(&impl->properties, (size_t)(numProp * numVertsNew));

  for (size_t oldIdx = 0; oldIdx < numVerts; oldIdx++) {
    if (keep.data[oldIdx] == 0) continue;
    for (int p = 0; p < numProp; p++) {
      impl->properties.data[propOld2New.data[oldIdx] * numProp + p] =
          oldProp.data[oldIdx * numProp + p];
    }
  }

  for (size_t i = 0; i < impl->halfedge.len; i++) {
    int pv = impl->halfedge.data[i].propVert;
    if (pv >= 0 && (size_t)pv < numVerts)
      impl->halfedge.data[i].propVert = propOld2New.data[pv];
  }

  vec_int_free(&keep);
  vec_int_free(&propOld2New);
  vec_double_free(&oldProp);
}

void manifold_impl_get_face_box_morton(const ManifoldImpl *impl,
                                        ManifoldVecBox *faceBox,
                                        ManifoldVecU32 *faceMorton) {
  size_t numTri = manifold_impl_num_tri(impl);
  ManifoldBox emptyBox = manifold_box_empty();
  *faceBox = vec_box_create_fill(numTri, emptyBox);
  *faceMorton = vec_u32_create_n(numTri);

  for (size_t face = 0; face < numTri; face++) {
    if (impl->halfedge.data[3 * face].pairedHalfedge < 0) {
      faceMorton->data[face] = 0xFFFFFFFFu;
      continue;
    }
    ManifoldVec3 center = manifold_vec3(0, 0, 0);
    bool valid = true;
    for (int i = 0; i < 3; i++) {
      int sv = impl->halfedge.data[3 * face + i].startVert;
      if (sv < 0 || (size_t)sv >= impl->vertPos.len) { valid = false; break; }
      ManifoldVec3 pos = impl->vertPos.data[sv];
      center = vec3_add(center, pos);
      manifold_box_union_point(&faceBox->data[face], pos);
    }
    if (!valid) { faceMorton->data[face] = 0xFFFFFFFFu; continue; }
    center = vec3_scale(center, 1.0 / 3.0);
    faceMorton->data[face] = manifold_morton_code(center, impl->bBox);
  }
}

void manifold_impl_sort_faces(ManifoldImpl *impl, ManifoldVecBox *faceBox,
                               ManifoldVecU32 *faceMorton) {
  size_t numTri = manifold_impl_num_tri(impl);
  ManifoldVecInt faceNew2Old = vec_int_create_n(numTri);
  vec_int_sequence(&faceNew2Old);

  sort_int_by_u32(faceNew2Old.data, numTri, faceMorton->data);

  // Remove deleted tris (pairedHalfedge=-1, sorted to end)
  size_t newNumTri = numTri;
  for (size_t i = 0; i < numTri; i++) {
    if (faceMorton->data[faceNew2Old.data[i]] == 0xFFFFFFFFu) {
      newNumTri = i;
      break;
    }
  }
  vec_int_resize(&faceNew2Old, newNumTri);

  vec_u32_permute(faceMorton, &faceNew2Old);
  vec_box_permute(faceBox, &faceNew2Old);
  manifold_impl_gather_faces(impl, &faceNew2Old);

  vec_int_free(&faceNew2Old);
}

void manifold_impl_gather_faces(ManifoldImpl *impl,
                                 const ManifoldVecInt *faceNew2Old) {
  size_t numTri = faceNew2Old->len;
  size_t oldNumTri = manifold_impl_num_tri(impl);

  if (impl->meshRelation.triRef.len == oldNumTri) {
    vec_triref_permute(&impl->meshRelation.triRef, faceNew2Old);
  }
  if (impl->faceNormal.len == oldNumTri) {
    vec_vec3_permute(&impl->faceNormal, faceNew2Old);
  }

  ManifoldVecHalfedge oldHalfedge = impl->halfedge;
  impl->halfedge = vec_halfedge_create_n(3 * numTri);

  ManifoldVecVec4 oldTangent = impl->halfedgeTangent;
  if (oldTangent.len != 0) {
    impl->halfedgeTangent = vec_vec4_create_n(3 * numTri);
  }

  ManifoldVecInt faceOld2New = vec_int_create_n(oldNumTri);
  for (size_t i = 0; i < numTri; i++) {
    faceOld2New.data[faceNew2Old->data[i]] = (int)i;
  }

  for (size_t newFace = 0; newFace < numTri; newFace++) {
    int oldFace = faceNew2Old->data[newFace];
    for (int i = 0; i < 3; i++) {
      int oldEdge = 3 * oldFace + i;
      ManifoldHalfedge edge = oldHalfedge.data[oldEdge];
      int pairedFace = edge.pairedHalfedge / 3;
      int offset = edge.pairedHalfedge - 3 * pairedFace;
      edge.pairedHalfedge = 3 * faceOld2New.data[pairedFace] + offset;
      int newEdge = 3 * (int)newFace + i;
      impl->halfedge.data[newEdge] = edge;
      if (oldTangent.len != 0) {
        impl->halfedgeTangent.data[newEdge] = oldTangent.data[oldEdge];
      }
    }
  }

  vec_halfedge_free(&oldHalfedge);
  vec_vec4_free(&oldTangent);
  vec_int_free(&faceOld2New);
}

void manifold_impl_reorder_halfedges(ManifoldImpl *impl) {
  size_t numTri = manifold_impl_num_tri(impl);

  // Step 1: reorder within same face so smallest startVert is first
  for (size_t tri = 0; tri < numTri; tri++) {
    ManifoldHalfedge face[3];
    face[0] = impl->halfedge.data[tri * 3];
    face[1] = impl->halfedge.data[tri * 3 + 1];
    face[2] = impl->halfedge.data[tri * 3 + 2];
    if (face[0].startVert < 0) continue;
    int index = 0;
    for (int i = 1; i < 3; i++) {
      if (face[i].startVert < face[index].startVert) index = i;
    }
    for (int i = 0; i < 3; i++) {
      impl->halfedge.data[tri * 3 + i] = face[(index + i) % 3];
    }
  }

  // Step 2: fix paired halfedge references
  for (size_t tri = 0; tri < numTri; tri++) {
    for (int i = 0; i < 3; i++) {
      ManifoldHalfedge *curr = &impl->halfedge.data[tri * 3 + i];
      if (curr->startVert < 0) continue;
      if (curr->pairedHalfedge < 0) continue;
      int oppositeFace = curr->pairedHalfedge / 3;
      if ((size_t)oppositeFace >= numTri) { curr->pairedHalfedge = -1; continue; }
      int index = -1;
      for (int j = 0; j < 3; j++) {
        if (curr->startVert == impl->halfedge.data[oppositeFace * 3 + j].endVert &&
            curr->endVert == impl->halfedge.data[oppositeFace * 3 + j].startVert)
          index = j;
      }
      if (index < 0) { curr->pairedHalfedge = -1; continue; }
      curr->pairedHalfedge = oppositeFace * 3 + index;
    }
  }
}

void manifold_impl_sort_geometry(ManifoldImpl *impl) {
  if (impl->halfedge.len == 0) return;

  manifold_impl_sort_verts(impl);

  ManifoldVecBox faceBox;
  ManifoldVecU32 faceMorton;
  manifold_impl_get_face_box_morton(impl, &faceBox, &faceMorton);
  manifold_impl_sort_faces(impl, &faceBox, &faceMorton);

  if (impl->halfedge.len == 0) {
    vec_box_free(&faceBox);
    vec_u32_free(&faceMorton);
    return;
  }

  // Build collider
  if (impl->colliderBuilt) manifold_collider_free(&impl->collider);
  impl->collider = manifold_collider_create(
      faceBox.data, faceMorton.data, faceBox.len);
  impl->colliderBuilt = true;

  manifold_impl_compact_props(impl);

  vec_box_free(&faceBox);
  vec_u32_free(&faceMorton);
}

void manifold_impl_calculate_vert_normals(ManifoldImpl *impl) {
  size_t numVert = manifold_impl_num_vert(impl);
  ManifoldVec3 zero = manifold_vec3(0, 0, 0);
  vec_vec3_free(&impl->vertNormal);
  impl->vertNormal = vec_vec3_create_fill(numVert, zero);

  // Build vertHalfedgeMap: for each vert, store the minimum halfedge index
  int *vertHalfedgeMap = (int *)malloc(numVert * sizeof(int));
  for (size_t i = 0; i < numVert; i++)
    vertHalfedgeMap[i] = INT_MAX;
  for (size_t i = 0; i < impl->halfedge.len; i++) {
    int v = impl->halfedge.data[i].startVert;
    if (v >= 0 && (size_t)v < numVert && (int)i < vertHalfedgeMap[v])
      vertHalfedgeMap[v] = (int)i;
  }

  for (size_t vert = 0; vert < numVert; vert++) {
    int firstEdge = vertHalfedgeMap[vert];
    if (firstEdge == INT_MAX) {
      impl->vertNormal.data[vert] = zero;
      continue;
    }
    ManifoldVec3 normal = zero;
    // ForVert loop: iterate over all halfedges starting at this vertex
    int edge = firstEdge;
    int iterations = 0;
    do {
      int startV = impl->halfedge.data[edge].startVert;
      int endV = impl->halfedge.data[edge].endVert;
      int nextEdge = manifold_next_halfedge(edge);
      int thirdV = impl->halfedge.data[nextEdge].endVert;

      ManifoldVec3 currEdge = vec3_normalize(
          vec3_sub(impl->vertPos.data[endV], impl->vertPos.data[startV]));
      ManifoldVec3 prevEdge = vec3_normalize(
          vec3_sub(impl->vertPos.data[startV], impl->vertPos.data[thirdV]));

      if (isfinite(currEdge.x) && isfinite(prevEdge.x)) {
        double dot = -vec3_dot(prevEdge, currEdge);
        double phi = dot >= 1.0 ? 0.0
                   : (dot <= -1.0 ? MANIFOLD_PI : acos(dot));
        normal = vec3_add(normal,
            vec3_scale(impl->faceNormal.data[edge / 3], phi));
      }

      edge = manifold_next_halfedge(
          impl->halfedge.data[edge].pairedHalfedge);
      if (++iterations > (int)impl->halfedge.len) break;
    } while (edge != firstEdge);

    impl->vertNormal.data[vert] = manifold_safe_normalize(normal);
  }

  free(vertHalfedgeMap);
}

void manifold_impl_set_normals_and_coplanar(ManifoldImpl *impl) {
  size_t numTri = manifold_impl_num_tri(impl);
  ManifoldVec3 zero = manifold_vec3(0, 0, 0);
  impl->faceNormal = vec_vec3_create_fill(numTri, zero);

  // Compute face normals and area priorities
  typedef struct { double area2; int tri; } TriPriority;
  TriPriority *triPriority = (TriPriority *)malloc(numTri * sizeof(TriPriority));

  for (size_t tri = 0; tri < numTri; tri++) {
    if (impl->meshRelation.triRef.len > tri)
      impl->meshRelation.triRef.data[tri].coplanarID = -1;
    int sv0 = impl->halfedge.data[3 * tri].startVert;
    int sv1 = impl->halfedge.data[3 * tri + 1].startVert;
    int sv2 = impl->halfedge.data[3 * tri + 2].startVert;
    if (sv0 < 0 || (size_t)sv0 >= impl->vertPos.len ||
        sv1 < 0 || (size_t)sv1 >= impl->vertPos.len ||
        sv2 < 0 || (size_t)sv2 >= impl->vertPos.len) {
      triPriority[tri] = (TriPriority){0, (int)tri};
      continue;
    }
    ManifoldVec3 v0 = impl->vertPos.data[sv0];
    ManifoldVec3 v1 = impl->vertPos.data[sv1];
    ManifoldVec3 v2 = impl->vertPos.data[sv2];
    ManifoldVec3 n = vec3_cross(vec3_sub(v1, v0), vec3_sub(v2, v0));
    impl->faceNormal.data[tri] = manifold_safe_normalize(n);
    triPriority[tri] = (TriPriority){vec3_dot(n, n), (int)tri};
  }

  // Sort by area (largest first) — O(n log n) merge sort for stability
  {
    TriPriority *tmp = (TriPriority *)malloc(numTri * sizeof(TriPriority));
    for (size_t width = 1; width < numTri; width *= 2) {
      for (size_t i = 0; i < numTri; i += 2 * width) {
        size_t mid = i + width < numTri ? i + width : numTri;
        size_t end = i + 2 * width < numTri ? i + 2 * width : numTri;
        size_t l = i, r = mid, k = i;
        while (l < mid && r < end) {
          if (triPriority[l].area2 > triPriority[r].area2 ||
              (triPriority[l].area2 == triPriority[r].area2 &&
               triPriority[l].tri <= triPriority[r].tri))
            tmp[k++] = triPriority[l++];
          else
            tmp[k++] = triPriority[r++];
        }
        while (l < mid) tmp[k++] = triPriority[l++];
        while (r < end) tmp[k++] = triPriority[r++];
      }
      memcpy(triPriority, tmp, numTri * sizeof(TriPriority));
    }
    free(tmp);
  }

  // BFS to assign coplanar groups
  int *stack = (int *)malloc((numTri * 3 + 1) * sizeof(int));
  for (size_t idx = 0; idx < numTri; idx++) {
    int tp_tri = triPriority[idx].tri;
    if (impl->meshRelation.triRef.len <= (size_t)tp_tri) continue;
    if (impl->meshRelation.triRef.data[tp_tri].coplanarID >= 0) continue;

    impl->meshRelation.triRef.data[tp_tri].coplanarID = tp_tri;
    if (impl->halfedge.data[3 * tp_tri].startVert < 0) continue;

    ManifoldVec3 base = impl->vertPos.data[impl->halfedge.data[3 * tp_tri].startVert];
    ManifoldVec3 normal = impl->faceNormal.data[tp_tri];
    int stackLen = 0;
    stack[stackLen++] = 3 * tp_tri;
    stack[stackLen++] = 3 * tp_tri + 1;
    stack[stackLen++] = 3 * tp_tri + 2;

    while (stackLen > 0) {
      int he = impl->halfedge.data[stack[stackLen - 1]].pairedHalfedge;
      stackLen--;
      if (he < 0) continue;
      int h = manifold_next_halfedge(he);
      size_t tri = (size_t)h / 3;
      if (tri >= numTri) continue;
      if (impl->meshRelation.triRef.data[tri].coplanarID >= 0) continue;

      int ev = impl->halfedge.data[h].endVert;
      if (ev < 0 || (size_t)ev >= impl->vertPos.len) continue;
      ManifoldVec3 v = impl->vertPos.data[ev];
      if (fabs(vec3_dot(vec3_sub(v, base), normal)) < impl->tolerance) {
        impl->meshRelation.triRef.data[tri].coplanarID = tp_tri;
        impl->faceNormal.data[tri] = normal;

        if (stackLen == 0 ||
            h != impl->halfedge.data[stack[stackLen - 1]].pairedHalfedge) {
          stack[stackLen++] = h;
        } else {
          stackLen--;
        }
        int hNext = manifold_next_halfedge(h);
        stack[stackLen++] = hNext;
      }
    }
  }

  free(triPriority);
  free(stack);
  manifold_impl_calculate_vert_normals(impl);
}

// Sort context for CreateHalfedges edge sorting
typedef struct {
  const uint64_t *keys;
} EdgeSortCtx;

static int edge_sort_compare(const void *a, const void *b, void *ctx) {
  EdgeSortCtx *c = (EdgeSortCtx *)ctx;
  int ia = *(const int *)a;
  int ib = *(const int *)b;
  uint64_t ka = c->keys[ia];
  uint64_t kb = c->keys[ib];
  if (ka < kb) return -1;
  if (ka > kb) return 1;
  // Stable sort tiebreaker: preserve original order
  return ia < ib ? -1 : (ia > ib ? 1 : 0);
}

// Portable mergesort for stable sorting (qsort is not guaranteed stable)
static void merge_sort_ids(int *arr, int *tmp, size_t n, const uint64_t *keys) {
  if (n <= 1) return;
  size_t mid = n / 2;
  merge_sort_ids(arr, tmp, mid, keys);
  merge_sort_ids(arr + mid, tmp + mid, n - mid, keys);
  // Merge
  size_t i = 0, j = mid, k = 0;
  while (i < mid && j < n) {
    uint64_t ki = keys[arr[i]], kj = keys[arr[j]];
    if (ki < kj || (ki == kj && arr[i] < arr[j]))
      tmp[k++] = arr[i++];
    else
      tmp[k++] = arr[j++];
  }
  while (i < mid) tmp[k++] = arr[i++];
  while (j < n) tmp[k++] = arr[j++];
  memcpy(arr, tmp, n * sizeof(int));
}

void manifold_impl_create_halfedges(ManifoldImpl *impl,
                                     const ManifoldVecIVec3 *triProp,
                                     const ManifoldVecIVec3 *triVert) {
  const size_t numTri = triProp->len;
  const size_t numHalfedge = numTri * 3;

  vec_halfedge_free(&impl->halfedge);
  impl->halfedge = vec_halfedge_create_n(numHalfedge);

  const ManifoldVecIVec3 *vertSource =
      (triVert && triVert->len > 0) ? triVert : triProp;

  // Build halfedges and edge keys for sorting
  uint64_t *edgeKeys = (uint64_t *)malloc(numHalfedge * sizeof(uint64_t));
  int *ids = (int *)malloc(numHalfedge * sizeof(int));
  int *tmpBuf = (int *)malloc(numHalfedge * sizeof(int));

  for (size_t tri = 0; tri < numTri; tri++) {
    ManifoldIVec3 verts = vertSource->data[tri];
    ManifoldIVec3 props = triProp->data[tri];
    for (int i = 0; i < 3; i++) {
      int j = (i + 1) % 3;
      int e = (int)(3 * tri + i);
      int v0 = ivec3_get(verts, i);
      int v1 = ivec3_get(verts, j);
      impl->halfedge.data[e].startVert = v0;
      impl->halfedge.data[e].endVert = v1;
      impl->halfedge.data[e].pairedHalfedge = -1;
      impl->halfedge.data[e].propVert = ivec3_get(props, i);
      // Key: min(v0,v1) in high 32 bits, max(v0,v1) in low 32 bits
      int mn = v0 < v1 ? v0 : v1;
      int mx = v0 < v1 ? v1 : v0;
      edgeKeys[e] = ((uint64_t)mn << 32) | (uint64_t)mx;
    }
  }

  for (int i = 0; i < (int)numHalfedge; i++) ids[i] = i;
  merge_sort_ids(ids, tmpBuf, numHalfedge, edgeKeys);
  free(tmpBuf);
  free(edgeKeys);

  // Group by undirected edge, detect opposite faces, then pair reverses
  // kRemovedHalfedge = -1 (already default for pairedHalfedge)
  // We use -2 to mark halfedges of removed opposite face pairs
  const int kRemovedHalfedge = -2;

  size_t i = 0;
  while (i < numHalfedge) {
    // Find group end
    size_t groupStart = i;
    ManifoldHalfedge *h0 = &impl->halfedge.data[ids[i]];
    int mn0 = h0->startVert < h0->endVert ? h0->startVert : h0->endVert;
    int mx0 = h0->startVert < h0->endVert ? h0->endVert : h0->startVert;
    i++;
    while (i < numHalfedge) {
      ManifoldHalfedge *h = &impl->halfedge.data[ids[i]];
      int mn = h->startVert < h->endVert ? h->startVert : h->endVert;
      int mx = h->startVert < h->endVert ? h->endVert : h->startVert;
      if (mn != mn0 || mx != mx0) break;
      i++;
    }
    size_t groupEnd = i;

    // Detect opposite face pairs: two reverse halfedges whose triangles share
    // all 3 vertices (the third vertex of each triangle is the same).
    for (size_t a = groupStart; a < groupEnd; a++) {
      int ea = ids[a];
      if (impl->halfedge.data[ea].pairedHalfedge == kRemovedHalfedge) continue;
      for (size_t b = a + 1; b < groupEnd; b++) {
        int eb = ids[b];
        if (impl->halfedge.data[eb].pairedHalfedge == kRemovedHalfedge)
          continue;
        if (impl->halfedge.data[ea].startVert ==
                impl->halfedge.data[eb].endVert &&
            impl->halfedge.data[ea].endVert ==
                impl->halfedge.data[eb].startVert) {
          // Check if third vertices match (opposite face)
          int nextA = 3 * (ea / 3) + (ea % 3 + 1) % 3;
          int nextB = 3 * (eb / 3) + (eb % 3 + 1) % 3;
          if (impl->halfedge.data[nextA].endVert ==
              impl->halfedge.data[nextB].endVert) {
            // Mark all 3 halfedges of both triangles as removed
            int triA = ea / 3, triB = eb / 3;
            for (int k = 0; k < 3; k++) {
              impl->halfedge.data[3 * triA + k].pairedHalfedge =
                  kRemovedHalfedge;
              impl->halfedge.data[3 * triB + k].pairedHalfedge =
                  kRemovedHalfedge;
            }
            break;
          }
        }
      }
    }

    // Pair reverses within the group (skip removed halfedges)
    for (size_t a = groupStart; a < groupEnd; a++) {
      int ea = ids[a];
      if (impl->halfedge.data[ea].pairedHalfedge != -1) continue;
      for (size_t b = a + 1; b < groupEnd; b++) {
        int eb = ids[b];
        if (impl->halfedge.data[eb].pairedHalfedge != -1) continue;
        if (impl->halfedge.data[ea].startVert ==
                impl->halfedge.data[eb].endVert &&
            impl->halfedge.data[ea].endVert ==
                impl->halfedge.data[eb].startVert) {
          impl->halfedge.data[ea].pairedHalfedge = eb;
          impl->halfedge.data[eb].pairedHalfedge = ea;
          break;
        }
      }
    }
  }

  // Set removed halfedges to {-1, -1, -1}
  for (size_t e = 0; e < numHalfedge; e++) {
    if (impl->halfedge.data[e].pairedHalfedge == kRemovedHalfedge) {
      impl->halfedge.data[e].startVert = -1;
      impl->halfedge.data[e].endVert = -1;
      impl->halfedge.data[e].pairedHalfedge = -1;
    }
  }

  free(ids);
}

void manifold_impl_remove_unreferenced_verts(ManifoldImpl *impl) {
  size_t numVert = manifold_impl_num_vert(impl);
  if (numVert == 0) return;

  ManifoldVecBool referenced = vec_bool_create_fill(numVert, false);
  for (size_t i = 0; i < impl->halfedge.len; i++) {
    int sv = impl->halfedge.data[i].startVert;
    if (sv >= 0 && (size_t)sv < numVert) referenced.data[sv] = true;
    int ev = impl->halfedge.data[i].endVert;
    if (ev >= 0 && (size_t)ev < numVert) referenced.data[ev] = true;
  }

  // Build compaction map
  ManifoldVecInt old2new = vec_int_create_fill(numVert, -1);
  int newCount = 0;
  for (size_t i = 0; i < numVert; i++) {
    if (referenced.data[i]) {
      old2new.data[i] = newCount++;
    }
  }

  if ((size_t)newCount == numVert) {
    vec_bool_free(&referenced);
    vec_int_free(&old2new);
    return;
  }

  // Compact vertex positions
  ManifoldVecVec3 newPos = vec_vec3_create_n((size_t)newCount);
  for (size_t i = 0; i < numVert; i++) {
    if (old2new.data[i] >= 0) {
      newPos.data[old2new.data[i]] = impl->vertPos.data[i];
    }
  }
  vec_vec3_free(&impl->vertPos);
  impl->vertPos = newPos;

  // Compact vertex normals if present
  if (impl->vertNormal.len == numVert) {
    ManifoldVecVec3 newNorm = vec_vec3_create_n((size_t)newCount);
    for (size_t i = 0; i < numVert; i++) {
      if (old2new.data[i] >= 0) {
        newNorm.data[old2new.data[i]] = impl->vertNormal.data[i];
      }
    }
    vec_vec3_free(&impl->vertNormal);
    impl->vertNormal = newNorm;
  }

  // Update halfedge references
  bool hasProp = impl->numProp > 0;
  for (size_t i = 0; i < impl->halfedge.len; i++) {
    ManifoldHalfedge *he = &impl->halfedge.data[i];
    if (he->startVert >= 0 && (size_t)he->startVert < old2new.len)
      he->startVert = old2new.data[he->startVert];
    if (he->endVert >= 0 && (size_t)he->endVert < old2new.len)
      he->endVert = old2new.data[he->endVert];
    if (!hasProp) he->propVert = he->startVert;
  }

  vec_bool_free(&referenced);
  vec_int_free(&old2new);
}

void manifold_impl_dedupe_prop_verts(ManifoldImpl *impl) {
  // Simple version: just ensure consistency
  if (impl->numProp == 0) return;
  // More sophisticated deduplication can be added later
}

void manifold_impl_initialize_original(ManifoldImpl *impl) {
  size_t numTri = manifold_impl_num_tri(impl);
  uint32_t meshID = manifold_reserve_ids(1);

  ManifoldRelation rel;
  rel.originalID = (int)meshID;
  rel.transform = mat3x4_identity();
  rel.backSide = false;

  impl->meshRelation.originalID = (int)meshID;
  manifold_meshrelation_insert(&impl->meshRelation, (int)meshID, rel);

  vec_triref_resize(&impl->meshRelation.triRef, numTri);
  for (size_t tri = 0; tri < numTri; tri++) {
    ManifoldTriRef *ref = &impl->meshRelation.triRef.data[tri];
    ref->meshID = (int)meshID;
    ref->originalID = (int)meshID;
    ref->faceID = -1;
    ref->coplanarID = (int)tri;
  }
}

// ============== Edge operations (ported from edge_op.cpp) ==============

static inline ManifoldIVec3 tri_of(int edge) {
  ManifoldIVec3 triEdge;
  triEdge.x = edge;
  triEdge.y = manifold_next_halfedge(triEdge.x);
  triEdge.z = manifold_next_halfedge(triEdge.y);
  return triEdge;
}

static inline bool is01_longest(ManifoldVec2 v0, ManifoldVec2 v1,
                                 ManifoldVec2 v2) {
  ManifoldVec2 e0 = vec2_sub(v1, v0);
  ManifoldVec2 e1 = vec2_sub(v2, v1);
  ManifoldVec2 e2 = vec2_sub(v0, v2);
  double l0 = vec2_dot(e0, e0), l1 = vec2_dot(e1, e1), l2 = vec2_dot(e2, e2);
  return l0 > l1 && l0 > l2;
}

// Orbit CW around startVert of halfedge, calling fn(edge) for each edge
// that starts at the same vertex. Returns to start.
typedef void (*ForVertFn)(int edge, void *ctx);

MANIFOLD_UNUSED static void impl_for_vert(const ManifoldImpl *impl, int startEdge,
                           ForVertFn fn, void *ctx) {
  int current = startEdge;
  do {
    fn(current, ctx);
    current = manifold_next_halfedge(impl->halfedge.data[current].pairedHalfedge);
  } while (current != startEdge);
}

static void impl_pair_up(ManifoldImpl *impl, int edge0, int edge1) {
  impl->halfedge.data[edge0].pairedHalfedge = edge1;
  impl->halfedge.data[edge1].pairedHalfedge = edge0;
}

static void impl_update_vert(ManifoldImpl *impl, int vert,
                              int startEdge, int endEdge) {
  int current = startEdge;
  int maxIter = (int)impl->halfedge.len + 2;
  int iter = 0;
  while (current != endEdge && iter++ < maxIter) {
    if (current < 0 || (size_t)current >= impl->halfedge.len) break;
    impl->halfedge.data[current].endVert = vert;
    current = manifold_next_halfedge(current);
    if ((size_t)current >= impl->halfedge.len) break;
    impl->halfedge.data[current].startVert = vert;
    int paired = impl->halfedge.data[current].pairedHalfedge;
    if (paired < 0 || (size_t)paired >= impl->halfedge.len) break;
    current = paired;
    if (current == startEdge) break;  // infinite loop detection
  }
}

static void impl_collapse_tri(ManifoldImpl *impl, ManifoldIVec3 triEdge) {
  if (impl->halfedge.data[triEdge.y].pairedHalfedge == -1) return;
  int pair1 = impl->halfedge.data[triEdge.y].pairedHalfedge;
  int pair2 = impl->halfedge.data[triEdge.z].pairedHalfedge;
  impl->halfedge.data[pair1].pairedHalfedge = pair2;
  impl->halfedge.data[pair2].pairedHalfedge = pair1;
  int edges[3] = {triEdge.x, triEdge.y, triEdge.z};
  for (int i = 0; i < 3; i++) {
    int propV = impl->halfedge.data[edges[i]].propVert;
    impl->halfedge.data[edges[i]] = (ManifoldHalfedge){-1, -1, -1, propV};
  }
}

static void impl_remove_if_folded(ManifoldImpl *impl, int edge) {
  if (edge < 0 || (size_t)edge >= impl->halfedge.len) return;
  int paired = impl->halfedge.data[edge].pairedHalfedge;
  if (paired < 0 || (size_t)paired >= impl->halfedge.len) return;
  ManifoldIVec3 tri0edge = tri_of(edge);
  ManifoldIVec3 tri1edge = tri_of(paired);
  if (impl->halfedge.data[tri0edge.y].pairedHalfedge == -1) return;
  if (impl->halfedge.data[tri0edge.y].endVert ==
      impl->halfedge.data[tri1edge.y].endVert) {
    if (impl->halfedge.data[tri0edge.y].pairedHalfedge == tri1edge.z) {
      if (impl->halfedge.data[tri0edge.z].pairedHalfedge == tri1edge.y) {
        int e3[3] = {tri0edge.x, tri0edge.y, tri0edge.z};
        for (int i = 0; i < 3; i++) {
          ManifoldVec3 nan = {NAN, NAN, NAN};
          impl->vertPos.data[impl->halfedge.data[e3[i]].startVert] = nan;
        }
      } else {
        ManifoldVec3 nan = {NAN, NAN, NAN};
        impl->vertPos.data[impl->halfedge.data[tri0edge.y].startVert] = nan;
      }
    } else {
      if (impl->halfedge.data[tri0edge.z].pairedHalfedge == tri1edge.y) {
        ManifoldVec3 nan = {NAN, NAN, NAN};
        impl->vertPos.data[impl->halfedge.data[tri1edge.y].startVert] = nan;
      }
    }
    impl_pair_up(impl,
                 impl->halfedge.data[tri0edge.y].pairedHalfedge,
                 impl->halfedge.data[tri1edge.z].pairedHalfedge);
    impl_pair_up(impl,
                 impl->halfedge.data[tri0edge.z].pairedHalfedge,
                 impl->halfedge.data[tri1edge.y].pairedHalfedge);
    int e0[3] = {tri0edge.x, tri0edge.y, tri0edge.z};
    int e1[3] = {tri1edge.x, tri1edge.y, tri1edge.z};
    for (int i = 0; i < 3; i++) {
      impl->halfedge.data[e0[i]] = (ManifoldHalfedge){-1, -1, -1, -1};
      impl->halfedge.data[e1[i]] = (ManifoldHalfedge){-1, -1, -1, -1};
    }
  }
}

static void impl_form_loop(ManifoldImpl *impl, int current, int end) {
  if (current < 0 || (size_t)current >= impl->halfedge.len) return;
  if (end < 0 || (size_t)end >= impl->halfedge.len) return;
  int sv = impl->halfedge.data[current].startVert;
  int ev = impl->halfedge.data[current].endVert;
  if (sv < 0 || (size_t)sv >= impl->vertPos.len) return;
  if (ev < 0 || (size_t)ev >= impl->vertPos.len) return;

  int startVert = (int)impl->vertPos.len;
  vec_vec3_push(&impl->vertPos, impl->vertPos.data[sv]);
  int endVert = (int)impl->vertPos.len;
  vec_vec3_push(&impl->vertPos, impl->vertPos.data[ev]);

  int oldMatch = impl->halfedge.data[current].pairedHalfedge;
  int newMatch = impl->halfedge.data[end].pairedHalfedge;
  if (oldMatch < 0 || (size_t)oldMatch >= impl->halfedge.len) return;
  if (newMatch < 0 || (size_t)newMatch >= impl->halfedge.len) return;

  impl_update_vert(impl, startVert, oldMatch, newMatch);
  impl_update_vert(impl, endVert, end, current);

  impl->halfedge.data[current].pairedHalfedge = newMatch;
  impl->halfedge.data[newMatch].pairedHalfedge = current;
  impl->halfedge.data[end].pairedHalfedge = oldMatch;
  impl->halfedge.data[oldMatch].pairedHalfedge = end;

  impl_remove_if_folded(impl, end);
}

static bool impl_collapse_edge(ManifoldImpl *impl, int edge,
                                ManifoldVecInt *edges) {
  ManifoldHalfedge toRemove = impl->halfedge.data[edge];
  if (toRemove.pairedHalfedge < 0) return false;
  if (toRemove.startVert == -1 || toRemove.endVert == -1) return false;

  const int endVert = toRemove.endVert;
  if ((size_t)endVert >= impl->vertPos.len) return false;
  if ((size_t)toRemove.startVert >= impl->vertPos.len) return false;
  ManifoldIVec3 tri0edge = tri_of(edge);
  ManifoldIVec3 tri1edge = tri_of(toRemove.pairedHalfedge);

  ManifoldVec3 pNew = impl->vertPos.data[endVert];
  ManifoldVec3 pOld = impl->vertPos.data[toRemove.startVert];
  ManifoldVec3 delta = vec3_sub(pNew, pOld);
  bool shortEdge = vec3_dot(delta, delta) < impl->epsilon * impl->epsilon;

  int start = impl->halfedge.data[tri1edge.y].pairedHalfedge;
  if (start < 0 || (size_t)start >= impl->halfedge.len) return false;
  // Verify tri1edge.y's pairing is consistent
  if (impl->halfedge.data[start].pairedHalfedge != tri1edge.y) return false;
  int current = tri1edge.z;
  if (!shortEdge) {
    current = start;
    if (current < 0 || (size_t)current >= impl->halfedge.len) return false;
    ManifoldTriRef refCheck =
        impl->meshRelation.triRef.data[toRemove.pairedHalfedge / 3];
    ManifoldVec3 pLast =
        impl->vertPos.data[impl->halfedge.data[tri1edge.y].endVert];
    int safeN = 0;
    while (current != tri1edge.x && safeN++ < (int)impl->halfedge.len) {
      current = manifold_next_halfedge(current);
      if ((size_t)current >= impl->halfedge.len) return false;
      int ev = impl->halfedge.data[current].endVert;
      if (ev < 0 || (size_t)ev >= impl->vertPos.len) return false;
      ManifoldVec3 pNext = impl->vertPos.data[ev];
      int tri = current / 3;
      if ((size_t)tri >= impl->meshRelation.triRef.len) return false;
      ManifoldTriRef ref = impl->meshRelation.triRef.data[tri];
      ManifoldMat2x3 projection =
          manifold_get_axis_aligned_projection(impl->faceNormal.data[tri]);
      if (!manifold_triref_same_face(&ref, &refCheck)) {
        ManifoldTriRef oldRef = refCheck;
        refCheck = impl->meshRelation.triRef.data[edge / 3];
        if (!manifold_triref_same_face(&ref, &refCheck)) {
          return false;
        }
        if (ref.meshID != oldRef.meshID || ref.faceID != oldRef.faceID ||
            vec3_dot(impl->faceNormal.data[toRemove.pairedHalfedge / 3],
                     impl->faceNormal.data[tri]) < -0.5) {
          if (manifold_ccw(mat2x3_mul_vec3(projection, pLast),
                           mat2x3_mul_vec3(projection, pOld),
                           mat2x3_mul_vec3(projection, pNew),
                           impl->epsilon) != 0)
            return false;
        }
      }
      if (manifold_ccw(mat2x3_mul_vec3(projection, pNext),
                       mat2x3_mul_vec3(projection, pLast),
                       mat2x3_mul_vec3(projection, pNew),
                       impl->epsilon) < 0)
        return false;

      pLast = pNext;
      int paired = impl->halfedge.data[current].pairedHalfedge;
      if (paired < 0 || (size_t)paired >= impl->halfedge.len) return false;
      current = paired;
    }
    if (safeN >= (int)impl->halfedge.len) return false;
  }

  // Verify endVert fan: from pair(tri0edge.y) to tri1edge.z
  {
    int cur = impl->halfedge.data[tri0edge.y].pairedHalfedge;
    if (cur < 0 || (size_t)cur >= impl->halfedge.len) return false;
    int safeN = 0;
    while (cur != tri1edge.z && safeN++ < (int)impl->halfedge.len) {
      cur = manifold_next_halfedge(cur);
      if ((size_t)cur >= impl->halfedge.len) break;
      int p = impl->halfedge.data[cur].pairedHalfedge;
      if (p < 0 || (size_t)p >= impl->halfedge.len) break;
      cur = p;
    }
    if (cur != tri1edge.z) return false;
  }

  // Orbit endVert
  {
    int cur = impl->halfedge.data[tri0edge.y].pairedHalfedge;
    if (cur < 0 || (size_t)cur >= impl->halfedge.len) return false;
    int safeN = 0;
    while (cur != tri1edge.z && safeN++ < (int)impl->halfedge.len) {
      cur = manifold_next_halfedge(cur);
      vec_int_push(edges, cur);
      int paired = impl->halfedge.data[cur].pairedHalfedge;
      if (paired < 0 || (size_t)paired >= impl->halfedge.len) {
        edges->len = 0;
        return false;
      }
      cur = paired;
    }
    if (safeN >= (int)impl->halfedge.len) {
      edges->len = 0;
      return false;
    }
  }

  // Verify startVert fan: from start to tri0edge.z
  {
    int cur = start;
    int safeN = 0;
    int safeLimit = 2 * (int)impl->halfedge.len;
    while (cur != tri0edge.z && safeN++ < safeLimit) {
      cur = manifold_next_halfedge(cur);
      if ((size_t)cur >= impl->halfedge.len) break;
      if (impl->halfedge.data[cur].startVert < 0) break;
      int nxt = impl->halfedge.data[cur].pairedHalfedge;
      if (nxt < 0 || (size_t)nxt >= impl->halfedge.len) break;
      cur = nxt;
    }
    if (cur != tri0edge.z) {
      edges->len = 0;
      return false;
    }
  }

  // Remove toRemove.startVert
  ManifoldVec3 nan = {NAN, NAN, NAN};
  impl->vertPos.data[toRemove.startVert] = nan;
  impl_collapse_tri(impl, tri1edge);

  // Orbit startVert
  const int tri0 = edge / 3;
  const int tri1 = toRemove.pairedHalfedge / 3;
  current = start;
  int safeN = 0;
  int safeLimit = 2 * (int)impl->halfedge.len;
  while (current != tri0edge.z && safeN++ < safeLimit) {
    current = manifold_next_halfedge(current);

    if (manifold_impl_num_prop(impl) > 0) {
      int tri = current / 3;
      ManifoldTriRef ref = impl->meshRelation.triRef.data[tri];
      ManifoldTriRef ref0 = impl->meshRelation.triRef.data[tri0];
      ManifoldTriRef ref1 = impl->meshRelation.triRef.data[tri1];
      if (manifold_triref_same_face(&ref, &ref0)) {
        impl->halfedge.data[current].propVert =
            impl->halfedge.data[manifold_next_halfedge(edge)].propVert;
      } else if (manifold_triref_same_face(&ref, &ref1)) {
        impl->halfedge.data[current].propVert =
            impl->halfedge.data[toRemove.pairedHalfedge].propVert;
      }
    }

    const int vert = impl->halfedge.data[current].endVert;
    const int next = impl->halfedge.data[current].pairedHalfedge;
    for (size_t i = 0; i < edges->len; i++) {
      if (vert == impl->halfedge.data[edges->data[i]].endVert) {
        impl_form_loop(impl, edges->data[i], current);
        start = next;
        edges->len = i;
        break;
      }
    }
    current = next;
  }
  if (safeN >= safeLimit) return false;

  impl_update_vert(impl, endVert, start, tri0edge.z);
  impl_collapse_tri(impl, tri0edge);
  impl_remove_if_folded(impl, start);
  return true;
}

static void impl_recursive_edge_swap(ManifoldImpl *impl, int edge,
                                      int *tag, int *visited,
                                      ManifoldVecInt *edgeSwapStack,
                                      ManifoldVecInt *edges) {
  if (edge < 0) return;
  int pair = impl->halfedge.data[edge].pairedHalfedge;
  if (pair < 0) return;
  if (visited[edge] == *tag && visited[pair] == *tag) return;

  ManifoldIVec3 tri0edge = tri_of(edge);
  ManifoldIVec3 tri1edge = tri_of(pair);

  ManifoldMat2x3 projection =
      manifold_get_axis_aligned_projection(impl->faceNormal.data[edge / 3]);
  ManifoldVec2 v[4];
  int t0[3] = {tri0edge.x, tri0edge.y, tri0edge.z};
  for (int i = 0; i < 3; i++)
    v[i] = mat2x3_mul_vec3(projection,
               impl->vertPos.data[impl->halfedge.data[t0[i]].startVert]);
  if (manifold_ccw(v[0], v[1], v[2], impl->tolerance) > 0 ||
      !is01_longest(v[0], v[1], v[2]))
    return;

  // Switch to neighbor's projection
  projection =
      manifold_get_axis_aligned_projection(impl->faceNormal.data[pair / 3]);
  for (int i = 0; i < 3; i++)
    v[i] = mat2x3_mul_vec3(projection,
               impl->vertPos.data[impl->halfedge.data[t0[i]].startVert]);
  v[3] = mat2x3_mul_vec3(projection,
             impl->vertPos.data[impl->halfedge.data[tri1edge.z].startVert]);

  // Lambda equivalent for SwapEdge
  #define DO_SWAP_EDGE() do { \
    int v0_ = impl->halfedge.data[tri0edge.z].startVert; \
    int v1_ = impl->halfedge.data[tri1edge.z].startVert; \
    impl->halfedge.data[tri0edge.x].startVert = v1_; \
    impl->halfedge.data[tri0edge.z].endVert = v1_; \
    impl->halfedge.data[tri1edge.x].startVert = v0_; \
    impl->halfedge.data[tri1edge.z].endVert = v0_; \
    impl_pair_up(impl, tri0edge.x, \
                 impl->halfedge.data[tri1edge.z].pairedHalfedge); \
    impl_pair_up(impl, tri1edge.x, \
                 impl->halfedge.data[tri0edge.z].pairedHalfedge); \
    impl_pair_up(impl, tri0edge.z, tri1edge.z); \
    int st0 = tri0edge.x / 3; \
    int st1 = tri1edge.x / 3; \
    impl->faceNormal.data[st0] = impl->faceNormal.data[st1]; \
    impl->meshRelation.triRef.data[st0] = impl->meshRelation.triRef.data[st1]; \
    double l01_ = vec3_length(vec3_sub( \
        impl->vertPos.data[impl->halfedge.data[tri0edge.x].endVert], \
        impl->vertPos.data[impl->halfedge.data[tri0edge.x].startVert])); \
    double l02_ = vec3_length(vec3_sub( \
        impl->vertPos.data[impl->halfedge.data[tri0edge.y].endVert], \
        impl->vertPos.data[impl->halfedge.data[tri0edge.x].startVert])); \
    double a_ = l01_ > 0 ? fmin(1.0, fmax(0.0, l02_ / l01_)) : 0; \
    if (impl->properties.len > 0) { \
      impl->halfedge.data[tri0edge.y].propVert = \
          impl->halfedge.data[tri1edge.x].propVert; \
      impl->halfedge.data[tri0edge.x].propVert = \
          impl->halfedge.data[tri1edge.z].propVert; \
      impl->halfedge.data[tri0edge.z].propVert = \
          impl->halfedge.data[tri1edge.z].propVert; \
      int numP = impl->numProp; \
      int newProp = (int)(impl->properties.len / numP); \
      int propIdx0 = impl->halfedge.data[tri1edge.x].propVert; \
      int propIdx1 = impl->halfedge.data[tri1edge.y].propVert; \
      for (int pp = 0; pp < numP; pp++) { \
        double val = a_ * impl->properties.data[numP * propIdx0 + pp] + \
                     (1 - a_) * impl->properties.data[numP * propIdx1 + pp]; \
        vec_double_push(&impl->properties, val); \
      } \
      impl->halfedge.data[tri1edge.x].propVert = newProp; \
      impl->halfedge.data[tri0edge.z].propVert = newProp; \
    } \
    { \
      int cur_ = impl->halfedge.data[tri1edge.x].pairedHalfedge; \
      int endV_ = impl->halfedge.data[tri1edge.y].endVert; \
      int safeN_ = 0; \
      while (cur_ != tri0edge.y && safeN_++ < (int)impl->halfedge.len) { \
        cur_ = manifold_next_halfedge(cur_); \
        if (cur_ < 0 || (size_t)cur_ >= impl->halfedge.len) break; \
        if (impl->halfedge.data[cur_].endVert == endV_) { \
          impl_form_loop(impl, tri0edge.z, cur_); \
          impl_remove_if_folded(impl, tri0edge.z); \
          break; \
        } \
        int paired_ = impl->halfedge.data[cur_].pairedHalfedge; \
        if (paired_ < 0 || (size_t)paired_ >= impl->halfedge.len) break; \
        cur_ = paired_; \
      } \
    } \
  } while(0)

  if (manifold_ccw(v[1], v[0], v[3], impl->tolerance) <= 0) {
    if (!is01_longest(v[1], v[0], v[3])) return;
    DO_SWAP_EDGE();
    ManifoldVec2 e23 = vec2_sub(v[3], v[2]);
    if (vec2_dot(e23, e23) < impl->tolerance * impl->tolerance) {
      (*tag)++;
      impl_collapse_edge(impl, tri0edge.z, edges);
      edges->len = 0;
    } else {
      visited[edge] = *tag;
      visited[pair] = *tag;
      vec_int_push(edgeSwapStack, tri1edge.y);
      vec_int_push(edgeSwapStack, tri1edge.x);
      vec_int_push(edgeSwapStack, tri0edge.y);
      vec_int_push(edgeSwapStack, tri0edge.x);
    }
    return;
  } else if (manifold_ccw(v[0], v[3], v[2], impl->tolerance) <= 0 ||
             manifold_ccw(v[1], v[2], v[3], impl->tolerance) <= 0) {
    return;
  }
  // Normal path
  DO_SWAP_EDGE();
  visited[edge] = *tag;
  visited[pair] = *tag;
  vec_int_push(edgeSwapStack,
               impl->halfedge.data[tri1edge.x].pairedHalfedge);
  vec_int_push(edgeSwapStack,
               impl->halfedge.data[tri0edge.y].pairedHalfedge);
  #undef DO_SWAP_EDGE
}

// ---------- DedupeEdge ----------

static void impl_dedupe_edge(ManifoldImpl *impl, int edge) {
  int nbEdges = (int)impl->halfedge.len;
  if (edge < 0 || edge >= nbEdges) return;
  int startVert = impl->halfedge.data[edge].startVert;
  int endVert = impl->halfedge.data[edge].endVert;
  if (startVert < 0 || endVert < 0) return;
  int nextEdge = manifold_next_halfedge(edge);
  if (nextEdge < 0 || nextEdge >= nbEdges) return;
  int endProp = impl->halfedge.data[nextEdge].propVert;
  int pairedNext = impl->halfedge.data[nextEdge].pairedHalfedge;
  if (pairedNext < 0 || pairedNext >= nbEdges) return;
  int current = pairedNext;
  int safe = 0;
  while (current != edge && safe++ < nbEdges) {
    if (current < 0 || current >= nbEdges) return;
    int vert = impl->halfedge.data[current].startVert;
    if (vert == startVert) {
      // Single topological unit - needs 2 faces to split
      int newVert = (int)impl->vertPos.len;
      vec_vec3_push(&impl->vertPos, impl->vertPos.data[endVert]);
      if (impl->vertNormal.len > 0)
        vec_vec3_push(&impl->vertNormal, impl->vertNormal.data[endVert]);
      int nc = manifold_next_halfedge(current);
      if (nc < 0 || nc >= (int)impl->halfedge.len) return;
      int pcur = impl->halfedge.data[nc].pairedHalfedge;
      if (pcur < 0 || pcur >= (int)impl->halfedge.len) return;
      current = pcur;
      int opposite = impl->halfedge.data[nextEdge].pairedHalfedge;
      if (opposite < 0 || opposite >= (int)impl->halfedge.len) return;

      impl_update_vert(impl, newVert, current, opposite);

      int newHalfedge = (int)impl->halfedge.len;
      int oldFace = current / 3;
      int outsideVert = impl->halfedge.data[current].startVert;
      ManifoldHalfedge h1 = {endVert, newVert, -1, endProp};
      ManifoldHalfedge h2 = {newVert, outsideVert, -1, endProp};
      ManifoldHalfedge h3 = {outsideVert, endVert, -1,
                              impl->halfedge.data[current].propVert};
      vec_halfedge_push(&impl->halfedge, h1);
      vec_halfedge_push(&impl->halfedge, h2);
      vec_halfedge_push(&impl->halfedge, h3);
      impl_pair_up(impl, newHalfedge + 2,
                   impl->halfedge.data[current].pairedHalfedge);
      impl_pair_up(impl, newHalfedge + 1, current);
      if (impl->meshRelation.triRef.len > 0)
        vec_triref_push(&impl->meshRelation.triRef,
                        impl->meshRelation.triRef.data[oldFace]);
      if (impl->faceNormal.len > 0)
        vec_vec3_push(&impl->faceNormal, impl->faceNormal.data[oldFace]);

      newHalfedge += 3;
      oldFace = opposite / 3;
      outsideVert = impl->halfedge.data[opposite].startVert;
      ManifoldHalfedge h4 = {newVert, endVert, -1, endProp};
      ManifoldHalfedge h5 = {endVert, outsideVert, -1, endProp};
      ManifoldHalfedge h6 = {outsideVert, newVert, -1,
                              impl->halfedge.data[opposite].propVert};
      vec_halfedge_push(&impl->halfedge, h4);
      vec_halfedge_push(&impl->halfedge, h5);
      vec_halfedge_push(&impl->halfedge, h6);
      impl_pair_up(impl, newHalfedge + 2,
                   impl->halfedge.data[opposite].pairedHalfedge);
      impl_pair_up(impl, newHalfedge + 1, opposite);
      impl_pair_up(impl, newHalfedge, newHalfedge - 3);
      if (impl->meshRelation.triRef.len > 0)
        vec_triref_push(&impl->meshRelation.triRef,
                        impl->meshRelation.triRef.data[oldFace]);
      if (impl->faceNormal.len > 0)
        vec_vec3_push(&impl->faceNormal, impl->faceNormal.data[oldFace]);
      break;
    }
    int nc2 = manifold_next_halfedge(current);
    if (nc2 < 0 || nc2 >= (int)impl->halfedge.len) return;
    int p2 = impl->halfedge.data[nc2].pairedHalfedge;
    if (p2 < 0 || p2 >= (int)impl->halfedge.len) return;
    current = p2;
  }

  if (current == edge) {
    // Separate topological unit
    int newVert = (int)impl->vertPos.len;
    vec_vec3_push(&impl->vertPos, impl->vertPos.data[endVert]);
    if (impl->vertNormal.len > 0)
      vec_vec3_push(&impl->vertNormal, impl->vertNormal.data[endVert]);

    int cur = manifold_next_halfedge(current);
    if (cur < 0 || cur >= (int)impl->halfedge.len) { return; }
    int startE = cur;
    int safe2 = 0;
    do {
      int p = impl->halfedge.data[cur].pairedHalfedge;
      if (p < 0 || p >= (int)impl->halfedge.len) break;
      impl->halfedge.data[cur].startVert = newVert;
      impl->halfedge.data[p].endVert = newVert;
      cur = manifold_next_halfedge(p);
      if (cur < 0 || cur >= (int)impl->halfedge.len) break;
    } while (cur != startE && safe2++ < nbEdges);
  }

  // Orbit startVert to check for pinched vert
  int pair = impl->halfedge.data[edge].pairedHalfedge;
  if (pair < 0 || pair >= (int)impl->halfedge.len) return;
  int np = manifold_next_halfedge(pair);
  if (np < 0 || np >= (int)impl->halfedge.len) return;
  current = impl->halfedge.data[np].pairedHalfedge;
  if (current < 0 || current >= (int)impl->halfedge.len) return;
  int safe3 = 0;
  while (current != pair && safe3++ < nbEdges) {
    int vert = impl->halfedge.data[current].startVert;
    if (vert == endVert) break;
    int nc3 = manifold_next_halfedge(current);
    if (nc3 < 0 || nc3 >= (int)impl->halfedge.len) return;
    int p3 = impl->halfedge.data[nc3].pairedHalfedge;
    if (p3 < 0 || p3 >= (int)impl->halfedge.len) return;
    current = p3;
  }

  if (current == pair) {
    // Split pinched vert
    int newVert = (int)impl->vertPos.len;
    vec_vec3_push(&impl->vertPos, impl->vertPos.data[endVert]);
    if (impl->vertNormal.len > 0)
      vec_vec3_push(&impl->vertNormal, impl->vertNormal.data[endVert]);

    int cur = manifold_next_halfedge(current);
    if (cur < 0 || cur >= (int)impl->halfedge.len) { return; }
    int startE = cur;
    int safe4 = 0;
    do {
      int p = impl->halfedge.data[cur].pairedHalfedge;
      if (p < 0 || p >= (int)impl->halfedge.len) break;
      impl->halfedge.data[cur].startVert = newVert;
      impl->halfedge.data[p].endVert = newVert;
      cur = manifold_next_halfedge(p);
      if (cur < 0 || cur >= (int)impl->halfedge.len) break;
    } while (cur != startE && safe4++ < nbEdges);
  }
}

// ---------- SplitPinchedVerts ----------

void manifold_impl_split_pinched_verts(ManifoldImpl *impl) {
  size_t nbEdges = impl->halfedge.len;
  if (nbEdges == 0) return;

  size_t numVert = manifold_impl_num_vert(impl);
  bool *vertProcessed = (bool *)calloc(numVert, sizeof(bool));
  bool *halfedgeProcessed = (bool *)calloc(nbEdges, sizeof(bool));
  if (!vertProcessed || !halfedgeProcessed) {
    free(vertProcessed);
    free(halfedgeProcessed);
    return;
  }

  for (size_t i = 0; i < nbEdges; i++) {
    if (halfedgeProcessed[i]) continue;
    int vert = impl->halfedge.data[i].startVert;
    if (vert < 0 || (size_t)vert >= numVert) continue;
    if (vertProcessed[vert]) {
      vec_vec3_push(&impl->vertPos, impl->vertPos.data[vert]);
      int newV = (int)manifold_impl_num_vert(impl) - 1;
      // ForVert loop with safety guard
      int cur = (int)i;
      int safe = 0;
      do {
        halfedgeProcessed[cur] = true;
        int paired = impl->halfedge.data[cur].pairedHalfedge;
        if (paired < 0 || (size_t)paired >= nbEdges) break;
        impl->halfedge.data[cur].startVert = newV;
        impl->halfedge.data[paired].endVert = newV;
        cur = manifold_next_halfedge(paired);
        if (cur < 0 || (size_t)cur >= nbEdges) break;
      } while (cur != (int)i && safe++ < (int)nbEdges);
    } else {
      vertProcessed[vert] = true;
      int cur = (int)i;
      int safe = 0;
      do {
        halfedgeProcessed[cur] = true;
        int paired = impl->halfedge.data[cur].pairedHalfedge;
        if (paired < 0 || (size_t)paired >= nbEdges) break;
        cur = manifold_next_halfedge(paired);
        if (cur < 0 || (size_t)cur >= nbEdges) break;
      } while (cur != (int)i && safe++ < (int)nbEdges);
    }
  }

  free(vertProcessed);
  free(halfedgeProcessed);
}

// ---------- DedupeEdges ----------

static void impl_dedupe_edges(ManifoldImpl *impl) {
  size_t prevFlagged = SIZE_MAX;
  int maxRounds = 10;
  while (maxRounds-- > 0) {
    size_t nbEdges = impl->halfedge.len;
    ManifoldVecInt duplicates = {0};
    bool *local = (bool *)calloc(nbEdges, sizeof(bool));
    if (!local) return;

    for (size_t i = 0; i < nbEdges; i++) {
      if (local[i] || impl->halfedge.data[i].startVert == -1 ||
          impl->halfedge.data[i].endVert == -1)
        continue;

      // Simple approach: track endVerts of edges around this vert
      typedef struct { int endV; int minEdge; } EvPair;
      EvPair endVerts[64];
      int numEv = 0;

      // First pass: find minimal edge per endVert
      int cur = (int)i;
      int safeCount1 = 0;
      do {
        local[cur] = true;
        if (impl->halfedge.data[cur].startVert != -1 &&
            impl->halfedge.data[cur].endVert != -1) {
          int endV = impl->halfedge.data[cur].endVert;
          bool found = false;
          for (int k = 0; k < numEv; k++) {
            if (endVerts[k].endV == endV) {
              if (cur < endVerts[k].minEdge) endVerts[k].minEdge = cur;
              found = true;
              break;
            }
          }
          if (!found && numEv < 64) {
            endVerts[numEv].endV = endV;
            endVerts[numEv].minEdge = cur;
            numEv++;
          }
        }
        int paired = impl->halfedge.data[cur].pairedHalfedge;
        if (paired < 0 || (size_t)paired >= nbEdges) break;
        cur = manifold_next_halfedge(paired);
        if (cur < 0 || (size_t)cur >= nbEdges) break;
      } while (cur != (int)i && safeCount1++ < (int)nbEdges);

      // Second pass: flag duplicates
      cur = (int)i;
      int safeCount2 = 0;
      do {
        if (impl->halfedge.data[cur].startVert != -1 &&
            impl->halfedge.data[cur].endVert != -1) {
          int endV = impl->halfedge.data[cur].endVert;
          for (int k = 0; k < numEv; k++) {
            if (endVerts[k].endV == endV) {
              if (endVerts[k].minEdge != cur) {
                vec_int_push(&duplicates, cur);
              }
              break;
            }
          }
        }
        int paired = impl->halfedge.data[cur].pairedHalfedge;
        if (paired < 0 || (size_t)paired >= nbEdges) break;
        cur = manifold_next_halfedge(paired);
        if (cur < 0 || (size_t)cur >= nbEdges) break;
      } while (cur != (int)i && safeCount2++ < (int)nbEdges);
    }

    free(local);

    size_t numFlagged = 0;
    for (size_t di = 0; di < duplicates.len; di++) {
      impl_dedupe_edge(impl, duplicates.data[di]);
      numFlagged++;
    }
    vec_int_free(&duplicates);

    if (numFlagged == 0) break;
    // If we're not making progress, stop trying
    if (numFlagged >= prevFlagged) break;
    prevFlagged = numFlagged;
  }
}

// ---------- Edge collapse/swap helpers ----------

// ShortEdge predicate: flag edges shorter than epsilon where at least one
// endpoint is a new vert (>= firstNewVert)
static bool is_short_edge(const ManifoldImpl *impl, int edge,
                          int firstNewVert) {
  const ManifoldHalfedge *half = &impl->halfedge.data[edge];
  if (half->pairedHalfedge < 0) return false;
  if (half->startVert < 0 || half->endVert < 0) return false;
  if (half->startVert < firstNewVert && half->endVert < firstNewVert)
    return false;
  if ((size_t)half->startVert >= impl->vertPos.len ||
      (size_t)half->endVert >= impl->vertPos.len) return false;
  ManifoldVec3 delta = vec3_sub(impl->vertPos.data[half->endVert],
                                impl->vertPos.data[half->startVert]);
  return vec3_dot(delta, delta) < impl->epsilon * impl->epsilon;
}

// FlagEdge predicate: flag edges where startVert is surrounded by only two
// original triangles (colinear edge)
static bool is_flag_edge(const ManifoldImpl *impl, int edge,
                         int firstNewVert) {
  const ManifoldHalfedge *half = &impl->halfedge.data[edge];
  if (half->pairedHalfedge < 0) return false;
  if (half->startVert < 0 || half->endVert < 0) return false;
  if (half->startVert < firstNewVert) return false;

  ManifoldTriRef ref0 = impl->meshRelation.triRef.data[edge / 3];
  int current = manifold_next_halfedge(half->pairedHalfedge);
  ManifoldTriRef ref1 = impl->meshRelation.triRef.data[current / 3];
  bool ref1Updated = !manifold_triref_same_face(&ref0, &ref1);

  int maxIter = (int)impl->halfedge.len + 2;
  int iter = 0;
  while (current != edge && iter++ < maxIter) {
    int paired = impl->halfedge.data[current].pairedHalfedge;
    if (paired < 0 || (size_t)paired >= impl->halfedge.len) return false;
    current = manifold_next_halfedge(paired);
    if (current < 0 || (size_t)current >= impl->halfedge.len) return false;
    int tri = current / 3;
    ManifoldTriRef ref = impl->meshRelation.triRef.data[tri];
    if (!manifold_triref_same_face(&ref, &ref0) &&
        !manifold_triref_same_face(&ref, &ref1)) {
      if (!ref1Updated) {
        ref1 = ref;
        ref1Updated = true;
      } else {
        return false;
      }
    }
  }
  return true;
}

// SwappableEdge predicate: flag degenerate edges that should be swapped
static bool is_swappable_edge(const ManifoldImpl *impl, int edge,
                              int firstNewVert) {
  const ManifoldHalfedge *half = &impl->halfedge.data[edge];
  if (half->pairedHalfedge < 0) return false;
  if (half->startVert < 0 || half->endVert < 0) return false;
  int pair = half->pairedHalfedge;

  // Skip if all 4 verts are old
  int v2 = impl->halfedge.data[manifold_next_halfedge(edge)].endVert;
  int v3 = impl->halfedge.data[manifold_next_halfedge(pair)].endVert;
  if (half->startVert < firstNewVert && half->endVert < firstNewVert &&
      v2 < firstNewVert && v3 < firstNewVert)
    return false;

  if ((size_t)half->startVert >= impl->vertPos.len ||
      (size_t)half->endVert >= impl->vertPos.len ||
      (size_t)v2 >= impl->vertPos.len ||
      (size_t)v3 >= impl->vertPos.len)
    return false;

  int tri = edge / 3;
  ManifoldIVec3 triEdge = tri_of(edge);
  ManifoldMat2x3 projection =
      manifold_get_axis_aligned_projection(impl->faceNormal.data[tri]);
  ManifoldVec2 v[3];
  int te[3] = {triEdge.x, triEdge.y, triEdge.z};
  for (int i = 0; i < 3; i++)
    v[i] = mat2x3_mul_vec3(projection,
               impl->vertPos.data[impl->halfedge.data[te[i]].startVert]);
  if (manifold_ccw(v[0], v[1], v[2], impl->tolerance) > 0 ||
      !is01_longest(v[0], v[1], v[2]))
    return false;

  // Check neighbor
  tri = pair / 3;
  triEdge = tri_of(pair);
  projection =
      manifold_get_axis_aligned_projection(impl->faceNormal.data[tri]);
  te[0] = triEdge.x; te[1] = triEdge.y; te[2] = triEdge.z;
  for (int i = 0; i < 3; i++)
    v[i] = mat2x3_mul_vec3(projection,
               impl->vertPos.data[impl->halfedge.data[te[i]].startVert]);
  return manifold_ccw(v[0], v[1], v[2], impl->tolerance) > 0 ||
         is01_longest(v[0], v[1], v[2]);
}

static void impl_collapse_short_edges(ManifoldImpl *impl, int firstNewVert) {
  size_t nbEdges = impl->halfedge.len;
  ManifoldVecInt scratchBuffer = vec_int_create(0);

  for (size_t i = 0; i < nbEdges; i++) {
    if (is_short_edge(impl, (int)i, firstNewVert)) {
      scratchBuffer.len = 0;
      impl_collapse_edge(impl, (int)i, &scratchBuffer);
    }
  }

  vec_int_free(&scratchBuffer);
}

static void impl_collapse_colinear_edges(ManifoldImpl *impl,
                                         int firstNewVert) {
  ManifoldVecInt flagged = vec_int_create(0);
  ManifoldVecInt scratchBuffer = vec_int_create(0);

  int maxIterations = 20;
  while (maxIterations-- > 0) {
    size_t nbEdges = impl->halfedge.len;
    flagged.len = 0;
    for (size_t i = 0; i < nbEdges; i++) {
      if (is_flag_edge(impl, (int)i, firstNewVert))
        vec_int_push(&flagged, (int)i);
    }
    if (flagged.len == 0) break;
    size_t numCollapsed = 0;
    for (size_t i = 0; i < flagged.len; i++) {
      scratchBuffer.len = 0;
      int e = flagged.data[i];
      bool collapsed = impl_collapse_edge(impl, e, &scratchBuffer);
      if (collapsed) numCollapsed++;
    }
    if (numCollapsed == 0) break;
  }

  vec_int_free(&flagged);
  vec_int_free(&scratchBuffer);
}

static void impl_swap_degenerates(ManifoldImpl *impl, int firstNewVert) {
  size_t nbEdges = impl->halfedge.len;
  ManifoldVecInt flagged = vec_int_create(0);
  ManifoldVecInt edgeSwapStack = vec_int_create(0);
  ManifoldVecInt scratchBuffer = vec_int_create(0);

  for (size_t i = 0; i < nbEdges; i++) {
    if (is_swappable_edge(impl, (int)i, firstNewVert))
      vec_int_push(&flagged, (int)i);
  }

  int *visited = (int *)calloc(nbEdges > 0 ? nbEdges : 1, sizeof(int));
  int tag = 0;
  for (size_t i = 0; i < flagged.len; i++) {
    tag++;
    edgeSwapStack.len = 0;
    impl_recursive_edge_swap(impl, flagged.data[i], &tag, visited,
                             &edgeSwapStack, &scratchBuffer);
    int maxSwaps = (int)(nbEdges * 4);
    int swapCount = 0;
    while (edgeSwapStack.len > 0 && swapCount++ < maxSwaps) {
      int last = edgeSwapStack.data[edgeSwapStack.len - 1];
      edgeSwapStack.len--;
      impl_recursive_edge_swap(impl, last, &tag, visited,
                               &edgeSwapStack, &scratchBuffer);
    }
  }

  free(visited);
  vec_int_free(&flagged);
  vec_int_free(&edgeSwapStack);
  vec_int_free(&scratchBuffer);
}

// ---------- Public entry points ----------

void manifold_impl_cleanup_topology(ManifoldImpl *impl) {
  if (!impl->halfedge.len) return;
  manifold_impl_split_pinched_verts(impl);
  impl_dedupe_edges(impl);
}

void manifold_impl_simplify_topology(ManifoldImpl *impl, int firstNewVert) {
  if (!impl->halfedge.len) return;

  manifold_impl_cleanup_topology(impl);

  impl_collapse_short_edges(impl, firstNewVert);

  impl_collapse_colinear_edges(impl, firstNewVert);
  impl_swap_degenerates(impl, firstNewVert);
  manifold_impl_calculate_vert_normals(impl);
}

void manifold_impl_remove_degenerates(ManifoldImpl *impl, int firstNewVert) {
  if (!impl->halfedge.len) return;
  manifold_impl_cleanup_topology(impl);
  impl_collapse_short_edges(impl, firstNewVert);
  impl_swap_degenerates(impl, firstNewVert);
  manifold_impl_calculate_vert_normals(impl);
}
