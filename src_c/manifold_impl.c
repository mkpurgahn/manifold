// Copyright 2021 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Core implementation for the C11 Manifold port.

#include "manifold_impl.h"
#include <stdio.h>

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
    keep.data[impl->halfedge.data[i].propVert] = 1;
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
    impl->halfedge.data[i].propVert =
        propOld2New.data[impl->halfedge.data[i].propVert];
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
    for (int i = 0; i < 3; i++) {
      ManifoldVec3 pos =
          impl->vertPos.data[impl->halfedge.data[3 * face + i].startVert];
      center = vec3_add(center, pos);
      manifold_box_union_point(&faceBox->data[face], pos);
    }
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
      int oppositeFace = curr->pairedHalfedge / 3;
      int index = -1;
      for (int j = 0; j < 3; j++) {
        if (curr->startVert == impl->halfedge.data[oppositeFace * 3 + j].endVert)
          index = j;
      }
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
  size_t numTri = manifold_impl_num_tri(impl);
  ManifoldVec3 zero = manifold_vec3(0, 0, 0);
  impl->vertNormal = vec_vec3_create_fill(numVert, zero);

  for (size_t tri = 0; tri < numTri; tri++) {
    ManifoldVec3 normal = impl->faceNormal.data[tri];
    for (int i = 0; i < 3; i++) {
      int v = impl->halfedge.data[3 * tri + i].startVert;
      impl->vertNormal.data[v] = vec3_add(impl->vertNormal.data[v], normal);
    }
  }

  for (size_t i = 0; i < numVert; i++) {
    impl->vertNormal.data[i] = manifold_safe_normalize(impl->vertNormal.data[i]);
  }
}

void manifold_impl_set_normals_and_coplanar(ManifoldImpl *impl) {
  size_t numTri = manifold_impl_num_tri(impl);
  ManifoldVec3 zero = manifold_vec3(0, 0, 0);
  impl->faceNormal = vec_vec3_create_fill(numTri, zero);

  for (size_t tri = 0; tri < numTri; tri++) {
    ManifoldVec3 v0 = impl->vertPos.data[impl->halfedge.data[3 * tri].startVert];
    ManifoldVec3 v1 = impl->vertPos.data[impl->halfedge.data[3 * tri + 1].startVert];
    ManifoldVec3 v2 = impl->vertPos.data[impl->halfedge.data[3 * tri + 2].startVert];
    impl->faceNormal.data[tri] = manifold_safe_normalize(
        vec3_cross(vec3_sub(v1, v0), vec3_sub(v2, v0)));
  }

  manifold_impl_calculate_vert_normals(impl);

  // Set coplanar IDs in triRef
  if (impl->meshRelation.triRef.len == numTri) {
    for (size_t tri = 0; tri < numTri; tri++) {
      impl->meshRelation.triRef.data[tri].coplanarID = (int)tri;
    }
  }
}

void manifold_impl_create_halfedges(ManifoldImpl *impl,
                                     const ManifoldVecIVec3 *triProp,
                                     const ManifoldVecIVec3 *triVert) {
  size_t numTri = triProp->len;
  size_t numHalfedge = numTri * 3;

  impl->halfedge = vec_halfedge_create_n(numHalfedge);

  // For each triangle, create 3 halfedges
  const ManifoldVecIVec3 *vertSource = (triVert && triVert->len > 0) ? triVert : triProp;

  for (size_t tri = 0; tri < numTri; tri++) {
    ManifoldIVec3 verts = vertSource->data[tri];
    ManifoldIVec3 props = triProp->data[tri];
    for (int i = 0; i < 3; i++) {
      int j = (i + 1) % 3;
      ManifoldHalfedge *he = &impl->halfedge.data[3 * tri + i];
      he->startVert = ivec3_get(verts, i);
      he->endVert = ivec3_get(verts, j);
      he->pairedHalfedge = -1;
      he->propVert = ivec3_get(props, i);
    }
  }

  // Pair halfedges: find matching start/end pairs
  // Build a hash map of (startVert, endVert) -> halfedge index
  // Simple O(n^2) pairing for now (can optimize later)
  for (size_t i = 0; i < numHalfedge; i++) {
    if (impl->halfedge.data[i].pairedHalfedge >= 0) continue;
    int sv = impl->halfedge.data[i].startVert;
    int ev = impl->halfedge.data[i].endVert;
    for (size_t j = i + 1; j < numHalfedge; j++) {
      if (impl->halfedge.data[j].pairedHalfedge >= 0) continue;
      if (impl->halfedge.data[j].startVert == ev &&
          impl->halfedge.data[j].endVert == sv) {
        impl->halfedge.data[i].pairedHalfedge = (int)j;
        impl->halfedge.data[j].pairedHalfedge = (int)i;
        break;
      }
    }
  }
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
    if (he->startVert >= 0) he->startVert = old2new.data[he->startVert];
    if (he->endVert >= 0) he->endVert = old2new.data[he->endVert];
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
    ref->faceID = (int)tri;
    ref->coplanarID = (int)tri;
  }
}

// Stub implementations for edge operations
void manifold_impl_cleanup_topology(ManifoldImpl *impl) {
  manifold_impl_split_pinched_verts(impl);
}

void manifold_impl_split_pinched_verts(ManifoldImpl *impl) {
  // Detect and split pinched vertices
  // For now, simple stub - full implementation needed for complex meshes
  (void)impl;
}

void manifold_impl_remove_degenerates(ManifoldImpl *impl, int firstNewVert) {
  // Remove degenerate triangles (zero-area)
  size_t numTri = manifold_impl_num_tri(impl);
  (void)firstNewVert;
  for (size_t tri = 0; tri < numTri; tri++) {
    ManifoldHalfedge *he = &impl->halfedge.data[3 * tri];
    // Check for duplicate vertices in triangle
    if (he[0].startVert == he[1].startVert ||
        he[1].startVert == he[2].startVert ||
        he[2].startVert == he[0].startVert) {
      // Mark for removal
      he[0].pairedHalfedge = -1;
      he[1].pairedHalfedge = -1;
      he[2].pairedHalfedge = -1;
    }
  }
}

void manifold_impl_simplify_topology(ManifoldImpl *impl, int firstNewVert) {
  (void)impl;
  (void)firstNewVert;
  // TODO: full topology simplification
}
