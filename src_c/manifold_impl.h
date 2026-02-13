// Copyright 2021 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Core mesh implementation structure for the C11 Manifold port.

#ifndef MANIFOLD_IMPL_H
#define MANIFOLD_IMPL_H

#include "manifold_collider.h"
#include "manifold_disjoint_sets.h"
#include "manifold_hashtable.h"
#include "manifold_vec.h"
#include "manifold_vec_math.h"

// Map from meshID → Relation (simple growable array for now)
typedef struct {
  int originalID;
  ManifoldMat3x4 transform;
  bool backSide;
} ManifoldRelation;

typedef struct {
  int key; // meshID
  ManifoldRelation value;
} ManifoldMeshIDEntry;

MANIFOLD_VEC_STRUCT(ManifoldMeshIDEntry, ManifoldVecMeshIDEntry)
MANIFOLD_VEC_FUNCS(ManifoldMeshIDEntry, ManifoldVecMeshIDEntry, vec_meshid)

typedef struct {
  int originalID;
  ManifoldVecMeshIDEntry meshIDtransform;
  ManifoldVecTriRef triRef;
} ManifoldMeshRelation;

typedef struct {
  int tri, start4, end4;
} ManifoldBaryIndices;

// The core internal mesh representation
typedef struct {
  ManifoldBox bBox;
  double epsilon;
  double tolerance;
  int numProp;
  ManifoldError status;

  ManifoldVecVec3 vertPos;
  ManifoldVecHalfedge halfedge;
  ManifoldVecDouble properties;
  ManifoldVecVec3 vertNormal;
  ManifoldVecVec3 faceNormal;
  ManifoldVecVec4 halfedgeTangent;
  ManifoldMeshRelation meshRelation;

  // Collider (built lazily)
  ManifoldCollider collider;
  bool colliderBuilt;
} ManifoldImpl;

// Global mesh ID counter
extern uint32_t manifold_mesh_id_counter;
uint32_t manifold_reserve_ids(uint32_t n);

// Init/destroy
void manifold_impl_init(ManifoldImpl *impl);
void manifold_impl_free(ManifoldImpl *impl);
void manifold_impl_make_empty(ManifoldImpl *impl, ManifoldError status);

// Basic queries
static inline bool manifold_impl_is_empty(const ManifoldImpl *impl) {
  return impl->halfedge.len / 3 == 0;
}
static inline size_t manifold_impl_num_vert(const ManifoldImpl *impl) {
  return impl->vertPos.len;
}
static inline size_t manifold_impl_num_edge(const ManifoldImpl *impl) {
  return impl->halfedge.len / 2;
}
static inline size_t manifold_impl_num_tri(const ManifoldImpl *impl) {
  return impl->halfedge.len / 3;
}
static inline size_t manifold_impl_num_prop(const ManifoldImpl *impl) {
  return (size_t)impl->numProp;
}
static inline size_t manifold_impl_num_prop_vert(const ManifoldImpl *impl) {
  return impl->numProp == 0 ? impl->vertPos.len
                            : impl->properties.len / impl->numProp;
}

// Mesh relation helpers
static inline ManifoldRelation *manifold_meshrelation_find(
    ManifoldMeshRelation *rel, int meshID) {
  for (size_t i = 0; i < rel->meshIDtransform.len; i++) {
    if (rel->meshIDtransform.data[i].key == meshID)
      return &rel->meshIDtransform.data[i].value;
  }
  return NULL;
}

static inline void manifold_meshrelation_insert(ManifoldMeshRelation *rel,
                                                 int meshID,
                                                 ManifoldRelation r) {
  // Check if exists
  for (size_t i = 0; i < rel->meshIDtransform.len; i++) {
    if (rel->meshIDtransform.data[i].key == meshID) {
      rel->meshIDtransform.data[i].value = r;
      return;
    }
  }
  ManifoldMeshIDEntry entry = {meshID, r};
  vec_meshid_push(&rel->meshIDtransform, entry);
}

// Forward declarations for impl methods
void manifold_impl_calculate_bbox(ManifoldImpl *impl);
void manifold_impl_set_epsilon(ManifoldImpl *impl, double minEpsilon,
                                bool useSingle);
bool manifold_impl_is_manifold(const ManifoldImpl *impl);
bool manifold_impl_is_2manifold(const ManifoldImpl *impl);
bool manifold_impl_is_finite(const ManifoldImpl *impl);
void manifold_impl_sort_geometry(ManifoldImpl *impl);
void manifold_impl_sort_verts(ManifoldImpl *impl);
void manifold_impl_reindex_verts(ManifoldImpl *impl,
                                  const ManifoldVecInt *vertNew2Old,
                                  size_t numOldVert);
void manifold_impl_compact_props(ManifoldImpl *impl);
void manifold_impl_get_face_box_morton(const ManifoldImpl *impl,
                                        ManifoldVecBox *faceBox,
                                        ManifoldVecU32 *faceMorton);
void manifold_impl_sort_faces(ManifoldImpl *impl, ManifoldVecBox *faceBox,
                               ManifoldVecU32 *faceMorton);
void manifold_impl_gather_faces(ManifoldImpl *impl,
                                 const ManifoldVecInt *faceNew2Old);
void manifold_impl_reorder_halfedges(ManifoldImpl *impl);

void manifold_impl_create_halfedges(ManifoldImpl *impl,
                                     const ManifoldVecIVec3 *triProp,
                                     const ManifoldVecIVec3 *triVert);
void manifold_impl_set_normals_and_coplanar(ManifoldImpl *impl);
void manifold_impl_dedupe_prop_verts(ManifoldImpl *impl);
void manifold_impl_remove_unreferenced_verts(ManifoldImpl *impl);
void manifold_impl_initialize_original(ManifoldImpl *impl);
void manifold_impl_calculate_vert_normals(ManifoldImpl *impl);

// edge_op.cpp equivalents
void manifold_impl_cleanup_topology(ManifoldImpl *impl);
void manifold_impl_simplify_topology(ManifoldImpl *impl, int firstNewVert);
void manifold_impl_remove_degenerates(ManifoldImpl *impl, int firstNewVert);
void manifold_impl_split_pinched_verts(ManifoldImpl *impl);

// properties.cpp equivalents
double manifold_impl_get_volume(const ManifoldImpl *impl);
double manifold_impl_get_surface_area(const ManifoldImpl *impl);
void manifold_impl_calculate_curvature(ManifoldImpl *impl, int gaussianIdx,
                                        int meanIdx);
void manifold_impl_set_properties(ManifoldImpl *impl, int numProp,
    void (*propFunc)(double *newProp, ManifoldVec3 pos, const double *oldProp, void *ctx),
    void *ctx);

// Decompose helper
int manifold_impl_decompose(const ManifoldImpl *impl, ManifoldImpl *components,
                             int maxComponents);

// constructors.cpp equivalents
void manifold_impl_tetrahedron(ManifoldImpl *impl);
void manifold_impl_cube(ManifoldImpl *impl, ManifoldMat3x4 transform);
void manifold_impl_octahedron(ManifoldImpl *impl, ManifoldMat3x4 transform);
void manifold_impl_extrude(ManifoldImpl *impl,
                           const ManifoldVec2 *polyVerts,
                           const int *polySizes, int nPolys,
                           double height, int nDivisions,
                           double twistDegrees, ManifoldVec2 scaleTop);
void manifold_impl_revolve(ManifoldImpl *impl,
                           const ManifoldVec2 *polyVerts,
                           const int *polySizes, int nPolys,
                           int circularSegments, double revolveDegrees);

// sdf.c equivalents
void manifold_impl_level_set(ManifoldImpl *impl,
                              double (*sdf)(double x, double y, double z, void *ctx),
                              void *ctx,
                              ManifoldBox bounds,
                              double edgeLength,
                              double level);

#endif // MANIFOLD_IMPL_H
