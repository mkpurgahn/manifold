// Copyright 2024 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Mesh fix utilities for the C11 Manifold port.
// C equivalent of src/mesh_fixes.h.

#ifndef MANIFOLD_MESH_FIXES_H_
#define MANIFOLD_MESH_FIXES_H_

#include "manifold_types.h"
#include "manifold_vec_math.h"

// Flip a halfedge index to the opposite corner of its triangle.
// Given halfedge h in triangle t (h = 3*t + v), returns 3*t + (2-v).
static inline int manifold_flip_halfedge(int halfedge) {
  int tri = halfedge / 3;
  int vert = 2 - (halfedge - 3 * tri);
  return 3 * tri + vert;
}

// Transform normals by a matrix, re-normalizing and handling NaN.
static inline ManifoldVec3 manifold_transform_normal(ManifoldMat3 transform,
                                                      ManifoldVec3 normal) {
  normal = mat3_mul_vec3(transform, normal);
  double len = vec3_length(normal);
  if (len > 0.0) {
    normal = vec3_scale(normal, 1.0 / len);
  } else {
    normal = (ManifoldVec3){{0.0, 0.0, 0.0}};
  }
  // Check for NaN
  if (normal.x != normal.x) normal = (ManifoldVec3){{0.0, 0.0, 0.0}};
  return normal;
}

// Transform tangent by a matrix.
static inline ManifoldVec4 manifold_transform_tangent(ManifoldMat3 transform,
                                                       ManifoldVec4 tangent) {
  ManifoldVec3 t3 = {{tangent.x, tangent.y, tangent.z}};
  t3 = mat3_mul_vec3(transform, t3);
  return (ManifoldVec4){{t3.x, t3.y, t3.z, tangent.w}};
}

// Flip triangle winding order: swap halfedges [3*tri] and [3*tri+2],
// then swap startVert/endVert and fix pairedHalfedge references.
static inline void manifold_flip_tri(ManifoldHalfedge *halfedge, int tri) {
  ManifoldHalfedge tmp = halfedge[3 * tri];
  halfedge[3 * tri] = halfedge[3 * tri + 2];
  halfedge[3 * tri + 2] = tmp;

  for (int i = 0; i < 3; i++) {
    ManifoldHalfedge *h = &halfedge[3 * tri + i];
    int sv = h->startVert;
    h->startVert = h->endVert;
    h->endVert = sv;
    h->pairedHalfedge = manifold_flip_halfedge(h->pairedHalfedge);
  }
}

#endif  // MANIFOLD_MESH_FIXES_H_
