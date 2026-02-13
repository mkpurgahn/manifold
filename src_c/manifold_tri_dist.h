// Copyright 2024 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Triangle-triangle distance for the C11 Manifold port.
// C equivalent of src/tri_dist.h.
//
// The triangle distance implementation (edge_edge_dist,
// distance_triangle_triangle_squared) is in manifold_properties.c.
// This header re-exports for structural parity.

#ifndef MANIFOLD_TRI_DIST_H_
#define MANIFOLD_TRI_DIST_H_

#include "manifold_vec_math.h"

// Edge-edge distance computation.
// Finds closest points on two line segments.
void manifold_edge_edge_dist(ManifoldVec3 *x, ManifoldVec3 *y,
                             ManifoldVec3 p, ManifoldVec3 a,
                             ManifoldVec3 q, ManifoldVec3 b);

// Triangle-triangle squared distance.
// Returns minimum squared distance between two triangles.
double manifold_distance_tri_tri_squared(const ManifoldVec3 p[3],
                                         const ManifoldVec3 q[3]);

#endif  // MANIFOLD_TRI_DIST_H_
