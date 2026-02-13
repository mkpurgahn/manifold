// Copyright 2021 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Boolean operations for the C11 Manifold port.
// Full implementation of the Manifold boolean algorithm with edge-face
// intersections and winding number computation.

#ifndef MANIFOLD_BOOLEAN_H
#define MANIFOLD_BOOLEAN_H

#include "manifold_impl.h"

// Intersection data: edges of A intersecting faces of B
typedef struct {
  ManifoldVecIntArr2 p1q2;  // pairs of (edge_idx, face_idx)
  ManifoldVecInt x12;       // winding number contributions
  ManifoldVecVec3 v12;      // intersection positions
} ManifoldIntersections;

// Boolean3 state
typedef struct {
  const ManifoldImpl *inP;
  const ManifoldImpl *inQ;
  bool expandP;
  ManifoldIntersections xv12;  // edges of P vs faces of Q
  ManifoldIntersections xv21;  // edges of Q vs faces of P
  ManifoldVecInt w03;          // winding numbers: verts of P in Q
  ManifoldVecInt w30;          // winding numbers: verts of Q in P
  bool valid;
} ManifoldBoolean3;

// Perform boolean operation between two manifolds
ManifoldError manifold_boolean_op(ManifoldImpl *result,
                                   const ManifoldImpl *p,
                                   const ManifoldImpl *q,
                                   ManifoldOpType op);

#endif // MANIFOLD_BOOLEAN_H
