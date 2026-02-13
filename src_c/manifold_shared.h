// Copyright 2021 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Shared types and utilities for the C11 Manifold port.
// C equivalent of src/shared.h.
//
// Core data structures (Halfedge, Barycentric, TriRef, Box) and utility
// functions (SafeNormalize, MaxEpsilon, NextHalfedge, GetBarycentric, CCW)
// are already defined in manifold_types.h and manifold_vec_math.h.
// This header re-exports them for consistency with the C++ structure.

#ifndef MANIFOLD_SHARED_H_
#define MANIFOLD_SHARED_H_

#include "manifold_types.h"
#include "manifold_vec_math.h"

// All shared types and functions are already provided by:
// - manifold_types.h: ManifoldHalfedge, ManifoldTriRef, ManifoldBox, etc.
// - manifold_vec_math.h: vec3_normalize, manifold_ccw, manifold_next_halfedge,
//   manifold_get_barycentric, manifold_safe_normalize, etc.

#endif  // MANIFOLD_SHARED_H_
