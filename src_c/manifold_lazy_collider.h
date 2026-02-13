// Copyright 2025 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Lazy collider for the C11 Manifold port.
// C equivalent of src/lazy_collider.cpp.
//
// In the C++ version, LazyCollider uses std::shared_ptr, std::mutex, and
// std::optional to defer BVH construction until needed. In this C11 port,
// colliders are built eagerly when needed (in manifold_impl.c sort_geometry).
// This header provides the equivalent data structures and operations.

#ifndef MANIFOLD_LAZY_COLLIDER_H_
#define MANIFOLD_LAZY_COLLIDER_H_

#include "manifold_collider.h"
#include <stdbool.h>

// In the C port, the collider is built eagerly and stored directly in
// ManifoldImpl. The "lazy" aspect is handled by a simple flag:
// impl->colliderBuilt indicates whether the collider has been constructed.
//
// The C++ LazyCollider supports three modes:
// 1. LeafData: raw boxes + morton codes, not yet built
// 2. Base + transform: reuse another collider with a transform
// 3. Base + updatedLeafBox: reuse another collider with updated boxes
//
// In the C port, mode 1 is the standard path (build from face boxes/morton).
// Modes 2 and 3 are not needed because we don't share colliders between
// manifold instances (each ManifoldImpl owns its collider).

// Check if a transform matrix is axis-aligned (only permutation + scale).
// Matching C++ LazyCollider::IsAxisAligned.
static inline bool manifold_collider_is_axis_aligned(const ManifoldMat3x4 *transform) {
  for (int row = 0; row < 3; row++) {
    int zeroCount = 0;
    for (int col = 0; col < 3; col++) {
      ManifoldVec3 c = transform->cols[col];
      double val = (row == 0) ? c.x : (row == 1) ? c.y : c.z;
      if (val == 0.0) zeroCount++;
    }
    if (zeroCount != 2) return false;
  }
  return true;
}

// Build or rebuild the collider for a ManifoldImpl.
// This is the eager equivalent of LazyCollider::EnsureBuilt().
// (Implementation is in manifold_impl.c as manifold_impl_build_collider)

#endif  // MANIFOLD_LAZY_COLLIDER_H_
