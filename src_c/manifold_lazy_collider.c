// Copyright 2025 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Lazy collider for the C11 Manifold port.
// C equivalent of src/lazy_collider.cpp.
//
// In the C++ version, LazyCollider defers BVH construction until first query.
// In this C11 port, colliders are built eagerly in sort_geometry (manifold_impl.c).
// This file provides the equivalent helper functions.

#include "manifold_lazy_collider.h"
#include <stdlib.h>

// The core collider build/update logic is in manifold_impl.c:
// - manifold_impl_sort_geometry builds the collider from face boxes/morton codes
// - The collider is stored directly in ManifoldImpl.collider
// - impl->colliderBuilt flag tracks whether the collider is valid
//
// This file exists for structural parity with the C++ source tree.
// The lazy collider pattern (defer construction, share via shared_ptr) is
// replaced by eager construction (build once, owned by impl).
