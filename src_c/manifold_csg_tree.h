// Copyright 2022 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// CSG tree operations for the C11 Manifold port.
// C equivalent of src/csg_tree.cpp.
//
// In the C++ version, CSG operations use a lazy evaluation tree with
// virtual dispatch (CsgNode/CsgLeafNode/CsgOpNode). In this C11 port,
// operations are evaluated eagerly since C lacks virtual methods and
// shared_ptr. The Compose and BatchBoolean functions are implemented
// directly, matching the C++ behavior.
//
// This file provides:
//   - manifold_compose: Efficient union of pairwise disjoint meshes
//   - manifold_batch_boolean_impl: Heap-based batch boolean (union/intersect)
//   - manifold_batch_union: Optimized batch union with disjoint composition

#ifndef MANIFOLD_CSG_TREE_H_
#define MANIFOLD_CSG_TREE_H_

#include "manifold_impl.h"

// Forward declarations (implemented in manifold_api.c / manifold_boolean.c)
// These are the core operations that the CSG tree delegates to.

// Compose: merge pairwise disjoint meshes without boolean operations.
// This is used by BatchUnion to efficiently combine non-overlapping meshes.
// Equivalent to CsgLeafNode::Compose in C++.
void manifold_compose_impls(const ManifoldImpl **impls, int count,
                            ManifoldImpl *out);

// SimpleBoolean: direct boolean operation between two meshes.
// Equivalent to SimpleBoolean in C++ csg_tree.cpp.
// (Already implemented as manifold_boolean3 in manifold_boolean.c)

#endif  // MANIFOLD_CSG_TREE_H_
