// Copyright 2025 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// 2D KD-tree for polygon triangulation.
// C equivalent of src/tree2d.cpp.
//
// The implementation is in manifold_tree2d.h (header-only, matching the C++
// pattern where tree2d.h contains the template query function).
// This file exists for structural parity with the C++ source tree.

#include "manifold_tree2d.h"

// All tree2d functions are implemented in manifold_tree2d.h:
// - manifold_build_2d_tree_impl: recursive median-split sort
// - manifold_build_2d_tree: entry point (skips for <=8 points)
// - manifold_query_2d_tree: range query with explicit stack
