// Copyright 2022 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Iterator utilities for the C11 Manifold port.
// C equivalent of src/iters.h.
//
// The C++ version provides template iterator adapters:
//   - TransformIterator: applies a function during iteration
//   - CountingIterator: generates sequential indices
//   - StridedRange: iterates with a stride
//
// In C11, these patterns are replaced by:
//   - Direct for-loops with inline transformation
//   - Index variables (size_t i)
//   - Explicit stride arithmetic (i * stride + offset)
//
// This header provides macros for common iteration patterns.

#ifndef MANIFOLD_ITERS_H_
#define MANIFOLD_ITERS_H_

#include <stddef.h>

// Iterate over a range [start, end) with a body.
#define MANIFOLD_FOR(var, start, end) \
  for (size_t var = (start); var < (end); var++)

// Iterate with stride: var takes values start, start+stride, ...
#define MANIFOLD_FOR_STRIDE(var, start, end, stride) \
  for (size_t var = (start); var < (end); var += (stride))

// Count elements satisfying a predicate in range [start, end).
#define MANIFOLD_COUNT_IF(result, start, end, predicate) \
  do {                                                    \
    size_t _count = 0;                                    \
    for (size_t _i = (start); _i < (end); _i++) {        \
      if (predicate(_i)) _count++;                        \
    }                                                     \
    (result) = _count;                                    \
  } while (0)

#endif  // MANIFOLD_ITERS_H_
