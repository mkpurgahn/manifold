// Copyright 2020 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// General utilities for the C11 Manifold port.
// C equivalent of src/utils.h.
//
// Most utilities (CCW, Mat4, hash64bit, AtomicAdd, etc.) are already
// implemented in manifold_vec_math.h. This header provides additional
// utilities that complement those.

#ifndef MANIFOLD_UTILS_H_
#define MANIFOLD_UTILS_H_

#include "manifold_vec_math.h"
#include <stdatomic.h>

// Next and Prev indices in a triangle (mod 3).
static inline int manifold_next3(int i) {
  static const int next[3] = {1, 2, 0};
  return next[i];
}

static inline int manifold_prev3(int i) {
  static const int prev[3] = {2, 0, 1};
  return prev[i];
}

// Atomic add for doubles (used in curvature computation).
// Uses compare-and-swap loop since C11 stdatomic doesn't support double.
static inline double manifold_atomic_add_double(double *target, double add) {
  // Note: In the C port, curvature and other accumulations are done
  // sequentially, so this is just a regular add. When parallelism is
  // enabled via pthreads, proper atomic operations would be needed.
  double old = *target;
  *target += add;
  return old;
}

// Atomic add for integers.
static inline int manifold_atomic_add_int(int *target, int add) {
  _Atomic int *at = (_Atomic int *)target;
  return atomic_fetch_add(at, add);
}

// Hash function for uint64_t (same as C++ version).
// Already defined in manifold_vec_math.h as manifold_hash64.

// Identity function (template replacement).
static inline double manifold_identity_double(double v) { return v; }
static inline int manifold_identity_int(int v) { return v; }

// Negate function (template replacement).
static inline double manifold_negate_double(double v) { return -v; }
static inline int manifold_negate_int(int v) { return -v; }

#endif  // MANIFOLD_UTILS_H_
