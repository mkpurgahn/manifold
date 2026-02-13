// Copyright 2021 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Boolean operations for the C11 Manifold port.
// This is a simplified implementation that handles basic cases.
// Full boolean algorithm (edge-face intersection, winding numbers)
// requires more porting work from boolean3.cpp and boolean_result.cpp.

#ifndef MANIFOLD_BOOLEAN_H
#define MANIFOLD_BOOLEAN_H

#include "manifold_impl.h"

// Perform boolean operation between two manifolds
// Returns a new ManifoldImpl with the result
// op: 0=Add(Union), 1=Subtract(Difference), 2=Intersect
ManifoldError manifold_boolean_op(ManifoldImpl *result,
                                   const ManifoldImpl *p,
                                   const ManifoldImpl *q,
                                   ManifoldOpType op);

#endif // MANIFOLD_BOOLEAN_H
