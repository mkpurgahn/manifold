// Copyright 2024 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Convex hull computation for the C11 Manifold port.
// Implements a simplified gift-wrapping / incremental approach.
// TODO: Port the full QuickHull algorithm from quickhull.cpp for better performance.

#ifndef MANIFOLD_HULL_H
#define MANIFOLD_HULL_H

#include "manifold_impl.h"

// Compute the convex hull of a set of 3D points.
// Result is stored in impl as a manifold mesh.
void manifold_convex_hull(ManifoldImpl *impl, const ManifoldVec3 *points, size_t numPoints);

#endif // MANIFOLD_HULL_H
