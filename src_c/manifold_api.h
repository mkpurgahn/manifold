// Copyright 2021 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Public API for the C11 Manifold port.

#ifndef MANIFOLD_API_H
#define MANIFOLD_API_H

#include "manifold_impl.h"

// The public Manifold object wraps an Impl
typedef struct {
  ManifoldImpl impl;
} Manifold;

// ---------- Lifecycle ----------
void manifold_create(Manifold *m);
void manifold_destroy(Manifold *m);
void manifold_copy(Manifold *dst, const Manifold *src);

// ---------- Constructors ----------
Manifold manifold_tetrahedron(void);
Manifold manifold_cube(ManifoldVec3 size, bool center);
Manifold manifold_sphere(double radius, int circularSegments);
Manifold manifold_cylinder(double height, double radiusLow, double radiusHigh,
                           int circularSegments, bool center);

// ---------- Information ----------
ManifoldError manifold_status(const Manifold *m);
bool manifold_is_empty(const Manifold *m);
size_t manifold_num_vert(const Manifold *m);
size_t manifold_num_edge(const Manifold *m);
size_t manifold_num_tri(const Manifold *m);
ManifoldBox manifold_bounding_box(const Manifold *m);

// ---------- Measurement ----------
double manifold_volume(const Manifold *m);
double manifold_surface_area(const Manifold *m);

// ---------- Transformations ----------
Manifold manifold_translate(const Manifold *m, ManifoldVec3 v);
Manifold manifold_scale(const Manifold *m, ManifoldVec3 v);
Manifold manifold_rotate(const Manifold *m, double xDeg, double yDeg,
                         double zDeg);
Manifold manifold_transform(const Manifold *m, ManifoldMat3x4 t);

// ---------- Boolean ----------
Manifold manifold_boolean(const Manifold *a, const Manifold *b,
                          ManifoldOpType op);
Manifold manifold_union(const Manifold *a, const Manifold *b);
Manifold manifold_difference(const Manifold *a, const Manifold *b);
Manifold manifold_intersection(const Manifold *a, const Manifold *b);

// ---------- Hull ----------
Manifold manifold_hull(const Manifold *m);
Manifold manifold_hull_points(const ManifoldVec3 *points, size_t numPoints);

// ---------- SDF Level Set ----------
Manifold manifold_level_set(double (*sdf)(double x, double y, double z, void *ctx),
                            void *ctx, ManifoldBox bounds, double edgeLength,
                            double level, double tolerance);

// ---------- Quality ----------
int manifold_get_circular_segments(double radius);
void manifold_set_circular_segments(int n);
void manifold_set_min_circular_angle(double angle);
void manifold_set_min_circular_edge_length(double length);
void manifold_quality_reset(void);

#endif // MANIFOLD_API_H
