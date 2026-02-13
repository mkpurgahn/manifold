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
Manifold manifold_empty(void);
Manifold manifold_tetrahedron(void);
Manifold manifold_cube(ManifoldVec3 size, bool center);
Manifold manifold_sphere(double radius, int circularSegments);
Manifold manifold_cylinder(double height, double radiusLow, double radiusHigh,
                           int circularSegments, bool center);

// Create a manifold from raw vertex positions and triangle indices.
Manifold manifold_from_mesh(const ManifoldVec3 *vertPos, size_t numVert,
                            const ManifoldIVec3 *triVerts, size_t numTri);

// Extrude a 2D polygon cross-section along the Z-axis.
// crossSection: array of simple polygons (array of vec2 arrays)
// nPolys: number of polygons
// polySizes: array of polygon vertex counts (length nPolys)
// polyVerts: flat array of vec2 vertices for all polygons
// height: extrusion height
// nDivisions: number of intermediate layers (0 = just top and bottom)
// twistDegrees: total twist from bottom to top
// scaleTop: scale factors at the top (1,1 = no scaling)
Manifold manifold_extrude(const ManifoldVec2 *polyVerts,
                          const int *polySizes, int nPolys,
                          double height, int nDivisions,
                          double twistDegrees, ManifoldVec2 scaleTop);

// Revolve a 2D polygon cross-section around the Y-axis.
Manifold manifold_revolve(const ManifoldVec2 *polyVerts,
                          const int *polySizes, int nPolys,
                          int circularSegments, double revolveDegrees);

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

// ---------- Deformations ----------
Manifold manifold_warp(const Manifold *m,
                       void (*warpFn)(double *x, double *y, double *z, void *ctx),
                       void *ctx);

// ---------- Quality ----------
int manifold_get_circular_segments(double radius);
void manifold_set_circular_segments(int n);
void manifold_set_min_circular_angle(double angle);
void manifold_set_min_circular_edge_length(double length);
void manifold_quality_reset(void);

// ---------- Additional Information ----------
int manifold_genus(const Manifold *m);
int manifold_original_id(const Manifold *m);
double manifold_get_epsilon(const Manifold *m);
double manifold_get_tolerance(const Manifold *m);
size_t manifold_num_prop(const Manifold *m);
size_t manifold_num_prop_vert(const Manifold *m);

// ---------- Mirror ----------
Manifold manifold_mirror(const Manifold *m, ManifoldVec3 normal);

// ---------- Split / Trim ----------
// Split by another manifold - returns intersection in *first and difference in *second
void manifold_split(const Manifold *m, const Manifold *cutter,
                    Manifold *first, Manifold *second);
// Split by plane - returns the portion in the direction of normal in *first
void manifold_split_by_plane(const Manifold *m, ManifoldVec3 normal,
                              double originOffset, Manifold *first, Manifold *second);
// Trim to the portion in the direction of normal
Manifold manifold_trim_by_plane(const Manifold *m, ManifoldVec3 normal,
                                 double originOffset);

// ---------- Decompose ----------
// Decompose into topologically disconnected components.
// Returns number of components. Caller provides array and frees each.
int manifold_decompose(const Manifold *m, Manifold **components, int maxComponents);

// ---------- Batch Boolean ----------
Manifold manifold_batch_boolean(const Manifold *manifolds, int count,
                                 ManifoldOpType op);

// ---------- SetProperties ----------
Manifold manifold_set_properties(const Manifold *m, int numProp,
    void (*propFunc)(double *newProp, ManifoldVec3 pos, const double *oldProp, void *ctx),
    void *ctx);

// ---------- CalculateCurvature ----------
Manifold manifold_calculate_curvature(const Manifold *m, int gaussianIdx,
                                       int meanIdx);

// ---------- AsOriginal ----------
Manifold manifold_as_original(const Manifold *m);
uint32_t manifold_reserve_ids_api(uint32_t n);

// ---------- Simplify ----------
Manifold manifold_simplify(const Manifold *m, double tolerance);
Manifold manifold_set_tolerance(const Manifold *m, double tolerance);

// ---------- Refine ----------
// Subdivide each triangle into n^2 sub-triangles (n>=2)
Manifold manifold_refine(const Manifold *m, int n);

// ---------- Mesh Data Access ----------
// Get raw vertex positions (read-only)
const ManifoldVec3 *manifold_get_vert_positions(const Manifold *m, size_t *count);
// Get triangle indices as triples of halfedge start vertices
void manifold_get_triangles(const Manifold *m, ManifoldIVec3 *out, size_t *count);
// Get mesh as MeshGL-like flat arrays
void manifold_get_mesh(const Manifold *m,
                       float **vertProps, size_t *numVert, size_t *numProp,
                       int **triVerts, size_t *numTri);
// Free mesh data returned by manifold_get_mesh
void manifold_free_mesh(float *vertProps, int *triVerts);

#endif // MANIFOLD_API_H
