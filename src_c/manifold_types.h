// Copyright 2021 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// C11 port of the Manifold geometry library.

#ifndef MANIFOLD_TYPES_H
#define MANIFOLD_TYPES_H

#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>
#include <float.h>

// ---------- Scalar constants ----------

#define MANIFOLD_PI      3.14159265358979323846264338327950288
#define MANIFOLD_TWO_PI  6.28318530717958647692528676655900576
#define MANIFOLD_HALF_PI 1.57079632679489661923132169163975144
#define MANIFOLD_PRECISION 1e-12

// ---------- Vector types ----------

typedef struct { double x, y; } ManifoldVec2;
typedef struct { double x, y, z; } ManifoldVec3;
typedef struct { double x, y, z, w; } ManifoldVec4;
typedef struct { int x, y; } ManifoldIVec2;
typedef struct { int x, y, z; } ManifoldIVec3;
typedef struct { int x, y, z, w; } ManifoldIVec4;
typedef struct { bool x, y, z, w; } ManifoldBVec4;

// ---------- Matrix types (column-major) ----------
// mat2 = 2x2, mat3 = 3x3, etc.
// mat3x4 = 3 rows, 4 columns (used as affine transform)

typedef struct { ManifoldVec2 cols[2]; } ManifoldMat2;
// matMxN = M rows, N columns, stored column-major: N columns of vecM
typedef struct { ManifoldVec3 cols[2]; } ManifoldMat3x2;
typedef struct { ManifoldVec4 cols[2]; } ManifoldMat4x2;
typedef struct { ManifoldVec2 cols[3]; } ManifoldMat2x3;
typedef struct { ManifoldVec3 cols[3]; } ManifoldMat3;
typedef struct { ManifoldVec3 cols[4]; } ManifoldMat3x4;
typedef struct { ManifoldVec4 cols[3]; } ManifoldMat4x3;
typedef struct { ManifoldVec4 cols[4]; } ManifoldMat4;

// ---------- Bounding types ----------

typedef struct {
  ManifoldVec3 min;
  ManifoldVec3 max;
} ManifoldBox;

typedef struct {
  ManifoldVec2 min;
  ManifoldVec2 max;
} ManifoldRect;

// ---------- Enums ----------

typedef enum {
  MANIFOLD_OP_ADD = 0,
  MANIFOLD_OP_SUBTRACT = 1,
  MANIFOLD_OP_INTERSECT = 2
} ManifoldOpType;

typedef enum {
  MANIFOLD_ERROR_NO_ERROR = 0,
  MANIFOLD_ERROR_NON_FINITE_VERTEX,
  MANIFOLD_ERROR_NOT_MANIFOLD,
  MANIFOLD_ERROR_VERTEX_OUT_OF_BOUNDS,
  MANIFOLD_ERROR_PROPERTIES_WRONG_LENGTH,
  MANIFOLD_ERROR_MISSING_POSITION_PROPERTIES,
  MANIFOLD_ERROR_MERGE_VECTORS_DIFFERENT_LENGTHS,
  MANIFOLD_ERROR_MERGE_INDEX_OUT_OF_BOUNDS,
  MANIFOLD_ERROR_TRANSFORM_WRONG_LENGTH,
  MANIFOLD_ERROR_RUN_INDEX_WRONG_LENGTH,
  MANIFOLD_ERROR_FACE_ID_WRONG_LENGTH,
  MANIFOLD_ERROR_INVALID_CONSTRUCTION,
  MANIFOLD_ERROR_RESULT_TOO_LARGE
} ManifoldError;

// ---------- Core data structures ----------

typedef struct {
  size_t halfedge;
  double smoothness;
} ManifoldSmoothness;

typedef struct {
  int startVert;
  int endVert;
  int pairedHalfedge;
  int propVert;
} ManifoldHalfedge;

typedef struct {
  int tri;
  ManifoldVec4 uvw;
} ManifoldBarycentric;

typedef struct {
  int meshID;
  int originalID;
  int faceID;
  int coplanarID;
} ManifoldTriRef;

typedef struct {
  int first;
  int second;
  int halfedgeIdx;
} ManifoldTmpEdge;

// ---------- Execution parameters ----------

typedef struct {
  bool intermediateChecks;
  bool selfIntersectionChecks;
  bool processOverlaps;
  bool suppressErrors;
  bool cleanupTriangles;
  int verbose;
} ManifoldExecParams;

// ---------- Quality settings ----------

typedef struct {
  int circularSegments;
  double minCircularAngle;
  double minCircularEdgeLength;
} ManifoldQuality;

// ---------- Default initializers ----------

static inline ManifoldVec2 manifold_vec2(double x, double y) {
  ManifoldVec2 v = {x, y};
  return v;
}

static inline ManifoldVec3 manifold_vec3(double x, double y, double z) {
  ManifoldVec3 v = {x, y, z};
  return v;
}

static inline ManifoldVec4 manifold_vec4(double x, double y, double z, double w) {
  ManifoldVec4 v = {x, y, z, w};
  return v;
}

static inline ManifoldIVec2 manifold_ivec2(int x, int y) {
  ManifoldIVec2 v = {x, y};
  return v;
}

static inline ManifoldIVec3 manifold_ivec3(int x, int y, int z) {
  ManifoldIVec3 v = {x, y, z};
  return v;
}

static inline ManifoldIVec4 manifold_ivec4(int x, int y, int z, int w) {
  ManifoldIVec4 v = {x, y, z, w};
  return v;
}

static inline ManifoldVec3 manifold_vec3_splat(double v) {
  return manifold_vec3(v, v, v);
}

static inline ManifoldVec2 manifold_vec2_splat(double v) {
  return manifold_vec2(v, v);
}

// Default box: empty (inverted)
static inline ManifoldBox manifold_box_empty(void) {
  ManifoldBox b;
  b.min = manifold_vec3_splat(INFINITY);
  b.max = manifold_vec3_splat(-INFINITY);
  return b;
}

static inline ManifoldRect manifold_rect_empty(void) {
  ManifoldRect r;
  r.min = manifold_vec2_splat(INFINITY);
  r.max = manifold_vec2_splat(-INFINITY);
  return r;
}

static inline ManifoldExecParams manifold_exec_params_default(void) {
  ManifoldExecParams p;
  p.intermediateChecks = false;
  p.selfIntersectionChecks = false;
  p.processOverlaps = true;
  p.suppressErrors = false;
  p.cleanupTriangles = true;
  p.verbose = 0;
  return p;
}

static inline bool manifold_halfedge_is_forward(const ManifoldHalfedge *h) {
  return h->startVert < h->endVert;
}

static inline bool manifold_triref_same_face(const ManifoldTriRef *a,
                                              const ManifoldTriRef *b) {
  return a->meshID == b->meshID && a->coplanarID == b->coplanarID &&
         a->faceID == b->faceID;
}

static inline ManifoldTmpEdge manifold_tmpedge(int start, int end, int idx) {
  ManifoldTmpEdge e;
  e.first = start < end ? start : end;
  e.second = start < end ? end : start;
  e.halfedgeIdx = idx;
  return e;
}

#endif // MANIFOLD_TYPES_H
