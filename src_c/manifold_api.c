// Copyright 2021 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Public API implementation for the C11 Manifold port.

#include "manifold_api.h"
#include <string.h>

// Quality settings (global)
static ManifoldQuality g_quality = {0, 10.0, 1.0};

int manifold_get_circular_segments(double radius) {
  if (g_quality.circularSegments > 0) return g_quality.circularSegments;
  int byAngle = (int)(360.0 / g_quality.minCircularAngle);
  int byLength = (int)(2.0 * MANIFOLD_PI * radius / g_quality.minCircularEdgeLength);
  int n = byAngle < byLength ? byAngle : byLength;
  n = ((n + 3) / 4) * 4;  // round up to multiple of 4
  if (n < 3) n = 3;
  return n;
}

void manifold_set_circular_segments(int n) { g_quality.circularSegments = n; }
void manifold_set_min_circular_angle(double a) { g_quality.minCircularAngle = a; }
void manifold_set_min_circular_edge_length(double l) { g_quality.minCircularEdgeLength = l; }
void manifold_quality_reset(void) {
  g_quality.circularSegments = 0;
  g_quality.minCircularAngle = 10.0;
  g_quality.minCircularEdgeLength = 1.0;
}

void manifold_create(Manifold *m) {
  manifold_impl_init(&m->impl);
}

void manifold_destroy(Manifold *m) {
  manifold_impl_free(&m->impl);
}

void manifold_copy(Manifold *dst, const Manifold *src) {
  manifold_impl_init(&dst->impl);
  dst->impl.bBox = src->impl.bBox;
  dst->impl.epsilon = src->impl.epsilon;
  dst->impl.tolerance = src->impl.tolerance;
  dst->impl.numProp = src->impl.numProp;
  dst->impl.status = src->impl.status;
  dst->impl.vertPos = vec_vec3_copy(&src->impl.vertPos);
  dst->impl.halfedge = vec_halfedge_copy(&src->impl.halfedge);
  dst->impl.properties = vec_double_copy(&src->impl.properties);
  dst->impl.vertNormal = vec_vec3_copy(&src->impl.vertNormal);
  dst->impl.faceNormal = vec_vec3_copy(&src->impl.faceNormal);
  dst->impl.halfedgeTangent = vec_vec4_copy(&src->impl.halfedgeTangent);
  dst->impl.meshRelation.originalID = src->impl.meshRelation.originalID;
  dst->impl.meshRelation.meshIDtransform = vec_meshid_copy(&src->impl.meshRelation.meshIDtransform);
  dst->impl.meshRelation.triRef = vec_triref_copy(&src->impl.meshRelation.triRef);
  // Don't copy collider - it will be rebuilt lazily
  dst->impl.colliderBuilt = false;
}

ManifoldError manifold_status(const Manifold *m) { return m->impl.status; }
bool manifold_is_empty(const Manifold *m) { return manifold_impl_is_empty(&m->impl); }
size_t manifold_num_vert(const Manifold *m) { return manifold_impl_num_vert(&m->impl); }
size_t manifold_num_edge(const Manifold *m) { return manifold_impl_num_edge(&m->impl); }
size_t manifold_num_tri(const Manifold *m) { return manifold_impl_num_tri(&m->impl); }
ManifoldBox manifold_bounding_box(const Manifold *m) { return m->impl.bBox; }
double manifold_volume(const Manifold *m) { return manifold_impl_get_volume(&m->impl); }
double manifold_surface_area(const Manifold *m) { return manifold_impl_get_surface_area(&m->impl); }

Manifold manifold_tetrahedron(void) {
  Manifold m;
  manifold_impl_tetrahedron(&m.impl);
  return m;
}

Manifold manifold_cube(ManifoldVec3 size, bool center) {
  Manifold m;
  if (size.x < 0 || size.y < 0 || size.z < 0 ||
      (size.x == 0 && size.y == 0 && size.z == 0)) {
    manifold_impl_init(&m.impl);
    manifold_impl_make_empty(&m.impl, MANIFOLD_ERROR_INVALID_CONSTRUCTION);
    return m;
  }
  ManifoldVec3 offset = center ? vec3_scale(size, -0.5) : manifold_vec3(0, 0, 0);
  ManifoldMat3x4 transform;
  transform.cols[0] = manifold_vec3(size.x, 0, 0);
  transform.cols[1] = manifold_vec3(0, size.y, 0);
  transform.cols[2] = manifold_vec3(0, 0, size.z);
  transform.cols[3] = offset;
  manifold_impl_cube(&m.impl, transform);
  return m;
}

Manifold manifold_translate(const Manifold *m, ManifoldVec3 v) {
  Manifold out;
  manifold_copy(&out, m);
  for (size_t i = 0; i < out.impl.vertPos.len; i++) {
    out.impl.vertPos.data[i] = vec3_add(out.impl.vertPos.data[i], v);
  }
  out.impl.bBox = manifold_box_shift(out.impl.bBox, v);
  return out;
}

Manifold manifold_scale(const Manifold *m, ManifoldVec3 v) {
  Manifold out;
  manifold_copy(&out, m);
  for (size_t i = 0; i < out.impl.vertPos.len; i++) {
    out.impl.vertPos.data[i] = vec3_mul(out.impl.vertPos.data[i], v);
  }
  manifold_impl_calculate_bbox(&out.impl);
  return out;
}

Manifold manifold_rotate(const Manifold *m, double xDeg, double yDeg,
                         double zDeg) {
  double cx = manifold_cosd(xDeg), sx = manifold_sind(xDeg);
  double cy = manifold_cosd(yDeg), sy = manifold_sind(yDeg);
  double cz = manifold_cosd(zDeg), sz = manifold_sind(zDeg);
  ManifoldMat3 rx, ry, rz;
  rx.cols[0] = manifold_vec3(1, 0, 0);
  rx.cols[1] = manifold_vec3(0, cx, sx);
  rx.cols[2] = manifold_vec3(0, -sx, cx);
  ry.cols[0] = manifold_vec3(cy, 0, -sy);
  ry.cols[1] = manifold_vec3(0, 1, 0);
  ry.cols[2] = manifold_vec3(sy, 0, cy);
  rz.cols[0] = manifold_vec3(cz, sz, 0);
  rz.cols[1] = manifold_vec3(-sz, cz, 0);
  rz.cols[2] = manifold_vec3(0, 0, 1);
  ManifoldMat3 rot = mat3_mul(rz, mat3_mul(ry, rx));

  Manifold out;
  manifold_copy(&out, m);
  for (size_t i = 0; i < out.impl.vertPos.len; i++) {
    out.impl.vertPos.data[i] = mat3_mul_vec3(rot, out.impl.vertPos.data[i]);
  }
  manifold_impl_calculate_bbox(&out.impl);
  return out;
}

Manifold manifold_transform(const Manifold *m, ManifoldMat3x4 t) {
  Manifold out;
  manifold_copy(&out, m);
  for (size_t i = 0; i < out.impl.vertPos.len; i++) {
    out.impl.vertPos.data[i] = mat3x4_transform_point(t, out.impl.vertPos.data[i]);
  }
  manifold_impl_calculate_bbox(&out.impl);
  return out;
}

// Boolean operations - stub implementation
// Full boolean implementation requires boolean3.c and boolean_result.c
Manifold manifold_boolean(const Manifold *a, const Manifold *b,
                          ManifoldOpType op) {
  (void)op;
  // TODO: implement boolean operations
  Manifold result;
  manifold_copy(&result, a);
  (void)b;
  return result;
}

Manifold manifold_union(const Manifold *a, const Manifold *b) {
  return manifold_boolean(a, b, MANIFOLD_OP_ADD);
}
Manifold manifold_difference(const Manifold *a, const Manifold *b) {
  return manifold_boolean(a, b, MANIFOLD_OP_SUBTRACT);
}
Manifold manifold_intersection(const Manifold *a, const Manifold *b) {
  return manifold_boolean(a, b, MANIFOLD_OP_INTERSECT);
}
