// Copyright 2021 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Vector and matrix math operations for the C11 Manifold port.

#ifndef MANIFOLD_VEC_MATH_H
#define MANIFOLD_VEC_MATH_H

#include "manifold_types.h"

// ============== Vec2 operations ==============

static inline ManifoldVec2 vec2_add(ManifoldVec2 a, ManifoldVec2 b) {
  return manifold_vec2(a.x + b.x, a.y + b.y);
}
static inline ManifoldVec2 vec2_sub(ManifoldVec2 a, ManifoldVec2 b) {
  return manifold_vec2(a.x - b.x, a.y - b.y);
}
static inline ManifoldVec2 vec2_mul(ManifoldVec2 a, ManifoldVec2 b) {
  return manifold_vec2(a.x * b.x, a.y * b.y);
}
static inline ManifoldVec2 vec2_div(ManifoldVec2 a, ManifoldVec2 b) {
  return manifold_vec2(a.x / b.x, a.y / b.y);
}
static inline ManifoldVec2 vec2_scale(ManifoldVec2 a, double s) {
  return manifold_vec2(a.x * s, a.y * s);
}
static inline ManifoldVec2 vec2_neg(ManifoldVec2 a) {
  return manifold_vec2(-a.x, -a.y);
}
static inline double vec2_dot(ManifoldVec2 a, ManifoldVec2 b) {
  return a.x * b.x + a.y * b.y;
}
static inline double vec2_cross(ManifoldVec2 a, ManifoldVec2 b) {
  return a.x * b.y - a.y * b.x;
}
static inline double vec2_length(ManifoldVec2 a) {
  return sqrt(vec2_dot(a, a));
}
static inline ManifoldVec2 vec2_normalize(ManifoldVec2 a) {
  double len = vec2_length(a);
  if (len == 0.0) return manifold_vec2(0.0, 0.0);
  return vec2_scale(a, 1.0 / len);
}
static inline ManifoldVec2 vec2_min(ManifoldVec2 a, ManifoldVec2 b) {
  return manifold_vec2(fmin(a.x, b.x), fmin(a.y, b.y));
}
static inline ManifoldVec2 vec2_max(ManifoldVec2 a, ManifoldVec2 b) {
  return manifold_vec2(fmax(a.x, b.x), fmax(a.y, b.y));
}
static inline ManifoldVec2 vec2_abs(ManifoldVec2 a) {
  return manifold_vec2(fabs(a.x), fabs(a.y));
}
static inline bool vec2_all_gequal(ManifoldVec2 a, ManifoldVec2 b) {
  return a.x >= b.x && a.y >= b.y;
}
static inline bool vec2_isfinite(ManifoldVec2 a) {
  return isfinite(a.x) && isfinite(a.y);
}
static inline double vec2_max_comp(ManifoldVec2 a) {
  return fmax(a.x, a.y);
}

// ============== Vec3 operations ==============

static inline ManifoldVec3 vec3_add(ManifoldVec3 a, ManifoldVec3 b) {
  return manifold_vec3(a.x + b.x, a.y + b.y, a.z + b.z);
}
static inline ManifoldVec3 vec3_sub(ManifoldVec3 a, ManifoldVec3 b) {
  return manifold_vec3(a.x - b.x, a.y - b.y, a.z - b.z);
}
static inline ManifoldVec3 vec3_mul(ManifoldVec3 a, ManifoldVec3 b) {
  return manifold_vec3(a.x * b.x, a.y * b.y, a.z * b.z);
}
static inline ManifoldVec3 vec3_div(ManifoldVec3 a, ManifoldVec3 b) {
  return manifold_vec3(a.x / b.x, a.y / b.y, a.z / b.z);
}
static inline ManifoldVec3 vec3_scale(ManifoldVec3 a, double s) {
  return manifold_vec3(a.x * s, a.y * s, a.z * s);
}
static inline ManifoldVec3 vec3_neg(ManifoldVec3 a) {
  return manifold_vec3(-a.x, -a.y, -a.z);
}
static inline double vec3_dot(ManifoldVec3 a, ManifoldVec3 b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}
static inline ManifoldVec3 vec3_cross(ManifoldVec3 a, ManifoldVec3 b) {
  return manifold_vec3(a.y * b.z - a.z * b.y,
                       a.z * b.x - a.x * b.z,
                       a.x * b.y - a.y * b.x);
}
static inline double vec3_length(ManifoldVec3 a) {
  return sqrt(vec3_dot(a, a));
}
static inline ManifoldVec3 vec3_normalize(ManifoldVec3 a) {
  double len = vec3_length(a);
  if (len == 0.0) return manifold_vec3(0.0, 0.0, 0.0);
  return vec3_scale(a, 1.0 / len);
}
static inline ManifoldVec3 vec3_min(ManifoldVec3 a, ManifoldVec3 b) {
  return manifold_vec3(fmin(a.x, b.x), fmin(a.y, b.y), fmin(a.z, b.z));
}
static inline ManifoldVec3 vec3_max(ManifoldVec3 a, ManifoldVec3 b) {
  return manifold_vec3(fmax(a.x, b.x), fmax(a.y, b.y), fmax(a.z, b.z));
}
static inline ManifoldVec3 vec3_abs(ManifoldVec3 a) {
  return manifold_vec3(fabs(a.x), fabs(a.y), fabs(a.z));
}
static inline bool vec3_equal(ManifoldVec3 a, ManifoldVec3 b) {
  return a.x == b.x && a.y == b.y && a.z == b.z;
}
static inline bool vec3_all_gequal(ManifoldVec3 a, ManifoldVec3 b) {
  return a.x >= b.x && a.y >= b.y && a.z >= b.z;
}
static inline bool vec3_isfinite(ManifoldVec3 a) {
  return isfinite(a.x) && isfinite(a.y) && isfinite(a.z);
}
static inline double vec3_max_comp(ManifoldVec3 a) {
  return fmax(a.x, fmax(a.y, a.z));
}
static inline double vec3_get(ManifoldVec3 v, int i) {
  return i == 0 ? v.x : (i == 1 ? v.y : v.z);
}
static inline ManifoldVec3 vec3_set(ManifoldVec3 v, int i, double val) {
  if (i == 0) v.x = val;
  else if (i == 1) v.y = val;
  else v.z = val;
  return v;
}

// ============== Vec4 operations ==============

static inline ManifoldVec4 vec4_add(ManifoldVec4 a, ManifoldVec4 b) {
  return manifold_vec4(a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w);
}
static inline ManifoldVec4 vec4_sub(ManifoldVec4 a, ManifoldVec4 b) {
  return manifold_vec4(a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w);
}
static inline ManifoldVec4 vec4_scale(ManifoldVec4 a, double s) {
  return manifold_vec4(a.x * s, a.y * s, a.z * s, a.w * s);
}
static inline ManifoldVec4 vec4_neg(ManifoldVec4 a) {
  return manifold_vec4(-a.x, -a.y, -a.z, -a.w);
}
static inline double vec4_dot(ManifoldVec4 a, ManifoldVec4 b) {
  return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
static inline double vec4_length(ManifoldVec4 a) {
  return sqrt(vec4_dot(a, a));
}
static inline ManifoldVec4 vec4_normalize(ManifoldVec4 a) {
  double len = vec4_length(a);
  if (len == 0.0) return manifold_vec4(0, 0, 0, 0);
  return vec4_scale(a, 1.0 / len);
}
static inline double vec4_get(ManifoldVec4 v, int i) {
  switch (i) {
    case 0: return v.x;
    case 1: return v.y;
    case 2: return v.z;
    default: return v.w;
  }
}
static inline ManifoldVec4 vec4_set(ManifoldVec4 v, int i, double val) {
  switch (i) {
    case 0: v.x = val; break;
    case 1: v.y = val; break;
    case 2: v.z = val; break;
    default: v.w = val; break;
  }
  return v;
}

// ============== IVec3 operations ==============

static inline ManifoldIVec3 ivec3_add(ManifoldIVec3 a, ManifoldIVec3 b) {
  return manifold_ivec3(a.x + b.x, a.y + b.y, a.z + b.z);
}
static inline int ivec3_get(ManifoldIVec3 v, int i) {
  return i == 0 ? v.x : (i == 1 ? v.y : v.z);
}

// ============== Vec3 ↔ Vec4 conversions ==============

static inline ManifoldVec4 vec3_to_vec4(ManifoldVec3 v, double w) {
  return manifold_vec4(v.x, v.y, v.z, w);
}
static inline ManifoldVec3 vec4_to_vec3(ManifoldVec4 v) {
  return manifold_vec3(v.x, v.y, v.z);
}
static inline ManifoldVec3 vec2_to_vec3(ManifoldVec2 v, double z) {
  return manifold_vec3(v.x, v.y, z);
}

// ============== Mat3 operations ==============

static inline ManifoldMat3 mat3_identity(void) {
  ManifoldMat3 m;
  m.cols[0] = manifold_vec3(1, 0, 0);
  m.cols[1] = manifold_vec3(0, 1, 0);
  m.cols[2] = manifold_vec3(0, 0, 1);
  return m;
}

static inline ManifoldMat3x4 mat3x4_identity(void) {
  ManifoldMat3x4 m;
  m.cols[0] = manifold_vec3(1, 0, 0);
  m.cols[1] = manifold_vec3(0, 1, 0);
  m.cols[2] = manifold_vec3(0, 0, 1);
  m.cols[3] = manifold_vec3(0, 0, 0);
  return m;
}

// mat3 * vec3
static inline ManifoldVec3 mat3_mul_vec3(ManifoldMat3 m, ManifoldVec3 v) {
  return vec3_add(vec3_add(vec3_scale(m.cols[0], v.x),
                           vec3_scale(m.cols[1], v.y)),
                  vec3_scale(m.cols[2], v.z));
}

// mat3x4 * vec4 → vec3 (affine transform)
static inline ManifoldVec3 mat3x4_mul_vec4(ManifoldMat3x4 m, ManifoldVec4 v) {
  return vec3_add(vec3_add(vec3_add(vec3_scale(m.cols[0], v.x),
                                    vec3_scale(m.cols[1], v.y)),
                           vec3_scale(m.cols[2], v.z)),
                  vec3_scale(m.cols[3], v.w));
}

// mat3x4 transform point (w=1)
static inline ManifoldVec3 mat3x4_transform_point(ManifoldMat3x4 m,
                                                    ManifoldVec3 p) {
  return mat3x4_mul_vec4(m, vec3_to_vec4(p, 1.0));
}

// mat3x4 transform vector (w=0)
static inline ManifoldVec3 mat3x4_transform_vec(ManifoldMat3x4 m,
                                                  ManifoldVec3 v) {
  return mat3x4_mul_vec4(m, vec3_to_vec4(v, 0.0));
}

// Extract the 3x3 rotation part from mat3x4
static inline ManifoldMat3 mat3x4_to_mat3(ManifoldMat3x4 m) {
  ManifoldMat3 r;
  r.cols[0] = m.cols[0];
  r.cols[1] = m.cols[1];
  r.cols[2] = m.cols[2];
  return r;
}

// mat3 transpose
static inline ManifoldMat3 mat3_transpose(ManifoldMat3 m) {
  ManifoldMat3 r;
  r.cols[0] = manifold_vec3(m.cols[0].x, m.cols[1].x, m.cols[2].x);
  r.cols[1] = manifold_vec3(m.cols[0].y, m.cols[1].y, m.cols[2].y);
  r.cols[2] = manifold_vec3(m.cols[0].z, m.cols[1].z, m.cols[2].z);
  return r;
}

// mat3 determinant
static inline double mat3_det(ManifoldMat3 m) {
  return m.cols[0].x * (m.cols[1].y * m.cols[2].z - m.cols[1].z * m.cols[2].y) -
         m.cols[1].x * (m.cols[0].y * m.cols[2].z - m.cols[0].z * m.cols[2].y) +
         m.cols[2].x * (m.cols[0].y * m.cols[1].z - m.cols[0].z * m.cols[1].y);
}

// mat3 inverse
static inline ManifoldMat3 mat3_inverse(ManifoldMat3 m) {
  double det = mat3_det(m);
  double inv_det = 1.0 / det;
  ManifoldMat3 r;
  r.cols[0] = vec3_scale(manifold_vec3(
    m.cols[1].y * m.cols[2].z - m.cols[2].y * m.cols[1].z,
    m.cols[2].y * m.cols[0].z - m.cols[0].y * m.cols[2].z,
    m.cols[0].y * m.cols[1].z - m.cols[1].y * m.cols[0].z), inv_det);
  r.cols[1] = vec3_scale(manifold_vec3(
    m.cols[2].x * m.cols[1].z - m.cols[1].x * m.cols[2].z,
    m.cols[0].x * m.cols[2].z - m.cols[2].x * m.cols[0].z,
    m.cols[1].x * m.cols[0].z - m.cols[0].x * m.cols[1].z), inv_det);
  r.cols[2] = vec3_scale(manifold_vec3(
    m.cols[1].x * m.cols[2].y - m.cols[2].x * m.cols[1].y,
    m.cols[2].x * m.cols[0].y - m.cols[0].x * m.cols[2].y,
    m.cols[0].x * m.cols[1].y - m.cols[1].x * m.cols[0].y), inv_det);
  return r;
}

// mat3 * mat3
static inline ManifoldMat3 mat3_mul(ManifoldMat3 a, ManifoldMat3 b) {
  ManifoldMat3 r;
  r.cols[0] = mat3_mul_vec3(a, b.cols[0]);
  r.cols[1] = mat3_mul_vec3(a, b.cols[1]);
  r.cols[2] = mat3_mul_vec3(a, b.cols[2]);
  return r;
}

// mat2x3 * vec3 → vec2
static inline ManifoldVec2 mat2x3_mul_vec3(ManifoldMat2x3 m, ManifoldVec3 v) {
  return vec2_add(vec2_add(vec2_scale(m.cols[0], v.x),
                           vec2_scale(m.cols[1], v.y)),
                  vec2_scale(m.cols[2], v.z));
}

// mat2 operations
static inline ManifoldMat2 mat2_identity(void) {
  ManifoldMat2 m;
  m.cols[0] = manifold_vec2(1, 0);
  m.cols[1] = manifold_vec2(0, 1);
  return m;
}

static inline double mat2_det(ManifoldMat2 m) {
  return m.cols[0].x * m.cols[1].y - m.cols[0].y * m.cols[1].x;
}

// ============== Box operations ==============

static inline ManifoldBox manifold_box(ManifoldVec3 p1, ManifoldVec3 p2) {
  ManifoldBox b;
  b.min = vec3_min(p1, p2);
  b.max = vec3_max(p1, p2);
  return b;
}

static inline ManifoldVec3 manifold_box_size(ManifoldBox b) {
  return vec3_sub(b.max, b.min);
}

static inline ManifoldVec3 manifold_box_center(ManifoldBox b) {
  return vec3_scale(vec3_add(b.max, b.min), 0.5);
}

static inline double manifold_box_scale(ManifoldBox b) {
  ManifoldVec3 absMax = vec3_max(vec3_abs(b.min), vec3_abs(b.max));
  return fmax(absMax.x, fmax(absMax.y, absMax.z));
}

static inline bool manifold_box_contains_point(ManifoldBox b, ManifoldVec3 p) {
  return vec3_all_gequal(p, b.min) && vec3_all_gequal(b.max, p);
}

static inline bool manifold_box_contains_box(ManifoldBox b, ManifoldBox box) {
  return vec3_all_gequal(box.min, b.min) && vec3_all_gequal(b.max, box.max);
}

static inline void manifold_box_union_point(ManifoldBox *b, ManifoldVec3 p) {
  b->min = vec3_min(b->min, p);
  b->max = vec3_max(b->max, p);
}

static inline ManifoldBox manifold_box_union(ManifoldBox a, ManifoldBox b) {
  ManifoldBox out;
  out.min = vec3_min(a.min, b.min);
  out.max = vec3_max(a.max, b.max);
  return out;
}

static inline ManifoldBox manifold_box_transform(ManifoldBox b,
                                                  ManifoldMat3x4 t) {
  ManifoldBox out;
  ManifoldVec3 minT = mat3x4_transform_point(t, b.min);
  ManifoldVec3 maxT = mat3x4_transform_point(t, b.max);
  out.min = vec3_min(minT, maxT);
  out.max = vec3_max(minT, maxT);
  return out;
}

static inline ManifoldBox manifold_box_shift(ManifoldBox b, ManifoldVec3 s) {
  ManifoldBox out;
  out.min = vec3_add(b.min, s);
  out.max = vec3_add(b.max, s);
  return out;
}

static inline ManifoldBox manifold_box_scale_vec(ManifoldBox b, ManifoldVec3 s) {
  ManifoldBox out;
  out.min = vec3_mul(b.min, s);
  out.max = vec3_mul(b.max, s);
  return out;
}

static inline bool manifold_box_overlaps(ManifoldBox a, ManifoldBox b) {
  return a.min.x <= b.max.x && a.min.y <= b.max.y && a.min.z <= b.max.z &&
         a.max.x >= b.min.x && a.max.y >= b.min.y && a.max.z >= b.min.z;
}

static inline bool manifold_box_overlaps_point(ManifoldBox b, ManifoldVec3 p) {
  return p.x <= b.max.x && p.x >= b.min.x && p.y <= b.max.y && p.y >= b.min.y;
}

static inline bool manifold_box_isfinite(ManifoldBox b) {
  return vec3_isfinite(b.min) && vec3_isfinite(b.max);
}

// ============== Rect operations ==============

static inline ManifoldRect manifold_rect(ManifoldVec2 a, ManifoldVec2 b) {
  ManifoldRect r;
  r.min = vec2_min(a, b);
  r.max = vec2_max(a, b);
  return r;
}

static inline ManifoldVec2 manifold_rect_size(ManifoldRect r) {
  return vec2_sub(r.max, r.min);
}

static inline double manifold_rect_area(ManifoldRect r) {
  ManifoldVec2 sz = manifold_rect_size(r);
  return sz.x * sz.y;
}

static inline double manifold_rect_scale(ManifoldRect r) {
  ManifoldVec2 absMax = vec2_max(vec2_abs(r.min), vec2_abs(r.max));
  return fmax(absMax.x, absMax.y);
}

static inline ManifoldVec2 manifold_rect_center(ManifoldRect r) {
  return vec2_scale(vec2_add(r.max, r.min), 0.5);
}

static inline bool manifold_rect_overlaps(ManifoldRect a, ManifoldRect b) {
  return a.min.x <= b.max.x && a.min.y <= b.max.y &&
         a.max.x >= b.min.x && a.max.y >= b.min.y;
}

static inline bool manifold_rect_isfinite(ManifoldRect r) {
  return vec2_isfinite(r.min) && vec2_isfinite(r.max);
}

static inline bool manifold_rect_is_empty(ManifoldRect r) {
  return r.max.y <= r.min.y || r.max.x <= r.min.x;
}

// vec3 outer product: result[col][row] = a[row] * b[col]
static inline ManifoldMat3 vec3_outerprod(ManifoldVec3 a, ManifoldVec3 b) {
  ManifoldMat3 m;
  m.cols[0] = vec3_scale(a, b.x);
  m.cols[1] = vec3_scale(a, b.y);
  m.cols[2] = vec3_scale(a, b.z);
  return m;
}

// mat3 subtract
static inline ManifoldMat3 mat3_sub(ManifoldMat3 a, ManifoldMat3 b) {
  ManifoldMat3 r;
  r.cols[0] = vec3_sub(a.cols[0], b.cols[0]);
  r.cols[1] = vec3_sub(a.cols[1], b.cols[1]);
  r.cols[2] = vec3_sub(a.cols[2], b.cols[2]);
  return r;
}

// mat3 scale (scalar * mat3)
static inline ManifoldMat3 mat3_scale_s(ManifoldMat3 m, double s) {
  ManifoldMat3 r;
  r.cols[0] = vec3_scale(m.cols[0], s);
  r.cols[1] = vec3_scale(m.cols[1], s);
  r.cols[2] = vec3_scale(m.cols[2], s);
  return r;
}

// mat3x4 from mat3 + translation
static inline ManifoldMat3x4 mat3x4_from_mat3_translate(ManifoldMat3 m,
                                                          ManifoldVec3 t) {
  ManifoldMat3x4 r;
  r.cols[0] = m.cols[0];
  r.cols[1] = m.cols[1];
  r.cols[2] = m.cols[2];
  r.cols[3] = t;
  return r;
}

// ============== Scalar helpers ==============

static inline double manifold_radians(double a) {
  return a * MANIFOLD_PI / 180.0;
}

static inline double manifold_degrees(double a) {
  return a * 180.0 / MANIFOLD_PI;
}

static inline double manifold_smoothstep(double edge0, double edge1, double a) {
  double t = (a - edge0) / (edge1 - edge0);
  if (t < 0.0) t = 0.0;
  if (t > 1.0) t = 1.0;
  return t * t * (3.0 - 2.0 * t);
}

static inline double manifold_sind(double x) {
  if (!isfinite(x)) return sin(x);
  if (x < 0.0) return -manifold_sind(-x);
  int quo;
  x = remquo(fabs(x), 90.0, &quo);
  switch (quo % 4) {
    case 0: return sin(manifold_radians(x));
    case 1: return cos(manifold_radians(x));
    case 2: return -sin(manifold_radians(x));
    case 3: return -cos(manifold_radians(x));
  }
  return 0.0;
}

static inline double manifold_cosd(double x) {
  return manifold_sind(x + 90.0);
}

// ============== Utility math from utils.h ==============

static inline ManifoldVec3 manifold_safe_normalize(ManifoldVec3 v) {
  v = vec3_normalize(v);
  return isfinite(v.x) ? v : manifold_vec3(0.0, 0.0, 0.0);
}

static inline double manifold_max_epsilon(double minEpsilon, ManifoldBox bBox) {
  double epsilon = fmax(minEpsilon, MANIFOLD_PRECISION * manifold_box_scale(bBox));
  return isfinite(epsilon) ? epsilon : -1.0;
}

static inline int manifold_next3(int i) {
  static const int next3[3] = {1, 2, 0};
  return next3[i];
}

static inline int manifold_prev3(int i) {
  static const int prev3[3] = {2, 0, 1};
  return prev3[i];
}

static inline int manifold_next_halfedge(int current) {
  ++current;
  if (current % 3 == 0) current -= 3;
  return current;
}

static inline ManifoldMat3 manifold_normal_transform(ManifoldMat3x4 transform) {
  ManifoldMat3 rot = mat3x4_to_mat3(transform);
  return mat3_inverse(mat3_transpose(rot));
}

static inline int manifold_ccw(ManifoldVec2 p0, ManifoldVec2 p1,
                                ManifoldVec2 p2, double tol) {
  ManifoldVec2 v1 = vec2_sub(p1, p0);
  ManifoldVec2 v2 = vec2_sub(p2, p0);
  double area = v1.x * v2.y - v1.y * v2.x;
  double base2 = fmax(vec2_dot(v1, v1), vec2_dot(v2, v2));
  if (area * area * 4.0 <= base2 * tol * tol)
    return 0;
  return area > 0 ? 1 : -1;
}

static inline uint64_t manifold_hash64bit(uint64_t x) {
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  x = x ^ (x >> 31);
  return x;
}

// ============== Axis-aligned projection ==============

static inline ManifoldMat2x3 manifold_get_axis_aligned_projection(
    ManifoldVec3 normal) {
  ManifoldVec3 absNormal = vec3_abs(normal);
  double xyzMax;
  // mat3x2 = 2 columns of vec3
  ManifoldMat3x2 projection;
  if (absNormal.z > absNormal.x && absNormal.z > absNormal.y) {
    projection.cols[0] = manifold_vec3(1.0, 0.0, 0.0);
    projection.cols[1] = manifold_vec3(0.0, 1.0, 0.0);
    xyzMax = normal.z;
  } else if (absNormal.y > absNormal.x) {
    projection.cols[0] = manifold_vec3(0.0, 0.0, 1.0);
    projection.cols[1] = manifold_vec3(1.0, 0.0, 0.0);
    xyzMax = normal.y;
  } else {
    projection.cols[0] = manifold_vec3(0.0, 1.0, 0.0);
    projection.cols[1] = manifold_vec3(0.0, 0.0, 1.0);
    xyzMax = normal.x;
  }
  if (xyzMax < 0) projection.cols[0] = vec3_scale(projection.cols[0], -1.0);

  // transpose mat3x2 (2 cols of vec3) → mat2x3 (3 cols of vec2)
  ManifoldMat2x3 result;
  result.cols[0] = manifold_vec2(projection.cols[0].x, projection.cols[1].x);
  result.cols[1] = manifold_vec2(projection.cols[0].y, projection.cols[1].y);
  result.cols[2] = manifold_vec2(projection.cols[0].z, projection.cols[1].z);
  return result;
}

// ============== Barycentric coordinates ==============

static inline ManifoldVec3 manifold_get_barycentric(ManifoldVec3 v,
    ManifoldVec3 triPos[3], double tolerance) {
  ManifoldVec3 edges[3];
  edges[0] = vec3_sub(triPos[2], triPos[1]);
  edges[1] = vec3_sub(triPos[0], triPos[2]);
  edges[2] = vec3_sub(triPos[1], triPos[0]);

  double d2[3];
  d2[0] = vec3_dot(edges[0], edges[0]);
  d2[1] = vec3_dot(edges[1], edges[1]);
  d2[2] = vec3_dot(edges[2], edges[2]);

  int longSide = (d2[0] > d2[1] && d2[0] > d2[2]) ? 0 : (d2[1] > d2[2] ? 1 : 2);
  ManifoldVec3 crossP = vec3_cross(edges[0], edges[1]);
  double area2 = vec3_dot(crossP, crossP);
  double tol2 = tolerance * tolerance;

  ManifoldVec3 uvw = manifold_vec3(0, 0, 0);
  for (int i = 0; i < 3; ++i) {
    ManifoldVec3 dv = vec3_sub(v, triPos[i]);
    if (vec3_dot(dv, dv) < tol2) {
      uvw = manifold_vec3(0, 0, 0);
      uvw = vec3_set(uvw, i, 1.0);
      return uvw;
    }
  }

  if (d2[longSide] < tol2) {
    return manifold_vec3(1, 0, 0);
  } else if (area2 > d2[longSide] * tol2) {
    for (int i = 0; i < 3; ++i) {
      int j = manifold_next3(i);
      ManifoldVec3 crossPv = vec3_cross(edges[i], vec3_sub(v, triPos[j]));
      double area2v = vec3_dot(crossPv, crossPv);
      double val = (area2v < d2[i] * tol2) ? 0.0 : vec3_dot(crossPv, crossP);
      uvw = vec3_set(uvw, i, val);
    }
    double sum = vec3_get(uvw, 0) + vec3_get(uvw, 1) + vec3_get(uvw, 2);
    uvw = vec3_scale(uvw, 1.0 / sum);
    return uvw;
  } else {
    int nextV = manifold_next3(longSide);
    double alpha = vec3_dot(vec3_sub(v, triPos[nextV]), edges[longSide]) /
                   d2[longSide];
    uvw = vec3_set(uvw, longSide, 0.0);
    uvw = vec3_set(uvw, nextV, 1.0 - alpha);
    int lastV = manifold_next3(nextV);
    uvw = vec3_set(uvw, lastV, alpha);
    return uvw;
  }
}

// ============== SVD (Singular Value Decomposition) ==============
// Ported from src/svd.h

#define SVD_GAMMA 5.82842712474619
#define SVD_CSTAR 0.9238795325112867
#define SVD_SSTAR 0.3826834323650898
#define SVD_EPSILON 1e-6
#define SVD_JACOBI_STEPS 12

typedef struct { double m_00, m_10, m_11, m_20, m_21, m_22; } Symmetric3x3;
typedef struct { double ch, sh; } Givens;
typedef struct { ManifoldMat3 Q, R; } QR;
typedef struct { ManifoldMat3 U, S, V; } SVDSet;

static inline void svd_cond_swap(bool c, double *X, double *Y) {
  double Z = *X;
  *X = c ? *Y : *X;
  *Y = c ? Z : *Y;
}

static inline void svd_cond_neg_swap(bool c, double *X, double *Y) {
  double Z = -*X;
  *X = c ? *Y : *X;
  *Y = c ? Z : *Y;
}

static inline double svd_dist2(ManifoldVec3 v) { return vec3_dot(v, v); }

static inline Givens svd_approx_givens(Symmetric3x3 *A) {
  Givens g = {2.0 * (A->m_00 - A->m_11), A->m_10};
  bool b = SVD_GAMMA * g.sh * g.sh < g.ch * g.ch;
  double w = 1.0 / hypot(g.ch, g.sh);
  if (!isfinite(w)) b = false;
  return (Givens){b ? w * g.ch : SVD_CSTAR, b ? w * g.sh : SVD_SSTAR};
}

static inline void svd_jacobi_conjugation(int x, int y, int z,
                                            Symmetric3x3 *S, double q[4]) {
  Givens g = svd_approx_givens(S);
  double scale = 1.0 / (g.ch * g.ch + g.sh * g.sh);
  double a = (g.ch * g.ch - g.sh * g.sh) * scale;
  double b = 2.0 * g.sh * g.ch * scale;
  Symmetric3x3 _S = *S;
  S->m_00 = a * (a * _S.m_00 + b * _S.m_10) + b * (a * _S.m_10 + b * _S.m_11);
  S->m_10 = a * (-b * _S.m_00 + a * _S.m_10) + b * (-b * _S.m_10 + a * _S.m_11);
  S->m_11 = -b * (-b * _S.m_00 + a * _S.m_10) + a * (-b * _S.m_10 + a * _S.m_11);
  S->m_20 = a * _S.m_20 + b * _S.m_21;
  S->m_21 = -b * _S.m_20 + a * _S.m_21;
  S->m_22 = _S.m_22;
  double tmp[3] = {g.sh * q[0], g.sh * q[1], g.sh * q[2]};
  double shw = g.sh * q[3];
  q[z] = q[z] * g.ch + shw;
  q[3] = q[3] * g.ch + (-tmp[z]);
  q[x] = q[x] * g.ch + tmp[y];
  q[y] = q[y] * g.ch + (-tmp[x]);
  _S.m_00 = S->m_11; _S.m_10 = S->m_21; _S.m_11 = S->m_22;
  _S.m_20 = S->m_10; _S.m_21 = S->m_20; _S.m_22 = S->m_00;
  *S = _S;
}

static inline ManifoldMat3 svd_jacobi_eigen(Symmetric3x3 S) {
  double q[4] = {0, 0, 0, 1};
  for (int i = 0; i < SVD_JACOBI_STEPS; i++) {
    svd_jacobi_conjugation(0, 1, 2, &S, q);
    svd_jacobi_conjugation(1, 2, 0, &S, q);
    svd_jacobi_conjugation(2, 0, 1, &S, q);
  }
  ManifoldMat3 V;
  V.cols[0] = manifold_vec3(1 - 2*(q[1]*q[1] + q[2]*q[2]),
                        2*(q[0]*q[1] + q[3]*q[2]),
                        2*(q[0]*q[2] - q[3]*q[1]));
  V.cols[1] = manifold_vec3(2*(q[0]*q[1] - q[3]*q[2]),
                        1 - 2*(q[0]*q[0] + q[2]*q[2]),
                        2*(q[1]*q[2] + q[3]*q[0]));
  V.cols[2] = manifold_vec3(2*(q[0]*q[2] + q[3]*q[1]),
                        2*(q[1]*q[2] - q[3]*q[0]),
                        1 - 2*(q[0]*q[0] + q[1]*q[1]));
  return V;
}

static inline void svd_sort_singular(ManifoldMat3 *B, ManifoldMat3 *V) {
  double rho1 = svd_dist2(B->cols[0]);
  double rho2 = svd_dist2(B->cols[1]);
  double rho3 = svd_dist2(B->cols[2]);
  bool c;
  c = rho1 < rho2;
  svd_cond_neg_swap(c, &B->cols[0].x, &B->cols[1].x);
  svd_cond_neg_swap(c, &V->cols[0].x, &V->cols[1].x);
  svd_cond_neg_swap(c, &B->cols[0].y, &B->cols[1].y);
  svd_cond_neg_swap(c, &V->cols[0].y, &V->cols[1].y);
  svd_cond_neg_swap(c, &B->cols[0].z, &B->cols[1].z);
  svd_cond_neg_swap(c, &V->cols[0].z, &V->cols[1].z);
  svd_cond_swap(c, &rho1, &rho2);
  c = rho1 < rho3;
  svd_cond_neg_swap(c, &B->cols[0].x, &B->cols[2].x);
  svd_cond_neg_swap(c, &V->cols[0].x, &V->cols[2].x);
  svd_cond_neg_swap(c, &B->cols[0].y, &B->cols[2].y);
  svd_cond_neg_swap(c, &V->cols[0].y, &V->cols[2].y);
  svd_cond_neg_swap(c, &B->cols[0].z, &B->cols[2].z);
  svd_cond_neg_swap(c, &V->cols[0].z, &V->cols[2].z);
  svd_cond_swap(c, &rho1, &rho3);
  c = rho2 < rho3;
  svd_cond_neg_swap(c, &B->cols[1].x, &B->cols[2].x);
  svd_cond_neg_swap(c, &V->cols[1].x, &V->cols[2].x);
  svd_cond_neg_swap(c, &B->cols[1].y, &B->cols[2].y);
  svd_cond_neg_swap(c, &V->cols[1].y, &V->cols[2].y);
  svd_cond_neg_swap(c, &B->cols[1].z, &B->cols[2].z);
  svd_cond_neg_swap(c, &V->cols[1].z, &V->cols[2].z);
}

static inline Givens svd_qr_givens(double a1, double a2) {
  double rho = hypot(a1, a2);
  Givens g = {fabs(a1) + fmax(rho, SVD_EPSILON), rho > SVD_EPSILON ? a2 : 0};
  bool b = a1 < 0.0;
  svd_cond_swap(b, &g.sh, &g.ch);
  double w = 1.0 / hypot(g.ch, g.sh);
  g.ch *= w; g.sh *= w;
  return g;
}

static inline QR svd_qr_decompose(ManifoldMat3 *B) {
  QR qr;
  ManifoldMat3 R;
  Givens g1 = svd_qr_givens(B->cols[0].x, B->cols[0].y);
  double a = -2.0 * g1.sh * g1.sh + 1.0;
  double b = 2.0 * g1.ch * g1.sh;
  R.cols[0].x = a*B->cols[0].x + b*B->cols[0].y; R.cols[1].x = a*B->cols[1].x + b*B->cols[1].y; R.cols[2].x = a*B->cols[2].x + b*B->cols[2].y;
  R.cols[0].y = -b*B->cols[0].x + a*B->cols[0].y; R.cols[1].y = -b*B->cols[1].x + a*B->cols[1].y; R.cols[2].y = -b*B->cols[2].x + a*B->cols[2].y;
  R.cols[0].z = B->cols[0].z; R.cols[1].z = B->cols[1].z; R.cols[2].z = B->cols[2].z;
  Givens g2 = svd_qr_givens(R.cols[0].x, R.cols[0].z);
  a = -2.0 * g2.sh * g2.sh + 1.0; b = 2.0 * g2.ch * g2.sh;
  B->cols[0].x = a*R.cols[0].x + b*R.cols[0].z; B->cols[1].x = a*R.cols[1].x + b*R.cols[1].z; B->cols[2].x = a*R.cols[2].x + b*R.cols[2].z;
  B->cols[0].y = R.cols[0].y; B->cols[1].y = R.cols[1].y; B->cols[2].y = R.cols[2].y;
  B->cols[0].z = -b*R.cols[0].x + a*R.cols[0].z; B->cols[1].z = -b*R.cols[1].x + a*R.cols[1].z; B->cols[2].z = -b*R.cols[2].x + a*R.cols[2].z;
  Givens g3 = svd_qr_givens(B->cols[1].y, B->cols[1].z);
  a = -2.0 * g3.sh * g3.sh + 1.0; b = 2.0 * g3.ch * g3.sh;
  R.cols[0].x = B->cols[0].x; R.cols[1].x = B->cols[1].x; R.cols[2].x = B->cols[2].x;
  R.cols[0].y = a*B->cols[0].y + b*B->cols[0].z; R.cols[1].y = a*B->cols[1].y + b*B->cols[1].z; R.cols[2].y = a*B->cols[2].y + b*B->cols[2].z;
  R.cols[0].z = -b*B->cols[0].y + a*B->cols[0].z; R.cols[1].z = -b*B->cols[1].y + a*B->cols[1].z; R.cols[2].z = -b*B->cols[2].y + a*B->cols[2].z;
  double sh12 = 2.0*(g1.sh*g1.sh - 0.5);
  double sh22 = 2.0*(g2.sh*g2.sh - 0.5);
  double sh32 = 2.0*(g3.sh*g3.sh - 0.5);
  ManifoldMat3 Q;
  Q.cols[0].x = sh12*sh22;
  Q.cols[0].y = 4*g2.ch*g3.ch*sh12*g2.sh*g3.sh + 2*g1.ch*g1.sh*sh32;
  Q.cols[0].z = 4*g1.ch*g3.ch*g1.sh*g3.sh + (-2)*g2.ch*sh12*g2.sh*sh32;
  Q.cols[1].x = -2*g1.ch*g1.sh*sh22;
  Q.cols[1].y = -8*g1.ch*g2.ch*g3.ch*g1.sh*g2.sh*g3.sh + sh12*sh32;
  Q.cols[1].z = -2*g3.ch*g3.sh + 4*g1.sh*(g3.ch*g1.sh*g3.sh + g1.ch*g2.ch*g2.sh*sh32);
  Q.cols[2].x = 2*g2.ch*g2.sh;
  Q.cols[2].y = -2*g3.ch*sh22*g3.sh;
  Q.cols[2].z = sh22*sh32;
  qr.Q = Q; qr.R = R;
  return qr;
}

static inline SVDSet manifold_svd(ManifoldMat3 A) {
  ManifoldMat3 AtA = mat3_mul(mat3_transpose(A), A);
  Symmetric3x3 sym = {AtA.cols[0].x, AtA.cols[0].y, AtA.cols[1].y,
                       AtA.cols[0].z, AtA.cols[1].z, AtA.cols[2].z};
  ManifoldMat3 V = svd_jacobi_eigen(sym);
  ManifoldMat3 B = mat3_mul(A, V);
  svd_sort_singular(&B, &V);
  QR qr = svd_qr_decompose(&B);
  return (SVDSet){qr.Q, qr.R, V};
}

static inline double manifold_spectral_norm(ManifoldMat3 A) {
  SVDSet usv = manifold_svd(A);
  return usv.S.cols[0].x;
}

#endif // MANIFOLD_VEC_MATH_H
