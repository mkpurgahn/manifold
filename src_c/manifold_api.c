// Copyright 2021 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Public API implementation for the C11 Manifold port.

#include "manifold_api.h"
#include "manifold_boolean.h"
#include "manifold_hull.h"
#include "manifold_collider.h"
#include "manifold_disjoint_sets.h"
#include <string.h>
#include <stdio.h>
#include <float.h>

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

// Global execution params
static ManifoldExecutionParams g_params = {false, false, true, false, 0};

ManifoldExecutionParams *manifold_get_params(void) { return &g_params; }

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

// Create from mesh
Manifold manifold_from_mesh(const ManifoldVec3 *vertPos, size_t numVert,
                            const ManifoldIVec3 *triVerts, size_t numTri) {
  Manifold m;
  manifold_impl_init(&m.impl);
  if (numVert == 0 || numTri == 0) return m;

  // Check for NaN/Inf vertices
  for (size_t i = 0; i < numVert; i++) {
    if (!isfinite(vertPos[i].x) || !isfinite(vertPos[i].y) ||
        !isfinite(vertPos[i].z)) {
      m.impl.status = MANIFOLD_ERROR_NON_FINITE_VERTEX;
      return m;
    }
  }

  // Check for out-of-bounds vertex indices
  for (size_t i = 0; i < numTri; i++) {
    if (triVerts[i].x < 0 || triVerts[i].x >= (int)numVert ||
        triVerts[i].y < 0 || triVerts[i].y >= (int)numVert ||
        triVerts[i].z < 0 || triVerts[i].z >= (int)numVert) {
      m.impl.status = MANIFOLD_ERROR_VERTEX_OUT_OF_BOUNDS;
      return m;
    }
  }

  m.impl.vertPos = vec_vec3_create_n(numVert);
  for (size_t i = 0; i < numVert; i++) {
    m.impl.vertPos.data[i] = vertPos[i];
  }

  ManifoldVecIVec3 tris = {0};
  for (size_t i = 0; i < numTri; i++) {
    vec_ivec3_push(&tris, triVerts[i]);
  }

  ManifoldVecIVec3 emptyTriVert = {0};
  manifold_impl_create_halfedges(&m.impl, &tris, &emptyTriVert);
  manifold_impl_initialize_original(&m.impl);
  manifold_impl_calculate_bbox(&m.impl);
  manifold_impl_set_epsilon(&m.impl, -1.0, false);
  // Only call CleanupTopology if mesh is manifold (all edges paired)
  if (manifold_impl_is_manifold(&m.impl)) {
    manifold_impl_cleanup_topology(&m.impl);
    manifold_impl_remove_unreferenced_verts(&m.impl);
  } else {
    m.impl.status = MANIFOLD_ERROR_NOT_MANIFOLD;
  }
  manifold_impl_sort_geometry(&m.impl);
  manifold_impl_set_normals_and_coplanar(&m.impl);

  vec_ivec3_free(&tris);
  vec_ivec3_free(&emptyTriVert);
  return m;
}

ManifoldError manifold_status(const Manifold *m) { return m->impl.status; }
bool manifold_is_empty(const Manifold *m) { return manifold_impl_is_empty(&m->impl); }
bool manifold_is_manifold(const Manifold *m) { return manifold_impl_is_manifold(&m->impl); }
bool manifold_is_2manifold(const Manifold *m) { return manifold_impl_is_2manifold(&m->impl); }
size_t manifold_num_vert(const Manifold *m) { return manifold_impl_num_vert(&m->impl); }
size_t manifold_num_edge(const Manifold *m) { return manifold_impl_num_edge(&m->impl); }
size_t manifold_num_tri(const Manifold *m) { return manifold_impl_num_tri(&m->impl); }
ManifoldBox manifold_bounding_box(const Manifold *m) { return m->impl.bBox; }
double manifold_volume(const Manifold *m) { return manifold_impl_get_volume(&m->impl); }
double manifold_surface_area(const Manifold *m) { return manifold_impl_get_surface_area(&m->impl); }

Manifold manifold_empty(void) {
  Manifold m;
  manifold_impl_init(&m.impl);
  return m;
}

Manifold manifold_invalid(void) {
  Manifold m;
  manifold_impl_init(&m.impl);
  m.impl.status = MANIFOLD_ERROR_INVALID_CONSTRUCTION;
  return m;
}

Manifold manifold_tetrahedron(void) {
  Manifold m;
  manifold_impl_tetrahedron(&m.impl);
  return m;
}

Manifold manifold_cube(ManifoldVec3 size, bool center) {
  Manifold m;
  if (size.x <= 0 || size.y <= 0 || size.z <= 0) {
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
  manifold_impl_calculate_bbox(&out.impl);
  manifold_impl_set_epsilon(&out.impl, out.impl.epsilon, false);

  // Update meshRelation transforms: compose translation into each relation
  for (size_t i = 0; i < out.impl.meshRelation.meshIDtransform.len; i++) {
    ManifoldMat3x4 *t = &out.impl.meshRelation.meshIDtransform.data[i].value.transform;
    t->cols[3] = vec3_add(t->cols[3], v);
  }

  // Reuse collider from source with translated boxes (matches C++ LazyCollider
  // transform behavior for axis-aligned translations)
  if (m->impl.colliderBuilt && m->impl.collider.internalChildren.len > 0) {
    out.impl.collider = manifold_collider_copy(&m->impl.collider);
    for (size_t i = 0; i < out.impl.collider.nodeBBox.len; i++) {
      out.impl.collider.nodeBBox.data[i].min =
          vec3_add(out.impl.collider.nodeBBox.data[i].min, v);
      out.impl.collider.nodeBBox.data[i].max =
          vec3_add(out.impl.collider.nodeBBox.data[i].max, v);
    }
    out.impl.colliderBuilt = true;
  }
  return out;
}

// Helper: FlipHalfedge remaps a halfedge index within its triangle (0↔2)
static inline int flip_halfedge(int e) {
  return 3 * (e / 3) + 2 - (e % 3);
}

// Transform tangent vectors by mat3. If invert, remap through flipped halfedges.
static void transform_tangents(ManifoldImpl *impl, ManifoldMat3 m3, bool invert,
                               const ManifoldVecVec4 *oldTangents,
                               const ManifoldVecHalfedge *oldHalfedge) {
  if (oldTangents->len == 0) return;
  vec_vec4_free(&impl->halfedgeTangent);
  vec_vec4_resize(&impl->halfedgeTangent, oldTangents->len);
  for (size_t edgeOut = 0; edgeOut < oldTangents->len; edgeOut++) {
    int edgeIn;
    if (invert) {
      edgeIn = oldHalfedge->data[flip_halfedge((int)edgeOut)].pairedHalfedge;
    } else {
      edgeIn = (int)edgeOut;
    }
    ManifoldVec3 t = vec4_to_vec3(oldTangents->data[edgeIn]);
    ManifoldVec3 transformed = mat3_mul_vec3(m3, t);
    impl->halfedgeTangent.data[edgeOut] =
        vec3_to_vec4(transformed, oldTangents->data[edgeIn].w);
  }
}

Manifold manifold_scale(const Manifold *m, ManifoldVec3 v) {
  Manifold out;
  manifold_copy(&out, m);
  for (size_t i = 0; i < out.impl.vertPos.len; i++) {
    out.impl.vertPos.data[i] = vec3_mul(out.impl.vertPos.data[i], v);
  }
  bool invert = v.x * v.y * v.z < 0;
  ManifoldMat3 scaleMat;
  scaleMat.cols[0] = manifold_vec3(v.x, 0, 0);
  scaleMat.cols[1] = manifold_vec3(0, v.y, 0);
  scaleMat.cols[2] = manifold_vec3(0, 0, v.z);

  // Update meshRelation transforms: compose scale into each relation
  ManifoldMat3x4 scaleT;
  scaleT.cols[0] = manifold_vec3(v.x, 0, 0);
  scaleT.cols[1] = manifold_vec3(0, v.y, 0);
  scaleT.cols[2] = manifold_vec3(0, 0, v.z);
  scaleT.cols[3] = manifold_vec3(0, 0, 0);
  for (size_t i = 0; i < out.impl.meshRelation.meshIDtransform.len; i++) {
    ManifoldRelation *rel = &out.impl.meshRelation.meshIDtransform.data[i].value;
    rel->transform = mat3x4_compose(scaleT, rel->transform);
  }

  // Transform tangents BEFORE flipping (uses old halfedge for mapping)
  if (m->impl.halfedgeTangent.len > 0) {
    transform_tangents(&out.impl, scaleMat, invert,
                       &m->impl.halfedgeTangent, &m->impl.halfedge);
  }

  // Flip triangle winding if negative determinant
  if (invert) {
    size_t numTri = out.impl.halfedge.len / 3;
    for (size_t tri = 0; tri < numTri; tri++) {
      ManifoldHalfedge tmp = out.impl.halfedge.data[3 * tri];
      out.impl.halfedge.data[3 * tri] = out.impl.halfedge.data[3 * tri + 2];
      out.impl.halfedge.data[3 * tri + 2] = tmp;
      for (int i = 0; i < 3; i++) {
        ManifoldHalfedge *he = &out.impl.halfedge.data[3 * tri + i];
        int sv = he->startVert;
        he->startVert = he->endVert;
        he->endVert = sv;
        if (he->pairedHalfedge >= 0) {
          he->pairedHalfedge = flip_halfedge(he->pairedHalfedge);
        }
      }
    }
  }
  // Transform normals with sign flips
  ManifoldVec3 signV = manifold_vec3(v.x < 0 ? -1 : 1, v.y < 0 ? -1 : 1,
                                     v.z < 0 ? -1 : 1);
  for (size_t i = 0; i < out.impl.faceNormal.len; i++) {
    out.impl.faceNormal.data[i] = vec3_normalize(
        vec3_mul(out.impl.faceNormal.data[i], signV));
  }
  for (size_t i = 0; i < out.impl.vertNormal.len; i++) {
    out.impl.vertNormal.data[i] = vec3_normalize(
        vec3_mul(out.impl.vertNormal.data[i], signV));
  }
  manifold_impl_calculate_bbox(&out.impl);
  out.impl.epsilon *= manifold_spectral_norm(scaleMat);
  manifold_impl_set_epsilon(&out.impl, out.impl.epsilon, false);
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

  ManifoldMat3x4 rotT;
  rotT.cols[0] = rot.cols[0]; rotT.cols[1] = rot.cols[1];
  rotT.cols[2] = rot.cols[2]; rotT.cols[3] = manifold_vec3(0, 0, 0);

  Manifold out;
  manifold_copy(&out, m);

  // Update meshRelation transforms
  for (size_t i = 0; i < out.impl.meshRelation.meshIDtransform.len; i++) {
    ManifoldRelation *rel = &out.impl.meshRelation.meshIDtransform.data[i].value;
    rel->transform = mat3x4_compose(rotT, rel->transform);
  }

  // Rotation has det=1, so no invert; tangent transform is just rotation
  if (m->impl.halfedgeTangent.len > 0) {
    transform_tangents(&out.impl, rot, false,
                       &m->impl.halfedgeTangent, &m->impl.halfedge);
  }
  for (size_t i = 0; i < out.impl.vertPos.len; i++) {
    out.impl.vertPos.data[i] = mat3_mul_vec3(rot, out.impl.vertPos.data[i]);
  }
  // Transform normals with rotation matrix
  ManifoldMat3 normalXf = mat3_inverse(mat3_transpose(rot));
  for (size_t i = 0; i < out.impl.faceNormal.len; i++) {
    out.impl.faceNormal.data[i] = vec3_normalize(
        mat3_mul_vec3(normalXf, out.impl.faceNormal.data[i]));
  }
  for (size_t i = 0; i < out.impl.vertNormal.len; i++) {
    out.impl.vertNormal.data[i] = vec3_normalize(
        mat3_mul_vec3(normalXf, out.impl.vertNormal.data[i]));
  }
  manifold_impl_calculate_bbox(&out.impl);
  manifold_impl_set_epsilon(&out.impl, out.impl.epsilon, false);
  return out;
}

Manifold manifold_transform(const Manifold *m, ManifoldMat3x4 t) {
  Manifold out;
  manifold_copy(&out, m);
  for (size_t i = 0; i < out.impl.vertPos.len; i++) {
    out.impl.vertPos.data[i] = mat3x4_transform_point(t, out.impl.vertPos.data[i]);
  }

  // Update meshRelation transforms
  for (size_t i = 0; i < out.impl.meshRelation.meshIDtransform.len; i++) {
    ManifoldRelation *rel = &out.impl.meshRelation.meshIDtransform.data[i].value;
    rel->transform = mat3x4_compose(t, rel->transform);
  }

  // Extract 3x3 rotation/scale part
  ManifoldMat3 m3;
  m3.cols[0] = manifold_vec3(t.cols[0].x, t.cols[0].y, t.cols[0].z);
  m3.cols[1] = manifold_vec3(t.cols[1].x, t.cols[1].y, t.cols[1].z);
  m3.cols[2] = manifold_vec3(t.cols[2].x, t.cols[2].y, t.cols[2].z);

  // Normal transform = inverse(transpose(m3))
  ManifoldMat3 normalXf = mat3_inverse(mat3_transpose(m3));

  bool invert = mat3_det(m3) < 0;
  // Transform tangents BEFORE flipping (uses old halfedge for mapping)
  if (m->impl.halfedgeTangent.len > 0) {
    transform_tangents(&out.impl, m3, invert,
                       &m->impl.halfedgeTangent, &m->impl.halfedge);
  }

  for (size_t i = 0; i < out.impl.faceNormal.len; i++) {
    out.impl.faceNormal.data[i] = vec3_normalize(
        mat3_mul_vec3(normalXf, out.impl.faceNormal.data[i]));
  }
  for (size_t i = 0; i < out.impl.vertNormal.len; i++) {
    out.impl.vertNormal.data[i] = vec3_normalize(
        mat3_mul_vec3(normalXf, out.impl.vertNormal.data[i]));
  }

  // Flip triangle winding if negative determinant
  if (invert) {
    size_t numTri = out.impl.halfedge.len / 3;
    for (size_t tri = 0; tri < numTri; tri++) {
      ManifoldHalfedge tmp = out.impl.halfedge.data[3 * tri];
      out.impl.halfedge.data[3 * tri] = out.impl.halfedge.data[3 * tri + 2];
      out.impl.halfedge.data[3 * tri + 2] = tmp;
      for (int i = 0; i < 3; i++) {
        ManifoldHalfedge *he = &out.impl.halfedge.data[3 * tri + i];
        int sv = he->startVert;
        he->startVert = he->endVert;
        he->endVert = sv;
        if (he->pairedHalfedge >= 0) {
          he->pairedHalfedge = flip_halfedge(he->pairedHalfedge);
        }
      }
    }
  }

  manifold_impl_calculate_bbox(&out.impl);
  // Scale epsilon by spectral norm of the 3x3 part
  double sn = manifold_spectral_norm(m3);
  out.impl.epsilon *= sn;
  manifold_impl_set_epsilon(&out.impl, out.impl.epsilon, false);
  return out;
}

// Boolean operations
Manifold manifold_boolean(const Manifold *a, const Manifold *b,
                          ManifoldOpType op) {
  Manifold result;
  manifold_impl_init(&result.impl);
  manifold_boolean_op(&result.impl, &a->impl, &b->impl, op);
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

// Sphere construction via icosphere subdivision
// Sphere construction helper (unused for now, reserved for subdivision)
#if 0
static void sphere_warp(ManifoldImpl *impl, double radius) {
  for (size_t i = 0; i < impl->vertPos.len; i++) {
    ManifoldVec3 v = impl->vertPos.data[i];
    // Apply cosine mapping then normalize (matching C++ code)
    v.x = cos(MANIFOLD_HALF_PI * (1.0 - v.x));
    v.y = cos(MANIFOLD_HALF_PI * (1.0 - v.y));
    v.z = cos(MANIFOLD_HALF_PI * (1.0 - v.z));
    v = vec3_normalize(v);
    if (isnan(v.x)) v = manifold_vec3(0, 0, 0);
    impl->vertPos.data[i] = vec3_scale(v, radius);
  }
}
#endif

// Edge divisions function for sphere: constant n-1
static int sphere_edge_divisions(ManifoldVec3 edgeVec, ManifoldVec4 t0,
                                 ManifoldVec4 t1, void *ctx) {
  (void)edgeVec; (void)t0; (void)t1;
  return *(int*)ctx;
}

Manifold manifold_sphere(double radius, int circularSegments) {
  Manifold m;
  if (radius <= 0.0) {
    manifold_impl_init(&m.impl);
    manifold_impl_make_empty(&m.impl, MANIFOLD_ERROR_INVALID_CONSTRUCTION);
    return m;
  }

  int n = circularSegments > 0 ? (circularSegments + 3) / 4
                               : manifold_get_circular_segments(radius) / 4;
  if (n < 1) n = 1;

  // Start with octahedron, matching C++: Impl::Shape::Octahedron
  manifold_impl_shape(&m.impl, 2, mat3x4_identity());  // 2 = Octahedron

  // Use proper Subdivide to split each edge into n segments (matching C++)
  if (n > 1) {
    int nMinusOne = n - 1;
    ManifoldVecBarycentric bary = manifold_impl_subdivide(
        &m.impl, sphere_edge_divisions, &nMinusOne, true);
    vec_bary_free(&bary);
  }

  // Apply cosine mapping (matching C++ code) then normalize
  for (size_t i = 0; i < m.impl.vertPos.len; i++) {
    ManifoldVec3 v = m.impl.vertPos.data[i];
    v.x = cos(MANIFOLD_HALF_PI * (1.0 - v.x));
    v.y = cos(MANIFOLD_HALF_PI * (1.0 - v.y));
    v.z = cos(MANIFOLD_HALF_PI * (1.0 - v.z));
    v = vec3_normalize(v);
    if (isnan(v.x)) v = manifold_vec3(0, 0, 0);
    m.impl.vertPos.data[i] = vec3_scale(v, radius);
  }

  // Re-initialize after vertex modification
  manifold_impl_initialize_original(&m.impl);
  manifold_impl_calculate_bbox(&m.impl);
  manifold_impl_set_epsilon(&m.impl, -1.0, false);
  manifold_impl_sort_geometry(&m.impl);
  manifold_impl_set_normals_and_coplanar(&m.impl);
  return m;
}

// Cylinder via extrusion of a circle polygon
Manifold manifold_cylinder(double height, double radiusLow, double radiusHigh,
                           int circularSegments, bool center) {
  Manifold m;
  if (height <= 0.0 || radiusLow <= 0.0) {
    manifold_impl_init(&m.impl);
    manifold_impl_make_empty(&m.impl, MANIFOLD_ERROR_INVALID_CONSTRUCTION);
    return m;
  }

  double scale = radiusHigh >= 0.0 ? radiusHigh / radiusLow : 1.0;
  double radius = fmax(radiusLow, radiusHigh >= 0.0 ? radiusHigh : radiusLow);
  int n = circularSegments > 2 ? circularSegments :
          manifold_get_circular_segments(radius);
  if (n < 3) n = 3;

  // Build circle polygon (matching C++ Cylinder which uses Extrude)
  ManifoldVec2 *circle = (ManifoldVec2 *)malloc(n * sizeof(ManifoldVec2));
  double dPhi = 360.0 / n;
  for (int i = 0; i < n; i++) {
    circle[i].x = radiusLow * manifold_cosd(dPhi * i);
    circle[i].y = radiusLow * manifold_sind(dPhi * i);
  }

  ManifoldVec2 scaleTop = {scale, scale};
  m = manifold_extrude(circle, &n, 1, height, 0, 0.0, scaleTop);
  free(circle);

  if (center) {
    Manifold translated = manifold_translate(&m, manifold_vec3(0, 0, -height / 2.0));
    manifold_destroy(&m);
    Manifold orig = manifold_as_original(&translated);
    manifold_destroy(&translated);
    return orig;
  }
  return m;
}

// Extrude
Manifold manifold_extrude(const ManifoldVec2 *polyVerts,
                          const int *polySizes, int nPolys,
                          double height, int nDivisions,
                          double twistDegrees, ManifoldVec2 scaleTop) {
  Manifold m;
  manifold_impl_extrude(&m.impl, polyVerts, polySizes, nPolys,
                        height, nDivisions, twistDegrees, scaleTop);
  return m;
}

// Revolve
Manifold manifold_revolve(const ManifoldVec2 *polyVerts,
                          const int *polySizes, int nPolys,
                          int circularSegments, double revolveDegrees) {
  Manifold m;
  manifold_impl_revolve(&m.impl, polyVerts, polySizes, nPolys,
                        circularSegments, revolveDegrees);
  return m;
}

// Warp
Manifold manifold_warp(const Manifold *m,
                       void (*warpFn)(double *x, double *y, double *z, void *ctx),
                       void *ctx) {
  Manifold out;
  manifold_copy(&out, m);
  for (size_t i = 0; i < out.impl.vertPos.len; i++) {
    warpFn(&out.impl.vertPos.data[i].x,
           &out.impl.vertPos.data[i].y,
           &out.impl.vertPos.data[i].z, ctx);
  }
  manifold_impl_calculate_bbox(&out.impl);
  manifold_impl_set_epsilon(&out.impl, -1.0, false);
  manifold_impl_sort_geometry(&out.impl);
  manifold_impl_set_normals_and_coplanar(&out.impl);
  out.impl.meshRelation.originalID = -1;
  return out;
}

Manifold manifold_warp_batch(const Manifold *m,
                            void (*warpFn)(ManifoldVec3 *verts, size_t numVerts, void *ctx),
                            void *ctx) {
  Manifold out;
  manifold_copy(&out, m);
  warpFn(out.impl.vertPos.data, out.impl.vertPos.len, ctx);
  manifold_impl_calculate_bbox(&out.impl);
  manifold_impl_set_epsilon(&out.impl, -1.0, false);
  manifold_impl_sort_geometry(&out.impl);
  manifold_impl_set_normals_and_coplanar(&out.impl);
  out.impl.meshRelation.originalID = -1;
  return out;
}

Manifold manifold_hull(const Manifold *m) {
  Manifold result;
  manifold_convex_hull(&result.impl, m->impl.vertPos.data, m->impl.vertPos.len);
  return result;
}

Manifold manifold_hull_points(const ManifoldVec3 *points, size_t numPoints) {
  Manifold result;
  manifold_convex_hull(&result.impl, points, numPoints);
  return result;
}

Manifold manifold_level_set(double (*sdf)(double x, double y, double z, void *ctx),
                            void *ctx, ManifoldBox bounds, double edgeLength,
                            double level, double tolerance) {
  Manifold m;
  manifold_impl_level_set(&m.impl, sdf, ctx, bounds, edgeLength, level, tolerance);
  return m;
}

// Mesh data access
const ManifoldVec3 *manifold_get_vert_positions(const Manifold *m, size_t *count) {
  if (count) *count = m->impl.vertPos.len;
  return m->impl.vertPos.data;
}

void manifold_get_triangles(const Manifold *m, ManifoldIVec3 *out, size_t *count) {
  size_t numTri = manifold_impl_num_tri(&m->impl);
  if (count) *count = numTri;
  if (!out) return;
  for (size_t i = 0; i < numTri; i++) {
    int base = (int)(i * 3);
    out[i].x = m->impl.halfedge.data[base].startVert;
    out[i].y = m->impl.halfedge.data[base + 1].startVert;
    out[i].z = m->impl.halfedge.data[base + 2].startVert;
  }
}

void manifold_get_mesh(const Manifold *m,
                       float **vertProps, size_t *numVert, size_t *numProp,
                       int **triVerts, size_t *numTri) {
  size_t nv = m->impl.vertPos.len;
  size_t nt = manifold_impl_num_tri(&m->impl);
  size_t np = 3 + m->impl.numProp; // xyz + custom properties

  if (numVert) *numVert = nv;
  if (numProp) *numProp = np;
  if (numTri) *numTri = nt;

  if (vertProps) {
    *vertProps = (float *)malloc(nv * np * sizeof(float));
    for (size_t i = 0; i < nv; i++) {
      (*vertProps)[i * np + 0] = (float)m->impl.vertPos.data[i].x;
      (*vertProps)[i * np + 1] = (float)m->impl.vertPos.data[i].y;
      (*vertProps)[i * np + 2] = (float)m->impl.vertPos.data[i].z;
      for (size_t p = 0; p < m->impl.numProp; p++) {
        size_t propIdx = i * m->impl.numProp + p;
        if (propIdx < m->impl.properties.len) {
          (*vertProps)[i * np + 3 + p] = (float)m->impl.properties.data[propIdx];
        } else {
          (*vertProps)[i * np + 3 + p] = 0.0f;
        }
      }
    }
  }

  if (triVerts) {
    *triVerts = (int *)malloc(nt * 3 * sizeof(int));
    for (size_t i = 0; i < nt; i++) {
      int base = (int)(i * 3);
      (*triVerts)[i * 3 + 0] = m->impl.halfedge.data[base].startVert;
      (*triVerts)[i * 3 + 1] = m->impl.halfedge.data[base + 1].startVert;
      (*triVerts)[i * 3 + 2] = m->impl.halfedge.data[base + 2].startVert;
    }
  }
}

void manifold_free_mesh(float *vertProps, int *triVerts) {
  free(vertProps);
  free(triVerts);
}

// ---------- Additional API functions ----------

int manifold_genus(const Manifold *m) {
  int chi = (int)manifold_num_vert(m) - (int)manifold_num_edge(m) +
            (int)manifold_num_tri(m);
  return 1 - chi / 2;
}

int manifold_original_id(const Manifold *m) {
  return m->impl.meshRelation.originalID;
}

double manifold_get_epsilon(const Manifold *m) { return m->impl.epsilon; }
double manifold_get_tolerance(const Manifold *m) { return m->impl.tolerance; }
size_t manifold_num_prop(const Manifold *m) {
  return manifold_impl_num_prop(&m->impl);
}
size_t manifold_num_prop_vert(const Manifold *m) {
  return manifold_impl_num_prop_vert(&m->impl);
}

bool manifold_matches_tri_normals(const Manifold *m) {
  const ManifoldImpl *impl = &m->impl;
  size_t numTri = impl->halfedge.len / 3;
  if (impl->halfedge.len == 0 || impl->faceNormal.len != numTri) return true;
  for (size_t face = 0; face < numTri; face++) {
    if (impl->halfedge.data[3 * face].pairedHalfedge < 0) continue;
    ManifoldMat2x3 projection = manifold_get_axis_aligned_projection(
        impl->faceNormal.data[face]);
    ManifoldVec2 v[3];
    double maxD = -1e308, minD = 1e308;
    for (int i = 0; i < 3; i++) {
      ManifoldVec3 p = impl->vertPos.data[impl->halfedge.data[3*face+i].startVert];
      v[i] = mat2x3_mul_vec3(projection, p);
      double d = vec3_dot(p, impl->faceNormal.data[face]);
      if (!isfinite(d)) continue;
      if (d > maxD) maxD = d;
      if (d < minD) minD = d;
    }
    if (maxD - minD > 2 * impl->tolerance) return false;
    if (manifold_ccw(v[0], v[1], v[2], impl->epsilon * 2) < 0) return false;
  }
  return true;
}

int manifold_num_degenerate_tris(const Manifold *m) {
  const ManifoldImpl *impl = &m->impl;
  size_t numTri = impl->halfedge.len / 3;
  if (impl->halfedge.len == 0 || impl->faceNormal.len != numTri) return 0;
  int count = 0;
  for (size_t face = 0; face < numTri; face++) {
    if (impl->halfedge.data[3 * face].pairedHalfedge < 0) { count++; continue; }
    ManifoldMat2x3 projection = manifold_get_axis_aligned_projection(
        impl->faceNormal.data[face]);
    ManifoldVec2 v[3];
    for (int i = 0; i < 3; i++)
      v[i] = mat2x3_mul_vec3(projection,
          impl->vertPos.data[impl->halfedge.data[3*face+i].startVert]);
    if (manifold_ccw(v[0], v[1], v[2], impl->epsilon / 2) == 0) count++;
  }
  return count;
}

static int refine_tol_cb(ManifoldVec3 edge, ManifoldVec4 t0, ManifoldVec4 t1, void *ctx) {
  double tolerance = *((double *)ctx);
  ManifoldVec3 edgeNorm = manifold_safe_normalize(edge);
  ManifoldVec3 tStart = vec4_to_vec3(t0);
  ManifoldVec3 tEnd = vec4_to_vec3(t1);
  ManifoldVec3 start = vec3_sub(tStart, vec3_scale(edgeNorm, vec3_dot(edgeNorm, tStart)));
  ManifoldVec3 end = vec3_sub(tEnd, vec3_scale(edgeNorm, vec3_dot(edgeNorm, tEnd)));
  double d = 0.5 * (vec3_length(start) + vec3_length(end)) +
             vec3_length(vec3_sub(start, end));
  return (int)sqrt(3.0 * d / (4.0 * tolerance));
}

Manifold manifold_refine_to_tolerance(const Manifold *m, double tolerance) {
  tolerance = fabs(tolerance);
  if (tolerance <= 0 || manifold_is_empty(m)) {
    Manifold out;
    manifold_copy(&out, m);
    return out;
  }

  Manifold out;
  manifold_copy(&out, m);
  if (out.impl.halfedgeTangent.len > 0) {
    manifold_impl_refine(&out.impl, refine_tol_cb, &tolerance, true);
  }
  return out;
}

Manifold manifold_mirror(const Manifold *m, ManifoldVec3 normal) {
  double len = vec3_length(normal);
  if (len == 0.0) {
    Manifold empty;
    manifold_create(&empty);
    return empty;
  }
  ManifoldVec3 n = vec3_scale(normal, 1.0 / len);
  // Mirror matrix: I - 2*n*n^T
  ManifoldMat3 mirror = mat3_sub(mat3_identity(),
                                  mat3_scale_s(vec3_outerprod(n, n), 2.0));

  Manifold out;
  manifold_copy(&out, m);
  // Apply mirror transform to all vertices
  for (size_t i = 0; i < out.impl.vertPos.len; i++) {
    out.impl.vertPos.data[i] = mat3_mul_vec3(mirror, out.impl.vertPos.data[i]);
  }
  // Mirror flips winding - reverse all triangle windings
  size_t numTri = manifold_impl_num_tri(&out.impl);
  for (size_t tri = 0; tri < numTri; tri++) {
    // Swap halfedges 1 and 2 of each triangle to flip winding
    ManifoldHalfedge tmp = out.impl.halfedge.data[3 * tri + 1];
    out.impl.halfedge.data[3 * tri + 1] = out.impl.halfedge.data[3 * tri + 2];
    out.impl.halfedge.data[3 * tri + 2] = tmp;
  }
  // Rebuild halfedges from scratch (re-pair)
  ManifoldVecIVec3 triVerts = {0};
  for (size_t tri = 0; tri < numTri; tri++) {
    ManifoldIVec3 tv;
    tv.x = out.impl.halfedge.data[3 * tri].startVert;
    tv.y = out.impl.halfedge.data[3 * tri + 1].startVert;
    tv.z = out.impl.halfedge.data[3 * tri + 2].startVert;
    vec_ivec3_push(&triVerts, tv);
  }
  // Clear and rebuild
  vec_halfedge_clear(&out.impl.halfedge);
  vec_vec3_clear(&out.impl.faceNormal);
  vec_vec3_clear(&out.impl.vertNormal);
  ManifoldVecIVec3 emptyTriProp = {0};
  manifold_impl_create_halfedges(&out.impl, &triVerts, &emptyTriProp);
  manifold_impl_initialize_original(&out.impl);
  manifold_impl_calculate_bbox(&out.impl);
  manifold_impl_set_epsilon(&out.impl, -1.0, false);
  manifold_impl_sort_geometry(&out.impl);
  manifold_impl_set_normals_and_coplanar(&out.impl);
  vec_ivec3_free(&triVerts);
  vec_ivec3_free(&emptyTriProp);
  return out;
}

// Helper: create a halfspace cutter
static Manifold manifold_halfspace(ManifoldBox bBox, ManifoldVec3 normal,
                                    double originOffset) {
  normal = vec3_normalize(normal);
  Manifold cutter = manifold_cube(manifold_vec3(2.0, 2.0, 2.0), true);
  cutter = manifold_translate(&cutter, manifold_vec3(1.0, 0.0, 0.0));

  double size = vec3_length(vec3_sub(manifold_box_center(bBox),
                                      vec3_scale(normal, originOffset))) +
                0.5 * vec3_length(manifold_box_size(bBox));
  ManifoldVec3 sv = manifold_vec3(size, size, size);
  Manifold scaled = manifold_scale(&cutter, sv);
  manifold_destroy(&cutter);

  Manifold translated = manifold_translate(&scaled,
                                            manifold_vec3(originOffset, 0, 0));
  manifold_destroy(&scaled);

  double yDeg = manifold_degrees(-asin(normal.z));
  double zDeg = manifold_degrees(atan2(normal.y, normal.x));
  Manifold rotated = manifold_rotate(&translated, 0.0, yDeg, zDeg);
  manifold_destroy(&translated);

  return rotated;
}

void manifold_split(const Manifold *m, const Manifold *cutter,
                    Manifold *first, Manifold *second) {
  *first = manifold_intersection(m, cutter);
  *second = manifold_difference(m, cutter);
}

void manifold_split_by_plane(const Manifold *m, ManifoldVec3 normal,
                              double originOffset, Manifold *first,
                              Manifold *second) {
  Manifold hs = manifold_halfspace(manifold_bounding_box(m), normal,
                                    originOffset);
  manifold_split(m, &hs, first, second);
  manifold_destroy(&hs);
}

Manifold manifold_trim_by_plane(const Manifold *m, ManifoldVec3 normal,
                                 double originOffset) {
  Manifold hs = manifold_halfspace(manifold_bounding_box(m), normal,
                                    originOffset);
  Manifold result = manifold_intersection(m, &hs);
  manifold_destroy(&hs);
  return result;
}

int manifold_decompose(const Manifold *m, Manifold **components,
                        int maxComponents) {
  if (manifold_is_empty(m) || maxComponents <= 0) return 0;

  // Allocate impl array
  ManifoldImpl *implArr = (ManifoldImpl *)malloc(
      (size_t)maxComponents * sizeof(ManifoldImpl));
  int n = manifold_impl_decompose(&m->impl, implArr, maxComponents);

  // Wrap in Manifold objects
  *components = (Manifold *)malloc((size_t)n * sizeof(Manifold));
  for (int i = 0; i < n; i++) {
    (*components)[i].impl = implArr[i];
  }
  free(implArr);
  return n;
}

Manifold manifold_batch_boolean(const Manifold *manifolds, int count,
                                 ManifoldOpType op) {
  if (count <= 0) {
    Manifold empty;
    manifold_create(&empty);
    return empty;
  }
  if (count == 1) {
    Manifold result;
    manifold_copy(&result, &manifolds[0]);
    return result;
  }
  if (count == 2) {
    return manifold_boolean(&manifolds[0], &manifolds[1], op);
  }
  // For Subtract: result = first - union(rest), matching C++ CSG tree semantics
  if (op == MANIFOLD_OP_SUBTRACT) {
    Manifold negUnion = manifold_batch_boolean(manifolds + 1, count - 1,
                                                MANIFOLD_OP_ADD);
    Manifold result = manifold_boolean(&manifolds[0], &negUnion,
                                        MANIFOLD_OP_SUBTRACT);
    manifold_destroy(&negUnion);
    return result;
  }
  // Heap-based batch boolean matching C++ behavior:
  // Always pair the two largest meshes (by vertex count) first.
  int n = count;
  Manifold *heap = (Manifold *)malloc((size_t)n * sizeof(Manifold));
  for (int i = 0; i < n; i++) {
    manifold_copy(&heap[i], &manifolds[i]);
  }

  // Max-heap by vertex count (largest at index 0)
  #define HEAP_VAL(i) ((int)heap[i].impl.vertPos.len)
  #define HEAP_SWAP(i, j) do { Manifold tmp_ = heap[i]; heap[i] = heap[j]; heap[j] = tmp_; } while(0)

  // Build max-heap
  for (int i = n / 2 - 1; i >= 0; i--) {
    int k = i;
    while (1) {
      int largest = k;
      int l = 2 * k + 1, r = 2 * k + 2;
      if (l < n && HEAP_VAL(l) > HEAP_VAL(largest)) largest = l;
      if (r < n && HEAP_VAL(r) > HEAP_VAL(largest)) largest = r;
      if (largest == k) break;
      HEAP_SWAP(k, largest);
      k = largest;
    }
  }

  while (n > 1) {
    Manifold *tmp = (Manifold *)malloc((size_t)n * sizeof(Manifold));
    int tmpLen = 0;
    // Pop pairs (up to 4 pairs per round, matching C++)
    for (int i = 0; i < 4 && n > 1; i++) {
      // Pop largest (heap[0])
      HEAP_SWAP(0, n - 1);
      Manifold a = heap[--n];
      // Sift down
      int k = 0;
      while (1) {
        int largest = k;
        int l = 2*k+1, r = 2*k+2;
        if (l < n && HEAP_VAL(l) > HEAP_VAL(largest)) largest = l;
        if (r < n && HEAP_VAL(r) > HEAP_VAL(largest)) largest = r;
        if (largest == k) break;
        HEAP_SWAP(k, largest);
        k = largest;
      }
      // Pop second largest
      HEAP_SWAP(0, n - 1);
      Manifold b = heap[--n];
      k = 0;
      while (1) {
        int largest = k;
        int l = 2*k+1, r = 2*k+2;
        if (l < n && HEAP_VAL(l) > HEAP_VAL(largest)) largest = l;
        if (r < n && HEAP_VAL(r) > HEAP_VAL(largest)) largest = r;
        if (largest == k) break;
        HEAP_SWAP(k, largest);
        k = largest;
      }
      tmp[tmpLen++] = manifold_boolean(&a, &b, op);
      manifold_destroy(&a);
      manifold_destroy(&b);
    }
    // Push results back into heap
    for (int i = 0; i < tmpLen; i++) {
      heap[n] = tmp[i];
      // Sift up
      int k = n++;
      while (k > 0) {
        int parent = (k - 1) / 2;
        if (HEAP_VAL(k) > HEAP_VAL(parent)) {
          HEAP_SWAP(k, parent);
          k = parent;
        } else break;
      }
    }
    free(tmp);
  }
  #undef HEAP_VAL
  #undef HEAP_SWAP

  Manifold result = heap[0];
  free(heap);
  return result;
}

Manifold manifold_set_properties(const Manifold *m, int numProp,
    void (*propFunc)(double *newProp, ManifoldVec3 pos, const double *oldProp, void *ctx),
    void *ctx) {
  Manifold out;
  manifold_copy(&out, m);
  manifold_impl_set_properties(&out.impl, numProp, propFunc, ctx);
  return out;
}

Manifold manifold_calculate_curvature(const Manifold *m, int gaussianIdx,
                                       int meanIdx) {
  Manifold out;
  manifold_copy(&out, m);
  manifold_impl_calculate_curvature(&out.impl, gaussianIdx, meanIdx);
  return out;
}

Manifold manifold_as_original(const Manifold *m) {
  Manifold out;
  manifold_copy(&out, m);
  manifold_impl_initialize_original(&out.impl);
  manifold_impl_set_normals_and_coplanar(&out.impl);
  return out;
}

uint32_t manifold_reserve_ids_api(uint32_t n) {
  return manifold_reserve_ids(n);
}

Manifold manifold_simplify(const Manifold *m, double tolerance) {
  Manifold out;
  manifold_copy(&out, m);
  double oldTolerance = out.impl.tolerance;
  if (tolerance == 0) tolerance = oldTolerance;
  if (tolerance > oldTolerance) {
    out.impl.tolerance = tolerance;
    manifold_impl_set_normals_and_coplanar(&out.impl);
  }
  manifold_impl_simplify_topology(&out.impl, 0);
  manifold_impl_sort_geometry(&out.impl);
  // Check for degenerate mesh (zero volume after simplification)
  if (out.impl.halfedge.len > 0) {
    double vol = fabs(manifold_volume(&out));
    double eps3 = out.impl.epsilon * out.impl.epsilon * out.impl.epsilon;
    if (vol < eps3) {
      manifold_destroy(&out);
      manifold_impl_init(&out.impl);
      return out;
    }
  }
  out.impl.tolerance = oldTolerance;
  return out;
}

Manifold manifold_set_tolerance(const Manifold *m, double tolerance) {
  Manifold out;
  manifold_copy(&out, m);
  if (tolerance > out.impl.tolerance) {
    out.impl.tolerance = tolerance;
    manifold_impl_set_normals_and_coplanar(&out.impl);
    manifold_impl_simplify_topology(&out.impl, 0);
    manifold_impl_sort_geometry(&out.impl);
  } else {
    out.impl.tolerance = fmax(out.impl.epsilon, tolerance);
  }
  return out;
}

static int refine_n_cb(ManifoldVec3 e, ManifoldVec4 t0, ManifoldVec4 t1, void *ctx) {
  (void)e; (void)t0; (void)t1;
  return *((int *)ctx) - 1;
}

Manifold manifold_refine(const Manifold *m, int n) {
  if (n < 2 || manifold_is_empty(m)) {
    Manifold out;
    manifold_copy(&out, m);
    return out;
  }

  Manifold out;
  manifold_copy(&out, m);
  int nMinus1 = n;
  manifold_impl_refine(&out.impl, refine_n_cb, &nMinus1, false);
  return out;
}

static int refine_length_cb(ManifoldVec3 edge, ManifoldVec4 t0, ManifoldVec4 t1, void *ctx) {
  (void)t0; (void)t1;
  double length = *((double *)ctx);
  return (int)(vec3_length(edge) / length);
}

Manifold manifold_refine_to_length(const Manifold *m, double length) {
  if (length <= 0 || manifold_is_empty(m)) {
    Manifold out;
    manifold_copy(&out, m);
    return out;
  }
  length = fabs(length);

  Manifold out;
  manifold_copy(&out, m);
  manifold_impl_refine(&out.impl, refine_length_cb, &length, false);
  return out;
}

Manifold manifold_smooth_by_normals(const Manifold *m, int normalIdx) {
  Manifold out;
  manifold_copy(&out, m);
  if (!manifold_is_empty(&out)) {
    manifold_impl_create_tangents_normals(&out.impl, normalIdx);
  }
  return out;
}

Manifold manifold_smooth_out(const Manifold *m, double minSharpAngle,
                              double minSmoothness) {
  Manifold out;
  manifold_copy(&out, m);
  if (!manifold_is_empty(&out)) {
    if (minSmoothness == 0) {
      int numProp = out.impl.numProp;
      ManifoldVecDouble savedProps = MANIFOLD_VEC_INIT;
      if (out.impl.properties.len > 0) {
        vec_double_resize(&savedProps, out.impl.properties.len);
        memcpy(savedProps.data, out.impl.properties.data,
               out.impl.properties.len * sizeof(double));
      }
      ManifoldVecHalfedge savedHE = MANIFOLD_VEC_INIT;
      vec_halfedge_resize(&savedHE, out.impl.halfedge.len);
      memcpy(savedHE.data, out.impl.halfedge.data,
             out.impl.halfedge.len * sizeof(ManifoldHalfedge));

      manifold_impl_set_normals_smooth(&out.impl, 0, minSharpAngle);
      manifold_impl_create_tangents_normals(&out.impl, 0);

      // Restore original properties and halfedge propVerts
      out.impl.numProp = numProp;
      vec_double_free(&out.impl.properties);
      out.impl.properties = savedProps;
      vec_halfedge_free(&out.impl.halfedge);
      out.impl.halfedge = savedHE;
    } else {
      ManifoldVecSmoothness sharp = manifold_impl_sharpen_edges(&out.impl,
          minSharpAngle, minSmoothness);
      manifold_impl_create_tangents_smooth(&out.impl,
          sharp.data, (int)sharp.len);
      vec_smooth_free(&sharp);
    }
  }
  return out;
}

Manifold manifold_smooth(const Manifold *m,
                          const ManifoldSmoothness *sharpenedEdges,
                          int numSharpened) {
  // For smooth with sharpened edges, we need the original triangle ordering.
  // When called from the API with an already-constructed Manifold, the original
  // ordering is the current halfedge ordering (since the input was already
  // a valid manifold). We set faceID=identity and go.
  Manifold out;
  manifold_copy(&out, m);
  if (manifold_is_empty(&out)) return out;

  size_t numTri = manifold_impl_num_tri(&out.impl);
  // Set faceID to identity mapping (current order IS the original order)
  for (size_t i = 0; i < numTri; i++) {
    out.impl.meshRelation.triRef.data[i].faceID = (int)i;
  }

  if (numSharpened > 0 && sharpenedEdges != NULL) {
    ManifoldVecSmoothness updated = manifold_impl_update_sharpened_edges(
        &out.impl, sharpenedEdges, numSharpened);
    manifold_impl_create_tangents_smooth(&out.impl,
        updated.data, (int)updated.len);
    vec_smooth_free(&updated);
  } else {
    manifold_impl_create_tangents_smooth(&out.impl, NULL, 0);
  }

  // Restore faceID
  numTri = manifold_impl_num_tri(&out.impl);
  for (size_t i = 0; i < numTri; i++) {
    out.impl.meshRelation.triRef.data[i].faceID = -1;
  }
  return out;
}

// Smooth from raw mesh data — matching C++ Manifold::Smooth(MeshGL64, edges)
// The sharpened edge halfedge indices refer to the original triangle ordering.
Manifold manifold_smooth_from_mesh(const ManifoldVec3 *vertPos, size_t numVert,
                                    const ManifoldIVec3 *triVerts, size_t numTri,
                                    const ManifoldSmoothness *sharpenedEdges,
                                    int numSharpened) {
  Manifold m;
  manifold_impl_init(&m.impl);
  if (numVert == 0 || numTri == 0) return m;

  fprintf(stderr, "smooth_from_mesh: numVert=%zu numTri=%zu\n", numVert, numTri);

  m.impl.vertPos = vec_vec3_create_n(numVert);
  for (size_t i = 0; i < numVert; i++)
    m.impl.vertPos.data[i] = vertPos[i];

  ManifoldVecIVec3 tris = {0};
  for (size_t i = 0; i < numTri; i++)
    vec_ivec3_push(&tris, triVerts[i]);

  ManifoldVecIVec3 emptyTriVert = {0};
  manifold_impl_create_halfedges(&m.impl, &tris, &emptyTriVert);
  vec_ivec3_free(&tris);
  vec_ivec3_free(&emptyTriVert);

  manifold_impl_initialize_original(&m.impl);

  // Set faceID BEFORE sorting (matching C++ SmoothImpl)
  // Must be after initialize_original which creates triRef entries
  size_t nTri = manifold_impl_num_tri(&m.impl);
  for (size_t i = 0; i < nTri; i++)
    m.impl.meshRelation.triRef.data[i].faceID = (int)i;

  manifold_impl_calculate_bbox(&m.impl);
  manifold_impl_set_epsilon(&m.impl, -1.0, false);
  if (manifold_impl_is_manifold(&m.impl)) {
    manifold_impl_cleanup_topology(&m.impl);
    manifold_impl_remove_unreferenced_verts(&m.impl);
  }
  manifold_impl_sort_geometry(&m.impl);
  manifold_impl_set_normals_and_coplanar(&m.impl);

  // Now UpdateSharpenedEdges maps original->sorted halfedge indices
  if (numSharpened > 0 && sharpenedEdges != NULL) {
    ManifoldVecSmoothness updated = manifold_impl_update_sharpened_edges(
        &m.impl, sharpenedEdges, numSharpened);
    manifold_impl_create_tangents_smooth(&m.impl,
        updated.data, (int)updated.len);
    vec_smooth_free(&updated);
  } else {
    manifold_impl_create_tangents_smooth(&m.impl, NULL, 0);
  }

  // Restore faceID to -1
  nTri = manifold_impl_num_tri(&m.impl);
  for (size_t i = 0; i < nTri; i++)
    m.impl.meshRelation.triRef.data[i].faceID = -1;

  return m;
}

Manifold manifold_smooth_from_meshgl(const ManifoldMeshGL *mesh,
                                      const ManifoldSmoothness *sharpenedEdges,
                                      int numSharpened) {
  Manifold out;
  manifold_impl_init(&out.impl);
  if (!mesh || mesh->vertLen == 0 || mesh->triLen == 0) return out;

  // Make a temporary copy with faceID set to identity
  ManifoldMeshGL tmp = *mesh;
  size_t numTri = mesh->triLen;
  int *faceID = (int *)malloc(numTri * sizeof(int));
  for (size_t i = 0; i < numTri; i++) faceID[i] = (int)i;
  tmp.faceID = faceID;
  tmp.faceIDLen = numTri;

  out = manifold_from_meshgl(&tmp);
  free(faceID);
  if (manifold_is_empty(&out)) return out;

  // Create tangents
  if (numSharpened > 0 && sharpenedEdges != NULL) {
    ManifoldVecSmoothness updated = manifold_impl_update_sharpened_edges(
        &out.impl, sharpenedEdges, numSharpened);
    manifold_impl_create_tangents_smooth(&out.impl,
        updated.data, (int)updated.len);
    vec_smooth_free(&updated);
  } else {
    manifold_impl_create_tangents_smooth(&out.impl, NULL, 0);
  }

  // Restore faceID: map through original MeshGL's faceID or set to -1
  size_t nTri = manifold_impl_num_tri(&out.impl);
  for (size_t i = 0; i < nTri; i++) {
    int origIdx = out.impl.meshRelation.triRef.data[i].faceID;
    if (mesh->faceIDLen == numTri && mesh->faceID &&
        origIdx >= 0 && (size_t)origIdx < numTri) {
      out.impl.meshRelation.triRef.data[i].faceID = mesh->faceID[origIdx];
    } else {
      out.impl.meshRelation.triRef.data[i].faceID = -1;
    }
  }
  return out;
}

bool manifold_is_convex(const Manifold *m) {
  return manifold_impl_is_convex(&m->impl);
}

static Manifold manifold_minkowski_impl(const Manifold *a, const Manifold *b,
                                         bool inset) {
  const ManifoldImpl *aImpl = &a->impl;
  const ManifoldImpl *bImpl = &b->impl;

  bool aConvex = manifold_impl_is_convex(aImpl);
  bool bConvex = manifold_impl_is_convex(bImpl);

  // If the convex manifold was supplied first, swap them
  if (aConvex && !bConvex) {
    const ManifoldImpl *tmp = aImpl;
    aImpl = bImpl;
    bImpl = tmp;
    bool tmpB = aConvex;
    aConvex = bConvex;
    bConvex = tmpB;
  }

  if (manifold_impl_is_empty(bImpl)) {
    Manifold result;
    manifold_copy(&result, (aImpl == &a->impl) ? a : b);
    return result;
  }
  if (manifold_impl_is_empty(aImpl)) {
    Manifold result;
    manifold_copy(&result, (bImpl == &b->impl) ? b : a);
    return result;
  }

  // Start with a copy of the base
  Manifold base;
  manifold_copy(&base, (aImpl == &a->impl) ? a : b);

  // Convex-Convex Minkowski Sum (not inset): hull of all pairwise sums
  if (!inset && aConvex && bConvex) {
    // Match C++: hull of all pairwise vertex sums
    size_t numA = aImpl->vertPos.len;
    size_t numB = bImpl->vertPos.len;
    size_t total = numA * numB;
    ManifoldVec3 *pts = (ManifoldVec3 *)malloc(total * sizeof(ManifoldVec3));
    size_t idx = 0;
    for (size_t i = 0; i < numA; i++)
      for (size_t j = 0; j < numB; j++)
        pts[idx++] = vec3_add(aImpl->vertPos.data[i], bImpl->vertPos.data[j]);
    Manifold hullResult = manifold_hull_points(pts, total);
    free(pts);
    manifold_destroy(&base);
    Manifold orig = manifold_as_original(&hullResult);
    manifold_destroy(&hullResult);
    return orig;
  }

  // Minkowski Difference (inset) with convex B: intersect translated copies of A
  // A ⊖ B = intersection of (A translated by -b_j) for all vertices b_j of B
  if (inset && bConvex) {
    size_t numBVerts = bImpl->vertPos.len;
    const Manifold *aManifold = (aImpl == &a->impl) ? a : b;
    Manifold *translated = (Manifold *)malloc(numBVerts * sizeof(Manifold));
    for (size_t j = 0; j < numBVerts; j++) {
      ManifoldVec3 bj = bImpl->vertPos.data[j];
      Manifold aCopy;
      manifold_copy(&aCopy, aManifold);
      translated[j] = manifold_translate(&aCopy, manifold_vec3(-bj.x, -bj.y, -bj.z));
      manifold_destroy(&aCopy);
    }
    Manifold result = manifold_batch_boolean(translated, (int)numBVerts,
                                              MANIFOLD_OP_INTERSECT);
    for (size_t j = 0; j < numBVerts; j++) manifold_destroy(&translated[j]);
    free(translated);
    manifold_destroy(&base);
    Manifold orig = manifold_as_original(&result);
    manifold_destroy(&result);
    return orig;
  }

  // Convex-NonConvex or NonConvex-Convex with inset
  if ((inset || !aConvex) && bConvex) {
    size_t numTri = manifold_impl_num_tri(aImpl);
    size_t numBVerts = bImpl->vertPos.len;
    size_t numBatches = (numTri + 999) / 1000;
    Manifold *composed = (Manifold *)malloc((1 + numBatches) * sizeof(Manifold));
    composed[0] = base;
    size_t composedLen = 1;

    for (size_t offset = 0; offset < numTri; offset += 1000) {
      size_t numIter = numTri - offset;
      if (numIter > 1000) numIter = 1000;
      Manifold *hulls = (Manifold *)malloc(numIter * sizeof(Manifold));
      for (size_t iter = 0; iter < numIter; iter++) {
        size_t triIdx = offset + iter;
        size_t numPts = 3 * numBVerts;
        ManifoldVec3 *pts = (ManifoldVec3 *)malloc(numPts * sizeof(ManifoldVec3));
        size_t pidx = 0;
        for (int i = 0; i < 3; i++) {
          ManifoldVec3 vert = aImpl->vertPos.data[
              aImpl->halfedge.data[triIdx * 3 + i].startVert];
          for (size_t j = 0; j < numBVerts; j++) {
            pts[pidx++] = vec3_add(vert, bImpl->vertPos.data[j]);
          }
        }
        hulls[iter] = manifold_hull_points(pts, numPts);
        free(pts);
      }
      composed[composedLen++] = manifold_batch_boolean(hulls, (int)numIter,
                                                        MANIFOLD_OP_ADD);
      for (size_t i = 0; i < numIter; i++) manifold_destroy(&hulls[i]);
      free(hulls);
    }

    ManifoldOpType op = inset ? MANIFOLD_OP_SUBTRACT : MANIFOLD_OP_ADD;
    Manifold result = manifold_batch_boolean(composed, (int)composedLen, op);
    for (size_t i = 0; i < composedLen; i++) manifold_destroy(&composed[i]);
    free(composed);
    Manifold orig = manifold_as_original(&result);
    manifold_destroy(&result);
    return orig;
  }

  // Non-Convex - Non-Convex
  {
    size_t numTriA = manifold_impl_num_tri(aImpl);
    size_t numTriB = manifold_impl_num_tri(bImpl);

    Manifold *accumulated = (Manifold *)malloc(200 * sizeof(Manifold));
    size_t accLen = 0;
    Manifold *allComposed = (Manifold *)malloc((2 + numTriA) * sizeof(Manifold));
    allComposed[0] = base;
    size_t allLen = 1;

    for (size_t aFace = 0; aFace < numTriA; aFace++) {
      ManifoldVec3 a1 = aImpl->vertPos.data[aImpl->halfedge.data[aFace * 3].startVert];
      ManifoldVec3 a2 = aImpl->vertPos.data[aImpl->halfedge.data[aFace * 3 + 1].startVert];
      ManifoldVec3 a3 = aImpl->vertPos.data[aImpl->halfedge.data[aFace * 3 + 2].startVert];
      ManifoldVec3 nA = aImpl->faceNormal.data[aFace];

      Manifold *validHulls = (Manifold *)malloc(numTriB * sizeof(Manifold));
      size_t validCount = 0;

      for (size_t bFace = 0; bFace < numTriB; bFace++) {
        ManifoldVec3 nB = bImpl->faceNormal.data[bFace];
        double dotSame = vec3_dot(nA, nB);
        double dotOpp = vec3_dot(nA, vec3_neg(nB));
        bool coplanar = (fabs(dotSame - 1.0) < 1e-12) ||
                        (fabs(dotOpp - 1.0) < 1e-12);
        if (coplanar) continue;

        ManifoldVec3 b1 = bImpl->vertPos.data[bImpl->halfedge.data[bFace * 3].startVert];
        ManifoldVec3 b2 = bImpl->vertPos.data[bImpl->halfedge.data[bFace * 3 + 1].startVert];
        ManifoldVec3 b3 = bImpl->vertPos.data[bImpl->halfedge.data[bFace * 3 + 2].startVert];

        ManifoldVec3 pts[9] = {
          vec3_add(a1, b1), vec3_add(a1, b2), vec3_add(a1, b3),
          vec3_add(a2, b1), vec3_add(a2, b2), vec3_add(a2, b3),
          vec3_add(a3, b1), vec3_add(a3, b2), vec3_add(a3, b3),
        };
        Manifold h = manifold_hull_points(pts, 9);
        if (!manifold_is_empty(&h)) {
          validHulls[validCount++] = h;
        } else {
          manifold_destroy(&h);
        }
      }

      if (validCount > 0) {
        accumulated[accLen++] = manifold_batch_boolean(validHulls, (int)validCount,
                                                        MANIFOLD_OP_ADD);
      }
      for (size_t i = 0; i < validCount; i++) manifold_destroy(&validHulls[i]);
      free(validHulls);

      if (accLen >= 200) {
        Manifold reduced = manifold_batch_boolean(accumulated, (int)accLen,
                                                   MANIFOLD_OP_ADD);
        for (size_t i = 0; i < accLen; i++) manifold_destroy(&accumulated[i]);
        accLen = 0;
        accumulated[accLen++] = reduced;
      }
    }

    if (accLen > 0) {
      allComposed[allLen++] = manifold_batch_boolean(accumulated, (int)accLen,
                                                      MANIFOLD_OP_ADD);
      for (size_t i = 0; i < accLen; i++) manifold_destroy(&accumulated[i]);
    }
    free(accumulated);

    ManifoldOpType op = inset ? MANIFOLD_OP_SUBTRACT : MANIFOLD_OP_ADD;
    Manifold result = manifold_batch_boolean(allComposed, (int)allLen, op);
    for (size_t i = 0; i < allLen; i++) manifold_destroy(&allComposed[i]);
    free(allComposed);

    Manifold orig = manifold_as_original(&result);
    manifold_destroy(&result);
    return orig;
  }
}

Manifold manifold_minkowski_sum(const Manifold *a, const Manifold *b) {
  return manifold_minkowski_impl(a, b, false);
}

Manifold manifold_minkowski_difference(const Manifold *a, const Manifold *b) {
  return manifold_minkowski_impl(a, b, true);
}

double manifold_min_gap(const Manifold *a, const Manifold *b,
                        double searchLength) {
  return manifold_impl_min_gap(&a->impl, &b->impl, searchLength);
}

Manifold manifold_calculate_normals(const Manifold *m, int normalIdx,
                                     double minSharpAngle) {
  Manifold result;
  manifold_copy(&result, m);
  manifold_impl_calculate_normals(&result.impl, normalIdx, minSharpAngle);
  return result;
}

// ---------- MeshGL Merge ----------
// Port of C++ MeshGL::Merge() from sort.cpp

// Edge pair for open-edge tracking (simulates std::multiset<pair<int,int>>)
typedef struct { int first, second; } MergeEdgePair;

// Collision callback context for merge
typedef struct {
  ManifoldDisjointSets *uf;
  int *openVerts;
} MergeCollisionCtx;

static void merge_collision_callback(int a, int b, void *ctx) {
  MergeCollisionCtx *mc = (MergeCollisionCtx *)ctx;
  manifold_disjoint_sets_unite(mc->uf, (uint32_t)mc->openVerts[a],
                               (uint32_t)mc->openVerts[b]);
}

// Morton sort comparator context
typedef struct {
  const uint32_t *morton;
} MergeMortonCtx;

static MergeMortonCtx g_merge_morton_ctx;

static int merge_morton_compare(const void *a, const void *b) {
  int ia = *(const int *)a;
  int ib = *(const int *)b;
  uint32_t ma = g_merge_morton_ctx.morton[ia];
  uint32_t mb = g_merge_morton_ctx.morton[ib];
  if (ma != mb) return ma < mb ? -1 : 1;
  // stable sort tie-breaking by index
  return ia < ib ? -1 : (ia > ib ? 1 : 0);
}

bool manifold_meshgl_merge(ManifoldMeshGL *mesh) {
  if (!mesh || mesh->vertLen == 0 || mesh->triLen == 0) return false;

  const size_t numVert = mesh->vertLen;
  const size_t numTri = mesh->triLen;

  // Build merge identity map, applying existing merge pairs
  int *merge = (int *)malloc(numVert * sizeof(int));
  for (size_t i = 0; i < numVert; i++) merge[i] = (int)i;
  for (size_t i = 0; i < mesh->mergeLen; i++) {
    merge[mesh->mergeFromVert[i]] = mesh->mergeToVert[i];
  }

  // Find open edges using sorted array (simulating multiset)
  // Each triangle has 3 edges; allocate worst case
  size_t maxEdges = numTri * 3;
  MergeEdgePair *edges = (MergeEdgePair *)malloc(maxEdges * sizeof(MergeEdgePair));
  size_t numEdges = 0;

  static const int next3[3] = {1, 2, 0};
  for (size_t tri = 0; tri < numTri; tri++) {
    for (int i = 0; i < 3; i++) {
      MergeEdgePair edge;
      edge.first = merge[mesh->triVerts[3 * tri + next3[i]]];
      edge.second = merge[mesh->triVerts[3 * tri + i]];
      // Search for this edge in openEdges (C++ does: find(edge))
      bool found = false;
      for (size_t j = 0; j < numEdges; j++) {
        if (edges[j].first == edge.first && edges[j].second == edge.second) {
          // Erase this entry by shifting
          memmove(&edges[j], &edges[j + 1], (numEdges - j - 1) * sizeof(MergeEdgePair));
          numEdges--;
          found = true;
          break;
        }
      }
      if (!found) {
        // Swap and insert (C++: swap(edge.first, edge.second); insert(edge))
        MergeEdgePair swapped = {edge.second, edge.first};
        edges[numEdges++] = swapped;
      }
    }
  }

  if (numEdges == 0) {
    free(merge);
    free(edges);
    return false;
  }

  // Collect open vertex indices
  int *openVerts = (int *)malloc(numEdges * sizeof(int));
  for (size_t i = 0; i < numEdges; i++) {
    openVerts[i] = edges[i].first;
  }
  free(edges);

  // Compute bounding box from all vertex positions
  ManifoldBox bBox = manifold_box_empty();
  for (size_t i = 0; i < numVert; i++) {
    ManifoldVec3 p = {
      (double)mesh->vertProperties[mesh->numProp * i + 0],
      (double)mesh->vertProperties[mesh->numProp * i + 1],
      (double)mesh->vertProperties[mesh->numProp * i + 2]
    };
    manifold_box_union_point(&bBox, p);
  }

  // Compute tolerance (float precision: use FLT_EPSILON)
  double tolerance = fmax((double)mesh->tolerance,
                          (double)FLT_EPSILON * manifold_box_scale(bBox));

  // Build vertex boxes and Morton codes for open vertices
  size_t numOpen = numEdges;
  ManifoldBox *vertBox = (ManifoldBox *)malloc(numOpen * sizeof(ManifoldBox));
  uint32_t *vertMorton = (uint32_t *)malloc(numOpen * sizeof(uint32_t));

  for (size_t i = 0; i < numOpen; i++) {
    int vert = openVerts[i];
    ManifoldVec3 center = {
      (double)mesh->vertProperties[mesh->numProp * vert + 0],
      (double)mesh->vertProperties[mesh->numProp * vert + 1],
      (double)mesh->vertProperties[mesh->numProp * vert + 2]
    };
    vertBox[i].min = vec3_sub(center, manifold_vec3_splat(tolerance / 2.0));
    vertBox[i].max = vec3_add(center, manifold_vec3_splat(tolerance / 2.0));
    vertMorton[i] = manifold_morton_code(center, bBox);
  }

  // Sort by Morton code (stable sort via qsort with index comparison)
  int *vertNew2Old = (int *)malloc(numOpen * sizeof(int));
  for (size_t i = 0; i < numOpen; i++) vertNew2Old[i] = (int)i;

  g_merge_morton_ctx.morton = vertMorton;
  qsort(vertNew2Old, numOpen, sizeof(int), merge_morton_compare);

  // Permute vertMorton, vertBox, openVerts by vertNew2Old
  {
    uint32_t *tmpMorton = (uint32_t *)malloc(numOpen * sizeof(uint32_t));
    ManifoldBox *tmpBox = (ManifoldBox *)malloc(numOpen * sizeof(ManifoldBox));
    int *tmpVerts = (int *)malloc(numOpen * sizeof(int));
    for (size_t i = 0; i < numOpen; i++) {
      tmpMorton[i] = vertMorton[vertNew2Old[i]];
      tmpBox[i] = vertBox[vertNew2Old[i]];
      tmpVerts[i] = openVerts[vertNew2Old[i]];
    }
    memcpy(vertMorton, tmpMorton, numOpen * sizeof(uint32_t));
    memcpy(vertBox, tmpBox, numOpen * sizeof(ManifoldBox));
    memcpy(openVerts, tmpVerts, numOpen * sizeof(int));
    free(tmpMorton);
    free(tmpBox);
    free(tmpVerts);
  }
  free(vertNew2Old);

  // Build collider and find collisions
  ManifoldCollider collider = manifold_collider_create(vertBox, vertMorton, numOpen);
  ManifoldDisjointSets uf = manifold_disjoint_sets_create((uint32_t)numVert);

  MergeCollisionCtx ctx = {&uf, openVerts};
  manifold_collider_collisions(&collider, vertBox, numOpen,
                               merge_collision_callback, &ctx, false);

  // Also unite existing merge pairs
  for (size_t i = 0; i < mesh->mergeLen; i++) {
    manifold_disjoint_sets_unite(&uf, (uint32_t)mesh->mergeFromVert[i],
                                 (uint32_t)mesh->mergeToVert[i]);
  }

  // Rebuild merge vectors
  free(mesh->mergeFromVert);
  free(mesh->mergeToVert);
  mesh->mergeFromVert = NULL;
  mesh->mergeToVert = NULL;
  mesh->mergeLen = 0;

  // Count first
  size_t mergeCount = 0;
  for (size_t v = 0; v < numVert; v++) {
    size_t mergeTo = manifold_disjoint_sets_find(&uf, (uint32_t)v);
    if (mergeTo != v) mergeCount++;
  }

  if (mergeCount > 0) {
    mesh->mergeFromVert = (int *)malloc(mergeCount * sizeof(int));
    mesh->mergeToVert = (int *)malloc(mergeCount * sizeof(int));
    mesh->mergeLen = mergeCount;
    size_t idx = 0;
    for (size_t v = 0; v < numVert; v++) {
      size_t mergeTo = manifold_disjoint_sets_find(&uf, (uint32_t)v);
      if (mergeTo != v) {
        mesh->mergeFromVert[idx] = (int)v;
        mesh->mergeToVert[idx] = (int)mergeTo;
        idx++;
      }
    }
  }

  // Cleanup
  manifold_disjoint_sets_free(&uf);
  manifold_collider_free(&collider);
  free(vertBox);
  free(vertMorton);
  free(openVerts);
  free(merge);

  return true;
}

Manifold manifold_from_meshgl(const ManifoldMeshGL *mesh) {
  if (!mesh || mesh->vertLen == 0 || mesh->triLen == 0 || mesh->numProp < 3) {
    Manifold m;
    manifold_impl_init(&m.impl);
    return m;
  }
  // Check for NaN/Inf vertices
  for (size_t i = 0; i < mesh->vertLen * (size_t)mesh->numProp; i++) {
    if (!isfinite(mesh->vertProperties[i])) {
      Manifold m;
      manifold_impl_init(&m.impl);
      return m;
    }
  }
  // Check for NaN/Inf transforms
  for (size_t i = 0; i < mesh->runTransformLen; i++) {
    if (!isfinite(mesh->runTransform[i])) {
      Manifold m;
      manifold_impl_init(&m.impl);
      m.impl.status = MANIFOLD_ERROR_INVALID_CONSTRUCTION;
      return m;
    }
  }

  // Validate runIndex length (C++ impl.h line 95-100)
  if (mesh->runOriginalIDLen > 0 && mesh->runIndexLen > 0 &&
      mesh->runOriginalIDLen + 1 != mesh->runIndexLen &&
      mesh->runOriginalIDLen != mesh->runIndexLen) {
    Manifold m;
    manifold_impl_init(&m.impl);
    m.impl.status = MANIFOLD_ERROR_RUN_INDEX_WRONG_LENGTH;
    return m;
  }

  const size_t numVert = mesh->vertLen;
  const size_t numTri = mesh->triLen;
  const int numProp = mesh->numProp - 3;  // custom properties (excluding xyz)

  // Build prop2vert merge map
  int *prop2vert = NULL;
  if (mesh->mergeLen > 0 && mesh->mergeFromVert && mesh->mergeToVert) {
    prop2vert = (int *)malloc(numVert * sizeof(int));
    for (size_t i = 0; i < numVert; i++) prop2vert[i] = (int)i;
    for (size_t i = 0; i < mesh->mergeLen; i++) {
      int from = mesh->mergeFromVert[i];
      int to = mesh->mergeToVert[i];
      if (from < 0 || (size_t)from >= numVert || to < 0 || (size_t)to >= numVert) {
        free(prop2vert);
        Manifold m;
        manifold_impl_init(&m.impl);
        m.impl.status = MANIFOLD_ERROR_MERGE_INDEX_OUT_OF_BOUNDS;
        return m;
      }
      prop2vert[from] = to;
    }
  }

  Manifold m;
  manifold_impl_init(&m.impl);
  m.impl.numProp = numProp;
  m.impl.tolerance = (double)mesh->tolerance;

  // Copy vertex positions and properties
  vec_vec3_resize(&m.impl.vertPos, numVert);
  for (size_t i = 0; i < numVert; i++) {
    m.impl.vertPos.data[i].x = (double)mesh->vertProperties[i * mesh->numProp + 0];
    m.impl.vertPos.data[i].y = (double)mesh->vertProperties[i * mesh->numProp + 1];
    m.impl.vertPos.data[i].z = (double)mesh->vertProperties[i * mesh->numProp + 2];
  }
  if (numProp > 0) {
    vec_double_resize(&m.impl.properties, numVert * numProp);
    for (size_t i = 0; i < numVert; i++)
      for (int j = 0; j < numProp; j++)
        m.impl.properties.data[i * numProp + j] =
            (double)mesh->vertProperties[i * mesh->numProp + 3 + j];
  }

  // Copy halfedge tangents
  if (mesh->halfedgeTangentLen > 0 && mesh->halfedgeTangent) {
    size_t nht = mesh->halfedgeTangentLen / 4;
    vec_vec4_resize(&m.impl.halfedgeTangent, nht);
    for (size_t i = 0; i < nht; i++) {
      m.impl.halfedgeTangent.data[i].x = (double)mesh->halfedgeTangent[4*i+0];
      m.impl.halfedgeTangent.data[i].y = (double)mesh->halfedgeTangent[4*i+1];
      m.impl.halfedgeTangent.data[i].z = (double)mesh->halfedgeTangent[4*i+2];
      m.impl.halfedgeTangent.data[i].w = (double)mesh->halfedgeTangent[4*i+3];
    }
  }

  // Build runIndex (ensure it has runOriginalIDLen + 1 entries)
  size_t numRuns = mesh->runOriginalIDLen;
  int *runIdx = NULL;
  size_t runIdxLen = 0;
  if (numRuns == 0) {
    runIdxLen = 2;
    runIdx = (int *)malloc(2 * sizeof(int));
    runIdx[0] = 0;
    runIdx[1] = (int)(numTri * 3);
  } else if (mesh->runIndexLen == numRuns) {
    runIdxLen = numRuns + 1;
    runIdx = (int *)malloc(runIdxLen * sizeof(int));
    memcpy(runIdx, mesh->runIndex, numRuns * sizeof(int));
    runIdx[numRuns] = (int)(numTri * 3);
  } else if (mesh->runIndexLen == numRuns + 1) {
    runIdxLen = numRuns + 1;
    runIdx = (int *)malloc(runIdxLen * sizeof(int));
    memcpy(runIdx, mesh->runIndex, runIdxLen * sizeof(int));
  } else {
    runIdxLen = 2;
    runIdx = (int *)malloc(2 * sizeof(int));
    runIdx[0] = 0;
    runIdx[1] = (int)(numTri * 3);
  }
  size_t actualRuns = runIdxLen - 1;

  // Build runOriginalID
  uint32_t startID = manifold_reserve_ids((uint32_t)(actualRuns > 0 ? actualRuns : 1));
  uint32_t *runOrigID = (uint32_t *)malloc(actualRuns * sizeof(uint32_t));
  if (numRuns == 0) {
    runOrigID[0] = startID;
    actualRuns = 1;
  } else {
    memcpy(runOrigID, mesh->runOriginalID, actualRuns * sizeof(uint32_t));
  }

  // Set up triRef and meshIDtransform
  vec_triref_resize(&m.impl.meshRelation.triRef, numTri);
  for (size_t r = 0; r < actualRuns; r++) {
    int meshID = (int)(startID + r);
    int origID = (int)runOrigID[r];
    for (size_t tri = (size_t)runIdx[r] / 3; tri < (size_t)runIdx[r+1] / 3; tri++) {
      if (tri >= numTri) break;
      ManifoldTriRef *ref = &m.impl.meshRelation.triRef.data[tri];
      ref->meshID = meshID;
      ref->originalID = origID;
      ref->faceID = (mesh->faceIDLen > 0 && mesh->faceID) ? mesh->faceID[tri] : -1;
      ref->coplanarID = (int)tri;
    }

    ManifoldRelation rel;
    rel.originalID = origID;
    rel.backSide = false;
    if (mesh->runTransformLen > 0 && mesh->runTransform) {
      const float *mt = mesh->runTransform + 12 * r;
      rel.transform.cols[0] = (ManifoldVec3){mt[0], mt[1], mt[2]};
      rel.transform.cols[1] = (ManifoldVec3){mt[3], mt[4], mt[5]};
      rel.transform.cols[2] = (ManifoldVec3){mt[6], mt[7], mt[8]};
      rel.transform.cols[3] = (ManifoldVec3){mt[9], mt[10], mt[11]};
    } else {
      rel.transform = mat3x4_identity();
    }
    manifold_meshrelation_insert(&m.impl.meshRelation, meshID, rel);
  }

  // Build triangles, applying merge map
  bool needsPropMap = (numProp > 0 && prop2vert != NULL);
  ManifoldVecIVec3 triProp = {0};
  ManifoldVecIVec3 triVert = {0};
  ManifoldVecTriRef filteredTriRef = {0};

  for (size_t i = 0; i < numTri; i++) {
    ManifoldIVec3 tp, tv;
    bool valid = true;
    for (int j = 0; j < 3; j++) {
      int vert = mesh->triVerts[3 * i + j];
      if (vert < 0 || (size_t)vert >= numVert) {
        valid = false;
        break;
      }
      ((int*)&tp)[j] = vert;
      ((int*)&tv)[j] = prop2vert ? prop2vert[vert] : vert;
    }
    if (!valid) continue;
    // Skip degenerate triangles
    if (tv.x == tv.y || tv.y == tv.z || tv.z == tv.x) continue;
    if (needsPropMap) {
      vec_ivec3_push(&triProp, tp);
      vec_ivec3_push(&triVert, tv);
    } else {
      vec_ivec3_push(&triProp, tv);
    }
    vec_triref_push(&filteredTriRef, m.impl.meshRelation.triRef.data[i]);
  }

  // Replace triRef with filtered version
  vec_triref_free(&m.impl.meshRelation.triRef);
  m.impl.meshRelation.triRef = filteredTriRef;

  ManifoldVecIVec3 emptyTriVert = {0};
  manifold_impl_create_halfedges(&m.impl, &triProp, needsPropMap ? &triVert : &emptyTriVert);

  if (!manifold_impl_is_manifold(&m.impl)) {
    manifold_impl_make_empty(&m.impl, MANIFOLD_ERROR_NOT_MANIFOLD);
    vec_ivec3_free(&emptyTriVert);
    goto cleanup;
  }

  manifold_impl_calculate_bbox(&m.impl);
  manifold_impl_set_epsilon(&m.impl, -1.0, true);  // float precision

  manifold_impl_cleanup_topology(&m.impl);
  manifold_impl_dedupe_prop_verts(&m.impl);
  manifold_impl_set_normals_and_coplanar(&m.impl);
  manifold_impl_remove_degenerates(&m.impl, 0);
  manifold_impl_remove_unreferenced_verts(&m.impl);
  manifold_impl_sort_geometry(&m.impl);

  if (!manifold_impl_is_finite(&m.impl)) {
    manifold_impl_make_empty(&m.impl, MANIFOLD_ERROR_NON_FINITE_VERTEX);
    vec_ivec3_free(&emptyTriVert);
    goto cleanup;
  }

  // A Manifold from input mesh is never an original
  m.impl.meshRelation.originalID = -1;

  vec_ivec3_free(&emptyTriVert);

cleanup:
  free(prop2vert);
  free(runIdx);
  free(runOrigID);
  vec_ivec3_free(&triProp);
  vec_ivec3_free(&triVert);
  return m;
}

// ============== GetMeshGL (mirrors C++ GetMeshGLImpl<float>) ==============

ManifoldMeshGL manifold_get_meshgl(const Manifold *m) {
  ManifoldMeshGL out = manifold_meshgl_empty();
  const ManifoldImpl *impl = &m->impl;
  const int numProp = impl->numProp;
  const int numVert = (int)manifold_impl_num_prop_vert(impl);
  const int numTri = (int)manifold_impl_num_tri(impl);

  if (numTri == 0) return out;

  const bool isOriginal = impl->meshRelation.originalID >= 0;

  out.numProp = 3 + numProp;
  out.tolerance = (float)impl->tolerance;
  // Ensure float-level tolerance
  {
    double bboxScale = 0;
    ManifoldVec3 sz = {impl->bBox.max.x - impl->bBox.min.x,
                       impl->bBox.max.y - impl->bBox.min.y,
                       impl->bBox.max.z - impl->bBox.min.z};
    if (sz.x > bboxScale) bboxScale = sz.x;
    if (sz.y > bboxScale) bboxScale = sz.y;
    if (sz.z > bboxScale) bboxScale = sz.z;
    double floatTol = (double)FLT_EPSILON * bboxScale;
    if (floatTol > out.tolerance) out.tolerance = (float)floatTol;
  }

  out.triVerts = (int *)malloc(3 * numTri * sizeof(int));
  out.triLen = numTri;

  // Copy halfedge tangents
  const int numHalfedge = (int)impl->halfedgeTangent.len;
  if (numHalfedge > 0) {
    out.halfedgeTangentLen = 4 * numHalfedge;
    out.halfedgeTangent = (float *)malloc(out.halfedgeTangentLen * sizeof(float));
    for (int i = 0; i < numHalfedge; i++) {
      ManifoldVec4 t = impl->halfedgeTangent.data[i];
      out.halfedgeTangent[4*i+0] = (float)t.x;
      out.halfedgeTangent[4*i+1] = (float)t.y;
      out.halfedgeTangent[4*i+2] = (float)t.z;
      out.halfedgeTangent[4*i+3] = (float)t.w;
    }
  }

  // Sort triangles into runs by originalID/meshID
  out.faceIDLen = numTri;
  out.faceID = (int *)malloc(numTri * sizeof(int));

  int *triNew2Old = (int *)malloc(numTri * sizeof(int));
  for (int i = 0; i < numTri; i++) triNew2Old[i] = i;

  const ManifoldTriRef *triRef = impl->meshRelation.triRef.data;

  if (!isOriginal && numTri > 0 && impl->meshRelation.triRef.len > 0) {
    // Stable sort by originalID, then meshID
    // Simple insertion sort (stable)
    for (int i = 1; i < numTri; i++) {
      int key = triNew2Old[i];
      const ManifoldTriRef *refKey = &triRef[key];
      int j = i - 1;
      while (j >= 0) {
        const ManifoldTriRef *refJ = &triRef[triNew2Old[j]];
        bool shouldSwap = (refJ->originalID > refKey->originalID) ||
                          (refJ->originalID == refKey->originalID &&
                           refJ->meshID > refKey->meshID);
        if (!shouldSwap) break;
        triNew2Old[j + 1] = triNew2Old[j];
        j--;
      }
      triNew2Old[j + 1] = key;
    }
  }

  // Build runs, faceID, triVerts
  // Temporary dynamic arrays for runs
  size_t runCap = 16;
  size_t runCount = 0;
  int *runIdxArr = (int *)malloc(runCap * sizeof(int));
  uint32_t *runOrigArr = (uint32_t *)malloc(runCap * sizeof(uint32_t));
  float *runTransArr = NULL;
  size_t runTransCap = 0;
  if (!isOriginal) {
    runTransCap = runCap * 12;
    runTransArr = (float *)malloc(runTransCap * sizeof(float));
  }

  // Track which meshIDtransform entries are used
  bool *meshIDused = (bool *)calloc(impl->meshRelation.meshIDtransform.len, sizeof(bool));

  int lastMeshID = -1;
  for (int tri = 0; tri < numTri; tri++) {
    int oldTri = triNew2Old[tri];
    const ManifoldTriRef *ref = &triRef[oldTri];
    int meshID = ref->meshID;

    out.faceID[tri] = (ref->faceID >= 0) ? ref->faceID : ref->coplanarID;
    for (int i = 0; i < 3; i++)
      out.triVerts[3 * tri + i] = impl->halfedge.data[3 * oldTri + i].startVert;

    if (meshID != lastMeshID) {
      if (runCount >= runCap) {
        runCap *= 2;
        runIdxArr = (int *)realloc(runIdxArr, runCap * sizeof(int));
        runOrigArr = (uint32_t *)realloc(runOrigArr, runCap * sizeof(uint32_t));
        if (!isOriginal) {
          runTransCap = runCap * 12;
          runTransArr = (float *)realloc(runTransArr, runTransCap * sizeof(float));
        }
      }
      runIdxArr[runCount] = 3 * tri;

      ManifoldRelation rel = {0, mat3x4_identity(), false};
      for (size_t k = 0; k < impl->meshRelation.meshIDtransform.len; k++) {
        if (impl->meshRelation.meshIDtransform.data[k].key == meshID) {
          rel = impl->meshRelation.meshIDtransform.data[k].value;
          meshIDused[k] = true;
          break;
        }
      }

      runOrigArr[runCount] = (uint32_t)rel.originalID;

      if (!isOriginal) {
        float *mt = runTransArr + 12 * runCount;
        mt[0] = (float)rel.transform.cols[0].x;
        mt[1] = (float)rel.transform.cols[0].y;
        mt[2] = (float)rel.transform.cols[0].z;
        mt[3] = (float)rel.transform.cols[1].x;
        mt[4] = (float)rel.transform.cols[1].y;
        mt[5] = (float)rel.transform.cols[1].z;
        mt[6] = (float)rel.transform.cols[2].x;
        mt[7] = (float)rel.transform.cols[2].y;
        mt[8] = (float)rel.transform.cols[2].z;
        mt[9] = (float)rel.transform.cols[3].x;
        mt[10] = (float)rel.transform.cols[3].y;
        mt[11] = (float)rel.transform.cols[3].z;
      }
      runCount++;
      lastMeshID = meshID;
    }
  }

  // Add runs for originals that did not contribute any faces
  for (size_t k = 0; k < impl->meshRelation.meshIDtransform.len; k++) {
    if (!meshIDused[k]) {
      if (runCount >= runCap) {
        runCap *= 2;
        runIdxArr = (int *)realloc(runIdxArr, runCap * sizeof(int));
        runOrigArr = (uint32_t *)realloc(runOrigArr, runCap * sizeof(uint32_t));
        if (!isOriginal) {
          runTransCap = runCap * 12;
          runTransArr = (float *)realloc(runTransArr, runTransCap * sizeof(float));
        }
      }
      ManifoldRelation *rel = &impl->meshRelation.meshIDtransform.data[k].value;
      runIdxArr[runCount] = 3 * numTri;
      runOrigArr[runCount] = (uint32_t)rel->originalID;
      if (!isOriginal) {
        float *mt = runTransArr + 12 * runCount;
        mt[0] = (float)rel->transform.cols[0].x;
        mt[1] = (float)rel->transform.cols[0].y;
        mt[2] = (float)rel->transform.cols[0].z;
        mt[3] = (float)rel->transform.cols[1].x;
        mt[4] = (float)rel->transform.cols[1].y;
        mt[5] = (float)rel->transform.cols[1].z;
        mt[6] = (float)rel->transform.cols[2].x;
        mt[7] = (float)rel->transform.cols[2].y;
        mt[8] = (float)rel->transform.cols[2].z;
        mt[9] = (float)rel->transform.cols[3].x;
        mt[10] = (float)rel->transform.cols[3].y;
        mt[11] = (float)rel->transform.cols[3].z;
      }
      runCount++;
    }
  }
  free(meshIDused);

  // Final runIndex: runCount + 1 entries
  out.runOriginalIDLen = runCount;
  out.runOriginalID = (uint32_t *)malloc(runCount * sizeof(uint32_t));
  memcpy(out.runOriginalID, runOrigArr, runCount * sizeof(uint32_t));

  out.runIndexLen = runCount + 1;
  out.runIndex = (int *)malloc((runCount + 1) * sizeof(int));
  memcpy(out.runIndex, runIdxArr, runCount * sizeof(int));
  out.runIndex[runCount] = 3 * numTri;

  if (!isOriginal && runTransArr) {
    out.runTransformLen = 12 * runCount;
    out.runTransform = (float *)malloc(out.runTransformLen * sizeof(float));
    memcpy(out.runTransform, runTransArr, out.runTransformLen * sizeof(float));
  }

  free(runIdxArr);
  free(runOrigArr);
  free(runTransArr);

  // Build vertex properties with deduplication for properties
  if (numProp == 0) {
    // No custom properties - simple case
    out.vertLen = impl->vertPos.len;
    out.vertProperties = (float *)malloc(3 * out.vertLen * sizeof(float));
    for (size_t i = 0; i < out.vertLen; i++) {
      ManifoldVec3 v = impl->vertPos.data[i];
      out.vertProperties[3*i+0] = (float)v.x;
      out.vertProperties[3*i+1] = (float)v.y;
      out.vertProperties[3*i+2] = (float)v.z;
    }
  } else {
    // Duplicate verts with different property indices
    int *vert2idx = (int *)malloc(impl->vertPos.len * sizeof(int));
    for (size_t i = 0; i < impl->vertPos.len; i++) vert2idx[i] = -1;

    // vertPropPair: for each geometric vertex, list of (propVert, outputIdx) pairs
    typedef struct { int prop; int idx; } PropPair;
    typedef struct { PropPair *data; size_t len; size_t cap; } PropPairVec;
    PropPairVec *vertPropPair = (PropPairVec *)calloc(impl->vertPos.len, sizeof(PropPairVec));

    // Dynamic array for output vert properties
    size_t vpCap = (size_t)numVert * out.numProp;
    size_t vpLen = 0;
    float *vpData = (float *)malloc(vpCap * sizeof(float));

    // Dynamic arrays for merge verts
    size_t mergeCap = 64;
    size_t mergeCount = 0;
    int *mergeFrom = (int *)malloc(mergeCap * sizeof(int));
    int *mergeTo = (int *)malloc(mergeCap * sizeof(int));

    for (size_t run = 0; run < out.runOriginalIDLen; run++) {
      for (size_t tri = (size_t)out.runIndex[run] / 3;
           tri < (size_t)out.runIndex[run + 1] / 3; tri++) {
        for (int i = 0; i < 3; i++) {
          int prop = impl->halfedge.data[3 * triNew2Old[tri] + i].propVert;
          int vert = out.triVerts[3 * tri + i];

          PropPairVec *bin = &vertPropPair[vert];
          bool found = false;
          for (size_t b = 0; b < bin->len; b++) {
            if (bin->data[b].prop == prop) {
              out.triVerts[3 * tri + i] = bin->data[b].idx;
              found = true;
              break;
            }
          }
          if (found) continue;

          int idx = (int)(vpLen / out.numProp);
          out.triVerts[3 * tri + i] = idx;

          if (bin->len >= bin->cap) {
            bin->cap = bin->cap ? bin->cap * 2 : 4;
            bin->data = (PropPair *)realloc(bin->data, bin->cap * sizeof(PropPair));
          }
          bin->data[bin->len++] = (PropPair){prop, idx};

          if (vpLen + out.numProp > vpCap) {
            vpCap *= 2;
            vpData = (float *)realloc(vpData, vpCap * sizeof(float));
          }
          ManifoldVec3 v = impl->vertPos.data[vert];
          vpData[vpLen++] = (float)v.x;
          vpData[vpLen++] = (float)v.y;
          vpData[vpLen++] = (float)v.z;
          for (int p = 0; p < numProp; p++) {
            vpData[vpLen++] = (float)impl->properties.data[prop * numProp + p];
          }

          if (vert2idx[vert] == -1) {
            vert2idx[vert] = idx;
          } else {
            if (mergeCount >= mergeCap) {
              mergeCap *= 2;
              mergeFrom = (int *)realloc(mergeFrom, mergeCap * sizeof(int));
              mergeTo = (int *)realloc(mergeTo, mergeCap * sizeof(int));
            }
            mergeFrom[mergeCount] = idx;
            mergeTo[mergeCount] = vert2idx[vert];
            mergeCount++;
          }
        }
      }
    }

    out.vertLen = vpLen / out.numProp;
    out.vertProperties = (float *)realloc(vpData, vpLen * sizeof(float));

    if (mergeCount > 0) {
      out.mergeLen = mergeCount;
      out.mergeFromVert = (int *)realloc(mergeFrom, mergeCount * sizeof(int));
      out.mergeToVert = (int *)realloc(mergeTo, mergeCount * sizeof(int));
    } else {
      free(mergeFrom);
      free(mergeTo);
    }

    for (size_t i = 0; i < impl->vertPos.len; i++) free(vertPropPair[i].data);
    free(vertPropPair);
    free(vert2idx);
  }

  free(triNew2Old);
  return out;
}

void manifold_free_meshgl(ManifoldMeshGL *mgl) {
  if (!mgl) return;
  free(mgl->vertProperties);
  free(mgl->triVerts);
  free(mgl->mergeFromVert);
  free(mgl->mergeToVert);
  free(mgl->runOriginalID);
  free(mgl->runIndex);
  free(mgl->runTransform);
  free(mgl->faceID);
  free(mgl->halfedgeTangent);
  *mgl = manifold_meshgl_empty();
}

// ---------- OBJ Import ----------

// ============== Project ==============

void manifold_polygons2d_free(ManifoldPolygons2D *p) {
  if (!p) return;
  if (p->polys) {
    free(p->polys);
    p->polys = NULL;
  }
  if (p->polySizes) {
    free(p->polySizes);
    p->polySizes = NULL;
  }
  p->numPolys = 0;
}

// Port of Manifold::Impl::Slice(double height)
// Slices the manifold at the given z-height, returning 2D cross-section polygons.
ManifoldPolygons2D manifold_slice(const Manifold *m, double height) {
  ManifoldPolygons2D result = {NULL, NULL, 0};
  const ManifoldImpl *impl = &m->impl;

  if (impl->halfedge.len == 0) return result;

  const size_t numTri = impl->halfedge.len / 3;
  const ManifoldHalfedge *he = impl->halfedge.data;
  const ManifoldVec3 *vp = impl->vertPos.data;

  // Find triangles that cross the slice plane: min z <= height && max z > height
  bool *inSet = (bool *)calloc(numTri, sizeof(bool));
  int numInSet = 0;
  for (size_t tri = 0; tri < numTri; tri++) {
    double zmin = 1e308, zmax = -1e308;
    for (int j = 0; j < 3; j++) {
      double z = vp[he[3 * tri + j].startVert].z;
      if (z < zmin) zmin = z;
      if (z > zmax) zmax = z;
    }
    if (zmin <= height && zmax > height) {
      inSet[tri] = true;
      numInSet++;
    }
  }

  if (numInSet == 0) {
    free(inSet);
    return result;
  }

  // Trace polygon loops
  int maxPolys = numInSet;
  ManifoldVec2 *polyBuf = (ManifoldVec2 *)malloc(numInSet * sizeof(ManifoldVec2));
  int *polySizes = (int *)malloc(maxPolys * sizeof(int));
  int numPolys = 0;
  int totalVerts = 0;

  while (numInSet > 0) {
    // Find first triangle still in set
    int startTri = -1;
    for (size_t tri = 0; tri < numTri; tri++) {
      if (inSet[tri]) { startTri = (int)tri; break; }
    }
    if (startTri < 0) break;

    // Find starting halfedge k: the one where startVert is above and endVert is at/below
    int k = 0;
    for (int j = 0; j < 3; j++) {
      int next = (j + 1) % 3;
      if (vp[he[3 * startTri + j].startVert].z > height &&
          vp[he[3 * startTri + next].startVert].z <= height) {
        k = next;
        break;
      }
    }

    int polyStart = totalVerts;
    int tri = startTri;
    do {
      inSet[tri] = false;
      numInSet--;

      // If endVert of current halfedge is at/below height, advance k
      if (vp[he[3 * tri + k].endVert].z <= height) {
        k = (k + 1) % 3;
      }

      // Halfedge k goes from below to above (crosses the plane upward)
      ManifoldHalfedge up = he[3 * tri + k];
      ManifoldVec3 below = vp[up.startVert];
      ManifoldVec3 above = vp[up.endVert];
      double a = (height - below.z) / (above.z - below.z);

      // Grow buffer if needed
      if (totalVerts >= numInSet + totalVerts + 100) {
        // shouldn't happen, but safe
      }
      polyBuf[totalVerts].x = below.x + a * (above.x - below.x);
      polyBuf[totalVerts].y = below.y + a * (above.y - below.y);
      totalVerts++;

      // Follow paired halfedge to next triangle
      int pair = up.pairedHalfedge;
      tri = pair / 3;
      k = (pair % 3 + 1) % 3;
    } while (tri != startTri);

    polySizes[numPolys] = totalVerts - polyStart;
    numPolys++;
  }

  // Build result
  result.numPolys = numPolys;
  result.polySizes = (int *)malloc(numPolys * sizeof(int));
  memcpy(result.polySizes, polySizes, numPolys * sizeof(int));
  result.polys = (ManifoldVec2 *)malloc(totalVerts * sizeof(ManifoldVec2));
  memcpy(result.polys, polyBuf, totalVerts * sizeof(ManifoldVec2));

  free(inSet);
  free(polyBuf);
  free(polySizes);
  return result;
}

ManifoldPolygons2D manifold_project(const Manifold *m) {
  ManifoldPolygons2D result = {NULL, NULL, 0};
  const ManifoldImpl *impl = &m->impl;

  if (impl->halfedge.len == 0 || impl->faceNormal.len == 0)
    return result;

  const size_t numHe = impl->halfedge.len;
  const ManifoldHalfedge *he = impl->halfedge.data;
  const ManifoldVec3 *fn = impl->faceNormal.data;

  // Find cusp halfedges where face of this edge has normal.z >= 0,
  // and face of paired edge has normal.z < 0 (silhouette edges from above).
  ManifoldHalfedge *cusps = (ManifoldHalfedge *)malloc(numHe * sizeof(ManifoldHalfedge));
  int numCusps = 0;

  for (size_t i = 0; i < numHe; i++) {
    int paired = he[i].pairedHalfedge;
    if (paired < 0 || (size_t)paired >= numHe) continue;
    int thisFace = he[paired].pairedHalfedge / 3;
    int pairedFace = paired / 3;
    if ((size_t)thisFace >= impl->faceNormal.len ||
        (size_t)pairedFace >= impl->faceNormal.len) continue;
    if (fn[thisFace].z >= 0 && fn[pairedFace].z < 0) {
      cusps[numCusps] = he[i];
      numCusps++;
    }
  }

  if (numCusps == 0) {
    free(cusps);
    return result;
  }

  // AssembleHalfedges: build closed loops from cusp edges
  bool *used = (bool *)calloc(numCusps, sizeof(bool));

  // Build loops
  int maxPolys = numCusps;
  int *loopBuf = (int *)malloc(numCusps * sizeof(int));  // halfedge indices
  int *loopStarts = (int *)malloc(maxPolys * sizeof(int));
  int *loopSizes = (int *)malloc(maxPolys * sizeof(int));
  int numLoops = 0;
  int totalVerts = 0;

  while (1) {
    // Find first unused cusp
    int startIdx = -1;
    for (int i = 0; i < numCusps; i++) {
      if (!used[i]) { startIdx = i; break; }
    }
    if (startIdx < 0) break;

    loopStarts[numLoops] = totalVerts;
    int loopLen = 0;
    int curIdx = startIdx;

    do {
      used[curIdx] = true;
      loopBuf[totalVerts++] = curIdx;
      loopLen++;

      // Find next cusp whose startVert == current endVert
      int endV = cusps[curIdx].endVert;
      int nextIdx = -1;
      for (int i = 0; i < numCusps; i++) {
        if (!used[i] && cusps[i].startVert == endV) {
          nextIdx = i;
          break;
        }
      }
      if (nextIdx < 0) {
        // Check if we've looped back to start
        if (cusps[startIdx].startVert == endV) break;
        // Broken loop
        break;
      }
      curIdx = nextIdx;
    } while (curIdx != startIdx);

    loopSizes[numLoops] = loopLen;
    numLoops++;
  }

  // Build polygon output: project vertices to 2D (just x,y for z-axis projection)
  result.numPolys = numLoops;
  result.polySizes = (int *)malloc(numLoops * sizeof(int));
  result.polys = (ManifoldVec2 *)malloc(totalVerts * sizeof(ManifoldVec2));
  int outIdx = 0;
  for (int p = 0; p < numLoops; p++) {
    result.polySizes[p] = loopSizes[p];
    for (int v = 0; v < loopSizes[p]; v++) {
      int ci = loopBuf[loopStarts[p] + v];
      int sv = cusps[ci].startVert;
      result.polys[outIdx].x = impl->vertPos.data[sv].x;
      result.polys[outIdx].y = impl->vertPos.data[sv].y;
      outIdx++;
    }
  }

  free(cusps);
  free(used);
  free(loopBuf);
  free(loopStarts);
  free(loopSizes);
  return result;
}

// Compute winding number at a point with respect to all polygon edges
static int winding_number_at_point(double px, double py,
                                    const ManifoldVec2 *allVerts,
                                    const int *polySizes, int numPolys) {
  int winding = 0;
  int offset = 0;
  for (int p = 0; p < numPolys; p++) {
    int n = polySizes[p];
    for (int i = 0; i < n; i++) {
      int j = (i + 1) % n;
      double ay = allVerts[offset + i].y, by = allVerts[offset + j].y;
      double ax = allVerts[offset + i].x, bx = allVerts[offset + j].x;
      if (ay <= py) {
        if (by > py) {
          double cross = (bx - ax) * (py - ay) - (px - ax) * (by - ay);
          if (cross > 0) winding++;
        }
      } else {
        if (by <= py) {
          double cross = (bx - ax) * (py - ay) - (px - ax) * (by - ay);
          if (cross < 0) winding--;
        }
      }
    }
    offset += n;
  }
  return winding;
}

// Compute area under positive fill rule using rasterization.
// Handles self-intersecting polygons correctly via winding number.
// Equivalent to CrossSection(polys, FillRule::Positive).Area() in C++.
double manifold_cross_section_area(const ManifoldPolygons2D *polys) {
  if (!polys || polys->numPolys == 0) return 0.0;

  int totalVerts = 0;
  for (int i = 0; i < polys->numPolys; i++)
    totalVerts += polys->polySizes[i];
  if (totalVerts == 0) return 0.0;

  double minX = polys->polys[0].x, maxX = minX;
  double minY = polys->polys[0].y, maxY = minY;
  for (int i = 1; i < totalVerts; i++) {
    if (polys->polys[i].x < minX) minX = polys->polys[i].x;
    if (polys->polys[i].x > maxX) maxX = polys->polys[i].x;
    if (polys->polys[i].y < minY) minY = polys->polys[i].y;
    if (polys->polys[i].y > maxY) maxY = polys->polys[i].y;
  }

  double rangeX = maxX - minX, rangeY = maxY - minY;
  if (rangeX < 1e-15 || rangeY < 1e-15) return 0.0;

  int gridN = 1000;
  double dx = rangeX / gridN;
  double dy = rangeY / gridN;

  int inside = 0;
  for (int iy = 0; iy < gridN; iy++) {
    double py = minY + (iy + 0.5) * dy;
    for (int ix = 0; ix < gridN; ix++) {
      double px = minX + (ix + 0.5) * dx;
      int w = winding_number_at_point(px, py, polys->polys,
                                       polys->polySizes, polys->numPolys);
      if (w > 0) inside++;
    }
  }

  return inside * dx * dy;
}

// Compute area of simple (non-self-intersecting) polygons using shoelace formula.
// Returns the exact area for simple polygon contours (e.g. from manifold_slice).
double manifold_polygons2d_area(const ManifoldPolygons2D *polys) {
  if (!polys || polys->numPolys == 0) return 0.0;

  double totalArea = 0.0;
  int offset = 0;
  for (int p = 0; p < polys->numPolys; p++) {
    int n = polys->polySizes[p];
    if (n < 3) { offset += n; continue; }
    double area = 0.0;
    for (int i = 0; i < n; i++) {
      int j = (i + 1) % n;
      area += polys->polys[offset + i].x * polys->polys[offset + j].y;
      area -= polys->polys[offset + j].x * polys->polys[offset + i].y;
    }
    totalArea += area;
    offset += n;
  }

  return fabs(totalArea) * 0.5;
}

Manifold manifold_read_obj(const char *path) {
  FILE *f = fopen(path, "r");
  if (!f) return manifold_invalid();

  // Dynamic arrays for vertex coords and face indices
  size_t vertCap = 256, vertLen = 0;
  double *verts = (double *)malloc(vertCap * sizeof(double));
  size_t triCap = 256, triLen = 0;
  int *tris = (int *)malloc(triCap * sizeof(int));
  double tolerance = 0.0;
  double epsilon = -1.0;

  char line[1024];
  while (fgets(line, sizeof(line), f)) {
    if (line[0] == '#') {
      // Check for "# tolerance = X" or "# epsilon = X"
      double val;
      if (sscanf(line, "# tolerance = %lf", &val) == 1) {
        tolerance = val;
      } else if (sscanf(line, "# epsilon = %lf", &val) == 1) {
        epsilon = val;
      }
    } else if (line[0] == 'v' && line[1] == ' ') {
      double x, y, z;
      if (sscanf(line + 2, "%lf %lf %lf", &x, &y, &z) == 3) {
        if (vertLen + 3 > vertCap) {
          vertCap *= 2;
          verts = (double *)realloc(verts, vertCap * sizeof(double));
        }
        verts[vertLen++] = x;
        verts[vertLen++] = y;
        verts[vertLen++] = z;
      }
    } else if (line[0] == 'f' && line[1] == ' ') {
      int a, b, c;
      if (sscanf(line + 2, "%d %d %d", &a, &b, &c) == 3) {
        if (triLen + 3 > triCap) {
          triCap *= 2;
          tris = (int *)realloc(tris, triCap * sizeof(int));
        }
        tris[triLen++] = a - 1;  // OBJ is 1-based
        tris[triLen++] = b - 1;
        tris[triLen++] = c - 1;
      }
    }
  }
  fclose(f);

  size_t numVert = vertLen / 3;
  size_t numTri = triLen / 3;

  // Build ManifoldVec3 array
  ManifoldVec3 *vertPos = (ManifoldVec3 *)malloc(numVert * sizeof(ManifoldVec3));
  for (size_t i = 0; i < numVert; i++) {
    vertPos[i].x = verts[3 * i + 0];
    vertPos[i].y = verts[3 * i + 1];
    vertPos[i].z = verts[3 * i + 2];
  }
  ManifoldIVec3 *triVerts = (ManifoldIVec3 *)malloc(numTri * sizeof(ManifoldIVec3));
  for (size_t i = 0; i < numTri; i++) {
    triVerts[i].x = tris[3 * i + 0];
    triVerts[i].y = tris[3 * i + 1];
    triVerts[i].z = tris[3 * i + 2];
  }

  Manifold result = manifold_from_mesh(vertPos, numVert, triVerts, numTri);

  // Apply tolerance if specified in OBJ comment
  if (tolerance > 0.0) {
    result.impl.tolerance = tolerance;
  }
  // Apply epsilon if specified in OBJ comment
  if (epsilon >= 0.0) {
    result.impl.epsilon = epsilon;
  }

  free(verts);
  free(tris);
  free(vertPos);
  free(triVerts);
  return result;
}

// ---------- GLB Import ----------
// Minimal JSON helpers for glTF parsing

static const char *glb_skip_ws(const char *p) {
  while (*p == ' ' || *p == '\t' || *p == '\n' || *p == '\r') p++;
  return p;
}

// Skip a JSON value (string, number, object, array, bool, null)
static const char *glb_skip_json(const char *p) {
  p = glb_skip_ws(p);
  if (*p == '"') {
    p++;
    while (*p) {
      if (*p == '\\') { p += 2; continue; }
      if (*p == '"') return p + 1;
      p++;
    }
    return p;
  }
  if (*p == '{' || *p == '[') {
    int depth = 0;
    while (*p) {
      if (*p == '\\') { p += 2; continue; }
      if (*p == '"') {
        p++;
        while (*p) {
          if (*p == '\\') { p += 2; continue; }
          if (*p == '"') break;
          p++;
        }
        if (*p) p++;
        continue;
      }
      if (*p == '{' || *p == '[') depth++;
      else if (*p == '}' || *p == ']') {
        depth--;
        if (depth == 0) return p + 1;
      }
      p++;
    }
    return p;
  }
  while (*p && *p != ',' && *p != '}' && *p != ']' &&
         *p != ' ' && *p != '\t' && *p != '\n' && *p != '\r') p++;
  return p;
}

// Find a key's value in a JSON object. p must point to '{'
static const char *glb_json_get(const char *p, const char *key) {
  p = glb_skip_ws(p);
  if (*p != '{') return NULL;
  size_t keylen = strlen(key);
  p = glb_skip_ws(p + 1);
  while (*p && *p != '}') {
    if (*p != '"') return NULL;
    p++;
    const char *kstart = p;
    while (*p && *p != '"') p++;
    size_t klen = (size_t)(p - kstart);
    p++; // skip closing "
    p = glb_skip_ws(p);
    if (*p == ':') p = glb_skip_ws(p + 1);
    if (klen == keylen && memcmp(kstart, key, klen) == 0) return p;
    p = glb_skip_json(p);
    p = glb_skip_ws(p);
    if (*p == ',') p = glb_skip_ws(p + 1);
  }
  return NULL;
}

// Get the nth element of a JSON array. p must point to '['
static const char *glb_json_idx(const char *p, int n) {
  p = glb_skip_ws(p);
  if (*p != '[') return NULL;
  p = glb_skip_ws(p + 1);
  for (int i = 0; i < n; i++) {
    if (*p == ']') return NULL;
    p = glb_skip_json(p);
    p = glb_skip_ws(p);
    if (*p == ',') p = glb_skip_ws(p + 1);
  }
  return (*p == ']') ? NULL : p;
}

static int glb_json_int(const char *p) {
  if (!p) return 0;
  p = glb_skip_ws(p);
  return atoi(p);
}

Manifold manifold_read_glb(const char *path) {
  FILE *f = fopen(path, "rb");
  if (!f) return manifold_invalid();

  // Read entire file
  fseek(f, 0, SEEK_END);
  long fsize = ftell(f);
  fseek(f, 0, SEEK_SET);
  unsigned char *buf = (unsigned char *)malloc(fsize);
  if (!buf || fread(buf, 1, fsize, f) != (size_t)fsize) {
    free(buf);
    fclose(f);
    return manifold_invalid();
  }
  fclose(f);

  // Parse GLB header: magic(4) + version(4) + length(4)
  if (fsize < 12) { free(buf); return manifold_invalid(); }
  uint32_t magic, version, total;
  memcpy(&magic, buf, 4);
  memcpy(&version, buf + 4, 4);
  memcpy(&total, buf + 8, 4);
  if (magic != 0x46546C67 || version != 2) { // "glTF"
    free(buf);
    return manifold_invalid();
  }

  // Chunk 0: JSON
  if (fsize < 20) { free(buf); return manifold_invalid(); }
  uint32_t json_len, json_type;
  memcpy(&json_len, buf + 12, 4);
  memcpy(&json_type, buf + 16, 4);
  if (json_type != 0x4E4F534A) { free(buf); return manifold_invalid(); } // "JSON"
  char *json = (char *)malloc(json_len + 1);
  memcpy(json, buf + 20, json_len);
  json[json_len] = '\0';

  // Chunk 1: BIN
  size_t bin_off = 20 + json_len;
  if ((long)bin_off + 8 > fsize) { free(json); free(buf); return manifold_invalid(); }
  uint32_t bin_len, bin_type;
  memcpy(&bin_len, buf + bin_off, 4);
  memcpy(&bin_type, buf + bin_off + 4, 4);
  if (bin_type != 0x004E4942) { free(json); free(buf); return manifold_invalid(); } // "BIN\0"
  const unsigned char *bin = buf + bin_off + 8;

  // Parse JSON: meshes[0].primitives[0].attributes.POSITION and .indices
  const char *meshes = glb_json_get(json, "meshes");
  const char *mesh0 = meshes ? glb_json_idx(meshes, 0) : NULL;
  const char *prims = mesh0 ? glb_json_get(mesh0, "primitives") : NULL;
  const char *prim0 = prims ? glb_json_idx(prims, 0) : NULL;
  const char *attrs = prim0 ? glb_json_get(prim0, "attributes") : NULL;
  const char *pos_val = attrs ? glb_json_get(attrs, "POSITION") : NULL;
  const char *idx_val = prim0 ? glb_json_get(prim0, "indices") : NULL;

  if (!pos_val || !idx_val) { free(json); free(buf); return manifold_invalid(); }
  int pos_acc = glb_json_int(pos_val);
  int idx_acc = glb_json_int(idx_val);

  // Parse accessors
  const char *accessors = glb_json_get(json, "accessors");
  const char *bufferViews = glb_json_get(json, "bufferViews");
  if (!accessors || !bufferViews) { free(json); free(buf); return manifold_invalid(); }

  // Position accessor
  const char *pa = glb_json_idx(accessors, pos_acc);
  int pa_bv = glb_json_int(glb_json_get(pa, "bufferView"));
  int pa_count = glb_json_int(glb_json_get(pa, "count"));
  const char *pa_bo_val = pa ? glb_json_get(pa, "byteOffset") : NULL;
  int pa_byte_offset = pa_bo_val ? glb_json_int(pa_bo_val) : 0;

  // Index accessor
  const char *ia = glb_json_idx(accessors, idx_acc);
  int ia_bv = glb_json_int(glb_json_get(ia, "bufferView"));
  int ia_count = glb_json_int(glb_json_get(ia, "count"));
  int ia_comp_type = glb_json_int(glb_json_get(ia, "componentType"));
  const char *ia_bo_val = ia ? glb_json_get(ia, "byteOffset") : NULL;
  int ia_byte_offset = ia_bo_val ? glb_json_int(ia_bo_val) : 0;

  // Position buffer view
  const char *pbv = glb_json_idx(bufferViews, pa_bv);
  const char *pbv_bo_val = pbv ? glb_json_get(pbv, "byteOffset") : NULL;
  int pbv_offset = pbv_bo_val ? glb_json_int(pbv_bo_val) : 0;

  // Index buffer view
  const char *ibv = glb_json_idx(bufferViews, ia_bv);
  const char *ibv_bo_val = ibv ? glb_json_get(ibv, "byteOffset") : NULL;
  int ibv_offset = ibv_bo_val ? glb_json_int(ibv_bo_val) : 0;

  size_t numVert = (size_t)pa_count;
  size_t numIdx = (size_t)ia_count;
  size_t numTri = numIdx / 3;

  // Build MeshGL (float precision, matching C++ ImportMesh)
  ManifoldMeshGL mgl = manifold_meshgl_empty();
  mgl.numProp = 3;
  mgl.vertLen = numVert;
  mgl.triLen = numTri;
  mgl.vertProperties = (float *)malloc(numVert * 3 * sizeof(float));
  mgl.triVerts = (int *)malloc(numTri * 3 * sizeof(int));

  // Read positions (float VEC3) with (z,x,y) coordinate swap
  const float *pos_data = (const float *)(bin + pbv_offset + pa_byte_offset);
  for (size_t i = 0; i < numVert; i++) {
    float x = pos_data[3 * i + 0];
    float y = pos_data[3 * i + 1];
    float z = pos_data[3 * i + 2];
    mgl.vertProperties[3 * i + 0] = z;
    mgl.vertProperties[3 * i + 1] = x;
    mgl.vertProperties[3 * i + 2] = y;
  }

  // Read indices
  const unsigned char *idx_data = bin + ibv_offset + ia_byte_offset;
  if (ia_comp_type == 5125) { // UNSIGNED_INT
    const uint32_t *idx32 = (const uint32_t *)idx_data;
    for (size_t i = 0; i < numIdx; i++)
      mgl.triVerts[i] = (int)idx32[i];
  } else if (ia_comp_type == 5123) { // UNSIGNED_SHORT
    const uint16_t *idx16 = (const uint16_t *)idx_data;
    for (size_t i = 0; i < numIdx; i++)
      mgl.triVerts[i] = (int)idx16[i];
  } else { // unsupported
    free(mgl.vertProperties);
    free(mgl.triVerts);
    free(json);
    free(buf);
    return manifold_invalid();
  }

  free(json);
  free(buf);

  Manifold result = manifold_from_meshgl(&mgl);
  free(mgl.vertProperties);
  free(mgl.triVerts);
  return result;
}
