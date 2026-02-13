// Copyright 2021 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// Public API implementation for the C11 Manifold port.

#include "manifold_api.h"
#include "manifold_boolean.h"
#include "manifold_hull.h"
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

// Create from mesh
Manifold manifold_from_mesh(const ManifoldVec3 *vertPos, size_t numVert,
                            const ManifoldIVec3 *triVerts, size_t numTri) {
  Manifold m;
  manifold_impl_init(&m.impl);
  if (numVert == 0 || numTri == 0) return m;

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
  // Preserve inherited epsilon; bbox may grow so use max of inherited and bbox-based
  manifold_impl_set_epsilon(&out.impl, out.impl.epsilon, false);
  return out;
}

Manifold manifold_scale(const Manifold *m, ManifoldVec3 v) {
  Manifold out;
  manifold_copy(&out, m);
  for (size_t i = 0; i < out.impl.vertPos.len; i++) {
    out.impl.vertPos.data[i] = vec3_mul(out.impl.vertPos.data[i], v);
  }
  // Negative determinant means we need to flip triangle winding
  if (v.x * v.y * v.z < 0) {
    size_t numTri = out.impl.halfedge.len / 3;
    for (size_t tri = 0; tri < numTri; tri++) {
      // Swap halfedges 0 and 2 within the triangle
      ManifoldHalfedge tmp = out.impl.halfedge.data[3 * tri];
      out.impl.halfedge.data[3 * tri] = out.impl.halfedge.data[3 * tri + 2];
      out.impl.halfedge.data[3 * tri + 2] = tmp;
      // Swap startVert/endVert and fix pairedHalfedge indices
      for (int i = 0; i < 3; i++) {
        ManifoldHalfedge *he = &out.impl.halfedge.data[3 * tri + i];
        int sv = he->startVert;
        he->startVert = he->endVert;
        he->endVert = sv;
        // FlipHalfedge: remap index within its triangle
        if (he->pairedHalfedge >= 0) {
          int pTri = he->pairedHalfedge / 3;
          int pVert = 2 - (he->pairedHalfedge - 3 * pTri);
          he->pairedHalfedge = 3 * pTri + pVert;
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
  // Scale epsilon by the spectral norm of the 3x3 scale matrix (matches C++)
  ManifoldMat3 scaleMat;
  scaleMat.cols[0] = manifold_vec3(v.x, 0, 0);
  scaleMat.cols[1] = manifold_vec3(0, v.y, 0);
  scaleMat.cols[2] = manifold_vec3(0, 0, v.z);
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

  Manifold out;
  manifold_copy(&out, m);
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

  // Extract 3x3 rotation/scale part
  ManifoldMat3 m3;
  m3.cols[0] = manifold_vec3(t.cols[0].x, t.cols[0].y, t.cols[0].z);
  m3.cols[1] = manifold_vec3(t.cols[1].x, t.cols[1].y, t.cols[1].z);
  m3.cols[2] = manifold_vec3(t.cols[2].x, t.cols[2].y, t.cols[2].z);

  // Normal transform = inverse(transpose(m3))
  ManifoldMat3 normalXf = mat3_inverse(mat3_transpose(m3));

  for (size_t i = 0; i < out.impl.faceNormal.len; i++) {
    out.impl.faceNormal.data[i] = vec3_normalize(
        mat3_mul_vec3(normalXf, out.impl.faceNormal.data[i]));
  }
  for (size_t i = 0; i < out.impl.vertNormal.len; i++) {
    out.impl.vertNormal.data[i] = vec3_normalize(
        mat3_mul_vec3(normalXf, out.impl.vertNormal.data[i]));
  }

  // Flip triangle winding if negative determinant
  double det = mat3_det(m3);
  if (det < 0) {
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
          int pTri = he->pairedHalfedge / 3;
          int pVert = 2 - (he->pairedHalfedge - 3 * pTri);
          he->pairedHalfedge = 3 * pTri + pVert;
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

  double rh = radiusHigh >= 0.0 ? radiusHigh : radiusLow;
  double maxR = fmax(radiusLow, rh);
  int n = circularSegments > 2 ? circularSegments :
          manifold_get_circular_segments(maxR);
  if (n < 3) n = 3;

  // Build cylinder from scratch: n verts on bottom, n on top, 2n triangles for sides, n-2 each for caps
  manifold_impl_init(&m.impl);

  // Bottom circle
  for (int i = 0; i < n; i++) {
    double angle = 360.0 * i / n;
    ManifoldVec3 p = manifold_vec3(
        radiusLow * manifold_cosd(angle),
        radiusLow * manifold_sind(angle),
        0.0);
    vec_vec3_push(&m.impl.vertPos, p);
  }
  // Top circle
  for (int i = 0; i < n; i++) {
    double angle = 360.0 * i / n;
    ManifoldVec3 p = manifold_vec3(
        rh * manifold_cosd(angle),
        rh * manifold_sind(angle),
        height);
    vec_vec3_push(&m.impl.vertPos, p);
  }

  ManifoldVecIVec3 triVerts = {0};
  // Side triangles
  for (int i = 0; i < n; i++) {
    int i2 = (i + 1) % n;
    vec_ivec3_push(&triVerts, manifold_ivec3(i, i2, n + i2));
    vec_ivec3_push(&triVerts, manifold_ivec3(i, n + i2, n + i));
  }
  // Bottom cap (fan from vertex 0)
  for (int i = 1; i < n - 1; i++) {
    vec_ivec3_push(&triVerts, manifold_ivec3(0, i + 1, i));
  }
  // Top cap (fan from vertex n)
  for (int i = 1; i < n - 1; i++) {
    vec_ivec3_push(&triVerts, manifold_ivec3(n, n + i, n + i + 1));
  }

  ManifoldVecIVec3 emptyTriVert = {0};
  manifold_impl_create_halfedges(&m.impl, &triVerts, &emptyTriVert);
  manifold_impl_initialize_original(&m.impl);
  manifold_impl_calculate_bbox(&m.impl);
  manifold_impl_set_epsilon(&m.impl, -1.0, false);
  manifold_impl_sort_geometry(&m.impl);
  manifold_impl_set_normals_and_coplanar(&m.impl);

  if (center) {
    ManifoldVec3 offset = manifold_vec3(0, 0, -height / 2.0);
    for (size_t i = 0; i < m.impl.vertPos.len; i++) {
      m.impl.vertPos.data[i] = vec3_add(m.impl.vertPos.data[i], offset);
    }
    m.impl.bBox = manifold_box_shift(m.impl.bBox, offset);
  }

  vec_ivec3_free(&triVerts);
  vec_ivec3_free(&emptyTriVert);
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
  manifold_impl_set_normals_and_coplanar(&out.impl);
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
  size_t np = 3; // just positions for now

  if (numVert) *numVert = nv;
  if (numProp) *numProp = np;
  if (numTri) *numTri = nt;

  if (vertProps) {
    *vertProps = (float *)malloc(nv * np * sizeof(float));
    for (size_t i = 0; i < nv; i++) {
      (*vertProps)[i * np + 0] = (float)m->impl.vertPos.data[i].x;
      (*vertProps)[i * np + 1] = (float)m->impl.vertPos.data[i].y;
      (*vertProps)[i * np + 2] = (float)m->impl.vertPos.data[i].z;
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
  // Apply operations sequentially
  Manifold result;
  manifold_copy(&result, &manifolds[0]);
  for (int i = 1; i < count; i++) {
    Manifold next = manifold_boolean(&result, &manifolds[i], op);
    manifold_destroy(&result);
    result = next;
  }
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
  Manifold out;
  manifold_copy(&out, m);
  if (!manifold_is_empty(&out)) {
    // Add faceID mapping for UpdateSharpenedEdges
    size_t numTri = manifold_impl_num_tri(&out.impl);
    for (size_t i = 0; i < numTri; i++) {
      out.impl.meshRelation.triRef.data[i].faceID = (int)i;
    }
    // Re-create from mesh to get proper sorting
    // Then create tangents with mapped sharpened edges
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
    for (size_t i = 0; i < numTri; i++) {
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
    // Ensure the smaller mesh is the one we iterate (bImpl)
    const ManifoldImpl *bigImpl = aImpl;
    const ManifoldImpl *smallImpl = bImpl;
    if (aImpl->vertPos.len < bImpl->vertPos.len) {
      bigImpl = bImpl;
      smallImpl = aImpl;
    }
    size_t numSmall = smallImpl->vertPos.len;
    size_t numBig = bigImpl->vertPos.len;
    if (numSmall * numBig <= 200) {
      size_t total = numSmall * numBig;
      ManifoldVec3 *pts = (ManifoldVec3 *)malloc(total * sizeof(ManifoldVec3));
      size_t idx = 0;
      for (size_t i = 0; i < numSmall; i++)
        for (size_t j = 0; j < numBig; j++)
          pts[idx++] = vec3_add(smallImpl->vertPos.data[i], bigImpl->vertPos.data[j]);
      Manifold hullResult = manifold_hull_points(pts, total);
      free(pts);
      manifold_destroy(&base);
      Manifold orig = manifold_as_original(&hullResult);
      manifold_destroy(&hullResult);
      return orig;
    }
    // Large point sets: per-face hulls + batch union
    size_t numSmallTri = manifold_impl_num_tri(smallImpl);
    Manifold *hulls = (Manifold *)malloc(numSmallTri * sizeof(Manifold));
    for (size_t tri = 0; tri < numSmallTri; tri++) {
      size_t numPts = 3 * numBig;
      ManifoldVec3 *pts = (ManifoldVec3 *)malloc(numPts * sizeof(ManifoldVec3));
      size_t pidx = 0;
      for (int i = 0; i < 3; i++) {
        ManifoldVec3 vert = smallImpl->vertPos.data[
            smallImpl->halfedge.data[tri * 3 + i].startVert];
        for (size_t j = 0; j < numBig; j++)
          pts[pidx++] = vec3_add(vert, bigImpl->vertPos.data[j]);
      }
      hulls[tri] = manifold_hull_points(pts, numPts);
      free(pts);
    }
    Manifold hullResult = manifold_batch_boolean(hulls, (int)numSmallTri, MANIFOLD_OP_ADD);
    for (size_t i = 0; i < numSmallTri; i++) manifold_destroy(&hulls[i]);
    free(hulls);
    manifold_destroy(&base);
    Manifold orig = manifold_as_original(&hullResult);
    manifold_destroy(&hullResult);
    return orig;
  }

  // Convex-Convex Minkowski Difference (inset): intersect translated copies of A
  if (inset && aConvex && bConvex) {
    // A ⊖ B = intersection of (A translated by -b_j) for all vertices b_j of B
    size_t numBVerts = bImpl->vertPos.len;
    Manifold result;
    manifold_copy(&result, (aImpl == &a->impl) ? a : b);
    // Translate A by -b_0
    ManifoldVec3 b0 = bImpl->vertPos.data[0];
    Manifold translated = manifold_translate(&result, manifold_vec3(-b0.x, -b0.y, -b0.z));
    manifold_destroy(&result);
    result = translated;
    for (size_t j = 1; j < numBVerts; j++) {
      ManifoldVec3 bj = bImpl->vertPos.data[j];
      Manifold aCopy;
      manifold_copy(&aCopy, (aImpl == &a->impl) ? a : b);
      Manifold tr = manifold_translate(&aCopy, manifold_vec3(-bj.x, -bj.y, -bj.z));
      manifold_destroy(&aCopy);
      Manifold inter = manifold_intersection(&result, &tr);
      manifold_destroy(&result);
      manifold_destroy(&tr);
      result = inter;
    }
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

Manifold manifold_from_meshgl(const ManifoldMeshGL *mesh) {
  if (!mesh || mesh->vertLen == 0 || mesh->triLen == 0 || mesh->numProp < 3) {
    Manifold m;
    manifold_impl_init(&m.impl);
    return m;
  }
  // Check for NaN vertices
  for (size_t i = 0; i < mesh->vertLen; i++) {
    for (int p = 0; p < 3; p++) {
      if (isnan(mesh->vertProperties[i * mesh->numProp + p]) ||
          isinf(mesh->vertProperties[i * mesh->numProp + p])) {
        Manifold m;
        manifold_impl_init(&m.impl);
        return m;
      }
    }
  }
  // Extract positions
  ManifoldVec3 *verts = (ManifoldVec3 *)malloc(mesh->vertLen * sizeof(ManifoldVec3));
  for (size_t i = 0; i < mesh->vertLen; i++) {
    verts[i].x = mesh->vertProperties[i * mesh->numProp + 0];
    verts[i].y = mesh->vertProperties[i * mesh->numProp + 1];
    verts[i].z = mesh->vertProperties[i * mesh->numProp + 2];
  }
  // Extract triangles
  ManifoldIVec3 *tris = (ManifoldIVec3 *)malloc(mesh->triLen * sizeof(ManifoldIVec3));
  for (size_t i = 0; i < mesh->triLen; i++) {
    tris[i].x = mesh->triVerts[i * 3 + 0];
    tris[i].y = mesh->triVerts[i * 3 + 1];
    tris[i].z = mesh->triVerts[i * 3 + 2];
  }
  Manifold result = manifold_from_mesh(verts, mesh->vertLen, tris, mesh->triLen);
  if (mesh->tolerance > 0) {
    Manifold tol = manifold_set_tolerance(&result, mesh->tolerance);
    manifold_destroy(&result);
    result = tol;
  }
  free(verts);
  free(tris);
  return result;
}
