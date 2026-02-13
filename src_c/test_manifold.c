// C11 test suite for the Manifold port.
// Simple assert-based tests, no external framework needed.

#include <assert.h>
#include <math.h>
#include <stdio.h>

#include "manifold_api.h"

#define ASSERT_NEAR(a, b, tol) \
  do { \
    double _a = (a), _b = (b), _t = (tol); \
    if (fabs(_a - _b) > _t) { \
      fprintf(stderr, "FAIL: %s:%d: |%g - %g| = %g > %g\n", \
              __FILE__, __LINE__, _a, _b, fabs(_a - _b), _t); \
      assert(0); \
    } \
  } while (0)

#define ASSERT_TRUE(x) \
  do { \
    if (!(x)) { \
      fprintf(stderr, "FAIL: %s:%d: %s\n", __FILE__, __LINE__, #x); \
      assert(0); \
    } \
  } while (0)

#define ASSERT_EQ(a, b) \
  do { \
    if ((a) != (b)) { \
      fprintf(stderr, "FAIL: %s:%d: %s != %s\n", __FILE__, __LINE__, #a, #b); \
      assert(0); \
    } \
  } while (0)

#define RUN_TEST(name) \
  do { \
    printf("  %-40s ", #name); \
    test_##name(); \
    printf("PASS\n"); \
  } while (0)

// ============== Vec Math Tests ==============

static void test_vec3_basics(void) {
  ManifoldVec3 a = manifold_vec3(1, 2, 3);
  ManifoldVec3 b = manifold_vec3(4, 5, 6);
  ManifoldVec3 sum = vec3_add(a, b);
  ASSERT_NEAR(sum.x, 5.0, 1e-10);
  ASSERT_NEAR(sum.y, 7.0, 1e-10);
  ASSERT_NEAR(sum.z, 9.0, 1e-10);

  ManifoldVec3 diff = vec3_sub(b, a);
  ASSERT_NEAR(diff.x, 3.0, 1e-10);
  ASSERT_NEAR(diff.y, 3.0, 1e-10);
  ASSERT_NEAR(diff.z, 3.0, 1e-10);

  double d = vec3_dot(a, b);
  ASSERT_NEAR(d, 32.0, 1e-10);

  ManifoldVec3 cross = vec3_cross(a, b);
  ASSERT_NEAR(cross.x, -3.0, 1e-10);
  ASSERT_NEAR(cross.y, 6.0, 1e-10);
  ASSERT_NEAR(cross.z, -3.0, 1e-10);
}

static void test_vec3_normalize(void) {
  ManifoldVec3 v = manifold_vec3(3, 4, 0);
  ManifoldVec3 n = vec3_normalize(v);
  ASSERT_NEAR(n.x, 0.6, 1e-10);
  ASSERT_NEAR(n.y, 0.8, 1e-10);
  ASSERT_NEAR(n.z, 0.0, 1e-10);
  ASSERT_NEAR(vec3_length(n), 1.0, 1e-10);
}

static void test_mat3_identity(void) {
  ManifoldMat3 I = mat3_identity();
  ManifoldVec3 v = manifold_vec3(1, 2, 3);
  ManifoldVec3 result = mat3_mul_vec3(I, v);
  ASSERT_NEAR(result.x, 1.0, 1e-10);
  ASSERT_NEAR(result.y, 2.0, 1e-10);
  ASSERT_NEAR(result.z, 3.0, 1e-10);
}

static void test_mat3_inverse(void) {
  ManifoldMat3 m;
  m.cols[0] = manifold_vec3(1, 0, 0);
  m.cols[1] = manifold_vec3(0, 2, 0);
  m.cols[2] = manifold_vec3(0, 0, 4);
  ManifoldMat3 inv = mat3_inverse(m);
  ASSERT_NEAR(inv.cols[0].x, 1.0, 1e-10);
  ASSERT_NEAR(inv.cols[1].y, 0.5, 1e-10);
  ASSERT_NEAR(inv.cols[2].z, 0.25, 1e-10);
}

static void test_box_operations(void) {
  ManifoldBox b = manifold_box_empty();
  manifold_box_union_point(&b, manifold_vec3(1, 2, 3));
  manifold_box_union_point(&b, manifold_vec3(-1, -2, -3));
  ASSERT_NEAR(b.min.x, -1.0, 1e-10);
  ASSERT_NEAR(b.max.z, 3.0, 1e-10);
  ASSERT_TRUE(manifold_box_contains_point(b, manifold_vec3(0, 0, 0)));
  ASSERT_TRUE(!manifold_box_contains_point(b, manifold_vec3(10, 0, 0)));
}

// ============== Vec (Dynamic Array) Tests ==============

static void test_vec_int(void) {
  ManifoldVecInt v = {0};
  vec_int_push(&v, 10);
  vec_int_push(&v, 20);
  vec_int_push(&v, 30);
  ASSERT_EQ(v.len, (size_t)3);
  ASSERT_EQ(v.data[0], 10);
  ASSERT_EQ(v.data[1], 20);
  ASSERT_EQ(v.data[2], 30);
  vec_int_pop(&v);
  ASSERT_EQ(v.len, (size_t)2);
  vec_int_free(&v);
}

// ============== Hashtable Tests ==============

static void test_hashtable(void) {
  ManifoldHashTable ht = manifold_hashtable_create(16, 1);
  manifold_hashtable_insert(&ht, 42, 100);
  manifold_hashtable_insert(&ht, 99, 200);
  ASSERT_EQ(manifold_hashtable_lookup(&ht, 42), 100);
  ASSERT_EQ(manifold_hashtable_lookup(&ht, 99), 200);
  ASSERT_EQ(manifold_hashtable_lookup(&ht, 1), -1);
  manifold_hashtable_free(&ht);
}

// ============== Tetrahedron Tests ==============

static void test_tetrahedron(void) {
  Manifold m = manifold_tetrahedron();
  ASSERT_EQ(manifold_status(&m), MANIFOLD_ERROR_NO_ERROR);
  ASSERT_TRUE(!manifold_is_empty(&m));
  ASSERT_EQ(manifold_num_vert(&m), (size_t)4);
  ASSERT_EQ(manifold_num_tri(&m), (size_t)4);
  ASSERT_EQ(manifold_num_edge(&m), (size_t)6);

  // Volume of tetrahedron with vertices at (±1,±1,±1) (alternating)
  double vol = manifold_volume(&m);
  ASSERT_NEAR(fabs(vol), 8.0 / 3.0, 1e-10);

  ManifoldBox bb = manifold_bounding_box(&m);
  ASSERT_TRUE(manifold_box_isfinite(bb));

  manifold_destroy(&m);
}

// ============== Cube Tests ==============

static void test_cube(void) {
  ManifoldVec3 size = manifold_vec3(1, 1, 1);
  Manifold m = manifold_cube(size, false);
  ASSERT_EQ(manifold_status(&m), MANIFOLD_ERROR_NO_ERROR);
  ASSERT_TRUE(!manifold_is_empty(&m));
  ASSERT_EQ(manifold_num_vert(&m), (size_t)8);
  ASSERT_EQ(manifold_num_tri(&m), (size_t)12);

  double vol = manifold_volume(&m);
  ASSERT_NEAR(vol, 1.0, 1e-10);

  double sa = manifold_surface_area(&m);
  ASSERT_NEAR(sa, 6.0, 1e-10);

  manifold_destroy(&m);
}

static void test_cube_centered(void) {
  ManifoldVec3 size = manifold_vec3(2, 2, 2);
  Manifold m = manifold_cube(size, true);
  ASSERT_EQ(manifold_status(&m), MANIFOLD_ERROR_NO_ERROR);

  double vol = manifold_volume(&m);
  ASSERT_NEAR(vol, 8.0, 1e-10);

  ManifoldBox bb = manifold_bounding_box(&m);
  ASSERT_NEAR(bb.min.x, -1.0, 1e-10);
  ASSERT_NEAR(bb.max.x, 1.0, 1e-10);

  manifold_destroy(&m);
}

static void test_cube_scaled(void) {
  ManifoldVec3 size = manifold_vec3(2, 3, 4);
  Manifold m = manifold_cube(size, false);

  double vol = manifold_volume(&m);
  ASSERT_NEAR(vol, 24.0, 1e-10);

  double sa = manifold_surface_area(&m);
  ASSERT_NEAR(sa, 2.0 * (2*3 + 3*4 + 2*4), 1e-10);

  manifold_destroy(&m);
}

// ============== Transform Tests ==============

static void test_translate(void) {
  Manifold m = manifold_cube(manifold_vec3(1, 1, 1), false);
  Manifold translated = manifold_translate(&m, manifold_vec3(5, 0, 0));

  ManifoldBox bb = manifold_bounding_box(&translated);
  ASSERT_NEAR(bb.min.x, 5.0, 1e-10);
  ASSERT_NEAR(bb.max.x, 6.0, 1e-10);

  // Volume should be preserved
  ASSERT_NEAR(manifold_volume(&translated), 1.0, 1e-10);

  manifold_destroy(&m);
  manifold_destroy(&translated);
}

static void test_scale(void) {
  Manifold m = manifold_cube(manifold_vec3(1, 1, 1), false);
  Manifold scaled = manifold_scale(&m, manifold_vec3(2, 3, 4));

  double vol = manifold_volume(&scaled);
  ASSERT_NEAR(vol, 24.0, 1e-10);

  manifold_destroy(&m);
  manifold_destroy(&scaled);
}

// ============== Quality Tests ==============

static void test_quality(void) {
  manifold_quality_reset();
  int seg = manifold_get_circular_segments(1.0);
  ASSERT_TRUE(seg >= 3);
  ASSERT_TRUE(seg % 4 == 0);

  manifold_set_circular_segments(32);
  ASSERT_EQ(manifold_get_circular_segments(1.0), 32);

  manifold_quality_reset();
}

// ============== Disjoint Sets Tests ==============

static void test_disjoint_sets(void) {
  ManifoldDisjointSets ds = manifold_disjoint_sets_create(10);
  manifold_disjoint_sets_unite(&ds, 0, 1);
  manifold_disjoint_sets_unite(&ds, 2, 3);
  manifold_disjoint_sets_unite(&ds, 0, 3);
  ASSERT_TRUE(manifold_disjoint_sets_same(&ds, 0, 1));
  ASSERT_TRUE(manifold_disjoint_sets_same(&ds, 0, 2));
  ASSERT_TRUE(manifold_disjoint_sets_same(&ds, 0, 3));
  ASSERT_TRUE(!manifold_disjoint_sets_same(&ds, 0, 4));
  manifold_disjoint_sets_free(&ds);
}

// ============== Manifold checks ==============

static void test_manifold_check(void) {
  Manifold m = manifold_tetrahedron();
  ASSERT_TRUE(manifold_impl_is_manifold(&m.impl));
  manifold_destroy(&m);

  Manifold c = manifold_cube(manifold_vec3(1,1,1), false);
  ASSERT_TRUE(manifold_impl_is_manifold(&c.impl));
  manifold_destroy(&c);
}

// ============== Sphere Tests ==============

static void test_sphere(void) {
  Manifold m = manifold_sphere(1.0, 0);
  ASSERT_EQ(manifold_status(&m), MANIFOLD_ERROR_NO_ERROR);
  ASSERT_TRUE(!manifold_is_empty(&m));
  ASSERT_TRUE(manifold_num_vert(&m) >= 6);
  ASSERT_TRUE(manifold_num_tri(&m) >= 8);

  ManifoldBox bb = manifold_bounding_box(&m);
  ASSERT_NEAR(bb.min.x, -1.0, 0.01);
  ASSERT_NEAR(bb.max.x, 1.0, 0.01);

  manifold_destroy(&m);
}

// ============== Cylinder Tests ==============

static void test_cylinder(void) {
  manifold_set_circular_segments(12);
  Manifold m = manifold_cylinder(2.0, 1.0, -1.0, 0, false);
  ASSERT_EQ(manifold_status(&m), MANIFOLD_ERROR_NO_ERROR);
  ASSERT_TRUE(!manifold_is_empty(&m));
  ASSERT_TRUE(manifold_num_vert(&m) >= 24);

  ManifoldBox bb = manifold_bounding_box(&m);
  ASSERT_NEAR(bb.max.z, 2.0, 0.01);
  ASSERT_NEAR(bb.min.z, 0.0, 0.01);

  manifold_destroy(&m);
  manifold_quality_reset();
}

static void test_cylinder_centered(void) {
  manifold_set_circular_segments(8);
  Manifold m = manifold_cylinder(4.0, 1.0, -1.0, 0, true);

  ManifoldBox bb = manifold_bounding_box(&m);
  ASSERT_NEAR(bb.min.z, -2.0, 0.01);
  ASSERT_NEAR(bb.max.z, 2.0, 0.01);

  manifold_destroy(&m);
  manifold_quality_reset();
}

// ============== SDF Level Set Tests ==============

static double sdf_sphere(double x, double y, double z, void *ctx) {
  double r = *(double *)ctx;
  return sqrt(x*x + y*y + z*z) - r;
}

static void test_sdf_sphere(void) {
  double radius = 1.0;
  ManifoldBox bounds = manifold_box(manifold_vec3(-1.5, -1.5, -1.5),
                                     manifold_vec3(1.5, 1.5, 1.5));
  Manifold m = manifold_level_set(sdf_sphere, &radius, bounds, 0.3, 0.0, -1.0);
  ASSERT_EQ(manifold_status(&m), MANIFOLD_ERROR_NO_ERROR);
  ASSERT_TRUE(!manifold_is_empty(&m));
  ASSERT_TRUE(manifold_num_tri(&m) > 0);

  // Volume of sphere = 4/3 * pi * r^3 ≈ 4.189
  // Marching cubes gives approximate result
  double vol = manifold_volume(&m);
  ASSERT_TRUE(vol > 3.0);   // should be close to 4.189
  ASSERT_TRUE(vol < 5.5);

  manifold_destroy(&m);
}

static double sdf_box(double x, double y, double z, void *ctx) {
  (void)ctx;
  // SDF of unit cube centered at origin
  double dx = fabs(x) - 0.5;
  double dy = fabs(y) - 0.5;
  double dz = fabs(z) - 0.5;
  double mx = fmax(dx, fmax(dy, dz));
  if (mx < 0) return mx;
  double ox = fmax(dx, 0), oy = fmax(dy, 0), oz = fmax(dz, 0);
  return sqrt(ox*ox + oy*oy + oz*oz);
}

static void test_sdf_box(void) {
  ManifoldBox bounds = manifold_box(manifold_vec3(-1, -1, -1),
                                     manifold_vec3(1, 1, 1));
  Manifold m = manifold_level_set(sdf_box, NULL, bounds, 0.2, 0.0, -1.0);
  ASSERT_EQ(manifold_status(&m), MANIFOLD_ERROR_NO_ERROR);
  ASSERT_TRUE(!manifold_is_empty(&m));

  manifold_destroy(&m);
}

// ============== Copy Test ==============

static void test_copy(void) {
  Manifold a = manifold_cube(manifold_vec3(1, 1, 1), false);
  Manifold b;
  manifold_copy(&b, &a);

  ASSERT_EQ(manifold_num_vert(&b), manifold_num_vert(&a));
  ASSERT_EQ(manifold_num_tri(&b), manifold_num_tri(&a));
  ASSERT_NEAR(manifold_volume(&b), manifold_volume(&a), 1e-10);

  manifold_destroy(&a);
  manifold_destroy(&b);
}

// ============== Rotation Test ==============

static void test_rotate(void) {
  Manifold m = manifold_cube(manifold_vec3(1, 1, 1), true);
  Manifold r = manifold_rotate(&m, 0, 0, 90);

  // After 90° rotation around Z, bbox should be the same
  ASSERT_NEAR(manifold_volume(&r), 1.0, 1e-10);
  ManifoldBox bb = manifold_bounding_box(&r);
  ASSERT_NEAR(bb.min.x, -0.5, 0.01);
  ASSERT_NEAR(bb.max.x, 0.5, 0.01);

  manifold_destroy(&m);
  manifold_destroy(&r);
}

// ============== Boolean Tests ==============

static void test_boolean_union_non_overlapping(void) {
  // Two cubes that don't overlap
  Manifold a = manifold_cube(manifold_vec3(1, 1, 1), false);
  Manifold b_base = manifold_cube(manifold_vec3(1, 1, 1), false);
  Manifold b = manifold_translate(&b_base, manifold_vec3(3, 0, 0));

  Manifold u = manifold_union(&a, &b);

  ASSERT_EQ(manifold_status(&u), MANIFOLD_ERROR_NO_ERROR);
  ASSERT_TRUE(!manifold_is_empty(&u));
  // Non-overlapping union should have combined vertices and tris
  ASSERT_EQ(manifold_num_vert(&u), manifold_num_vert(&a) + manifold_num_vert(&b));
  ASSERT_EQ(manifold_num_tri(&u), manifold_num_tri(&a) + manifold_num_tri(&b));
  // Volume should be sum of both cubes
  ASSERT_NEAR(manifold_volume(&u), 2.0, 0.01);

  manifold_destroy(&a);
  manifold_destroy(&b_base);
  manifold_destroy(&b);
  manifold_destroy(&u);
}

static void test_boolean_empty(void) {
  Manifold a = manifold_cube(manifold_vec3(1, 1, 1), false);
  Manifold empty;
  manifold_create(&empty);

  // Union with empty
  Manifold u = manifold_union(&a, &empty);
  ASSERT_EQ(manifold_num_vert(&u), manifold_num_vert(&a));
  ASSERT_NEAR(manifold_volume(&u), 1.0, 0.01);

  // Intersect with empty
  Manifold inter = manifold_intersection(&a, &empty);
  ASSERT_TRUE(manifold_is_empty(&inter));

  manifold_destroy(&a);
  manifold_destroy(&empty);
  manifold_destroy(&u);
  manifold_destroy(&inter);
}

static void test_boolean_subtract_non_overlapping(void) {
  // Two cubes that don't overlap - subtraction should return the first
  Manifold a = manifold_cube(manifold_vec3(1, 1, 1), false);
  Manifold b_base = manifold_cube(manifold_vec3(1, 1, 1), false);
  Manifold b = manifold_translate(&b_base, manifold_vec3(5, 0, 0));

  Manifold d = manifold_difference(&a, &b);

  ASSERT_EQ(manifold_status(&d), MANIFOLD_ERROR_NO_ERROR);
  ASSERT_TRUE(!manifold_is_empty(&d));
  ASSERT_EQ(manifold_num_vert(&d), manifold_num_vert(&a));
  ASSERT_NEAR(manifold_volume(&d), 1.0, 0.01);

  manifold_destroy(&a);
  manifold_destroy(&b_base);
  manifold_destroy(&b);
  manifold_destroy(&d);
}

static void test_boolean_intersect_non_overlapping(void) {
  // Two cubes that don't overlap - intersection should be empty
  Manifold a = manifold_cube(manifold_vec3(1, 1, 1), false);
  Manifold b_base = manifold_cube(manifold_vec3(1, 1, 1), false);
  Manifold b = manifold_translate(&b_base, manifold_vec3(5, 0, 0));

  Manifold inter = manifold_intersection(&a, &b);

  ASSERT_EQ(manifold_status(&inter), MANIFOLD_ERROR_NO_ERROR);
  ASSERT_TRUE(manifold_is_empty(&inter));

  manifold_destroy(&a);
  manifold_destroy(&b_base);
  manifold_destroy(&b);
  manifold_destroy(&inter);
}

// ============== Hull Tests ==============

static void test_hull_cube(void) {
  // Hull of a cube's vertices should give back the same cube
  Manifold cube = manifold_cube(manifold_vec3(1, 1, 1), false);
  Manifold hull = manifold_hull(&cube);

  ASSERT_EQ(manifold_status(&hull), MANIFOLD_ERROR_NO_ERROR);
  ASSERT_TRUE(!manifold_is_empty(&hull));
  ASSERT_EQ(manifold_num_vert(&hull), (size_t)8);
  // Hull of cube has 12 triangles (6 faces * 2 tris)
  ASSERT_EQ(manifold_num_tri(&hull), (size_t)12);
  ASSERT_NEAR(manifold_volume(&hull), 1.0, 0.01);

  manifold_destroy(&cube);
  manifold_destroy(&hull);
}

static void test_hull_points(void) {
  // Create a set of random points including some interior ones
  ManifoldVec3 pts[] = {
    {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1},
    {1, 1, 0}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1},
    {0.5, 0.5, 0.5},  // interior point
    {0.3, 0.3, 0.3},  // interior point
  };
  Manifold hull = manifold_hull_points(pts, 10);

  ASSERT_EQ(manifold_status(&hull), MANIFOLD_ERROR_NO_ERROR);
  ASSERT_TRUE(!manifold_is_empty(&hull));
  // Should have 8 vertices (the cube corners), not the interior points
  ASSERT_EQ(manifold_num_vert(&hull), (size_t)8);
  ASSERT_NEAR(manifold_volume(&hull), 1.0, 0.01);

  manifold_destroy(&hull);
}

// ============== Mesh Access Tests ==============

static void warp_double_x(double *x, double *y, double *z, void *ctx) {
  (void)y; (void)z; (void)ctx;
  *x *= 2.0;
}

static void test_warp(void) {
  Manifold cube = manifold_cube(manifold_vec3(1, 1, 1), true);

  Manifold warped = manifold_warp(&cube, warp_double_x, NULL);
  ASSERT_NEAR(manifold_volume(&warped), 2.0, 0.01);  // volume should double

  ManifoldBox bb = manifold_bounding_box(&warped);
  ASSERT_NEAR(bb.min.x, -1.0, 0.01);  // was -0.5, now -1.0
  ASSERT_NEAR(bb.max.x, 1.0, 0.01);

  manifold_destroy(&cube);
  manifold_destroy(&warped);
}

static void test_mesh_access(void) {
  Manifold m = manifold_cube(manifold_vec3(1, 1, 1), false);

  size_t count;
  const ManifoldVec3 *verts = manifold_get_vert_positions(&m, &count);
  ASSERT_EQ(count, (size_t)8);
  ASSERT_TRUE(verts != NULL);

  size_t numTri;
  manifold_get_triangles(&m, NULL, &numTri);
  ASSERT_EQ(numTri, (size_t)12);

  ManifoldIVec3 *tris = (ManifoldIVec3 *)malloc(numTri * sizeof(ManifoldIVec3));
  manifold_get_triangles(&m, tris, NULL);
  // All vertex indices should be valid
  for (size_t i = 0; i < numTri; i++) {
    ASSERT_TRUE(tris[i].x >= 0 && tris[i].x < (int)count);
    ASSERT_TRUE(tris[i].y >= 0 && tris[i].y < (int)count);
    ASSERT_TRUE(tris[i].z >= 0 && tris[i].z < (int)count);
  }
  free(tris);

  float *vertProps;
  int *triVerts;
  size_t nv, np, nt;
  manifold_get_mesh(&m, &vertProps, &nv, &np, &triVerts, &nt);
  ASSERT_EQ(nv, (size_t)8);
  ASSERT_EQ(np, (size_t)3);
  ASSERT_EQ(nt, (size_t)12);
  manifold_free_mesh(vertProps, triVerts);

  manifold_destroy(&m);
}

// ============== Stress Tests ==============

static void test_multiple_operations(void) {
  // Test a chain of operations: create, transform, copy, measure
  Manifold base = manifold_cube(manifold_vec3(2, 2, 2), true);
  Manifold moved = manifold_translate(&base, manifold_vec3(1, 0, 0));
  Manifold scaled = manifold_scale(&moved, manifold_vec3(1, 2, 1));
  Manifold rotated = manifold_rotate(&scaled, 0, 0, 45);

  // Volume should be 2*2*2 * 1*2*1 = 16
  ASSERT_NEAR(manifold_volume(&scaled), 16.0, 0.01);

  // Rotated volume should be the same
  ASSERT_NEAR(manifold_volume(&rotated), 16.0, 0.01);

  manifold_destroy(&base);
  manifold_destroy(&moved);
  manifold_destroy(&scaled);
  manifold_destroy(&rotated);
}

static void test_hull_tetrahedron(void) {
  // Hull of 4 non-coplanar points should be a tetrahedron
  ManifoldVec3 pts[] = {
    {0, 0, 0}, {1, 0, 0}, {0.5, 1, 0}, {0.5, 0.5, 1}
  };
  Manifold hull = manifold_hull_points(pts, 4);

  ASSERT_EQ(manifold_status(&hull), MANIFOLD_ERROR_NO_ERROR);
  ASSERT_TRUE(!manifold_is_empty(&hull));
  ASSERT_EQ(manifold_num_vert(&hull), (size_t)4);
  ASSERT_EQ(manifold_num_tri(&hull), (size_t)4);

  // Volume of tetrahedron with these vertices
  double vol = manifold_volume(&hull);
  ASSERT_TRUE(fabs(vol) > 0.01);

  manifold_destroy(&hull);
}

static void test_sdf_volume_accuracy(void) {
  // Test SDF sphere with finer resolution for better accuracy
  double radius = 1.0;
  ManifoldBox bounds = manifold_box(manifold_vec3(-1.5, -1.5, -1.5),
                                     manifold_vec3(1.5, 1.5, 1.5));
  Manifold m = manifold_level_set(sdf_sphere, &radius, bounds, 0.1, 0.0, -1.0);

  double vol = manifold_volume(&m);
  double expected = 4.0 / 3.0 * MANIFOLD_PI; // ~4.189
  // With edge length 0.1, should be within 5% of exact
  ASSERT_TRUE(fabs(vol - expected) / expected < 0.05);

  manifold_destroy(&m);
}

// ============== Overlapping Boolean Tests ==============

static void test_boolean_union_overlapping(void) {
  // Two overlapping cubes - union should produce a single solid
  Manifold a = manifold_cube(manifold_vec3(2.0, 2.0, 2.0), false);
  Manifold b_base = manifold_cube(manifold_vec3(2.0, 2.0, 2.0), false);
  Manifold b = manifold_translate(&b_base, manifold_vec3(1.0, 0.0, 0.0));

  Manifold u = manifold_boolean(&a, &b, MANIFOLD_OP_ADD);

  ASSERT_EQ(manifold_status(&u), MANIFOLD_ERROR_NO_ERROR);
  ASSERT_TRUE(!manifold_is_empty(&u));

  // Volume of union should be 2*2*2 + 2*2*2 - 1*2*2 = 12
  double vol = manifold_volume(&u);
  ASSERT_NEAR(vol, 12.0, 0.5);

  manifold_destroy(&a);
  manifold_destroy(&b_base);
  manifold_destroy(&b);
  manifold_destroy(&u);
}

static void test_boolean_subtract_overlapping(void) {
  // Subtract a smaller offset cube from a larger one
  Manifold a = manifold_cube(manifold_vec3(2.0, 2.0, 2.0), false);
  Manifold b_base = manifold_cube(manifold_vec3(2.0, 2.0, 2.0), false);
  Manifold b = manifold_translate(&b_base, manifold_vec3(1.0, 0.0, 0.0));

  Manifold d = manifold_boolean(&a, &b, MANIFOLD_OP_SUBTRACT);

  ASSERT_EQ(manifold_status(&d), MANIFOLD_ERROR_NO_ERROR);
  ASSERT_TRUE(!manifold_is_empty(&d));

  // Volume = 2*2*2 - 1*2*2 = 4
  double vol = manifold_volume(&d);
  ASSERT_NEAR(vol, 4.0, 0.5);

  manifold_destroy(&a);
  manifold_destroy(&b_base);
  manifold_destroy(&b);
  manifold_destroy(&d);
}

static void test_boolean_intersect_overlapping(void) {
  // Intersection of two offset cubes
  Manifold a = manifold_cube(manifold_vec3(2.0, 2.0, 2.0), false);
  Manifold b_base = manifold_cube(manifold_vec3(2.0, 2.0, 2.0), false);
  Manifold b = manifold_translate(&b_base, manifold_vec3(1.0, 0.0, 0.0));

  Manifold inter = manifold_boolean(&a, &b, MANIFOLD_OP_INTERSECT);

  ASSERT_EQ(manifold_status(&inter), MANIFOLD_ERROR_NO_ERROR);
  ASSERT_TRUE(!manifold_is_empty(&inter));

  // Volume = 1*2*2 = 4
  double vol = manifold_volume(&inter);
  ASSERT_NEAR(vol, 4.0, 0.5);

  manifold_destroy(&a);
  manifold_destroy(&b_base);
  manifold_destroy(&b);
  manifold_destroy(&inter);
}

static void test_boolean_3d_offset(void) {
  // Cubes offset in all 3 axes - tests more complex intersections
  Manifold a = manifold_cube(manifold_vec3(2.0, 2.0, 2.0), false);
  Manifold b_base = manifold_cube(manifold_vec3(2.0, 2.0, 2.0), false);
  Manifold b = manifold_translate(&b_base, manifold_vec3(1.0, 1.0, 1.0));

  Manifold u = manifold_boolean(&a, &b, MANIFOLD_OP_ADD);
  ASSERT_EQ(manifold_status(&u), MANIFOLD_ERROR_NO_ERROR);
  // Union volume = 2*8 - 1*1*1 = 15
  ASSERT_NEAR(manifold_volume(&u), 15.0, 0.5);

  Manifold inter = manifold_boolean(&a, &b, MANIFOLD_OP_INTERSECT);
  ASSERT_EQ(manifold_status(&inter), MANIFOLD_ERROR_NO_ERROR);
  // Intersection volume = 1*1*1 = 1
  ASSERT_NEAR(manifold_volume(&inter), 1.0, 0.5);

  manifold_destroy(&a);
  manifold_destroy(&b_base);
  manifold_destroy(&b);
  manifold_destroy(&u);
  manifold_destroy(&inter);
}

static void test_boolean_self_difference(void) {
  // A - A should be empty
  Manifold a = manifold_cube(manifold_vec3(2.0, 2.0, 2.0), false);
  Manifold b = manifold_cube(manifold_vec3(2.0, 2.0, 2.0), false);

  Manifold d = manifold_boolean(&a, &b, MANIFOLD_OP_SUBTRACT);
  ASSERT_EQ(manifold_status(&d), MANIFOLD_ERROR_NO_ERROR);
  // Should be empty or near-zero volume
  double vol = manifold_volume(&d);
  ASSERT_NEAR(vol, 0.0, 0.5);

  manifold_destroy(&a);
  manifold_destroy(&b);
  manifold_destroy(&d);
}

static void test_boolean_sequential(void) {
  // (A union B) - C, all cubes
  Manifold a = manifold_cube(manifold_vec3(2.0, 2.0, 2.0), false);
  Manifold b_base = manifold_cube(manifold_vec3(2.0, 2.0, 2.0), false);
  Manifold b = manifold_translate(&b_base, manifold_vec3(1.0, 0.0, 0.0));
  Manifold c_base = manifold_cube(manifold_vec3(1.0, 1.0, 1.0), false);
  Manifold c = manifold_translate(&c_base, manifold_vec3(0.5, 0.5, 0.5));

  Manifold u = manifold_boolean(&a, &b, MANIFOLD_OP_ADD);
  Manifold result = manifold_boolean(&u, &c, MANIFOLD_OP_SUBTRACT);

  ASSERT_EQ(manifold_status(&result), MANIFOLD_ERROR_NO_ERROR);
  // Volume = 12 - 1 = 11
  ASSERT_NEAR(manifold_volume(&result), 11.0, 0.5);

  manifold_destroy(&a);
  manifold_destroy(&b_base);
  manifold_destroy(&b);
  manifold_destroy(&c_base);
  manifold_destroy(&c);
  manifold_destroy(&u);
  manifold_destroy(&result);
}

static void test_boolean_sphere(void) {
  // Boolean on spheres
  Manifold a = manifold_sphere(1.0, 16);
  Manifold b_base = manifold_sphere(1.0, 16);
  Manifold b = manifold_translate(&b_base, manifold_vec3(1.0, 0.0, 0.0));

  Manifold u = manifold_boolean(&a, &b, MANIFOLD_OP_ADD);
  ASSERT_EQ(manifold_status(&u), MANIFOLD_ERROR_NO_ERROR);
  ASSERT_TRUE(!manifold_is_empty(&u));

  // Union of two low-poly spheres: should be non-empty and non-zero volume
  double vol = manifold_volume(&u);
  ASSERT_TRUE(vol > 2.0);  // should be noticeably larger than zero
  ASSERT_TRUE(vol < 10.0); // but bounded

  manifold_destroy(&a);
  manifold_destroy(&b_base);
  manifold_destroy(&b);
  manifold_destroy(&u);
}

static void test_boolean_face_union(void) {
  // Two cubes sharing a face: A + A.translate(1,0,0)
  // Should produce a single mesh with volume 2
  Manifold a = manifold_cube(manifold_vec3(1.0, 1.0, 1.0), false);
  Manifold b_base = manifold_cube(manifold_vec3(1.0, 1.0, 1.0), false);
  Manifold b = manifold_translate(&b_base, manifold_vec3(1.0, 0.0, 0.0));

  Manifold u = manifold_boolean(&a, &b, MANIFOLD_OP_ADD);
  ASSERT_EQ(manifold_status(&u), MANIFOLD_ERROR_NO_ERROR);
  ASSERT_NEAR(manifold_volume(&u), 2.0, 0.1);
  ASSERT_NEAR(manifold_surface_area(&u), 10.0, 1.0);

  manifold_destroy(&a);
  manifold_destroy(&b_base);
  manifold_destroy(&b);
  manifold_destroy(&u);
}

static void test_boolean_corner_union(void) {
  // Two cubes sharing only a corner: A + A.translate(1,1,1)
  // Should be two separate meshes (non-overlapping)
  Manifold a = manifold_cube(manifold_vec3(1.0, 1.0, 1.0), false);
  Manifold b_base = manifold_cube(manifold_vec3(1.0, 1.0, 1.0), false);
  Manifold b = manifold_translate(&b_base, manifold_vec3(1.0, 1.0, 1.0));

  Manifold u = manifold_boolean(&a, &b, MANIFOLD_OP_ADD);
  ASSERT_EQ(manifold_status(&u), MANIFOLD_ERROR_NO_ERROR);
  ASSERT_NEAR(manifold_volume(&u), 2.0, 0.1);

  manifold_destroy(&a);
  manifold_destroy(&b_base);
  manifold_destroy(&b);
  manifold_destroy(&u);
}

// NOTE: multi_coplanar test is disabled - crashes on 2nd sequential subtract
// with coplanar faces. This requires more robust halfedge handling in boolean
// result when the input mesh comes from a previous boolean operation.

static void test_boolean_multi_coplanar(void) {
  // Sequential subtracts with coplanar faces (previously crashed)
  // Matches C++ test: volume=0.18, surfaceArea=2.76
  Manifold c1 = manifold_cube(manifold_vec3(1.0, 1.0, 1.0), false);
  Manifold c2_base = manifold_cube(manifold_vec3(1.0, 1.0, 1.0), false);
  Manifold c2 = manifold_translate(&c2_base, manifold_vec3(0.3, 0.3, 0.0));
  Manifold first = manifold_boolean(&c1, &c2, MANIFOLD_OP_SUBTRACT);
  ASSERT_TRUE(manifold_num_vert(&first) > 0);

  Manifold c3_base = manifold_cube(manifold_vec3(1.0, 1.0, 1.0), false);
  Manifold c3 = manifold_translate(&c3_base, manifold_vec3(-0.3, -0.3, 0.0));
  Manifold result = manifold_boolean(&first, &c3, MANIFOLD_OP_SUBTRACT);
  ASSERT_TRUE(manifold_num_vert(&result) > 0);
  ASSERT_NEAR(manifold_volume(&result), 0.18, 0.05);
  ASSERT_NEAR(manifold_surface_area(&result), 2.76, 1.0);

  manifold_destroy(&c1);
  manifold_destroy(&c2_base);
  manifold_destroy(&c2);
  manifold_destroy(&first);
  manifold_destroy(&c3_base);
  manifold_destroy(&c3);
  manifold_destroy(&result);
}

static void test_extrude_square(void) {
  // Extrude a unit square to height 2
  ManifoldVec2 square[] = {{0,0}, {1,0}, {1,1}, {0,1}};
  int sizes[] = {4};
  ManifoldVec2 scaleTop = {1.0, 1.0};
  Manifold m = manifold_extrude(square, sizes, 1, 2.0, 0, 0.0, scaleTop);
  ASSERT_EQ(manifold_status(&m), MANIFOLD_ERROR_NO_ERROR);
  ASSERT_NEAR(manifold_volume(&m), 2.0, 0.01);
  ASSERT_TRUE(manifold_num_vert(&m) >= 8);
  manifold_destroy(&m);
}

static void test_extrude_triangle(void) {
  // Extrude a triangle to height 1
  ManifoldVec2 tri[] = {{0,0}, {1,0}, {0.5, 0.866}};
  int sizes[] = {3};
  ManifoldVec2 scaleTop = {1.0, 1.0};
  Manifold m = manifold_extrude(tri, sizes, 1, 1.0, 0, 0.0, scaleTop);
  ASSERT_EQ(manifold_status(&m), MANIFOLD_ERROR_NO_ERROR);
  double area = 0.5 * 1.0 * 0.866; // triangle area
  ASSERT_NEAR(manifold_volume(&m), area, 0.05);
  manifold_destroy(&m);
}

static void test_extrude_cone(void) {
  // Extrude a square to a cone (scale 0,0 at top)
  ManifoldVec2 square[] = {{0,0}, {1,0}, {1,1}, {0,1}};
  int sizes[] = {4};
  ManifoldVec2 scaleTop = {0.0, 0.0};
  Manifold m = manifold_extrude(square, sizes, 1, 3.0, 0, 0.0, scaleTop);
  ASSERT_EQ(manifold_status(&m), MANIFOLD_ERROR_NO_ERROR);
  // Volume of pyramid: base_area * height / 3 = 1 * 3 / 3 = 1
  ASSERT_NEAR(manifold_volume(&m), 1.0, 0.1);
  manifold_destroy(&m);
}

static void test_revolve_circle(void) {
  // Revolve a small rectangle around Y-axis to create a torus-like shape
  // Rectangle at x=[1,2], y=[0,0.5]
  ManifoldVec2 rect[] = {{1.0, 0.0}, {2.0, 0.0}, {2.0, 0.5}, {1.0, 0.5}};
  int sizes[] = {4};
  Manifold m = manifold_revolve(rect, sizes, 1, 8, 360.0);
  ASSERT_EQ(manifold_status(&m), MANIFOLD_ERROR_NO_ERROR);
  ASSERT_TRUE(manifold_num_vert(&m) > 0);
  // Volume of revolution: 2*pi*r_avg * cross_section_area
  // r_avg = 1.5, area = 1*0.5 = 0.5, so V = 2*pi*1.5*0.5 ≈ 4.71
  double vol = manifold_volume(&m);
  ASSERT_NEAR(vol, 4.71, 1.0);
  manifold_destroy(&m);
}

static void test_boolean_tetra(void) {
  // Simplest boolean: subtract translated tetrahedra (from C++ test suite)
  Manifold tetra = manifold_tetrahedron();
  Manifold tetra2_base = manifold_tetrahedron();
  Manifold tetra2 = manifold_translate(&tetra2_base, manifold_vec3(0.5, 0.5, 0.5));
  Manifold result = manifold_boolean(&tetra2, &tetra, MANIFOLD_OP_SUBTRACT);
  ASSERT_EQ(manifold_status(&result), MANIFOLD_ERROR_NO_ERROR);
  ASSERT_TRUE(manifold_num_vert(&result) >= 4);
  ASSERT_TRUE(manifold_volume(&result) > 0.0);
  manifold_destroy(&tetra);
  manifold_destroy(&tetra2_base);
  manifold_destroy(&tetra2);
  manifold_destroy(&result);
}

static void test_boolean_mirrored(void) {
  // Cube minus itself (same position) should be empty or near-zero
  Manifold cube = manifold_cube(manifold_vec3(2.0, 2.0, 2.0), true);
  Manifold cube2 = manifold_cube(manifold_vec3(2.0, 2.0, 2.0), true);
  Manifold result = manifold_boolean(&cube, &cube2, MANIFOLD_OP_SUBTRACT);
  ASSERT_EQ(manifold_status(&result), MANIFOLD_ERROR_NO_ERROR);
  // Self-subtract should give zero volume
  ASSERT_NEAR(manifold_volume(&result), 0.0, 0.01);
  manifold_destroy(&cube);
  manifold_destroy(&cube2);
  manifold_destroy(&result);
}

static void test_boolean_union_difference(void) {
  // UnionDifference test: volume conservation
  Manifold cube = manifold_cube(manifold_vec3(1.0, 1.0, 1.0), false);
  Manifold cube2_base = manifold_cube(manifold_vec3(1.0, 1.0, 1.0), false);
  Manifold cube2 = manifold_translate(&cube2_base, manifold_vec3(0.5, 0.5, 0.5));
  
  Manifold u = manifold_boolean(&cube, &cube2, MANIFOLD_OP_ADD);
  Manifold d = manifold_boolean(&cube, &cube2, MANIFOLD_OP_SUBTRACT);
  Manifold i = manifold_boolean(&cube, &cube2, MANIFOLD_OP_INTERSECT);
  
  double vol_a = manifold_volume(&cube);
  double vol_b = manifold_volume(&cube2);
  double vol_u = manifold_volume(&u);
  double vol_d = manifold_volume(&d);
  double vol_i = manifold_volume(&i);
  
  // Volume(A ∪ B) = Volume(A) + Volume(B) - Volume(A ∩ B)
  ASSERT_NEAR(vol_u, vol_a + vol_b - vol_i, 0.2);
  // Volume(A - B) = Volume(A) - Volume(A ∩ B)
  ASSERT_NEAR(vol_d, vol_a - vol_i, 0.2);
  
  manifold_destroy(&cube);
  manifold_destroy(&cube2_base);
  manifold_destroy(&cube2);
  manifold_destroy(&u);
  manifold_destroy(&d);
  manifold_destroy(&i);
}

static void test_boolean_split(void) {
  // Split: cube intersection + difference = whole cube volume
  Manifold cube = manifold_cube(manifold_vec3(2.0, 2.0, 2.0), true);
  Manifold sphere = manifold_sphere(1.0, 4);
  Manifold sphere_t = manifold_translate(&sphere, manifold_vec3(0.0, 0.0, 1.0));
  
  Manifold inside = manifold_boolean(&cube, &sphere_t, MANIFOLD_OP_INTERSECT);
  Manifold outside = manifold_boolean(&cube, &sphere_t, MANIFOLD_OP_SUBTRACT);
  
  double vol_cube = manifold_volume(&cube);
  double vol_in = manifold_volume(&inside);
  double vol_out = manifold_volume(&outside);
  
  ASSERT_NEAR(vol_in + vol_out, vol_cube, 1.0);
  
  manifold_destroy(&cube);
  manifold_destroy(&sphere);
  manifold_destroy(&sphere_t);
  manifold_destroy(&inside);
  manifold_destroy(&outside);
}

static void test_boolean_cylinder_subtract(void) {
  // Subtract a smaller cylinder from a larger one (pipe-like shape)
  Manifold outer = manifold_cylinder(2.0, 1.0, 1.0, 16, false);
  Manifold inner_base = manifold_cylinder(2.0, 0.5, 0.5, 16, false);
  Manifold inner = manifold_translate(&inner_base, manifold_vec3(0.0, 0.0, -0.1));
  
  Manifold result = manifold_boolean(&outer, &inner, MANIFOLD_OP_SUBTRACT);
  ASSERT_EQ(manifold_status(&result), MANIFOLD_ERROR_NO_ERROR);
  
  // Volume should be pi*(R^2 - r^2)*h ≈ pi*(1-0.25)*2 ≈ 4.71
  double vol = manifold_volume(&result);
  ASSERT_TRUE(vol > 3.0 && vol < 6.0);
  
  manifold_destroy(&outer);
  manifold_destroy(&inner_base);
  manifold_destroy(&inner);
  manifold_destroy(&result);
}

static void test_boolean_extrude_subtract(void) {
  // Extrude a square, then subtract another extruded square
  ManifoldVec2 sq1[] = {{0,0}, {2,0}, {2,2}, {0,2}};
  int sizes[] = {4};
  ManifoldVec2 scale1 = {1.0, 1.0};
  Manifold a = manifold_extrude(sq1, sizes, 1, 2.0, 0, 0.0, scale1);
  
  ManifoldVec2 sq2[] = {{0.5, 0.5}, {1.5, 0.5}, {1.5, 1.5}, {0.5, 1.5}};
  Manifold b = manifold_extrude(sq2, sizes, 1, 3.0, 0, 0.0, scale1);
  
  Manifold result = manifold_boolean(&a, &b, MANIFOLD_OP_SUBTRACT);
  ASSERT_EQ(manifold_status(&result), MANIFOLD_ERROR_NO_ERROR);
  
  // Volume: 2*2*2 - 1*1*2 = 8 - 2 = 6
  ASSERT_NEAR(manifold_volume(&result), 6.0, 0.5);
  
  manifold_destroy(&a);
  manifold_destroy(&b);
  manifold_destroy(&result);
}

// ============== Main ==============

int main(void) {
  printf("=== Manifold C11 Port Tests ===\n\n");
  printf("Vec Math:\n");
  RUN_TEST(vec3_basics);
  RUN_TEST(vec3_normalize);
  RUN_TEST(mat3_identity);
  RUN_TEST(mat3_inverse);
  RUN_TEST(box_operations);

  printf("\nDynamic Arrays:\n");
  RUN_TEST(vec_int);

  printf("\nHash Table:\n");
  RUN_TEST(hashtable);

  printf("\nDisjoint Sets:\n");
  RUN_TEST(disjoint_sets);

  printf("\nConstructors:\n");
  RUN_TEST(tetrahedron);
  RUN_TEST(cube);
  RUN_TEST(cube_centered);
  RUN_TEST(cube_scaled);

  printf("\nTransforms:\n");
  RUN_TEST(translate);
  RUN_TEST(scale);

  printf("\nQuality:\n");
  RUN_TEST(quality);

  printf("\nManifold Checks:\n");
  RUN_TEST(manifold_check);

  printf("\nSphere & Cylinder:\n");
  RUN_TEST(sphere);
  RUN_TEST(cylinder);
  RUN_TEST(cylinder_centered);

  printf("\nSDF Level Set:\n");
  RUN_TEST(sdf_sphere);
  RUN_TEST(sdf_box);

  printf("\nMisc:\n");
  RUN_TEST(copy);
  RUN_TEST(rotate);

  printf("\nBoolean Operations:\n");
  RUN_TEST(boolean_union_non_overlapping);
  RUN_TEST(boolean_empty);
  RUN_TEST(boolean_subtract_non_overlapping);
  RUN_TEST(boolean_intersect_non_overlapping);
  RUN_TEST(boolean_union_overlapping);
  RUN_TEST(boolean_subtract_overlapping);
  RUN_TEST(boolean_intersect_overlapping);
  RUN_TEST(boolean_3d_offset);
  RUN_TEST(boolean_self_difference);
  RUN_TEST(boolean_sequential);
  RUN_TEST(boolean_sphere);
  RUN_TEST(boolean_face_union);
  RUN_TEST(boolean_corner_union);
  RUN_TEST(boolean_multi_coplanar);

  printf("\nExtrude:\n");
  RUN_TEST(extrude_square);
  RUN_TEST(extrude_triangle);
  RUN_TEST(extrude_cone);
  RUN_TEST(revolve_circle);

  printf("\nAdvanced Boolean:\n");
  RUN_TEST(boolean_tetra);
  RUN_TEST(boolean_mirrored);
  RUN_TEST(boolean_union_difference);
  RUN_TEST(boolean_split);
  RUN_TEST(boolean_cylinder_subtract);
  RUN_TEST(boolean_extrude_subtract);

  printf("\nConvex Hull:\n");
  RUN_TEST(hull_cube);
  RUN_TEST(hull_points);

  printf("\nMesh Access:\n");
  RUN_TEST(warp);
  RUN_TEST(mesh_access);

  printf("\nStress Tests:\n");
  RUN_TEST(multiple_operations);
  RUN_TEST(hull_tetrahedron);
  RUN_TEST(sdf_volume_accuracy);

  printf("\n=== All %d tests passed! ===\n", 54);
  return 0;
}
