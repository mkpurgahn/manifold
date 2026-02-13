// C11 test suite for the Manifold port.
// Simple assert-based tests, no external framework needed.

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

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

#define ASSERT_GE(a, b) \
  do { \
    if ((a) < (b)) { \
      fprintf(stderr, "FAIL: %s:%d: %s < %s\n", __FILE__, __LINE__, #a, #b); \
      assert(0); \
    } \
  } while (0)

#define ASSERT_LE(a, b) \
  do { \
    if ((a) > (b)) { \
      fprintf(stderr, "FAIL: %s:%d: %s > %s\n", __FILE__, __LINE__, #a, #b); \
      assert(0); \
    } \
  } while (0)

#define ASSERT_LT(a, b) \
  do { \
    if ((a) >= (b)) { \
      fprintf(stderr, "FAIL: %s:%d: %s >= %s\n", __FILE__, __LINE__, #a, #b); \
      assert(0); \
    } \
  } while (0)

#define ASSERT_FALSE(x) ASSERT_TRUE(!(x))

#define RUN_TEST(name) \
  do { \
    printf("  %-40s ", #name); \
    test_##name(); \
    printf("PASS\n"); \
  } while (0)


// ============== Tests will be ported 1:1 from C++ ==============
// Each test file in test/ gets its tests converted here.
// Port order: boolean_test, boolean_complex_test, manifold_test,
// smooth_test, hull_test, sdf_test, properties_test, samples_test,
// polygon_test, manifoldc_test

int main(void) {
  int passed = 0, failed = 0;
  printf("\n=== Manifold C Port Test Suite ===\n\n");

  // Tests will be added here as they are ported from C++

  printf("\n=== %d tests passed ===\n", passed);
  return failed;
}
