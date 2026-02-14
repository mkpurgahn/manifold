# Manifold C Port — Test Completion Prompt

> **Status:** 161/161 ported tests **ALL PASS** (0 failures). 58 C++ tests remain unported. The core port is functionally complete — this prompt focuses on porting the remaining tests and fixing any failures they reveal.

> **⚠️ NEVER WEAKEN A TEST.** Do not increase tolerances, loosen assertions, or skip checks. The C++ tests define correct behavior. If the C port doesn't match, the C implementation is wrong.

> **⚠️ NEVER DISABLE OR COMMENT OUT A FAILING TEST.**

> **⚠️ DO NOT SPAWN SUB-AGENTS.** Do all work yourself, directly.

---

## Build & Run

```bash
cd /Users/admin/projects/manifold
make            # builds test_manifold
TEST=Boolean_Perturb1 ./test_manifold   # run single test
./test_manifold                          # run all tests (slow tests take 2-5 min each)
SKIP_SLOW=1 ./test_manifold              # skip tests marked slow
```

To run a test with timeout (macOS has no `timeout`, use perl):
```bash
TEST=Hull_Sphere perl -e 'alarm 240; exec @ARGV' ./test_manifold 2>&1
```

Note: slow tests (hull, sponge, scallop, minkowski) need **240+ seconds**. Don't assume a test hangs until 300s.

---

## Unported Tests — Port These (58 remaining)

### Priority 1: Boolean & BooleanComplex (15 tests)
These exercise the core boolean engine which is already ported. Most should just need test scaffolding.

| Test | C++ File | Notes |
|------|----------|-------|
| `Boolean_MeshGLRoundTrip` | `test/boolean_test.cpp` | MeshGL serialization round-trip |
| `Boolean_MixedProperties` | `test/boolean_test.cpp` | Property handling across booleans |
| `BooleanComplex_CraycloudBool` | `test/boolean_complex_test.cpp` | External mesh boolean |
| `BooleanComplex_GenericTwinBooleanTest7081` | `test/boolean_complex_test.cpp` | Regression test |
| `BooleanComplex_GenericTwinBooleanTest7863` | `test/boolean_complex_test.cpp` | Regression test |
| `BooleanComplex_Havocglass8Bool` | `test/boolean_complex_test.cpp` | External mesh boolean |
| `BooleanComplex_HullMask` | `test/boolean_complex_test.cpp` | Hull with mask |
| `BooleanComplex_InterpolatedNormals` | `test/boolean_complex_test.cpp` | Normal interpolation |
| `BooleanComplex_MeshRelation` | `test/boolean_complex_test.cpp` | Mesh relation tracking |
| `BooleanComplex_OffsetSelfIntersect` | `test/boolean_complex_test.cpp` | Offset self-intersection |
| `BooleanComplex_OffsetTriangulationFailure` | `test/boolean_complex_test.cpp` | Edge case |
| `BooleanComplex_Ring` | `test/boolean_complex_test.cpp` | Ring topology |
| `BooleanComplex_SimpleOffset` | `test/boolean_complex_test.cpp` | Simple offset |
| `BooleanComplex_Sphere` | `test/boolean_complex_test.cpp` | Sphere boolean |
| `BooleanComplex_Sweep` | `test/boolean_complex_test.cpp` | Sweep operation |

### Priority 2: Manifold Core (14 tests)
Core manifold operations — merge, slice, project, mesh relations.

| Test | C++ File | Notes |
|------|----------|-------|
| `Manifold_DecomposeProps` | `test/manifold_test.cpp` | Decompose with properties |
| `Manifold_FaceIDRoundTrip` | `test/manifold_test.cpp` | Face ID preservation |
| `Manifold_GetMeshGL` | `test/manifold_test.cpp` | MeshGL extraction |
| `Manifold_InvalidInput5` | `test/manifold_test.cpp` | Invalid input handling |
| `Manifold_InvalidInput7` | `test/manifold_test.cpp` | Invalid input handling |
| `Manifold_Merge` | `test/manifold_test.cpp` | Vertex merging |
| `Manifold_MergeEmpty` | `test/manifold_test.cpp` | Merge empty manifold |
| `Manifold_MergeRefine` | `test/manifold_test.cpp` | Merge + refine |
| `Manifold_MeshID` | `test/manifold_test.cpp` | Mesh ID tracking |
| `Manifold_MeshRelation` | `test/manifold_test.cpp` | Mesh relation |
| `Manifold_MeshRelationRefinePrecision` | `test/manifold_test.cpp` | Precision of mesh relation after refine |
| `Manifold_OpenscadCrash` | `test/manifold_test.cpp` | Regression crash fix |
| `Manifold_Project` | `test/manifold_test.cpp` | 2D projection |
| `Manifold_Slice` | `test/manifold_test.cpp` | Cross-section slicing |
| `Manifold_SliceEmptyObject` | `test/manifold_test.cpp` | Slice of empty |
| `Manifold_ValidInputOneRunIndex` | `test/manifold_test.cpp` | Edge case |
| `Manifold_WarpBatch` | `test/manifold_test.cpp` | Batch warp |
| `ManifoldFuzz_SimpleCube` | `test/manifold_test.cpp` | Fuzz regression |

### Priority 3: CrossSection (15 tests)
2D cross-section operations. May require `manifold_cross_section.c` if not yet ported.

| Test | C++ File | Notes |
|------|----------|-------|
| `CrossSection_BatchBoolean` | `test/cross_section_test.cpp` | 2D batch boolean |
| `CrossSection_BevelOffset` | `test/cross_section_test.cpp` | Bevel offset |
| `CrossSection_Decompose` | `test/cross_section_test.cpp` | Decompose |
| `CrossSection_Empty` | `test/cross_section_test.cpp` | Empty input |
| `CrossSection_FillRule` | `test/cross_section_test.cpp` | Fill rule |
| `CrossSection_Hull` | `test/cross_section_test.cpp` | 2D hull |
| `CrossSection_HullError` | `test/cross_section_test.cpp` | Hull error handling |
| `CrossSection_MirrorCheckAxis` | `test/cross_section_test.cpp` | Mirror axis check |
| `CrossSection_MirrorUnion` | `test/cross_section_test.cpp` | Mirror + union |
| `CrossSection_NegativeOffset` | `test/cross_section_test.cpp` | Negative offset |
| `CrossSection_Rect` | `test/cross_section_test.cpp` | Rectangle |
| `CrossSection_RoundOffset` | `test/cross_section_test.cpp` | Round offset |
| `CrossSection_Square` | `test/cross_section_test.cpp` | Square |
| `CrossSection_Transform` | `test/cross_section_test.cpp` | 2D transforms |
| `CrossSection_Warp` | `test/cross_section_test.cpp` | 2D warp |

### Priority 4: Smooth, Samples, Properties, Polygon (10 tests)

| Test | C++ File | Notes |
|------|----------|-------|
| `Smooth_Manual` | `test/smooth_test.cpp` | Manual smooth |
| `Smooth_Torus` | `test/smooth_test.cpp` | Torus smooth |
| `Samples_Bracelet` | `test/samples_test.cpp` | Complex sample |
| `Samples_CondensedMatter16` | `test/samples_test.cpp` | Condensed matter sim |
| `Samples_CondensedMatter64` | `test/samples_test.cpp` | Larger condensed matter |
| `Samples_GyroidModule` | `test/samples_test.cpp` | Gyroid surface |
| `Properties_Coplanar` | `test/properties_test.cpp` | Coplanar properties |
| `Properties_MingapStretchyBracelet` | `test/properties_test.cpp` | Complex property test |
| `PolygonFuzz_TriangulationNoCrash` | `test/polygon_test.cpp` | Fuzz — no crash |
| `PolygonFuzz_TriangulationNoCrashRounded` | `test/polygon_test.cpp` | Fuzz — rounded |

---

## Approach

1. **Port tests in priority order** — Priority 1 (Boolean) first since the implementation already exists.
2. **For each test:**
   a. Read the C++ test to understand inputs, operations, and expected values.
   b. Write the C equivalent in `src_c/test_manifold.c` using existing `manifold_*` API.
   c. Add `RUN_TEST(TestName);` to the appropriate section in `main()`.
   d. Run single test: `TEST=TestName ./test_manifold`
   e. If it fails, debug the C implementation (not the test).
   f. Run full suite to check regressions: `SKIP_SLOW=1 ./test_manifold`
   g. Commit: `test: port [TestName] from C++`
3. **Some tests may need new API functions** — if a C++ method isn't yet exposed in `manifold_api.h`, port it.
4. **CrossSection tests (Priority 3)** may need `manifold_cross_section.c` ported if it doesn't exist yet. Check before starting.

---

## C++ Reference Files

| C Port File | C++ Original | Purpose |
|-------------|-------------|---------|
| `src_c/manifold_boolean.c` | `src/boolean3.cpp`, `src/boolean_result.cpp` | Boolean operations |
| `src_c/manifold_smooth.c` | `src/smoothing.cpp` | Smooth subdivision |
| `src_c/manifold_impl.c` | `src/impl.cpp` | Core manifold operations |
| `src_c/manifold_csg_tree.c` | `src/csg_tree.cpp` | CSG tree evaluation |
| `src_c/manifold_hull.c` | `src/quickhull.cpp` | Convex hull |
| `src_c/manifold_properties.c` | `src/properties.cpp` | Mesh properties |
| `src_c/manifold_sdf.c` | `src/sdf.cpp` | Signed distance field |
| `src_c/manifold_api.c` | `bindings/c/manifoldc.cpp`, `src/minkowski.cpp` | C API + Minkowski |

C++ tests: `test/boolean_test.cpp`, `test/boolean_complex_test.cpp`, `test/smooth_test.cpp`, `test/hull_test.cpp`, `test/properties_test.cpp`, `test/samples_test.cpp`, `test/sdf_test.cpp`, `test/manifold_test.cpp`, `test/manifoldc_test.cpp`, `test/cross_section_test.cpp`, `test/polygon_test.cpp`.

---

## Rules

1. Port one test at a time. Do not scatter changes across multiple tests.
2. Always run the full test suite after porting to check for regressions.
3. Never weaken assertions. The C++ values are the ground truth.
4. Compare C vs C++ implementation line by line when a ported test fails.
5. Commit after each test or small batch: `test: port [TestName] from C++`
6. If a test needs a new C API function, port the function first, then the test.

---

## Flaky Test Investigation: `Boolean_NonConvexConvexMinkowskiDifference`

This test has been observed to **pass in isolation but fail in the full suite** — a classic sign of state leaking between tests or a race condition.

**Observed failure (in full suite):**
```
area: |16.818 - 16.70| = 0.117984 > 0.01
genus: 0 != 5
```

**Passes when run alone:** `TEST=Boolean_NonConvexConvexMinkowskiDifference ./test_manifold` → PASS

**Investigation steps:**
1. Run the full suite 5 times. Record pass/fail for this test each run.
2. If it fails, note which test ran immediately before it (test execution order matters).
3. Check for **global/static state** in the Minkowski code path:
   - `src_c/manifold_api.c` — `manifold_minkowski_impl()` uses `manifold_batch_boolean()` internally
   - `src_c/manifold_boolean.c` — any `static` variables or global caches?
   - `src_c/manifold_impl.c` — `manifold_impl_is_convex()` — does it cache results?
4. Check if `processOverlaps` param (set in C++ via `ManifoldParams().processOverlaps = true`) has a C equivalent and whether it persists between tests.
5. If state leak found: reset the state in the test or fix the implementation to not leak.
6. If timing-dependent: add `fprintf(stderr, ...)` to the Minkowski path to compare values between solo and suite runs.

**Verification (do this LAST, after all test porting is done):**
```bash
for i in $(seq 1 5); do
  echo -n "Run $i: "
  ./test_manifold 2>&1 | grep -E 'NonConvexConvexMinkowskiDiff|=== .* tests'
done
```
If it fails even once in 5 runs, the flakiness is confirmed and must be fixed.
