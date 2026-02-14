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

## Unported Tests — Port ONE at a Time (58 remaining)

**⚠️ FOCUS: Pick the NEXT unported test from the list below. Read its individual prompt file in `prompts/tests/`. Port it completely. If it needs a new C API function, port that function. If it needs CrossSection, port CrossSection. If it needs MeshGL, port MeshGL. Nothing is "blocked" — everything is achievable. Do not skip tests. Do not categorize tests as "not achievable." Port them in order.**

### Workflow: Ralph Loop Per Test

Each test has its own prompt file in `prompts/tests/NN-test-name.md`. The workflow is:

1. Start a ralph loop with the next test's prompt file (e.g., `prompts/tests/01-boolean-meshglroundtrip.md`).
2. The loop ports the test, makes it pass, and commits.
3. Once both the main loop and adversary agree the test passes with no regressions — move to the next test.
4. Kick off a new loop with the next prompt file (e.g., `prompts/tests/02-boolean-mixedproperties.md`).
5. Repeat until all 58 are done.

### Per-test prompt files

Each file contains: the exact C++ source code, any helper dependencies, step-by-step porting instructions, and build/run commands. Read the prompt file BEFORE starting work on that test.

```
prompts/tests/
├── 01-boolean-meshglroundtrip.md
├── 02-boolean-mixedproperties.md
├── 03-booleancomplex-sphere.md
├── ...
├── 56-samples-gyroidmodule.md
├── 57-polygonfuzz-triangulationnocrash.md
└── 58-polygonfuzz-triangulationnocrashrounded.md
```

### The list (port in this order):

```
 1. Boolean_MeshGLRoundTrip              test/boolean_test.cpp
 2. Boolean_MixedProperties              test/boolean_test.cpp
 3. BooleanComplex_Sphere                test/boolean_complex_test.cpp
 4. BooleanComplex_Ring                  test/boolean_complex_test.cpp
 5. BooleanComplex_MeshRelation          test/boolean_complex_test.cpp
 6. BooleanComplex_InterpolatedNormals   test/boolean_complex_test.cpp
 7. BooleanComplex_HullMask              test/boolean_complex_test.cpp
 8. BooleanComplex_SimpleOffset          test/boolean_complex_test.cpp
 9. BooleanComplex_OffsetSelfIntersect   test/boolean_complex_test.cpp
10. BooleanComplex_OffsetTriangulationFailure  test/boolean_complex_test.cpp
11. BooleanComplex_Sweep                 test/boolean_complex_test.cpp
12. BooleanComplex_CraycloudBool         test/boolean_complex_test.cpp
13. BooleanComplex_GenericTwinBooleanTest7081   test/boolean_complex_test.cpp
14. BooleanComplex_GenericTwinBooleanTest7863   test/boolean_complex_test.cpp
15. BooleanComplex_Havocglass8Bool       test/boolean_complex_test.cpp
16. Manifold_GetMeshGL                   test/manifold_test.cpp
17. Manifold_MeshID                      test/manifold_test.cpp
18. Manifold_MeshRelation                test/manifold_test.cpp
19. Manifold_MeshRelationRefinePrecision test/manifold_test.cpp
20. Manifold_FaceIDRoundTrip             test/manifold_test.cpp
21. Manifold_Merge                       test/manifold_test.cpp
22. Manifold_MergeEmpty                  test/manifold_test.cpp
23. Manifold_MergeRefine                 test/manifold_test.cpp
24. Manifold_DecomposeProps              test/manifold_test.cpp
25. Manifold_InvalidInput5               test/manifold_test.cpp
26. Manifold_InvalidInput7               test/manifold_test.cpp
27. Manifold_ValidInputOneRunIndex       test/manifold_test.cpp
28. Manifold_WarpBatch                   test/manifold_test.cpp
29. Manifold_OpenscadCrash               test/manifold_test.cpp
30. Manifold_Project                     test/manifold_test.cpp
31. Manifold_Slice                       test/manifold_test.cpp
32. Manifold_SliceEmptyObject            test/manifold_test.cpp
33. ManifoldFuzz_SimpleCube              test/manifold_test.cpp
34. Smooth_Manual                        test/smooth_test.cpp
35. Smooth_Torus                         test/smooth_test.cpp
36. Properties_Coplanar                  test/properties_test.cpp
37. Properties_MingapStretchyBracelet    test/properties_test.cpp
38. CrossSection_Empty                   test/cross_section_test.cpp
39. CrossSection_Rect                    test/cross_section_test.cpp
40. CrossSection_Square                  test/cross_section_test.cpp
41. CrossSection_Transform               test/cross_section_test.cpp
42. CrossSection_Hull                    test/cross_section_test.cpp
43. CrossSection_HullError               test/cross_section_test.cpp
44. CrossSection_MirrorCheckAxis         test/cross_section_test.cpp
45. CrossSection_MirrorUnion             test/cross_section_test.cpp
46. CrossSection_Warp                    test/cross_section_test.cpp
47. CrossSection_Decompose               test/cross_section_test.cpp
48. CrossSection_FillRule                 test/cross_section_test.cpp
49. CrossSection_BatchBoolean            test/cross_section_test.cpp
50. CrossSection_NegativeOffset          test/cross_section_test.cpp
51. CrossSection_RoundOffset             test/cross_section_test.cpp
52. CrossSection_BevelOffset             test/cross_section_test.cpp
53. Samples_Bracelet                     test/samples_test.cpp
54. Samples_CondensedMatter16            test/samples_test.cpp
55. Samples_CondensedMatter64            test/samples_test.cpp
56. Samples_GyroidModule                 test/samples_test.cpp
57. PolygonFuzz_TriangulationNoCrash          test/polygon_test.cpp
58. PolygonFuzz_TriangulationNoCrashRounded   test/polygon_test.cpp
```

---

## Approach

**Port ONE test. Make it pass. Commit. Repeat.**

If a test needs something that doesn't exist yet in the C port:
- **Missing API function?** Port it from the C++ source. Add to `manifold_api.h` and `manifold_api.c`.
- **Missing helper?** Port the helper (e.g., `WithPositionColors` from `test/test_main.cpp`).
- **Missing subsystem (CrossSection, MeshGL)?** Port the minimum needed for the test to work. You'll build it incrementally — each subsequent test adds more.
- **Needs external mesh data?** Check `test/meshIO/` for the `.glb` files the C++ tests load.

**Do NOT skip a test because it's hard.** If test #1 needs MeshGL round-trip support, then porting MeshGL round-trip IS the work for test #1.

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
5. If state leak found: fix the implementation to not leak.
6. If timing-dependent: add `fprintf(stderr, ...)` to the Minkowski path to compare values between solo and suite runs.

**Verification (do this LAST, after all test porting is done):**
```bash
for i in $(seq 1 5); do
  echo -n "Run $i: "
  ./test_manifold 2>&1 | grep -E 'NonConvexConvexMinkowskiDiff|=== .* tests'
done
```
If it fails even once in 5 runs, the flakiness is confirmed and must be fixed.
