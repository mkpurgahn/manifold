# Manifold C Port — Test Debugging Prompt

> **Status:** 151 tests ported. **129 PASS, 6 FAIL, 7 TIMEOUT (>300s), 2 TIMEOUT (>60s, may just be slow).** The porting phase is complete. This prompt focuses on debugging the remaining failures.

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
TEST=Hull_Sphere perl -e 'alarm 120; exec @ARGV' ./test_manifold 2>&1
```

---

## Failing Tests — Fix These

### 1. Boolean_CreatePropertiesSlow — `num_prop=0 != 3`
**Error:** `manifold_num_prop(&result)=0 != 3` — properties are not propagated through boolean operations.
**C++ ref:** `test/boolean_test.cpp` — search for `CreatePropertiesSlow`
**Root cause:** The boolean operation loses the property channels. Check `manifold_boolean.c` — the property merge/interpolation logic during face splitting. Compare with C++ `Boolean3::Result()` property handling.
**Category:** Properties propagation bug.

### 2. Boolean_AlmostCoplanar — `genus=-1 != 0`
**Error:** `manifold_genus(&result)=-1 != 0` — result is non-manifold (genus should be 0 for a valid solid).
**C++ ref:** `test/boolean_test.cpp` — search for `AlmostCoplanar`
**Root cause:** Coplanar face handling in the boolean algorithm produces a non-manifold mesh. The genus=-1 means edges aren't properly paired. Check `manifold_boolean.c` coplanar face detection and stitching.
**Category:** Boolean topology / coplanar handling.

### 3. Boolean_Perturb1 — `num_vert=25 != 24`, volume/area off
**Error:** Extra vertex (25 vs 24), volume off by 0.67, area off by 0.88.
**C++ ref:** `test/boolean_test.cpp` — search for `Perturb1`
**Root cause:** Boolean perturbation produces one extra vertex that shouldn't exist. Related to degenerate edge handling or vertex merging tolerance. Check the simplification/cleanup pass after boolean.
**Category:** Boolean topology / vertex merging.

### 4. Boolean_BatchBoolean — volume off by 88.57
**Error:** `|33245.17 - 33156.59| = 88.57 > 0.33` — large volume error in batch boolean.
**C++ ref:** `test/boolean_test.cpp` — search for `BatchBoolean`
**Root cause:** Accumulated error in batch boolean (multiple booleans chained). Could be the same coplanar issue as #2 compounding, or a CSG tree evaluation bug. Check `manifold_csg_tree.c`.
**Category:** Boolean accuracy / CSG tree.

### 5. Smooth_TruncatedCone — area/volume off
**Error:** area off by 21.4, volume off by 33.8 (expected ~1158/768, got ~1137/802).
**C++ ref:** `test/smooth_test.cpp` — search for `TruncatedCone`
**Root cause:** Smooth subdivision (Catmull-Clark or similar) produces slightly wrong geometry. Check `manifold_smooth.c` — the tangent/normal calculation for cone-like inputs, or the subdivision weight calculation.
**Category:** Smooth subdivision accuracy.

### 6. Smooth_Csaszar — area/volume way off
**Error:** area off by 3109 (82999 vs 79890), volume off by 383.5 (12333 vs 11950).
**C++ ref:** `test/smooth_test.cpp` — search for `Csaszar`
**Root cause:** Csaszar polyhedron is a complex non-convex surface. The smooth subdivision is significantly off — likely a fundamental issue in the smoothing algorithm for non-trivial topology. May share root cause with #5 but amplified.
**Category:** Smooth subdivision accuracy.

---

## Timeout Tests — Investigate These

These tests take >60s or >300s. Some may be legitimate slow tests that just need more time. Others may be infinite loops.

### Likely Slow (timed out at 60s — try 300s):
- **Boolean_Perturb3** — complex perturbation test, may just be slow
- **Hull_Sphere** — convex hull of many points, may just be slow
- **SDF_SphereShell** — SDF evaluation, may just be slow
- **Properties_ToleranceSphere** — tolerance computation, may just be slow

### Likely Hanging (timed out at 300s):
- **Hull_MengerSponge** — Menger sponge hull (extremely complex geometry). May be O(n³) in hull algorithm.
- **Samples_Sponge1** — Sponge sample with boolean ops. May loop infinitely in boolean.
- **Samples_Sponge4** — Same but deeper recursion.
- **Samples_Scallop** — Complex sample. Check for infinite loops.

**Debugging approach for hangs:** Run with `TEST=Hull_MengerSponge` and attach a debugger or add `fprintf(stderr, ...)` progress prints to find which function loops forever. The most likely hang points are:
- `manifold_hull.c` — gift wrapping or incremental hull with degenerate faces
- `manifold_boolean.c` — face splitting loop with coplanar faces
- `manifold_impl.c` — mesh cleanup/simplification

---

## Approach

Work on ONE failing test at a time, in this order (easiest → hardest):

1. **Boolean_CreatePropertiesSlow** — likely a simple property propagation miss
2. **Boolean_Perturb1** — vertex count off by 1, should be traceable
3. **Boolean_AlmostCoplanar** — coplanar handling bug
4. **Boolean_BatchBoolean** — may fix itself once #2/#3 are fixed
5. **Smooth_TruncatedCone** — subdivision accuracy
6. **Smooth_Csaszar** — may fix itself once #5 is fixed

For timeouts, try the 60s ones at 300s first — they may just be slow. For the 300s hangs, add progress logging to find the loop.

### For each test:

1. Read the C++ test in `test/` to understand what it expects
2. Read the C++ implementation of the failing function
3. Compare with the C port in `src_c/`
4. Find the divergence
5. Fix the C implementation
6. Run the single test: `TEST=Boolean_Perturb1 ./test_manifold`
7. Verify it passes, then run all tests to check for regressions
8. Commit: `fix: [test_name] — [what was wrong]`

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
| `src_c/manifold_api.c` | `bindings/c/manifoldc.cpp` | C API bindings |

C++ tests are in `test/boolean_test.cpp`, `test/smooth_test.cpp`, `test/hull_test.cpp`, `test/properties_test.cpp`, `test/samples_test.cpp`, `test/sdf_test.cpp`, `test/manifold_test.cpp`, `test/manifoldc_test.cpp`.

---

## Rules

1. Fix one test at a time. Do not scatter changes across multiple bugs.
2. Always run the full test suite after a fix to check for regressions.
3. Never weaken assertions. The C++ values are the ground truth.
4. Compare C vs C++ implementation line by line when debugging.
5. Commit after each fix: `fix: [test_name] — [root cause]`
