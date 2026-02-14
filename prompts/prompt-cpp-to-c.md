# Manifold C Port — Test Porting

> **Status:** 161/161 ported tests **ALL PASS** (0 failures). 58 C++ tests remain unported. The core port is functionally complete.

> **⚠️ NEVER WEAKEN A TEST.** The C++ tests define correct behavior.
> **⚠️ DO NOT SPAWN SUB-AGENTS.** Do all work yourself, directly.

---

## Build & Run

```bash
cd /Users/admin/projects/manifold
make                                    # builds test_manifold
TEST=<TestName> ./test_manifold         # run single test
SKIP_SLOW=1 ./test_manifold             # run fast tests, verify no regressions
./test_manifold                          # run ALL tests (slow ones take 2-5 min each)
```

Slow tests (hull, sponge, scallop, minkowski) need **240+ seconds**. Don't assume a test hangs until 300s.

---

## Phase Loop

Each unported test is its own phase. One prompt per ralph loop. The sequencer advances only when both loop and adversary agree the test passes.

```bash
./run-all-tests.sh                    # run all 58 phases
./run-all-tests.sh --start-phase 5    # resume from phase 5
```

Phase prompts: `prompts/tests/NN-test-name.md`
Shared rules: `prompts/tests/PREAMBLE.md`
Completion promise: **"Phase complete — test ported and verified"**

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

## Flaky Test: `Boolean_NonConvexConvexMinkowskiDifference`

Observed to **pass alone, fail in full suite** — possible state leak.

**Failure:** `area off by 0.12, genus=0 vs 5`

**Investigate after all 58 phases complete:**
```bash
for i in $(seq 1 5); do
  echo -n "Run $i: "
  ./test_manifold 2>&1 | grep -E 'NonConvexConvexMinkowskiDiff|=== .* tests'
done
```
Check for global/static state in `manifold_minkowski_impl()`, `manifold_batch_boolean()`, and `ManifoldParams` that persists between tests.
