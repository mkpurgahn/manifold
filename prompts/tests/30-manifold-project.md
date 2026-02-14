# Port Test: `Manifold_Project`

> **Read `prompts/tests/PREAMBLE.md` first.**

> **Test #30 of 58.** Port this single test to C. Make it pass. Commit. Then move to the next test.

> **⚠️ NEVER WEAKEN A TEST.** The C++ values are ground truth.
> **⚠️ DO NOT SKIP THIS TEST.** If it needs a new C API function, port it. If it needs a subsystem, port the minimum needed.
> **⚠️ DO NOT SPAWN SUB-AGENTS.** Do all work yourself, directly.

---

## Build & Run

```bash
cd /Users/admin/projects/manifold
make
TEST=Manifold_Project ./test_manifold          # run this test alone
SKIP_SLOW=1 ./test_manifold               # verify no regressions
```

---

## C++ Source (`test/manifold_test.cpp`)

```cpp
TEST(Manifold, Project) {
  MeshGL input;
  input.numProp = 3;
  input.vertProperties = {0,    0,       0,     //
                          -2,   -0.7,    -0.1,  //
                          -2,   -0.7,    0,     //
                          -1.9, -0.7,    -0.1,  //
                          -1.9, -0.6901, -0.1,  //
                          -1.9, -0.7,    0,     //
                          -1.9, -0.6901, 0,     //
                          -2,   -1,      3,     //
                          -1.9, -1,      3,     //
                          -2,   -1,      4,     //
                          -1.9, -1,      4,     //
                          -1.9, -0.6901, 3,     //
                          -1.9, -0.6901, 4,     //
                          -1.7, -0.6901, 3,     //
                          -1.7, -0.6901, 3.2,   //
                          -2,   0,       -0.1,  //
                          -2,   0,       0,     //
                          -2,   0,       3,     //
                          -2,   0,       4,     //
                          -1.7, 0,       3,     //
                          -1.7, 0,       3.2,   //
                          -1,   -0.6901, -0.1,  //
                          -1,   -0.6901, 0,     //
                          -1,   -0.6901, 3.2,   //
                          -1,   -0.6901, 4,     //
                          -1,   0,       -0.1,  //
                          -1,   0,       0,     //
                          -1,   0,       3.2,   //
                          -1,   0,       4};
  input.triVerts = {1,  3,  2,   //
                    1,  4,  3,   //
                    2,  3,  5,   //
                    5,  6,  2,   //
                    3,  4,  6,   //
                    5,  3,  6,   //
                    6,  4,  21,  //
                    26, 22, 25,  //
                    21, 25, 22,  //
                    25, 15, 26,  //
                    26, 6,  22,  //
                    21, 4,  25,  //
                    21, 22, 6,   //
                    16, 26, 15,  //
                    16, 6,  26,  //
                    4,  15, 25,  //
                    15, 1,  16,  //
                    16, 2,  6,   //
                    4,  1,  15,  //
                    1,  2,  16,  //
                    12, 14, 23,  //
                    12, 13, 14,  //
                    12, 11, 13,  //
                    18, 9,  12,  //
                    11, 7,  17,  //
                    7,  9,  18,  //
                    17, 7,  18,  //
                    13, 11, 19,  //
                    17, 18, 20,  //
                    19, 11, 17,  //
                    19, 17, 20,  //
                    14, 13, 20,  //
                    18, 12, 24,  //
                    20, 13, 19,  //
                    20, 18, 27,  //
                    12, 10, 11,  //
                    24, 12, 23,  //
                    9,  10, 12,  //
                    9,  8,  10,  //
                    8,  11, 10,  //
                    8,  7,  11,  //
                    8,  9,  7,   //
                    14, 20, 27,  //
                    24, 28, 18,  //
                    27, 18, 28,  //
                    23, 14, 27,  //
                    24, 23, 28,  //
                    28, 23, 27};
  Manifold in(input);
  CrossSection projected = in.Project();
  EXPECT_NEAR(projected.Area(), 0.72, 0.01);
}
```

## Dependencies (helpers/subsystems this test needs)

- **`CrossSection`** — find in `include/manifold/cross_section.h`. If not yet ported to C, port it first.
- **`MeshGL`** — find in `include/manifold/manifold.h`. If not yet ported to C, port it first.


## Steps

1. Read the C++ test above. Understand every line — what objects are created, what operations are performed, what values are asserted.
2. Check if all needed `manifold_*` API functions exist in `src_c/manifold_api.h`. If any are missing, port them from the C++ source first.
3. Write the C test function `test_Manifold_Project(void)` in `src_c/test_manifold.c`.
4. Add `RUN_TEST(Manifold_Project);` to the appropriate section in `main()`.
5. Build and run: `TEST=Manifold_Project ./test_manifold`
6. If it fails — the C **implementation** is wrong, not the test. Debug and fix the implementation.
7. Run full suite: `SKIP_SLOW=1 ./test_manifold` — no regressions.
8. Commit: `test: port Manifold_Project from C++`

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

---

## Rules

1. This prompt is about ONE test: `Manifold_Project`. Do not work on other tests.
2. Never weaken assertions. The C++ expected values are the ground truth.
3. If a C++ method is missing from the C API, port it before writing the test.
4. Run `SKIP_SLOW=1 ./test_manifold` after to verify no regressions.
5. Commit: `test: port Manifold_Project from C++`


---

## Completion

When this test passes with no regressions, output exactly:

**Phase complete — test ported and verified**