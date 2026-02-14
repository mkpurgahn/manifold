# Port Test: `Samples_GyroidModule`

> **Read `prompts/tests/PREAMBLE.md` first.**

> **Test #56 of 58.** Port this single test to C. Make it pass. Commit. Then move to the next test.

> **⚠️ NEVER WEAKEN A TEST.** The C++ values are ground truth.
> **⚠️ DO NOT SKIP THIS TEST.** If it needs a new C API function, port it. If it needs a subsystem, port the minimum needed.
> **⚠️ DO NOT SPAWN SUB-AGENTS.** Do all work yourself, directly.

---

## Build & Run

```bash
cd /Users/admin/projects/manifold
make
TEST=Samples_GyroidModule ./test_manifold          # run this test alone
SKIP_SLOW=1 ./test_manifold               # verify no regressions
```

---

## C++ Source (`test/samples_test.cpp`)

```cpp
TEST(Samples, GyroidModule) {
  const double size = 20;
  Manifold gyroid = GyroidModule(size);
  EXPECT_LE(gyroid.NumDegenerateTris(), 4);
  EXPECT_EQ(gyroid.Genus(), 15);
  CheckGL(gyroid);

  const Box bounds = gyroid.BoundingBox();
  const double epsilon = gyroid.GetEpsilon();
  EXPECT_NEAR(bounds.min.z, 0, epsilon);
  EXPECT_NEAR(bounds.max.z, size * std::sqrt(2.0), epsilon);

  CrossSection slice(gyroid.Slice(5));
  EXPECT_EQ(slice.NumContour(), 4);
  EXPECT_NEAR(slice.Area(), 121.9, 0.1);
  Manifold extrusion = Manifold::Extrude(slice.ToPolygons(), 1);
  EXPECT_EQ(extrusion.Genus(), -3);

#ifdef MANIFOLD_EXPORT
  if (options.exportModels)
    ExportMesh("gyroidModule.glb", gyroid.GetMeshGL(), {});
#endif
}
```

## Dependencies (helpers/subsystems this test needs)

- **`Gyroid`** — find in `samples/src/gyroid_module.cpp`. If not yet ported to C, port it first.
- **`CheckGL`** — find in `test/test_main.cpp`. If not yet ported to C, port it first.
- **`CrossSection`** — find in `include/manifold/cross_section.h`. If not yet ported to C, port it first.
- **`MeshGL`** — find in `include/manifold/manifold.h`. If not yet ported to C, port it first.


## Steps

1. Read the C++ test above. Understand every line — what objects are created, what operations are performed, what values are asserted.
2. Check if all needed `manifold_*` API functions exist in `src_c/manifold_api.h`. If any are missing, port them from the C++ source first.
3. Write the C test function `test_Samples_GyroidModule(void)` in `src_c/test_manifold.c`.
4. Add `RUN_TEST(Samples_GyroidModule);` to the appropriate section in `main()`.
5. Build and run: `TEST=Samples_GyroidModule ./test_manifold`
6. If it fails — the C **implementation** is wrong, not the test. Debug and fix the implementation.
7. Run full suite: `SKIP_SLOW=1 ./test_manifold` — no regressions.
8. Commit: `test: port Samples_GyroidModule from C++`

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

1. This prompt is about ONE test: `Samples_GyroidModule`. Do not work on other tests.
2. Never weaken assertions. The C++ expected values are the ground truth.
3. If a C++ method is missing from the C API, port it before writing the test.
4. Run `SKIP_SLOW=1 ./test_manifold` after to verify no regressions.
5. Commit: `test: port Samples_GyroidModule from C++`


---

## Completion

When this test passes with no regressions, output exactly:

**Phase complete — test ported and verified**