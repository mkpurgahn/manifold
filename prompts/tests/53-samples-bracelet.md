# Port Test: `Samples_Bracelet`

> **Test #53 of 58.** Port this single test to C. Make it pass. Commit. Then move to the next test.

> **⚠️ NEVER WEAKEN A TEST.** The C++ values are ground truth.
> **⚠️ DO NOT SKIP THIS TEST.** If it needs a new C API function, port it. If it needs a subsystem, port the minimum needed.
> **⚠️ DO NOT SPAWN SUB-AGENTS.** Do all work yourself, directly.

---

## Build & Run

```bash
cd /Users/admin/projects/manifold
make
TEST=Samples_Bracelet ./test_manifold          # run this test alone
SKIP_SLOW=1 ./test_manifold               # verify no regressions
```

---

## C++ Source (`test/samples_test.cpp`)

```cpp
TEST(Samples, Bracelet) {
  Manifold bracelet = StretchyBracelet();
  EXPECT_EQ(bracelet.NumDegenerateTris(), 0);
  EXPECT_EQ(bracelet.Genus(), 1);
  CheckGL(bracelet);

  CrossSection projection(bracelet.Project());
  projection = projection.Simplify(bracelet.BoundingBox().Scale() * 1e-8);
  Rect rect = projection.Bounds();
  Box box = bracelet.BoundingBox();
  EXPECT_FLOAT_EQ(rect.min.x, box.min.x);
  EXPECT_FLOAT_EQ(rect.min.y, box.min.y);
  EXPECT_FLOAT_EQ(rect.max.x, box.max.x);
  EXPECT_FLOAT_EQ(rect.max.y, box.max.y);
  EXPECT_NEAR(projection.Area(), 649, 1);
  EXPECT_EQ(projection.NumContour(), 2);
  Manifold extrusion = Manifold::Extrude(projection.ToPolygons(), 1);
  EXPECT_EQ(extrusion.NumDegenerateTris(), 0);
  EXPECT_EQ(extrusion.Genus(), 1);

  CrossSection slice(bracelet.Slice());
  EXPECT_EQ(slice.NumContour(), 2);
  EXPECT_NEAR(slice.Area(), 230.6, 0.1);
  extrusion = Manifold::Extrude(slice.ToPolygons(), 1);
  EXPECT_EQ(extrusion.Genus(), 1);

#ifdef MANIFOLD_EXPORT
  if (options.exportModels)
    ExportMesh("bracelet.glb", bracelet.GetMeshGL(), {});
#endif
}
```

## Dependencies (helpers/subsystems this test needs)

- **`StretchyBracelet`** — find in `samples/src/bracelet.cpp`. If not yet ported to C, port it first.
- **`CheckGL`** — find in `test/test_main.cpp`. If not yet ported to C, port it first.
- **`Bracelet`** — find in `samples/src/bracelet.cpp`. If not yet ported to C, port it first.
- **`CrossSection`** — find in `include/manifold/cross_section.h`. If not yet ported to C, port it first.
- **`MeshGL`** — find in `include/manifold/manifold.h`. If not yet ported to C, port it first.


## Steps

1. Read the C++ test above. Understand every line — what objects are created, what operations are performed, what values are asserted.
2. Check if all needed `manifold_*` API functions exist in `src_c/manifold_api.h`. If any are missing, port them from the C++ source first.
3. Write the C test function `test_Samples_Bracelet(void)` in `src_c/test_manifold.c`.
4. Add `RUN_TEST(Samples_Bracelet);` to the appropriate section in `main()`.
5. Build and run: `TEST=Samples_Bracelet ./test_manifold`
6. If it fails — the C **implementation** is wrong, not the test. Debug and fix the implementation.
7. Run full suite: `SKIP_SLOW=1 ./test_manifold` — no regressions.
8. Commit: `test: port Samples_Bracelet from C++`

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

1. This prompt is about ONE test: `Samples_Bracelet`. Do not work on other tests.
2. Never weaken assertions. The C++ expected values are the ground truth.
3. If a C++ method is missing from the C API, port it before writing the test.
4. Run `SKIP_SLOW=1 ./test_manifold` after to verify no regressions.
5. Commit: `test: port Samples_Bracelet from C++`
