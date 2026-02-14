# Port Test: `Smooth_Manual`

> **Read `prompts/tests/PREAMBLE.md` first.**

> **Test #34 of 58.** Port this single test to C. Make it pass. Commit. Then move to the next test.

> **⚠️ NEVER WEAKEN A TEST.** The C++ values are ground truth.
> **⚠️ DO NOT SKIP THIS TEST.** If it needs a new C API function, port it. If it needs a subsystem, port the minimum needed.
> **⚠️ DO NOT SPAWN SUB-AGENTS.** Do all work yourself, directly.

---

## Build & Run

```bash
cd /Users/admin/projects/manifold
make
TEST=Smooth_Manual ./test_manifold          # run this test alone
SKIP_SLOW=1 ./test_manifold               # verify no regressions
```

---

## C++ Source (`test/smooth_test.cpp`)

```cpp
TEST(Smooth, Manual) {
  // Unit Octahedron
  const auto oct = Manifold::Sphere(1, 4).GetMeshGL();
  MeshGL smooth = Manifold::Smooth(oct).GetMeshGL();
  // Sharpen the edge from vert 4 to 5
  smooth.halfedgeTangent[4 * 6 + 3] = 0;
  smooth.halfedgeTangent[4 * 22 + 3] = 0;
  smooth.halfedgeTangent[4 * 16 + 3] = 0;
  smooth.halfedgeTangent[4 * 18 + 3] = 0;
  Manifold interp(smooth);
  interp = interp.Refine(100);

  ExpectMeshes(interp, {{40002, 80000}});
  EXPECT_NEAR(interp.Volume(), 3.74, 0.01);
  EXPECT_NEAR(interp.SurfaceArea(), 11.78, 0.01);

#ifdef MANIFOLD_EXPORT
  if (options.exportModels) {
    interp = interp.CalculateCurvature(-1, 0).SetProperties(
        3, [](double* newProp, vec3 pos, const double* oldProp) {
          const vec3 red(1, 0, 0);
          const vec3 purple(1, 0, 1);
          vec3 color = la::lerp(purple, red, smoothstep(0.0, 2.0, oldProp[0]));
          for (const int i : {0, 1, 2}) newProp[i] = color[i];
        });
    const MeshGL out = interp.GetMeshGL();
    ExportOptions options;
    options.mat.roughness = 0.1;
    options.mat.colorIdx = 0;
    ExportMesh("manual.glb", out, options);
  }
#endif
}
```

## Dependencies (helpers/subsystems this test needs)

- **`MeshGL`** — find in `include/manifold/manifold.h`. If not yet ported to C, port it first.


## Steps

1. Read the C++ test above. Understand every line — what objects are created, what operations are performed, what values are asserted.
2. Check if all needed `manifold_*` API functions exist in `src_c/manifold_api.h`. If any are missing, port them from the C++ source first.
3. Write the C test function `test_Smooth_Manual(void)` in `src_c/test_manifold.c`.
4. Add `RUN_TEST(Smooth_Manual);` to the appropriate section in `main()`.
5. Build and run: `TEST=Smooth_Manual ./test_manifold`
6. If it fails — the C **implementation** is wrong, not the test. Debug and fix the implementation.
7. Run full suite: `SKIP_SLOW=1 ./test_manifold` — no regressions.
8. Commit: `test: port Smooth_Manual from C++`

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

1. This prompt is about ONE test: `Smooth_Manual`. Do not work on other tests.
2. Never weaken assertions. The C++ expected values are the ground truth.
3. If a C++ method is missing from the C API, port it before writing the test.
4. Run `SKIP_SLOW=1 ./test_manifold` after to verify no regressions.
5. Commit: `test: port Smooth_Manual from C++`


---

## Completion

When this test passes with no regressions, output exactly:

**Phase complete — test ported and verified**