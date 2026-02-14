# Port Test: `Smooth_Torus`

> **Read `prompts/tests/PREAMBLE.md` first.**

> **Test #35 of 58.** Port this single test to C. Make it pass. Commit. Then move to the next test.

> **⚠️ NEVER WEAKEN A TEST.** The C++ values are ground truth.
> **⚠️ DO NOT SKIP THIS TEST.** If it needs a new C API function, port it. If it needs a subsystem, port the minimum needed.
> **⚠️ DO NOT SPAWN SUB-AGENTS.** Do all work yourself, directly.

---

## Build & Run

```bash
cd /Users/admin/projects/manifold
make
TEST=Smooth_Torus ./test_manifold          # run this test alone
SKIP_SLOW=1 ./test_manifold               # verify no regressions
```

---

## C++ Source (`test/smooth_test.cpp`)

```cpp
TEST(Smooth, Torus) {
  MeshGL64 torusMesh =
      Manifold::Revolve(
          CrossSection::Circle(1, 8).Translate({2, 0}).ToPolygons(), 6)
          .GetMeshGL64();
  const int numTri = torusMesh.NumTri();

  // Create correct toroidal halfedge tangents - SmoothOut() is too generic to
  // do this perfectly.
  torusMesh.halfedgeTangent.resize(4 * 3 * numTri);
  for (int tri = 0; tri < numTri; ++tri) {
    const auto triVerts = torusMesh.GetTriVerts(tri);
    for (const int i : {0, 1, 2}) {
      vec4 tangent;
      const vec3 v = torusMesh.GetVertPos(triVerts[i]);
      const vec3 v1 = torusMesh.GetVertPos(triVerts[(i + 1) % 3]);
      const vec3 edge = v1 - v;
      if (edge.z == 0) {
        vec3 tan(v.y, -v.x, 0);
        tan *= la::dot(tan, edge) < 0 ? -1.0 : 1.0;
        tangent = CircularTangent(tan, edge);
      } else if (std::abs(la::determinant(mat2(vec2(v), vec2(edge)))) < 1e-5) {
        const double theta = std::asin(v.z);
        vec2 xy(v);
        const double r = la::length(xy);
        xy = xy / r * v.z * (r > 2 ? -1.0 : 1.0);
        vec3 tan(xy.x, xy.y, std::cos(theta));
        tan *= la::dot(tan, edge) < 0 ? -1.0 : 1.0;
        tangent = CircularTangent(tan, edge);
      } else {
        tangent = {0, 0, 0, -1};
      }
      const int e = 3 * tri + i;
      for (const int j : {0, 1, 2, 3})
        torusMesh.halfedgeTangent[4 * e + j] = tangent[j];
    }
  }

  Manifold smooth = Manifold(torusMesh)
                        .RefineToLength(0.1)
                        .CalculateCurvature(-1, 0)
                        .CalculateNormals(1);
  MeshGL out = smooth.GetMeshGL();
  float maxMeanCurvature = 0;
  for (size_t i = 0; i < out.vertProperties.size(); i += 7) {
    vec3 v(out.vertProperties[i], out.vertProperties[i + 1],
           out.vertProperties[i + 2]);
    vec3 p(v.x, v.y, 0);
    p = la::normalize(p) * 2.0;
    double r = la::length(v - p);
    ASSERT_NEAR(r, 1, 0.006);
    maxMeanCurvature =
        std::max(maxMeanCurvature, std::abs(out.vertProperties[i + 3]));
  }
  EXPECT_NEAR(maxMeanCurvature, 1.63, 0.01);

#ifdef MANIFOLD_EXPORT
  ExportOptions options2;
  options2.faceted = false;
  options2.mat.normalIdx = 1;
  options2.mat.roughness = 0;
  if (options.exportModels) ExportMesh("smoothTorus.glb", out, options2);
#endif
}
```

## Dependencies (helpers/subsystems this test needs)

- **`CrossSection`** — find in `include/manifold/cross_section.h`. If not yet ported to C, port it first.
- **`MeshGL`** — find in `include/manifold/manifold.h`. If not yet ported to C, port it first.


## Steps

1. Read the C++ test above. Understand every line — what objects are created, what operations are performed, what values are asserted.
2. Check if all needed `manifold_*` API functions exist in `src_c/manifold_api.h`. If any are missing, port them from the C++ source first.
3. Write the C test function `test_Smooth_Torus(void)` in `src_c/test_manifold.c`.
4. Add `RUN_TEST(Smooth_Torus);` to the appropriate section in `main()`.
5. Build and run: `TEST=Smooth_Torus ./test_manifold`
6. If it fails — the C **implementation** is wrong, not the test. Debug and fix the implementation.
7. Run full suite: `SKIP_SLOW=1 ./test_manifold` — no regressions.
8. Commit: `test: port Smooth_Torus from C++`

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

1. This prompt is about ONE test: `Smooth_Torus`. Do not work on other tests.
2. Never weaken assertions. The C++ expected values are the ground truth.
3. If a C++ method is missing from the C API, port it before writing the test.
4. Run `SKIP_SLOW=1 ./test_manifold` after to verify no regressions.
5. Commit: `test: port Smooth_Torus from C++`


---

## Completion

When this test passes with no regressions, output exactly:

**Phase complete — test ported and verified**