# Port Test: `BooleanComplex_SimpleOffset`

> **Test #8 of 58.** Port this single test to C. Make it pass. Commit. Then move to the next test.

> **⚠️ NEVER WEAKEN A TEST.** The C++ values are ground truth.
> **⚠️ DO NOT SKIP THIS TEST.** If it needs a new C API function, port it. If it needs a subsystem, port the minimum needed.
> **⚠️ DO NOT SPAWN SUB-AGENTS.** Do all work yourself, directly.

---

## Build & Run

```bash
cd /Users/admin/projects/manifold
make
TEST=BooleanComplex_SimpleOffset ./test_manifold          # run this test alone
SKIP_SLOW=1 ./test_manifold               # verify no regressions
```

---

## C++ Source (`test/boolean_complex_test.cpp`)

```cpp
TEST(BooleanComplex, SimpleOffset) {
  std::string file = __FILE__;
  std::string dir = file.substr(0, file.rfind('/'));
  MeshGL seeds = ImportMesh(dir + "/models/" + "Generic_Twin_91.1.t0.glb");
  EXPECT_TRUE(seeds.NumTri() > 10);
  EXPECT_TRUE(seeds.NumVert() > 10);
  // Unique edges
  std::vector<std::pair<int, int>> edges;
  for (size_t i = 0; i < seeds.NumTri(); i++) {
    const int k[3] = {1, 2, 0};
    for (const int j : {0, 1, 2}) {
      int v1 = seeds.triVerts[i * 3 + j];
      int v2 = seeds.triVerts[i * 3 + k[j]];
      if (v2 > v1) edges.push_back(std::make_pair(v1, v2));
    }
  }
  manifold::Manifold c;
  // Vertex Spheres
  Manifold sph = Manifold::Sphere(1, 8);
  for (size_t i = 0; i < seeds.NumVert(); i++) {
    vec3 vpos(seeds.vertProperties[3 * i + 0], seeds.vertProperties[3 * i + 1],
              seeds.vertProperties[3 * i + 2]);
    Manifold vsph = sph.Translate(vpos);
    c += vsph;
  }
  // Edge Cylinders
  for (size_t i = 0; i < edges.size(); i++) {
    vec3 ev1 = vec3(seeds.vertProperties[3 * edges[i].first + 0],
                    seeds.vertProperties[3 * edges[i].first + 1],
                    seeds.vertProperties[3 * edges[i].first + 2]);
    vec3 ev2 = vec3(seeds.vertProperties[3 * edges[i].second + 0],
                    seeds.vertProperties[3 * edges[i].second + 1],
                    seeds.vertProperties[3 * edges[i].second + 2]);
    vec3 edge = ev2 - ev1;
    double len = la::length(edge);
    if (len < std::numeric_limits<float>::min()) continue;
    manifold::Manifold origin_cyl = manifold::Manifold::Cylinder(len, 1, 1, 8);
    vec3 evec(-1 * edge.x, -1 * edge.y, edge.z);
    quat q = rotation_quat(normalize(evec), vec3(0, 0, 1));
    manifold::Manifold right = origin_cyl.Transform({la::qmat(q), ev1});
    c += right;
  }
  // Triangle Volumes
  for (size_t i = 0; i < seeds.NumTri(); i++) {
    int eind[3];
    for (int j = 0; j < 3; j++) eind[j] = seeds.triVerts[i * 3 + j];
    std::vector<vec3> ev;
    for (int j = 0; j < 3; j++) {
      ev.push_back(vec3(seeds.vertProperties[3 * eind[j] + 0],
                        seeds.vertProperties[3 * eind[j] + 1],
                        seeds.vertProperties[3 * eind[j] + 2]));
    }
    vec3 a = ev[0] - ev[2];
    vec3 b = ev[1] - ev[2];
    vec3 n = la::normalize(la::cross(a, b));
    if (!all(isfinite(n))) continue;
    // Extrude the points above and below the plane of the triangle
    vec3 pnts[6];
    for (int j = 0; j < 3; j++) pnts[j] = ev[j] + n;
    for (int j = 3; j < 6; j++) pnts[j] = ev[j - 3] - n;
    // Construct the points and faces of the new manifold
    double pts[3 * 6] = {pnts[4].x, pnts[4].y, pnts[4].z, pnts[3].x, pnts[3].y,
                         pnts[3].z, pnts[0].x, pnts[0].y, pnts[0].z, pnts[1].x,
                         pnts[1].y, pnts[1].z, pnts[5].x, pnts[5].y, pnts[5].z,
                         pnts[2].x, pnts[2].y, pnts[2].z};
    int faces[24] = {
        faces[0] = 0,  faces[1] = 1,  faces[2] = 4,   // 1 2 5
        faces[3] = 2,  faces[4] = 3,  faces[5] = 5,   // 3 4 6
        faces[6] = 1,  faces[7] = 0,  faces[8] = 3,   // 2 1 4
        faces[9] = 3,  faces[10] = 2, faces[11] = 1,  // 4 3 2
        faces[12] = 3, faces[13] = 0, faces[14] = 4,  // 4 1 5
        faces[15] = 4, faces[16] = 5, faces[17] = 3,  // 5 6 4
        faces[18] = 5, faces[19] = 4, faces[20] = 1,  // 6 5 2
        faces[21] = 1, faces[22] = 2, faces[23] = 5   // 2 3 6
    };
    manifold::MeshGL64 tri_m;
    for (int j = 0; j < 18; j++)
      tri_m.vertProperties.insert(tri_m.vertProperties.end(), pts[j]);
    for (int j = 0; j < 24; j++)
      tri_m.triVerts.insert(tri_m.triVerts.end(), faces[j]);
    manifold::Manifold right(tri_m);
    c += right;
  }
  EXPECT_EQ(c.Status(), Manifold::Error::NoError);
}
```

## Dependencies (helpers/subsystems this test needs)

- **`MeshGL`** — find in `include/manifold/manifold.h`. If not yet ported to C, port it first.


## Steps

1. Read the C++ test above. Understand every line — what objects are created, what operations are performed, what values are asserted.
2. Check if all needed `manifold_*` API functions exist in `src_c/manifold_api.h`. If any are missing, port them from the C++ source first.
3. Write the C test function `test_BooleanComplex_SimpleOffset(void)` in `src_c/test_manifold.c`.
4. Add `RUN_TEST(BooleanComplex_SimpleOffset);` to the appropriate section in `main()`.
5. Build and run: `TEST=BooleanComplex_SimpleOffset ./test_manifold`
6. If it fails — the C **implementation** is wrong, not the test. Debug and fix the implementation.
7. Run full suite: `SKIP_SLOW=1 ./test_manifold` — no regressions.
8. Commit: `test: port BooleanComplex_SimpleOffset from C++`

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

1. This prompt is about ONE test: `BooleanComplex_SimpleOffset`. Do not work on other tests.
2. Never weaken assertions. The C++ expected values are the ground truth.
3. If a C++ method is missing from the C API, port it before writing the test.
4. Run `SKIP_SLOW=1 ./test_manifold` after to verify no regressions.
5. Commit: `test: port BooleanComplex_SimpleOffset from C++`
