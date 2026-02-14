# Port Test: `Properties_Coplanar`

> **Read `prompts/tests/PREAMBLE.md` first.**

> **Test #36 of 58.** Port this single test to C. Make it pass. Commit. Then move to the next test.

> **⚠️ NEVER WEAKEN A TEST.** The C++ values are ground truth.
> **⚠️ DO NOT SKIP THIS TEST.** If it needs a new C API function, port it. If it needs a subsystem, port the minimum needed.
> **⚠️ DO NOT SPAWN SUB-AGENTS.** Do all work yourself, directly.

---

## Build & Run

```bash
cd /Users/admin/projects/manifold
make
TEST=Properties_Coplanar ./test_manifold          # run this test alone
SKIP_SLOW=1 ./test_manifold               # verify no regressions
```

---

## C++ Source (`test/properties_test.cpp`)

```cpp
TEST(Properties, Coplanar) {
  Manifold peg = Manifold::Cube({1, 1, 2}).Translate({1, 1, 0}).AsOriginal();
  const size_t pegID = peg.OriginalID();
  Manifold hole = Manifold::Cube({3, 3, 1}) - peg;
  hole = hole.AsOriginal();
  std::vector<MeshGL> input = {hole.GetMeshGL(), peg.GetMeshGL()};
  EXPECT_EQ(peg.Genus(), 0);
  EXPECT_EQ(hole.Genus(), 1);

  Manifold result = hole + peg;
  EXPECT_EQ(result.Genus(), 0);
  RelatedGL(result, input);

  MeshGL resultGL = result.GetMeshGL();
  float minPegZ = std::numeric_limits<float>::max();
  for (size_t run = 0; run < resultGL.runOriginalID.size(); run++) {
    if (resultGL.runOriginalID[run] == pegID) {
      for (size_t t3 = resultGL.runIndex[run]; t3 < resultGL.runIndex[run + 1];
           t3++) {
        const size_t v = resultGL.triVerts[t3];
        minPegZ = std::min(minPegZ,
                           resultGL.vertProperties[v * resultGL.numProp + 2]);
      }
    }
  }
  EXPECT_EQ(minPegZ, 0);

#ifdef MANIFOLD_EXPORT
  if (options.exportModels) ExportMesh("coplanar.glb", resultGL, {});
#endif
}
```

## Dependencies (helpers/subsystems this test needs)

- **`RelatedGL`** — find in `test/test_main.cpp`. If not yet ported to C, port it first.
- **`MeshGL`** — find in `include/manifold/manifold.h`. If not yet ported to C, port it first.


## Steps

1. Read the C++ test above. Understand every line — what objects are created, what operations are performed, what values are asserted.
2. Check if all needed `manifold_*` API functions exist in `src_c/manifold_api.h`. If any are missing, port them from the C++ source first.
3. Write the C test function `test_Properties_Coplanar(void)` in `src_c/test_manifold.c`.
4. Add `RUN_TEST(Properties_Coplanar);` to the appropriate section in `main()`.
5. Build and run: `TEST=Properties_Coplanar ./test_manifold`
6. If it fails — the C **implementation** is wrong, not the test. Debug and fix the implementation.
7. Run full suite: `SKIP_SLOW=1 ./test_manifold` — no regressions.
8. Commit: `test: port Properties_Coplanar from C++`

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

1. This prompt is about ONE test: `Properties_Coplanar`. Do not work on other tests.
2. Never weaken assertions. The C++ expected values are the ground truth.
3. If a C++ method is missing from the C API, port it before writing the test.
4. Run `SKIP_SLOW=1 ./test_manifold` after to verify no regressions.
5. Commit: `test: port Properties_Coplanar from C++`


---

## Completion

When this test passes with no regressions, output exactly:

**Phase complete — test ported and verified**