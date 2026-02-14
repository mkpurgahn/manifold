# Port Test: `PolygonFuzz_TriangulationNoCrash`

> **Read `prompts/tests/PREAMBLE.md` first.**

> **Test #57 of 58.** Port this single test to C. Make it pass. Commit. Then move to the next test.

> **⚠️ NEVER WEAKEN A TEST.** The C++ values are ground truth.
> **⚠️ DO NOT SKIP THIS TEST.**
> **⚠️ DO NOT SPAWN SUB-AGENTS.** Do all work yourself, directly.

---

## Build & Run

```bash
cd /Users/admin/projects/manifold
make
TEST=PolygonFuzz_TriangulationNoCrash ./test_manifold
SKIP_SLOW=1 ./test_manifold
```

---

## C++ Source (`test/polygon_fuzz.cpp`)

This is a **fuzz test**. It feeds random polygons into `Triangulate()` and asserts it doesn't crash.

```cpp
void TriangulationNoCrash(
    std::vector<std::vector<std::pair<float, float>>> input, float precision) {
  if (precision < 0) precision = -1;
  manifold::ManifoldParams().intermediateChecks = true;
  manifold::Polygons polys;
  for (const auto& simplePoly : input) {
    polys.emplace_back();
    for (const auto& p : simplePoly) {
      polys.back().emplace_back(p.first, p.second);
    }
  }
  manifold::Triangulate(polys, precision);
  // must not crash
}
FUZZ_TEST(PolygonFuzz, TriangulationNoCrash)
    .WithDomains(PolygonDomain, InRange<float>(-1.0, 0.1))
    .WithSeeds(SeedProvider);  // loads from test/polygons/polygon_corpus.txt
```

## Porting Strategy

The C++ test loads seed data from `test/polygons/polygon_corpus.txt`. Port as a deterministic test:

1. Load the seed corpus from `test/polygons/polygon_corpus.txt` (or hardcode a subset of 5-10 representative polygon sets).
2. For each: build polygon array, call `manifold_triangulate()`, assert it returns without error.
3. This tests that the triangulation doesn't crash on various polygon inputs.

---

## Steps

1. Check if `manifold_triangulate` exists in `manifold_api.h`. If not, port it.
2. Write `test_PolygonFuzz_TriangulationNoCrash(void)` using seed polygon data.
3. Add `RUN_TEST(PolygonFuzz_TriangulationNoCrash);` to `main()`.
4. Build and run: `TEST=PolygonFuzz_TriangulationNoCrash ./test_manifold`
5. Run full suite: `SKIP_SLOW=1 ./test_manifold`
6. Commit: `test: port PolygonFuzz_TriangulationNoCrash from C++`

---

## Rules

1. This prompt is about ONE test: `PolygonFuzz_TriangulationNoCrash`. Do not work on other tests.
2. Never weaken assertions.
3. Run `SKIP_SLOW=1 ./test_manifold` after to verify no regressions.
4. Commit: `test: port PolygonFuzz_TriangulationNoCrash from C++`


---

## Completion

When this test passes with no regressions, output exactly:

**Phase complete — test ported and verified**