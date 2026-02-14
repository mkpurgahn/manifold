# Port Test: `PolygonFuzz_TriangulationNoCrashRounded`

> **Read `prompts/tests/PREAMBLE.md` first.**

> **Test #58 of 58.** Port this single test to C. Make it pass. Commit. Then move to the next test.

> **⚠️ NEVER WEAKEN A TEST.** The C++ values are ground truth.
> **⚠️ DO NOT SKIP THIS TEST.**
> **⚠️ DO NOT SPAWN SUB-AGENTS.** Do all work yourself, directly.

---

## Build & Run

```bash
cd /Users/admin/projects/manifold
make
TEST=PolygonFuzz_TriangulationNoCrashRounded ./test_manifold
SKIP_SLOW=1 ./test_manifold
```

---

## C++ Source (`test/polygon_fuzz.cpp`)

Same as `TriangulationNoCrash` but rounds all input coordinates to integers first:

```cpp
void TriangulationNoCrashRounded(
    std::vector<std::vector<std::pair<float, float>>> input, float precision) {
  // rounds all coords to nearest integer, then calls TriangulationNoCrash
  TriangulationNoCrash(std::move(input), precision);
}
```

## Porting Strategy

Same as test #57 but round all polygon vertex coordinates to nearest integer before triangulating.

---

## Steps

1. Write `test_PolygonFuzz_TriangulationNoCrashRounded(void)` — same as #57 but `roundf()` all coords.
2. Add `RUN_TEST(PolygonFuzz_TriangulationNoCrashRounded);` to `main()`.
3. Build and run: `TEST=PolygonFuzz_TriangulationNoCrashRounded ./test_manifold`
4. Run full suite: `SKIP_SLOW=1 ./test_manifold`
5. Commit: `test: port PolygonFuzz_TriangulationNoCrashRounded from C++`

---

## Rules

1. This prompt is about ONE test: `PolygonFuzz_TriangulationNoCrashRounded`. Do not work on other tests.
2. Never weaken assertions.
3. Run `SKIP_SLOW=1 ./test_manifold` after to verify no regressions.
4. Commit: `test: port PolygonFuzz_TriangulationNoCrashRounded from C++`


---

## Completion

When this test passes with no regressions, output exactly:

**Phase complete — test ported and verified**