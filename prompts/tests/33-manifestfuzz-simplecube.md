# Port Test: `ManifoldFuzz_SimpleCube`

> **Test #33 of 58.** Port this single test to C. Make it pass. Commit. Then move to the next test.

> **⚠️ NEVER WEAKEN A TEST.** The C++ values are ground truth.
> **⚠️ DO NOT SKIP THIS TEST.** If it needs a new C API function, port it. If it needs a subsystem, port the minimum needed.
> **⚠️ DO NOT SPAWN SUB-AGENTS.** Do all work yourself, directly.

---

## Build & Run

```bash
cd /Users/admin/projects/manifold
make
TEST=ManifoldFuzz_SimpleCube ./test_manifold          # run this test alone
SKIP_SLOW=1 ./test_manifold               # verify no regressions
```

---

## C++ Source (`test/manifold_fuzz.cpp`)

This is a **fuzz test** using `fuzztest/fuzztest.h`. It generates random sequences of cube CSG operations (translate, rotate, union, subtract) and asserts `result.Status() == NoError`.

```cpp
void SimpleCube(const std::vector<CubeOp>& inputs) {
  ManifoldParams().intermediateChecks = true;
  ManifoldParams().processOverlaps = false;
  Manifold result;
  for (const auto& input : inputs) {
    auto cube = Manifold::Cube();
    for (const auto& transform : input.transforms) {
      switch (transform.ty) {
        case TransformType::Translate:
          cube = cube.Translate({...});
          break;
        case TransformType::Rotate:
          cube = cube.Rotate(...);
          break;
        case TransformType::Scale:
          cube = cube.Scale({...});
          break;
      }
    }
    if (input.isUnion) result += cube;
    else               result -= cube;
    EXPECT_EQ(result.Status(), Manifold::Error::NoError);
  }
}
FUZZ_TEST(ManifoldFuzz, SimpleCube).WithDomains(CsgDomain);
```

## Porting Strategy

Since we can't use `fuzztest` in C, create a **deterministic version** with a fixed set of representative inputs:

1. Create 5-10 hardcoded `CubeOp` sequences that exercise: translate+union, rotate+subtract, scale+union, chained transforms, empty result.
2. For each sequence: create cube, apply transforms, boolean with result, assert `manifold_status() == MANIFOLD_NO_ERROR`.
3. This tests that the boolean engine doesn't crash or produce invalid results on arbitrary CSG trees.

---

## Steps

1. Write `test_ManifoldFuzz_SimpleCube(void)` in `src_c/test_manifold.c`.
2. Create 5-10 deterministic CSG sequences covering translate/rotate/scale + union/subtract.
3. Assert `manifold_status(&result) == MANIFOLD_NO_ERROR` after each boolean.
4. Add `RUN_TEST(ManifoldFuzz_SimpleCube);` to `main()`.
5. Build and run: `TEST=ManifoldFuzz_SimpleCube ./test_manifold`
6. Run full suite: `SKIP_SLOW=1 ./test_manifold`
7. Commit: `test: port ManifoldFuzz_SimpleCube from C++`

---

## Rules

1. This prompt is about ONE test: `ManifoldFuzz_SimpleCube`. Do not work on other tests.
2. Never weaken assertions. The C++ expected values are the ground truth.
3. If a C++ method is missing from the C API, port it before writing the test.
4. Run `SKIP_SLOW=1 ./test_manifold` after to verify no regressions.
5. Commit: `test: port ManifoldFuzz_SimpleCube from C++`
