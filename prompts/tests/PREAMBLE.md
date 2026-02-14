# Manifold C Port — Test Porting Preamble

> **Include this at the start of every test phase.**

---

## ⚠️ SCOPE: ONE TEST ONLY

You are porting **one specific C++ test** to C. Do not work on any other test. Do not refactor unrelated code. Focus entirely on making this single test pass.

---

## 🔄 WRITE DRY CODE

- If you port a helper function, put it where other tests can reuse it
- If you add a new `manifold_*` API function, add it to `manifold_api.h` AND `manifold_api.c`
- Check if similar helpers already exist before writing new ones

---

## 🐢 SLOW DOWN — BE THOROUGH

1. Read the C++ test completely — understand every operation and assertion
2. Check which `manifold_*` API functions it needs — are they all ported?
3. If any API is missing, port it from the C++ source FIRST
4. Write the C test function
5. Build and run the single test
6. If it fails — the C **implementation** is wrong, not the test. Debug and fix.
7. Run `SKIP_SLOW=1 ./test_manifold` to verify no regressions
8. Commit

---

## ⚠️ MANDATORY: Build-Test Protocol

```bash
cd /Users/admin/projects/manifold
make                                    # must compile clean
TEST=<TestName> ./test_manifold         # must PASS
SKIP_SLOW=1 ./test_manifold             # must show 0 failed
```

**Do NOT claim completion until ALL THREE pass.**

---

## ⛔ NEVER

- **NEVER weaken a test.** Do not increase tolerances, loosen assertions, or skip checks.
- **NEVER disable or comment out a failing test.**
- **NEVER skip this test** because it needs a missing API function — port the function.
- **NEVER spawn sub-agents.** Do all work yourself, directly.

---

## Version Control

Commit after the test passes:
```bash
git add -A && git commit -m "test: port <TestName> from C++"
```

---

## Completion Promise

When the test passes with no regressions, output exactly:

**Phase complete — test ported and verified**
