# Manifold C++ → Pure C Port

> **Goal:** Convert this C++ library to pure C11. The result must compile with `cc -std=c11` with no C++ dependency. This is a **1:1 complete port** — every source file, every function, all functionality.

> **⚠️ FIRST STEP OF EVERY ITERATION:** Build the current state. If it fails, fix it before moving on.

> **⚠️ HARD PROBLEMS ARE WORTH FIXING.** Research and solve each one.

> **⚠️ NEVER DISABLE, SKIP, OR COMMENT OUT A FAILING TEST.** A failing test means the implementation is wrong. Fix the implementation, not the test. There is no such thing as "pre-existing" or "known issue" — if it crashes, debug it and fix the root cause.

---

## The Task

Convert every `.cpp` and `.h` file in `src/` and `include/` to pure C11. No modules skipped. No functions stubbed. All functionality integrated.

### How to Know What to Port

1. List every `.cpp` and `.h` file under `src/` and `include/` — that is the scope.
2. For each file, read it fully and translate every function to C11.
3. If a file exists in the original and has no corresponding C port in `src_c/`, port it.
4. If a function exists in the original and has no corresponding C implementation, port it.

**Do not cherry-pick. Do not skip files you think are unimportant. Port everything.**

### Translation Rules

**This is a mechanical translation, not a redesign.** Translate every C++ construct to its C equivalent. Do NOT change the architecture. Do NOT decide a module is "not needed" and skip it. Do NOT replace lazy evaluation with eager evaluation. Do NOT merge files or remove abstraction layers. If the C++ has it, the C port has it.

1. **Classes → structs + functions.** `Foo::Bar()` becomes `foo_bar(Foo *self)`.
2. **Constructors/destructors → init/destroy functions.**
3. **std::vector → dynamic arrays.** `{T *data; size_t len; size_t cap;}` with helpers.
4. **Templates → concrete types or macros.**
5. **Lambdas → function pointers + void* context.**
6. **RAII → explicit cleanup.** Every allocation has a matching free.
7. **Operator overloading → named functions.**
8. **Exceptions → error codes.**
9. **TBB parallelism → pthreads.** Full parallel implementation, not serial stubs.
10. **namespace → prefix.** Everything gets a `manifold_` prefix.

### No Exceptions

- **Nothing is dropped.** Every file gets ported.
- **Nothing is simplified.** Every function is translated 1:1 to C11.
- TBB parallelism → pthreads. Not "serial for now" — actual parallel implementation.
- External dependencies → C equivalents with the same behavior.
- Bindings, build system, I/O, 2D modules — all of it. Full conversion.

### Tests

- Convert every test case from the C++ test suite (`test/`) to C.
- Count the original `TEST_F` / `TEST_P` / `TEST` cases. The C port must have the same number or more.
- **All tests must pass. If a test crashes, fix the code, not the test. Never comment out or `#if 0` a failing test.**

## Approach

Work file by file, **in order**:

1. Read the C++ source completely
2. Create the equivalent `.c` and `.h` in `src_c/`
3. Translate all logic to C11
4. Ensure it compiles: `cc -std=c11 -c src_c/filename.c -I src_c/`
5. **Port all tests for that file.** Find the corresponding test cases in `test/` and convert them to C. Run them. All must pass.
6. Commit: `port: convert filename.cpp → filename.c`
7. **Only then move to the next file.**

**Do NOT skip to an easier file.** If a file is hard, that's where you stay until it's done. Do not move on leaving broken or incomplete work behind. One file at a time, fully ported, fully tested, then next.

Use the existing C binding API (`bindings/c/`) as reference for naming conventions.

## Validation

The port is done when:

1. Every `.cpp` and `.h` file under `src/` has a corresponding C file in `src_c/`
2. Every function in the original has a C equivalent
3. Every test case from `test/` has been ported and passes
4. The whole thing compiles with `cc -std=c11` and zero C++ dependencies

**Verify by diffing the file lists.** `ls src/*.cpp src/*.h` vs `ls src_c/*.c src_c/*.h`. If anything is missing, port it.

## Commit Format

```
port: convert [filename].cpp → [filename].c

[brief notes on translation decisions]
```

---

## Recursive Improvement

Each iteration: find the next unported file, translate it, compile it, test it, commit it. Keep going until every file is ported and all tests pass.
