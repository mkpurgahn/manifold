# Manifold C++ → Pure C Port

> **Goal:** Convert the `elalish/manifold` geometry library from C++ to pure C11. The result should compile with `cc -std=c11` with no C++ dependency. Preserve all core functionality: mesh boolean operations (union, intersection, difference), manifold construction, SDF level set, and mesh repair.

> **⚠️ FIRST STEP OF EVERY ITERATION:** Build the current state with `make -j8` (or the build system in place). If it fails, fix it before moving on.

> **⚠️ HARD PROBLEMS ARE WORTH FIXING.** This is a large mechanical translation with tricky parts. Research and solve each one.

---

## The Task

The `src/` directory contains ~17K lines of C++ and `include/` has ~4K lines of headers. Convert all of it to C11.

### Translation Rules

1. **Classes → structs + functions.** `Manifold::Method()` becomes `manifold_method(Manifold *self)`.
2. **Constructors/destructors → init/destroy functions.** `Manifold()` becomes `manifold_init()`, destructor becomes `manifold_destroy()`.
3. **std::vector → dynamic arrays.** Use a simple `{T *data; size_t len; size_t cap;}` pattern with `_push`, `_resize`, `_free` helpers.
4. **Templates → concrete types or macros.** If a template is used with 2-3 types, write concrete versions. If many types, use macros.
5. **Lambdas → function pointers + void* context.**
6. **RAII → explicit cleanup.** Every allocation must have a matching free on every code path.
7. **Operator overloading → named functions.** `vec3 + vec3` becomes `vec3_add(a, b)`.
8. **Exceptions → error codes.** Return status enum or use an error field on the struct.
9. **TBB parallelism → simple for-loops initially.** Can add pthreads/OpenMP later.
10. **namespace → prefix.** Everything gets a `manifold_` prefix.

### What to Keep

- All mesh boolean algorithms (the core value)
- SDF level set meshing
- Mesh repair and manifold verification
- The existing C binding API (`bindings/c/`) as reference — the port should be API-compatible
- All test cases — convert from C++ gtest to simple assert-based C tests

### What to Drop

- Python/JS/WASM bindings (can be re-added later)
- Assimp mesh I/O (we have our own GLB export)
- CMake build system (replace with a simple Makefile)
- Cross-section 2D support (not needed for our use case)

## File-by-File Approach

Work through the source files one at a time. For each file:

1. Read the C++ source completely
2. Create the equivalent `.c` and `.h` files in a `src_c/` directory
3. Translate all classes, methods, and logic to C11
4. Ensure it compiles: `cc -std=c11 -c src_c/filename.c -I src_c/`
5. Commit: `port: convert filename.cpp → filename.c`

### Suggested Order (dependencies first)

1. `src/vec.h` → math types (vec2, vec3, mat3, etc.)
2. `src/utils.h` → utility functions
3. `src/hashtable.h` → hash table implementation
4. `src/parallel.h` → parallel primitives (start with serial)
5. `src/sort.cpp` → sorting
6. `src/collider.h` / `src/lazy_collider.*` → spatial acceleration
7. `src/tree2d.*` → 2D tree
8. `src/impl.*` → core mesh data structure
9. `src/mesh_fixes.h` → mesh repair
10. `src/constructors.cpp` → primitive constructors (sphere, cube, cylinder, etc.)
11. `src/boolean3.*` → boolean algorithm core
12. `src/boolean_result.cpp` → boolean result processing
13. `src/csg_tree.*` → CSG tree evaluation
14. `src/manifold.cpp` → main API
15. `src/properties.cpp` → vertex properties
16. `src/sdf.cpp` → SDF level set
17. `src/smoothing.cpp` → mesh smoothing
18. `src/subdivision.cpp` → subdivision
19. `src/polygon.cpp` → polygon triangulation
20. `src/quickhull.*` → convex hull
21. `src/edge_op.cpp` / `src/face_op.cpp` → mesh operations
22. `src/minkowski.cpp` → Minkowski sum

## Validation

After all files are ported:

```bash
make test
```

Every test from the original C++ test suite must pass in the C version. The boolean operations must produce identical output (bit-for-bit manifold meshes).

## Commit Format

```
port: convert [filename].cpp → [filename].c

[brief notes on translation decisions]
```

---

## Recursive Improvement

Each iteration: pick the next file, translate it, compile it, commit it. Keep going until the entire library compiles as pure C11 and all tests pass.
