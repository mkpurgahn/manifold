## Adversary Finding (Iteration 1)

**Requirement:** The spec requires covering **scale** operations. The C++ source has `case TransformType::Scale: cube = cube.Scale({...}); break;`. The porting strategy explicitly says: "Create 5-10 hardcoded CubeOp sequences that exercise: translate+union, rotate+subtract, **scale+union**, chained transforms, empty result." Step 2 also says: "Create 5-10 deterministic CSG sequences covering **translate/rotate/scale** + union/subtract."

**Actual:** The test only defines two transform types: `ty: 0=translate, 1=rotate` (line 2535). The transform dispatch (lines 2571-2574) only handles translate and rotate with no scale branch. None of the 10 ops in the array use a scale transform. The `manifold_scale` API exists and is used elsewhere in the test file (e.g., lines 365, 376, 779).

**Fix:** 
1. Add `ty=2` for scale in the Xform comment and dispatch logic:
   ```c
   if (xf->ty == 0) {
       tmp = manifold_translate(&cur, ...);
   } else if (xf->ty == 1) {
       tmp = manifold_rotate(&cur, ...);
   } else {
       tmp = manifold_scale(&cur, (ManifoldVec3){xf->v[0], xf->v[1], xf->v[2]});
   }
   ```
2. Add at least one op using scale, e.g. replace one of the existing ops or add a new one:
   ```c
   // scale + union
   { .xforms = {{2, {2.0, 0.5, 1.5}}}, .nxforms = 1, .isUnion = 1 },
   ```
3. Add a chained op with scale, e.g. translate+scale+union to cover the "chained transforms" requirement with scale included.
