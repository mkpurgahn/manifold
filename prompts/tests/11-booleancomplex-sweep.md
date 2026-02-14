# Port Test: `BooleanComplex_Sweep`

> **Test #11 of 58.** Port this single test to C. Make it pass. Commit. Then move to the next test.

> **⚠️ NEVER WEAKEN A TEST.** The C++ values are ground truth.
> **⚠️ DO NOT SKIP THIS TEST.** If it needs a new C API function, port it. If it needs a subsystem, port the minimum needed.
> **⚠️ DO NOT SPAWN SUB-AGENTS.** Do all work yourself, directly.

---

## Build & Run

```bash
cd /Users/admin/projects/manifold
make
TEST=BooleanComplex_Sweep ./test_manifold          # run this test alone
SKIP_SLOW=1 ./test_manifold               # verify no regressions
```

---

## C++ Source (`test/boolean_complex_test.cpp`)

```cpp
TEST(BooleanComplex, Sweep) {
  ManifoldParamGuard guard;
  ManifoldParams().processOverlaps = true;

  // generate the minimum equivalent positive angle
  auto minPosAngle = [](double angle) {
    double div = angle / kTwoPi;
    double wholeDiv = floor(div);
    return angle - wholeDiv * kTwoPi;
  };

  // calculate determinant
  auto det = [](vec2 v1, vec2 v2) { return v1.x * v2.y - v1.y * v2.x; };

  // generate sweep profile
  auto generateProfile = []() {
    double filletRadius = 2.5;
    double filletWidth = 5;
    int numberOfArcPoints = 10;
    vec2 arcCenterPoint = vec2(filletWidth - filletRadius, filletRadius);
    std::vector<vec2> arcPoints;

    for (int i = 0; i < numberOfArcPoints; i++) {
      double angle = i * kPi / numberOfArcPoints;
      double y = arcCenterPoint.y - cos(angle) * filletRadius;
      double x = arcCenterPoint.x + sin(angle) * filletRadius;
      arcPoints.push_back(vec2(x, y));
    }

    std::vector<vec2> profile;
    profile.push_back(vec2(0, 0));
    profile.push_back(vec2(filletWidth - filletRadius, 0));
    for (int i = 0; i < numberOfArcPoints; i++) {
      profile.push_back(arcPoints[i]);
    }
    profile.push_back(vec2(0, filletWidth));

    CrossSection profileCrossSection = CrossSection(profile);
    return profileCrossSection;
  };

  CrossSection profile = generateProfile();

  auto partialRevolve = [minPosAngle, profile](double startAngle,
                                               double endAngle,
                                               int nSegmentsPerRotation) {
    double posEndAngle = minPosAngle(endAngle);
    double totalAngle = 0;
    if (startAngle < 0 && endAngle < 0 && startAngle < endAngle) {
      totalAngle = endAngle - startAngle;
    } else {
      totalAngle = posEndAngle - startAngle;
    }

    int nSegments = ceil(totalAngle / kTwoPi * nSegmentsPerRotation + 1);
    if (nSegments < 2) {
      nSegments = 2;
    }

    double angleStep = totalAngle / (nSegments - 1);
    auto warpFunc = [nSegments, angleStep, startAngle](vec3& vertex) {
      double zIndex = nSegments - 1 - vertex.z;
      double angle = zIndex * angleStep + startAngle;

      // transform
      vertex.z = vertex.y;
      vertex.y = vertex.x * sin(angle);
      vertex.x = vertex.x * cos(angle);
    };

    return Manifold::Extrude(profile.ToPolygons(), nSegments - 1, nSegments - 2)
        .Warp(warpFunc);
  };

  auto cutterPrimitives = [det, partialRevolve, profile](vec2 p1, vec2 p2,
                                                         vec2 p3) {
    vec2 diff = p2 - p1;
    vec2 vec1 = p1 - p2;
    vec2 vec2 = p3 - p2;
    double determinant = det(vec1, vec2);

    double startAngle = atan2(vec1.x, -vec1.y);
    double endAngle = atan2(-vec2.x, vec2.y);

    Manifold round =
        partialRevolve(startAngle, endAngle, 20).Translate(vec3(p2.x, p2.y, 0));

    double distance = sqrt(diff.x * diff.x + diff.y * diff.y);
    double angle = atan2(diff.y, diff.x);
    Manifold extrusionPrimitive =
        Manifold::Extrude(profile.ToPolygons(), distance)
            .Rotate(90, 0, -90)
            .Translate(vec3(distance, 0, 0))
            .Rotate(0, 0, angle * 180 / kPi)
            .Translate(vec3(p1.x, p1.y, 0));

    std::vector<Manifold> result;

    if (determinant < 0) {
      result.push_back(round);
      result.push_back(extrusionPrimitive);
    } else {
      result.push_back(extrusionPrimitive);
    }

    return result;
  };

  auto scalePath = [](std::vector<vec2> path, double scale) {
    std::vector<vec2> newPath;
    for (vec2 point : path) {
      newPath.push_back(scale * point);
    }
    return newPath;
  };

  std::vector<vec2> pathPoints = {{-21.707751473606564, 10.04202769267855},
                                  {-21.840846948218307, 9.535474475521578},
                                  {-21.940954413815387, 9.048287386171369},
                                  {-22.005569458385835, 8.587741145234093},
                                  {-22.032187669917704, 8.16111047331591},
                                  {-22.022356960178296, 7.755456475810721},
                                  {-21.9823319178086, 7.356408291345673},
                                  {-21.91208498286602, 6.964505631629036},
                                  {-21.811437268778267, 6.579251589515578},
                                  {-21.68020988897306, 6.200149257860059},
                                  {-21.51822395687812, 5.82670172951726},
                                  {-21.254086890521585, 5.336709200579579},
                                  {-21.01963533308061, 4.974523796623895},
                                  {-20.658228140926262, 4.497743844638198},
                                  {-20.350337020134603, 4.144115181723373},
                                  {-19.9542029967, 3.7276501717684054},
                                  {-20.6969129296381, 3.110639833377638},
                                  {-21.026318197401537, 2.793796378245609},
                                  {-21.454710558515973, 2.3418076758544806},
                                  {-21.735944543382722, 2.014266362004704},
                                  {-21.958999535447845, 1.7205197644485681},
                                  {-22.170169612837164, 1.3912359628761894},
                                  {-22.376940405634056, 1.0213515348242117},
                                  {-22.62545385249271, 0.507889651991388},
                                  {-22.77620002102207, 0.13973666928102288},
                                  {-22.8689989640578, -0.135962138067232},
                                  {-22.974385239894364, -0.5322784681448909},
                                  {-23.05966775687304, -0.9551466941218276},
                                  {-23.102914137841445, -1.2774406685179822},
                                  {-23.14134824916783, -1.8152432718003662},
                                  {-23.152085124298473, -2.241104719188421},
                                  {-23.121576743285054, -2.976332948223073},
                                  {-23.020491352156856, -3.6736813934577914},
                                  {-22.843552165110886, -4.364810769710428},
                                  {-22.60334013490563, -5.033012850282157},
                                  {-22.305015243491663, -5.67461444847819},
                                  {-21.942709324216615, -6.330962778427178},
                                  {-21.648491707764062, -6.799117771996025},
                                  {-21.15330508818782, -7.496539096945377},
                                  {-21.10687739725184, -7.656798276710632},
                                  {-21.01253055778545, -8.364144493707382},
                                  {-20.923211927856293, -8.782280691344269},
                                  {-20.771325204062215, -9.258087073404687},
                                  {-20.554404009259198, -9.72613360625344},
                                  {-20.384050989017144, -9.985885743112847},
                                  {-20.134404839253612, -10.263023004626703},
                                  {-19.756998832033442, -10.613109670467736},
                                  {-18.83161393127597, -15.68768837402245},
                                  {-19.155593463785983, -17.65410871259763},
                                  {-17.930304365744544, -19.005810988385562},
                                  {-16.893408103100064, -19.50558228186199},
                                  {-16.27514960757635, -19.8288501942628},
                                  {-15.183033464853374, -20.47781203017123},
                                  {-14.906850387751492, -20.693472553142833},
                                  {-14.585198957236713, -21.015257964547136},
                                  {-11.013839210807205, -34.70394287828328},
                                  {-8.79778020674896, -36.17434400175442},
                                  {-7.850491148257242, -36.48835987119041},
                                  {-6.982497182376991, -36.74546968896842},
                                  {-6.6361688522576, -36.81653354539242},
                                  {-6.0701080598244035, -36.964332993204},
                                  {-5.472439187922815, -37.08824838436714},
                                  {-4.802871164820756, -37.20127157090685},
                                  {-3.6605994233344745, -37.34427653957914},
                                  {-1.7314396363710867, -37.46415201430501},
                                  {-0.7021130485987349, -37.5},
                                  {0.01918509410483974, -37.49359541901704},
                                  {1.2107837650065625, -37.45093992812552},
                                  {3.375529069920302, 32.21823383780513},
                                  {1.9041980552754056, 32.89839543047101},
                                  {1.4107184651094313, 33.16556804736585},
                                  {1.1315552947605065, 33.34344755450097},
                                  {0.8882931135353977, 33.52377699790175},
                                  {0.6775397019893341, 33.708817857198056},
                                  {0.49590284067753837, 33.900831612019715},
                                  {0.2291596803839543, 34.27380625039597},
                                  {0.03901816126171688, 34.66402375075138},
                                  {-0.02952797094655369, 34.8933309389416},
                                  {-0.0561772851849209, 35.044928843125824},
                                  {-0.067490756643705, 35.27129875796868},
                                  {-0.05587453990569748, 35.42204271802184},
                                  {0.013497378362074697, 35.72471438137191},
                                  {0.07132375113026912, 35.877348797053145},
                                  {0.18708820875448923, 36.108917464873215},
                                  {0.39580614140195136, 36.424415957998825},
                                  {0.8433687814267005, 36.964365016108914},
                                  {0.7078417131710703, 37.172455373435916},
                                  {0.5992848016685662, 37.27482757003058},
                                  {0.40594743344375905, 37.36664006036318},
                                  {0.1397973410299913, 37.434752779117005}};

  int numPoints = pathPoints.size();
  pathPoints = scalePath(pathPoints, 0.9);

  std::vector<Manifold> result;

  for (int i = 0; i < numPoints; i++) {
    std::vector<Manifold> primitives =
        cutterPrimitives(pathPoints[i], pathPoints[(i + 1) % numPoints],
                         pathPoints[(i + 2) % numPoints]);

    for (Manifold primitive : primitives) {
      result.push_back(primitive);
    }
  }

  // all primitives should be valid
  for (Manifold primitive : result) {
    if (primitive.Volume() < 0) {
      std::cerr << "INVALID PRIMITIVE" << std::endl;
    }
  }

  Manifold shape = Manifold::BatchBoolean(result, OpType::Add);

  EXPECT_NEAR(shape.Volume(), 3757, 1);
#ifdef MANIFOLD_EXPORT
  if (options.exportModels) ExportMesh("unionError.glb", shape.GetMeshGL(), {});
#endif
}
```

## Dependencies (helpers/subsystems this test needs)

- **`ManifoldParams`** — find in `include/manifold/manifold.h`. If not yet ported to C, port it first.
- **`CrossSection`** — find in `include/manifold/cross_section.h`. If not yet ported to C, port it first.
- **`MeshGL`** — find in `include/manifold/manifold.h`. If not yet ported to C, port it first.


## Steps

1. Read the C++ test above. Understand every line — what objects are created, what operations are performed, what values are asserted.
2. Check if all needed `manifold_*` API functions exist in `src_c/manifold_api.h`. If any are missing, port them from the C++ source first.
3. Write the C test function `test_BooleanComplex_Sweep(void)` in `src_c/test_manifold.c`.
4. Add `RUN_TEST(BooleanComplex_Sweep);` to the appropriate section in `main()`.
5. Build and run: `TEST=BooleanComplex_Sweep ./test_manifold`
6. If it fails — the C **implementation** is wrong, not the test. Debug and fix the implementation.
7. Run full suite: `SKIP_SLOW=1 ./test_manifold` — no regressions.
8. Commit: `test: port BooleanComplex_Sweep from C++`

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

1. This prompt is about ONE test: `BooleanComplex_Sweep`. Do not work on other tests.
2. Never weaken assertions. The C++ expected values are the ground truth.
3. If a C++ method is missing from the C API, port it before writing the test.
4. Run `SKIP_SLOW=1 ./test_manifold` after to verify no regressions.
5. Commit: `test: port BooleanComplex_Sweep from C++`
