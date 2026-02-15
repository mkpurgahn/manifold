#include "clipper2_c.h"
#include "clipper2/clipper.h"
#include <cstdlib>

using namespace Clipper2Lib;

extern "C" {

int clipper2_union(const double *verts, const int *contourSizes, int numContours,
                   int fillRule, Clipper2Paths *result) {
  if (!result) return -1;
  result->verts = nullptr;
  result->contourSizes = nullptr;
  result->numContours = 0;
  result->totalVerts = 0;

  PathsD subjects;
  subjects.reserve(numContours);
  int offset = 0;
  for (int c = 0; c < numContours; c++) {
    int n = contourSizes[c];
    PathD path;
    path.reserve(n);
    for (int i = 0; i < n; i++) {
      path.push_back(PointD(verts[2 * (offset + i)], verts[2 * (offset + i) + 1]));
    }
    subjects.push_back(std::move(path));
    offset += n;
  }

  FillRule fr;
  switch (fillRule) {
    case CLIPPER2_FILL_EVEN_ODD: fr = FillRule::EvenOdd; break;
    case CLIPPER2_FILL_NON_ZERO: fr = FillRule::NonZero; break;
    case CLIPPER2_FILL_POSITIVE: fr = FillRule::Positive; break;
    case CLIPPER2_FILL_NEGATIVE: fr = FillRule::Negative; break;
    default: fr = FillRule::Positive; break;
  }

  // Use precision=8 to match C++ manifold's Clipper2 usage
  PathsD solution = Union(subjects, fr, 8);

  int nc = (int)solution.size();
  if (nc == 0) return 0;

  int totalV = 0;
  for (auto &p : solution) totalV += (int)p.size();

  result->numContours = nc;
  result->totalVerts = totalV;
  result->contourSizes = (int *)malloc(nc * sizeof(int));
  result->verts = (double *)malloc(totalV * 2 * sizeof(double));

  int vi = 0;
  for (int c = 0; c < nc; c++) {
    int n = (int)solution[c].size();
    result->contourSizes[c] = n;
    for (int i = 0; i < n; i++) {
      result->verts[2 * vi] = solution[c][i].x;
      result->verts[2 * vi + 1] = solution[c][i].y;
      vi++;
    }
  }

  return 0;
}

void clipper2_free_paths(Clipper2Paths *p) {
  if (!p) return;
  free(p->verts);
  free(p->contourSizes);
  p->verts = nullptr;
  p->contourSizes = nullptr;
  p->numContours = 0;
  p->totalVerts = 0;
}

}  // extern "C"
