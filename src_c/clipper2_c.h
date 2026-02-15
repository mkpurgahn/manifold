#ifndef CLIPPER2_C_H
#define CLIPPER2_C_H

#ifdef __cplusplus
extern "C" {
#endif

// Fill rules matching Clipper2
#define CLIPPER2_FILL_EVEN_ODD 0
#define CLIPPER2_FILL_NON_ZERO 1
#define CLIPPER2_FILL_POSITIVE 2
#define CLIPPER2_FILL_NEGATIVE 3

typedef struct {
  double *verts;     // x,y pairs: [x0,y0,x1,y1,...]
  int *contourSizes; // number of vertices per contour
  int numContours;
  int totalVerts;
} Clipper2Paths;

// Perform a polygon union with the given fill rule.
// Input: contours defined by verts (x,y pairs), contourSizes, numContours.
// Output: resolved contours in result (caller must free with clipper2_free_paths).
int clipper2_union(const double *verts, const int *contourSizes, int numContours,
                   int fillRule, Clipper2Paths *result);

void clipper2_free_paths(Clipper2Paths *p);

#ifdef __cplusplus
}
#endif

#endif // CLIPPER2_C_H
