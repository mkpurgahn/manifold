// Copyright 2024 The Manifold Authors.
// SPDX-License-Identifier: Apache-2.0
//
// QuickHull convex hull algorithm - ported from quickhull.cpp.
// Derived from the public domain work of Antti Kuukka.

#include "manifold_hull.h"
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>

// ─── Dynamic int array ──────────────────────────────────────────────────────
typedef struct { int *data; size_t len, cap; } IntVec;
static void iv_init(IntVec *v) { v->data = NULL; v->len = v->cap = 0; }
static void iv_free(IntVec *v) { free(v->data); v->data = NULL; v->len = v->cap = 0; }
static void iv_push(IntVec *v, int val) {
  if (v->len >= v->cap) { v->cap = v->cap ? v->cap * 2 : 16; v->data = (int*)realloc(v->data, v->cap * sizeof(int)); }
  v->data[v->len++] = val;
}
static void iv_clear(IntVec *v) { v->len = 0; }

// ─── Dynamic size_t array ───────────────────────────────────────────────────
typedef struct { size_t *data; size_t len, cap; } SzVec;
static void sv_init(SzVec *v) { v->data = NULL; v->len = v->cap = 0; }
static void sv_free(SzVec *v) { free(v->data); v->data = NULL; v->len = v->cap = 0; }
static void sv_push(SzVec *v, size_t val) {
  if (v->len >= v->cap) { v->cap = v->cap ? v->cap * 2 : 16; v->data = (size_t*)realloc(v->data, v->cap * sizeof(size_t)); }
  v->data[v->len++] = val;
}
static void sv_clear(SzVec *v) { v->len = 0; }

// ─── Plane ──────────────────────────────────────────────────────────────────
typedef struct { ManifoldVec3 N; double D; double sqrNLength; } QHPlane;

static QHPlane qh_plane(ManifoldVec3 N, ManifoldVec3 P) {
  QHPlane p;
  p.N = N;
  p.D = -vec3_dot(N, P);
  p.sqrNLength = vec3_dot(N, N);
  return p;
}

static double qh_signed_dist(ManifoldVec3 v, const QHPlane *p) {
  return vec3_dot(p->N, v) + p->D;
}

// ─── Triangle normal (normalized) ───────────────────────────────────────────
static ManifoldVec3 qh_tri_normal(ManifoldVec3 a, ManifoldVec3 b, ManifoldVec3 c) {
  ManifoldVec3 n = vec3_cross(vec3_sub(a, c), vec3_sub(b, c));
  double len = vec3_length(n);
  if (len > 1e-30) n = vec3_scale(n, 1.0 / len);
  return n;
}

// ─── MeshBuilder Face ───────────────────────────────────────────────────────
typedef struct {
  int he;                    // index of first halfedge
  QHPlane P;
  double mostDistantPointDist;
  size_t mostDistantPoint;
  size_t visibilityCheckedOnIteration;
  uint8_t isVisibleFaceOnCurrentIteration;
  uint8_t inFaceStack;
  uint8_t horizonEdgesOnCurrentIteration;
  // Points on positive side: indices into original vertex array
  SzVec pointsOnPositiveSide;
  bool hasPoints; // whether pointsOnPositiveSide is active
} QHFace;

static void qhface_init(QHFace *f, int he) {
  f->he = he;
  memset(&f->P, 0, sizeof(f->P));
  f->mostDistantPointDist = 0;
  f->mostDistantPoint = 0;
  f->visibilityCheckedOnIteration = 0;
  f->isVisibleFaceOnCurrentIteration = 0;
  f->inFaceStack = 0;
  f->horizonEdgesOnCurrentIteration = 0;
  sv_init(&f->pointsOnPositiveSide);
  f->hasPoints = false;
}

static bool qhface_is_disabled(const QHFace *f) { return f->he == -1; }
static void qhface_disable(QHFace *f) { f->he = -1; }

// ─── MeshBuilder ────────────────────────────────────────────────────────────
// Halfedge: {startVert, endVert, pairedHalfedge}
typedef struct { int startVert, endVert, pairedHalfedge; } QHHalfedge;

typedef struct {
  QHFace *faces;       size_t faceLen, faceCap;
  QHHalfedge *he;      size_t heLen, heCap;
  int *heToFace;       // halfedge -> face index
  int *heNext;         // halfedge -> next halfedge
  SzVec disabledFaces;
  SzVec disabledHE;
} QHMesh;

static void qhm_init(QHMesh *m) {
  m->faces = NULL; m->faceLen = m->faceCap = 0;
  m->he = NULL; m->heLen = m->heCap = 0;
  m->heToFace = NULL; m->heNext = NULL;
  sv_init(&m->disabledFaces);
  sv_init(&m->disabledHE);
}

static void qhm_free(QHMesh *m) {
  for (size_t i = 0; i < m->faceLen; i++) sv_free(&m->faces[i].pointsOnPositiveSide);
  free(m->faces); free(m->he); free(m->heToFace); free(m->heNext);
  sv_free(&m->disabledFaces); sv_free(&m->disabledHE);
}

static void qhm_ensure_he(QHMesh *m, size_t needed) {
  if (needed > m->heCap) {
    size_t nc = m->heCap ? m->heCap : 16;
    while (nc < needed) nc *= 2;
    m->he = (QHHalfedge*)realloc(m->he, nc * sizeof(QHHalfedge));
    m->heToFace = (int*)realloc(m->heToFace, nc * sizeof(int));
    m->heNext = (int*)realloc(m->heNext, nc * sizeof(int));
    m->heCap = nc;
  }
}

static void qhm_ensure_face(QHMesh *m, size_t needed) {
  if (needed > m->faceCap) {
    size_t nc = m->faceCap ? m->faceCap : 8;
    while (nc < needed) nc *= 2;
    m->faces = (QHFace*)realloc(m->faces, nc * sizeof(QHFace));
    m->faceCap = nc;
  }
}

static size_t qhm_add_face(QHMesh *m) {
  if (m->disabledFaces.len > 0) {
    size_t idx = m->disabledFaces.data[--m->disabledFaces.len];
    QHFace *f = &m->faces[idx];
    f->mostDistantPointDist = 0;
    f->mostDistantPoint = 0;
    f->isVisibleFaceOnCurrentIteration = 0;
    f->inFaceStack = 0;
    f->horizonEdgesOnCurrentIteration = 0;
    f->visibilityCheckedOnIteration = 0;
    return idx;
  }
  qhm_ensure_face(m, m->faceLen + 1);
  qhface_init(&m->faces[m->faceLen], -1);
  return m->faceLen++;
}

static size_t qhm_add_he(QHMesh *m) {
  if (m->disabledHE.len > 0) {
    size_t idx = m->disabledHE.data[--m->disabledHE.len];
    return idx;
  }
  qhm_ensure_he(m, m->heLen + 1);
  memset(&m->he[m->heLen], 0, sizeof(QHHalfedge));
  m->heToFace[m->heLen] = 0;
  m->heNext[m->heLen] = 0;
  return m->heLen++;
}

static void qhm_disable_face(QHMesh *m, size_t fi) {
  qhface_disable(&m->faces[fi]);
  sv_push(&m->disabledFaces, fi);
}

static void qhm_disable_he(QHMesh *m, size_t hi) {
  m->he[hi].pairedHalfedge = -1;
  sv_push(&m->disabledHE, hi);
}

static void qhm_get_face_he(const QHMesh *m, const QHFace *f, int out[3]) {
  out[0] = f->he;
  out[1] = m->heNext[out[0]];
  out[2] = m->heNext[out[1]];
}

static void qhm_get_face_verts(const QHMesh *m, const QHFace *f, int out[3]) {
  int hes[3]; qhm_get_face_he(m, f, hes);
  out[0] = m->he[hes[0]].endVert;
  out[1] = m->he[hes[1]].endVert;
  out[2] = m->he[hes[2]].endVert;
}

// Setup initial tetrahedron ABCD. dot(AB, normal(ABC)) should be negative.
static void qhm_setup(QHMesh *m, int a, int b, int c, int d) {
  m->faceLen = 0; m->heLen = 0;
  m->disabledFaces.len = 0; m->disabledHE.len = 0;

  qhm_ensure_he(m, 12);
  qhm_ensure_face(m, 4);

  // 12 halfedges for a tetrahedron
  #define ADDHE(sv, ev, pair, face, next) do { \
    m->he[m->heLen] = (QHHalfedge){0, ev, pair}; \
    m->heToFace[m->heLen] = face; \
    m->heNext[m->heLen] = next; \
    m->heLen++; \
  } while(0)

  // Face 0: ABC  (halfedges 0,1,2)
  ADDHE(0, b, 6,  0, 1);   // 0: AB paired=6(BA)
  ADDHE(0, c, 9,  0, 2);   // 1: BC paired=9(CB)
  ADDHE(0, a, 3,  0, 0);   // 2: CA paired=3(AC)
  // Face 1: ACD  (halfedges 3,4,5)
  ADDHE(0, c, 2,  1, 4);   // 3: AC paired=2(CA)
  ADDHE(0, d, 11, 1, 5);   // 4: CD paired=11(DC)
  ADDHE(0, a, 7,  1, 3);   // 5: DA paired=7(AD)
  // Face 2: ADB  (halfedges 6,7,8)
  ADDHE(0, a, 0,  2, 7);   // 6: BA paired=0(AB)
  ADDHE(0, d, 5,  2, 8);   // 7: AD paired=5(DA)
  ADDHE(0, b, 10, 2, 6);   // 8: DB paired=10(BD)
  // Face 3: BDC  (halfedges 9,10,11)
  ADDHE(0, b, 1,  3, 10);  // 9: CB paired=1(BC)
  ADDHE(0, d, 8,  3, 11);  // 10: BD paired=8(DB)
  ADDHE(0, c, 4,  3, 9);   // 11: DC paired=4(CD)

  #undef ADDHE

  for (int i = 0; i < 4; i++) {
    qhface_init(&m->faces[i], i * 3);
    m->faceLen++;
  }
}

// ─── QuickHull algorithm state ──────────────────────────────────────────────
typedef struct {
  int faceIndex;
  int enteredFromHalfedge;
} QHFaceData;

typedef struct {
  double epsilon, epsilonSq, scale;
  bool planar;
  const ManifoldVec3 *verts;
  size_t vertCount;
  QHMesh mesh;
  size_t extremes[6];

  // Temporaries
  SzVec newFaceIndices, newHEIndices, visibleFaces, horizonEdges;
  QHFaceData *pvfStack; size_t pvfLen, pvfCap;
  IntVec faceList;
} QHState;

static void qhs_init(QHState *s, const ManifoldVec3 *pts, size_t n) {
  memset(s, 0, sizeof(*s));
  s->verts = pts;
  s->vertCount = n;
  qhm_init(&s->mesh);
  sv_init(&s->newFaceIndices); sv_init(&s->newHEIndices);
  sv_init(&s->visibleFaces); sv_init(&s->horizonEdges);
  s->pvfStack = NULL; s->pvfLen = s->pvfCap = 0;
  iv_init(&s->faceList);
}

static void qhs_free(QHState *s) {
  qhm_free(&s->mesh);
  sv_free(&s->newFaceIndices); sv_free(&s->newHEIndices);
  sv_free(&s->visibleFaces); sv_free(&s->horizonEdges);
  free(s->pvfStack);
  iv_free(&s->faceList);
}

static void pvf_push(QHState *s, QHFaceData fd) {
  if (s->pvfLen >= s->pvfCap) {
    s->pvfCap = s->pvfCap ? s->pvfCap * 2 : 32;
    s->pvfStack = (QHFaceData*)realloc(s->pvfStack, s->pvfCap * sizeof(QHFaceData));
  }
  s->pvfStack[s->pvfLen++] = fd;
}

static void qhs_get_extremes(QHState *s) {
  for (int i = 0; i < 6; i++) s->extremes[i] = 0;
  double ev[6];
  ev[0] = ev[1] = s->verts[0].x;
  ev[2] = ev[3] = s->verts[0].y;
  ev[4] = ev[5] = s->verts[0].z;
  for (size_t i = 1; i < s->vertCount; i++) {
    ManifoldVec3 p = s->verts[i];
    if (p.x > ev[0]) { ev[0] = p.x; s->extremes[0] = i; }
    else if (p.x < ev[1]) { ev[1] = p.x; s->extremes[1] = i; }
    if (p.y > ev[2]) { ev[2] = p.y; s->extremes[2] = i; }
    else if (p.y < ev[3]) { ev[3] = p.y; s->extremes[3] = i; }
    if (p.z > ev[4]) { ev[4] = p.z; s->extremes[4] = i; }
    else if (p.z < ev[5]) { ev[5] = p.z; s->extremes[5] = i; }
  }
}

static double qhs_get_scale(QHState *s) {
  double sc = 0;
  for (int i = 0; i < 6; i++) {
    const double *v = (const double *)&s->verts[s->extremes[i]];
    double a = fabs(v[i / 2]);
    if (a > sc) sc = a;
  }
  return sc;
}

static bool qhs_add_point_to_face(QHState *s, QHFace *f, size_t ptIdx) {
  double D = qh_signed_dist(s->verts[ptIdx], &f->P);
  if (D > 0 && D * D > s->epsilonSq * f->P.sqrNLength) {
    if (!f->hasPoints) {
      sv_clear(&f->pointsOnPositiveSide);
      f->hasPoints = true;
    }
    sv_push(&f->pointsOnPositiveSide, ptIdx);
    if (D > f->mostDistantPointDist) {
      f->mostDistantPointDist = D;
      f->mostDistantPoint = ptIdx;
    }
    return true;
  }
  return false;
}

static bool qhs_reorder_horizon(QHState *s) {
  SzVec *he = &s->horizonEdges;
  size_t n = he->len;
  for (size_t i = 0; i + 1 < n; i++) {
    int endV = s->mesh.he[he->data[i]].endVert;
    bool found = false;
    for (size_t j = i + 1; j < n; j++) {
      int beginV = s->mesh.he[s->mesh.he[he->data[j]].pairedHalfedge].endVert;
      if (beginV == endV) {
        size_t tmp = he->data[i + 1]; he->data[i + 1] = he->data[j]; he->data[j] = tmp;
        found = true;
        break;
      }
    }
    if (!found) return false;
  }
  return true;
}

static void qhs_setup_initial_tet(QHState *s) {
  size_t vc = s->vertCount;

  // Handle <= 4 points: create degenerate tetrahedron (matching C++)
  if (vc <= 4) {
    if (vc < 4) {
      // Pad to 4 by duplicating the last point
      s->planar = true;
      return;
    }
    size_t v[4] = {0, 1, 2, 3};
    ManifoldVec3 N = qh_tri_normal(s->verts[v[0]], s->verts[v[1]], s->verts[v[2]]);
    QHPlane tp = qh_plane(N, s->verts[v[0]]);
    double d = fabs(qh_signed_dist(s->verts[v[3]], &tp));
    if (d < s->epsilon || tp.sqrNLength < 1e-30) {
      s->planar = true; return;
    }
    if (qh_signed_dist(s->verts[v[3]], &tp) >= 0) {
      size_t tmp = v[0]; v[0] = v[1]; v[1] = tmp;
    }
    qhm_setup(&s->mesh, (int)v[0], (int)v[1], (int)v[2], (int)v[3]);
    return;
  }

  // Find two most distant extreme points
  double maxD = s->epsilonSq;
  size_t sp0 = 0, sp1 = 0;
  for (int i = 0; i < 6; i++) {
    for (int j = i + 1; j < 6; j++) {
      ManifoldVec3 diff = vec3_sub(s->verts[s->extremes[i]], s->verts[s->extremes[j]]);
      double d = vec3_dot(diff, diff);
      if (d > maxD) { maxD = d; sp0 = s->extremes[i]; sp1 = s->extremes[j]; }
    }
  }
  if (maxD == s->epsilonSq) {
    qhm_setup(&s->mesh, 0, 1 % (int)vc, 2 % (int)vc, 3 % (int)vc);
    return;
  }

  // Find point farthest from line sp0-sp1
  ManifoldVec3 rayV = vec3_sub(s->verts[sp1], s->verts[sp0]);
  double rayInvLen2 = 1.0 / vec3_dot(rayV, rayV);
  maxD = s->epsilonSq;
  size_t maxI = SIZE_MAX;
  for (size_t i = 0; i < vc; i++) {
    ManifoldVec3 sv = vec3_sub(s->verts[i], s->verts[sp0]);
    double t = vec3_dot(sv, rayV);
    double d2 = vec3_dot(sv, sv) - t * t * rayInvLen2;
    if (d2 > maxD) { maxD = d2; maxI = i; }
  }
  if (maxI == SIZE_MAX) {
    qhm_setup(&s->mesh, (int)sp0, (int)sp1, 0, 1);
    return;
  }

  size_t baseTri[3] = {sp0, sp1, maxI};

  // Find 4th vertex farthest from base triangle plane
  ManifoldVec3 N = qh_tri_normal(s->verts[baseTri[0]], s->verts[baseTri[1]], s->verts[baseTri[2]]);
  QHPlane triPlane = qh_plane(N, s->verts[baseTri[0]]);
  maxD = s->epsilon;
  maxI = 0;
  for (size_t i = 0; i < vc; i++) {
    double d = fabs(qh_signed_dist(s->verts[i], &triPlane));
    if (d > maxD) { maxD = d; maxI = i; }
  }
  if (maxD == s->epsilon) {
    // All points are coplanar - set planar flag
    s->planar = true;
    return;
  }

  // Enforce CCW orientation
  if (qh_signed_dist(s->verts[maxI], &triPlane) >= 0) {
    size_t tmp = baseTri[0]; baseTri[0] = baseTri[1]; baseTri[1] = tmp;
  }

  qhm_setup(&s->mesh, (int)baseTri[0], (int)baseTri[1], (int)baseTri[2], (int)maxI);

  // Compute planes for initial faces
  for (size_t fi = 0; fi < 4; fi++) {
    int v[3]; qhm_get_face_verts(&s->mesh, &s->mesh.faces[fi], v);
    ManifoldVec3 fn = qh_tri_normal(s->verts[v[0]], s->verts[v[1]], s->verts[v[2]]);
    s->mesh.faces[fi].P = qh_plane(fn, s->verts[v[0]]);
  }

  // Assign each point to first visible face
  for (size_t i = 0; i < vc; i++) {
    for (size_t fi = 0; fi < 4; fi++) {
      if (qhs_add_point_to_face(s, &s->mesh.faces[fi], i)) break;
    }
  }
}

static void qhs_build(QHState *s) {
  s->planar = false;
  qhs_setup_initial_tet(s);
  if (s->planar) return;

  // Init face stack
  iv_clear(&s->faceList);
  for (size_t i = 0; i < 4; i++) {
    QHFace *f = &s->mesh.faces[i];
    if (f->hasPoints && f->pointsOnPositiveSide.len > 0) {
      iv_push(&s->faceList, (int)i);
      f->inFaceStack = 1;
    }
  }

  size_t iter = 0;
  size_t faceListFront = 0;
  while (faceListFront < s->faceList.len) {
    iter++;
    if (iter == SIZE_MAX) iter = 0;

    int topFaceIndex = s->faceList.data[faceListFront++];
    if (faceListFront > 1024 && faceListFront > s->faceList.len / 2) {
      size_t remaining = s->faceList.len - faceListFront;
      memmove(s->faceList.data, s->faceList.data + faceListFront,
              remaining * sizeof(int));
      s->faceList.len = remaining;
      faceListFront = 0;
    }

    QHFace *tf = &s->mesh.faces[topFaceIndex];
    tf->inFaceStack = 0;

    if (!tf->hasPoints || tf->pointsOnPositiveSide.len == 0 || qhface_is_disabled(tf))
      continue;

    ManifoldVec3 activePoint = s->verts[tf->mostDistantPoint];
    size_t activePointIndex = tf->mostDistantPoint;

    // Find visible faces via DFS from topFace
    sv_clear(&s->horizonEdges);
    s->pvfLen = 0;
    sv_clear(&s->visibleFaces);

    pvf_push(s, (QHFaceData){topFaceIndex, -1});

    while (s->pvfLen > 0) {
      QHFaceData fd = s->pvfStack[--s->pvfLen];
      QHFace *pvf = &s->mesh.faces[fd.faceIndex];

      if (pvf->visibilityCheckedOnIteration == iter) {
        if (pvf->isVisibleFaceOnCurrentIteration) continue;
      } else {
        pvf->visibilityCheckedOnIteration = iter;
        double d = vec3_dot(pvf->P.N, activePoint) + pvf->P.D;
        if (d > 0) {
          pvf->isVisibleFaceOnCurrentIteration = 1;
          pvf->horizonEdgesOnCurrentIteration = 0;
          sv_push(&s->visibleFaces, (size_t)fd.faceIndex);
          int hes[3]; qhm_get_face_he(&s->mesh, pvf, hes);
          for (int ei = 0; ei < 3; ei++) {
            if (s->mesh.he[hes[ei]].pairedHalfedge != fd.enteredFromHalfedge) {
              int paired = s->mesh.he[hes[ei]].pairedHalfedge;
              pvf_push(s, (QHFaceData){s->mesh.heToFace[paired], hes[ei]});
            }
          }
          continue;
        }
      }

      pvf->isVisibleFaceOnCurrentIteration = 0;
      if (fd.enteredFromHalfedge >= 0) {
        sv_push(&s->horizonEdges, (size_t)fd.enteredFromHalfedge);
        // Mark which halfedge of the visible face is a horizon edge
        int faceOfHE = s->mesh.heToFace[fd.enteredFromHalfedge];
        int hes[3]; qhm_get_face_he(&s->mesh, &s->mesh.faces[faceOfHE], hes);
        int8_t ind = (hes[0] == fd.enteredFromHalfedge) ? 0 : (hes[1] == fd.enteredFromHalfedge ? 1 : 2);
        s->mesh.faces[faceOfHE].horizonEdgesOnCurrentIteration |= (1 << ind);
      }
    }

    size_t horizonEdgeCount = s->horizonEdges.len;

    // Try to reorder horizon edges into a loop
    if (!qhs_reorder_horizon(s)) {
      // Failed: remove active point and continue
      SzVec *pts = &tf->pointsOnPositiveSide;
      bool found = false;
      for (size_t i = 0; i < pts->len; i++) {
        if (pts->data[i] == activePointIndex) found = true;
        if (found && i + 1 < pts->len) pts->data[i] = pts->data[i + 1];
      }
      if (found) pts->len--;
      if (pts->len == 0) tf->hasPoints = false;
      continue;
    }

    // Disable visible faces and collect their point lists
    // Collect disabled face point vectors
    SzVec *savedPointVecs = (SzVec*)malloc(s->visibleFaces.len * sizeof(SzVec));
    bool *savedHasPoints = (bool*)calloc(s->visibleFaces.len, sizeof(bool));
    size_t savedCount = 0;

    sv_clear(&s->newFaceIndices);
    sv_clear(&s->newHEIndices);
    size_t disableCounter = 0;

    for (size_t vi = 0; vi < s->visibleFaces.len; vi++) {
      size_t faceIndex = s->visibleFaces.data[vi];
      QHFace *df = &s->mesh.faces[faceIndex];
      int hes[3]; qhm_get_face_he(&s->mesh, df, hes);
      for (int j = 0; j < 3; j++) {
        if ((df->horizonEdgesOnCurrentIteration & (1 << j)) == 0) {
          if (disableCounter < horizonEdgeCount * 2) {
            sv_push(&s->newHEIndices, (size_t)hes[j]);
            disableCounter++;
          } else {
            qhm_disable_he(&s->mesh, (size_t)hes[j]);
          }
        }
      }
      // Save points
      if (df->hasPoints && df->pointsOnPositiveSide.len > 0) {
        savedPointVecs[savedCount] = df->pointsOnPositiveSide;
        savedHasPoints[savedCount] = true;
        savedCount++;
        // Reset without freeing (we took ownership)
        sv_init(&df->pointsOnPositiveSide);
        df->hasPoints = false;
      }
      qhm_disable_face(&s->mesh, faceIndex);
    }

    // Add more halfedges if needed
    while (disableCounter < horizonEdgeCount * 2) {
      sv_push(&s->newHEIndices, qhm_add_he(&s->mesh));
      disableCounter++;
    }

    // Create new faces connecting active point to horizon edges
    for (size_t i = 0; i < horizonEdgeCount; i++) {
      size_t AB = s->horizonEdges.data[i];
      int A = s->mesh.he[s->mesh.he[AB].pairedHalfedge].endVert;
      int B = s->mesh.he[AB].endVert;
      int C = (int)activePointIndex;

      size_t newFI = qhm_add_face(&s->mesh);
      sv_push(&s->newFaceIndices, newFI);

      size_t CA = s->newHEIndices.data[2 * i + 0];
      size_t BC = s->newHEIndices.data[2 * i + 1];

      s->mesh.heNext[AB] = (int)BC;
      s->mesh.heNext[BC] = (int)CA;
      s->mesh.heNext[CA] = (int)AB;

      s->mesh.heToFace[BC] = (int)newFI;
      s->mesh.heToFace[CA] = (int)newFI;
      s->mesh.heToFace[AB] = (int)newFI;

      s->mesh.he[CA].endVert = A;
      s->mesh.he[BC].endVert = C;

      QHFace *nf = &s->mesh.faces[newFI];
      ManifoldVec3 planeN = qh_tri_normal(s->verts[A], s->verts[B], activePoint);
      nf->P = qh_plane(planeN, activePoint);
      nf->he = (int)AB;

      // Pair CA with previous face's BC, and BC with next face's CA
      s->mesh.he[CA].pairedHalfedge =
          (int)s->newHEIndices.data[i > 0 ? i * 2 - 1 : 2 * horizonEdgeCount - 1];
      s->mesh.he[BC].pairedHalfedge =
          (int)s->newHEIndices.data[((i + 1) * 2) % (horizonEdgeCount * 2)];
    }

    // Reassign points from disabled faces to new faces
    for (size_t si = 0; si < savedCount; si++) {
      if (!savedHasPoints[si]) continue;
      SzVec *pts = &savedPointVecs[si];
      for (size_t pi = 0; pi < pts->len; pi++) {
        size_t ptIdx = pts->data[pi];
        if (ptIdx == activePointIndex) continue;
        for (size_t j = 0; j < horizonEdgeCount; j++) {
          if (qhs_add_point_to_face(s, &s->mesh.faces[s->newFaceIndices.data[j]], ptIdx))
            break;
        }
      }
      sv_free(pts);
    }
    free(savedPointVecs);
    free(savedHasPoints);

    // Push new faces with points to face stack
    for (size_t i = 0; i < s->newFaceIndices.len; i++) {
      QHFace *nf = &s->mesh.faces[s->newFaceIndices.data[i]];
      if (nf->hasPoints && nf->pointsOnPositiveSide.len > 0 && !nf->inFaceStack) {
        iv_push(&s->faceList, (int)s->newFaceIndices.data[i]);
        nf->inFaceStack = 1;
      }
    }
  }
}

// ─── Build the manifold mesh from QuickHull result ──────────────────────────
void manifold_convex_hull(ManifoldImpl *impl, const ManifoldVec3 *points,
                           size_t numPoints) {
  manifold_impl_init(impl);

  if (numPoints == 0) {
    manifold_impl_make_empty(impl, MANIFOLD_ERROR_INVALID_CONSTRUCTION);
    return;
  }

  // Pad to at least 4 points (matching C++ behavior for degenerate cases)
  ManifoldVec3 *paddedPoints = NULL;
  size_t paddedCount = numPoints;
  if (numPoints < 4) {
    paddedCount = 4;
    paddedPoints = (ManifoldVec3*)malloc(4 * sizeof(ManifoldVec3));
    for (size_t i = 0; i < numPoints; i++) paddedPoints[i] = points[i];
    for (size_t i = numPoints; i < 4; i++) paddedPoints[i] = points[numPoints - 1];
    points = paddedPoints;
    numPoints = 4;
  }

  QHState state;
  qhs_init(&state, points, numPoints);

  state.extremes[0] = 0;
  qhs_get_extremes(&state);
  state.scale = qhs_get_scale(&state);

  double eps = 0.0000001;
  state.epsilon = eps * state.scale;
  state.epsilonSq = state.epsilon * state.epsilon;

  qhs_build(&state);

  if (state.planar) {
    // Create a degenerate tetrahedron for planar/collinear/few-point cases
    // Pad points to at least 4 by duplicating
    size_t padCount = numPoints < 4 ? 4 : numPoints;
    ManifoldVec3 *padPoints = (ManifoldVec3*)malloc(padCount * sizeof(ManifoldVec3));
    for (size_t i = 0; i < numPoints; i++) padPoints[i] = points[i];
    for (size_t i = numPoints; i < 4; i++) padPoints[i] = points[numPoints - 1];

    // Try adding a point offset by the triangle normal for planar case
    if (numPoints >= 3) {
      ManifoldVec3 N = qh_tri_normal(padPoints[0], padPoints[1], padPoints[2]);
      double nlen = sqrt(N.x*N.x + N.y*N.y + N.z*N.z);
      if (nlen > 1e-30) {
        // Add extra point offset by normal
        padPoints = (ManifoldVec3*)realloc(padPoints, (padCount + 1) * sizeof(ManifoldVec3));
        padPoints[padCount] = (ManifoldVec3){
          padPoints[0].x + N.x, padPoints[0].y + N.y, padPoints[0].z + N.z};
        
        // Build hull with extra point
        qhs_free(&state);
        memset(&state, 0, sizeof(state));
        state.verts = padPoints;
        state.vertCount = padCount + 1;
        state.extremes[0] = 0;
        qhs_get_extremes(&state);
        state.scale = qhs_get_scale(&state);
        state.epsilon = eps * state.scale;
        state.epsilonSq = state.epsilon * state.epsilon;
        qhs_build(&state);

        if (!state.planar) {
          // Keep ALL faces including those with extra point, but remap
          // extra point's position to vertex[0]. This creates degenerate
          // faces that are topologically valid (all edges paired).
          // The degenerate faces get removed later by Simplify().
          size_t extraIdx = padCount;
          ManifoldVecIVec3 triVerts = {0};
          
          for (size_t fi = 0; fi < state.mesh.faceLen; fi++) {
            if (qhface_is_disabled(&state.mesh.faces[fi])) continue;
            int hes[3]; qhm_get_face_he(&state.mesh, &state.mesh.faces[fi], hes);
            if (state.mesh.he[hes[0]].pairedHalfedge < 0) continue;
            int sv0 = state.mesh.he[state.mesh.heNext[state.mesh.heNext[hes[0]]]].endVert;
            int sv1 = state.mesh.he[hes[0]].endVert;
            int sv2 = state.mesh.he[state.mesh.heNext[hes[0]]].endVert;
            // Clamp padded duplicate vertices to their originals
            if ((size_t)sv0 >= numPoints && (size_t)sv0 != extraIdx) sv0 = (int)(numPoints - 1);
            if ((size_t)sv1 >= numPoints && (size_t)sv1 != extraIdx) sv1 = (int)(numPoints - 1);
            if ((size_t)sv2 >= numPoints && (size_t)sv2 != extraIdx) sv2 = (int)(numPoints - 1);
            // Skip degenerate (two same vertices)
            if (sv0 == sv1 || sv1 == sv2 || sv2 == sv0) continue;
            vec_ivec3_push(&triVerts, manifold_ivec3(sv0, sv1, sv2));
          }
          
          if (triVerts.len > 0) {
            // Build vertex map for used vertices (including extra point)
            int *vertUsed = (int*)calloc(padCount + 1, sizeof(int));
            for (size_t t = 0; t < triVerts.len; t++) {
              vertUsed[triVerts.data[t].x] = 1;
              vertUsed[triVerts.data[t].y] = 1;
              vertUsed[triVerts.data[t].z] = 1;
            }
            int vertCount = 0;
            int *vertMap = (int*)malloc((padCount + 1) * sizeof(int));
            for (size_t i = 0; i < padCount + 1; i++) {
              vertMap[i] = vertUsed[i] ? vertCount++ : -1;
            }
            
            // Remap triangle indices
            for (size_t t = 0; t < triVerts.len; t++) {
              triVerts.data[t].x = vertMap[triVerts.data[t].x];
              triVerts.data[t].y = vertMap[triVerts.data[t].y];
              triVerts.data[t].z = vertMap[triVerts.data[t].z];
            }
            
            impl->vertPos = vec_vec3_create_n((size_t)vertCount);
            for (size_t i = 0; i < numPoints; i++) {
              if (vertUsed[i]) impl->vertPos.data[vertMap[i]] = points[i];
            }
            // Set extra point position to same as vertex[0]
            if (vertUsed[extraIdx]) {
              impl->vertPos.data[vertMap[extraIdx]] = points[0];
            }
            
            ManifoldVecIVec3 emptyTriVert = {0};
            manifold_impl_create_halfedges(impl, &triVerts, &emptyTriVert);
            manifold_impl_initialize_original(impl);
            manifold_impl_calculate_bbox(impl);
            manifold_impl_set_epsilon(impl, -1.0, false);
            manifold_impl_sort_geometry(impl);
            manifold_impl_set_normals_and_coplanar(impl);
            
            vec_ivec3_free(&triVerts);
            vec_ivec3_free(&emptyTriVert);
            free(vertMap);
            free(vertUsed);
          } else {
            vec_ivec3_free(&triVerts);
          }
          
          free(padPoints);
          qhs_free(&state);
          free(paddedPoints);
          return;
        }
      }
    }
    
    // Truly degenerate (collinear or single point) - create a degenerate 
    // tetrahedron with valid topology but zero volume
    if (numPoints >= 2) {
      // Find most distant pair
      double maxDist = 0;
      size_t p0 = 0, p1 = 1;
      for (size_t i = 0; i < numPoints; i++) {
        for (size_t j = i + 1; j < numPoints; j++) {
          double dx = points[i].x - points[j].x;
          double dy = points[i].y - points[j].y;
          double dz = points[i].z - points[j].z;
          double d = dx*dx + dy*dy + dz*dz;
          if (d > maxDist) { maxDist = d; p0 = i; p1 = j; }
        }
      }
      
      // Create a tetrahedron with 4 vertices (v2,v3 at same position as v0)
      impl->vertPos = vec_vec3_create_n(4);
      impl->vertPos.data[0] = points[p0];
      impl->vertPos.data[1] = points[p1];
      impl->vertPos.data[2] = points[p0];
      impl->vertPos.data[3] = points[p0];
      
      ManifoldVecIVec3 triVerts = {0};
      vec_ivec3_push(&triVerts, manifold_ivec3(0, 1, 2));
      vec_ivec3_push(&triVerts, manifold_ivec3(0, 2, 3));
      vec_ivec3_push(&triVerts, manifold_ivec3(0, 3, 1));
      vec_ivec3_push(&triVerts, manifold_ivec3(1, 3, 2));
      ManifoldVecIVec3 emptyTriVert = {0};
      manifold_impl_create_halfedges(impl, &triVerts, &emptyTriVert);
      manifold_impl_initialize_original(impl);
      manifold_impl_calculate_bbox(impl);
      manifold_impl_set_epsilon(impl, -1.0, false);
      // Don't call sort_geometry — it may collapse degenerate faces
      
      vec_ivec3_free(&triVerts);
      vec_ivec3_free(&emptyTriVert);
    }
    
    free(padPoints);
    qhs_free(&state);
    free(paddedPoints);
    return;
  }

  // Collect enabled faces and build triangle mesh
  // First pass: count faces and map vertices
  int *vertUsed = (int*)calloc(numPoints, sizeof(int));
  size_t triCount = 0;
  for (size_t fi = 0; fi < state.mesh.faceLen; fi++) {
    if (qhface_is_disabled(&state.mesh.faces[fi])) continue;
    int hes[3]; qhm_get_face_he(&state.mesh, &state.mesh.faces[fi], hes);
    // Validate halfedges
    if (state.mesh.he[hes[0]].pairedHalfedge < 0) continue;
    triCount++;
    int v[3]; qhm_get_face_verts(&state.mesh, &state.mesh.faces[fi], v);
    for (int j = 0; j < 3; j++) {
      if (v[j] >= 0 && (size_t)v[j] < numPoints) vertUsed[v[j]] = 1;
    }
  }

  // Count unique vertices used
  int usedVertCount = 0;
  for (size_t i = 0; i < numPoints; i++) {
    if (vertUsed[i]) usedVertCount++;
  }
  
  if (triCount < 4 || usedVertCount < 4) {
    // For degenerate cases (collinear/coplanar with < 4 faces),
    // create a degenerate tetrahedron to match C++ behavior
    if (numPoints >= 2) {
      // Find most distant pair
      double maxDist = 0;
      size_t dp0 = 0, dp1 = 1;
      for (size_t i = 0; i < numPoints; i++) {
        for (size_t j = i + 1; j < numPoints; j++) {
          double dx = points[i].x - points[j].x;
          double dy = points[i].y - points[j].y;
          double dz = points[i].z - points[j].z;
          double d = dx*dx + dy*dy + dz*dz;
          if (d > maxDist) { maxDist = d; dp0 = i; dp1 = j; }
        }
      }
      if (maxDist > 0) {
        impl->vertPos = vec_vec3_create_n(4);
        impl->vertPos.data[0] = points[dp0];
        impl->vertPos.data[1] = points[dp1];
        impl->vertPos.data[2] = points[dp0];
        impl->vertPos.data[3] = points[dp0];
        ManifoldVecIVec3 tv = {0};
        vec_ivec3_push(&tv, manifold_ivec3(0, 1, 2));
        vec_ivec3_push(&tv, manifold_ivec3(0, 2, 3));
        vec_ivec3_push(&tv, manifold_ivec3(0, 3, 1));
        vec_ivec3_push(&tv, manifold_ivec3(1, 3, 2));
        ManifoldVecIVec3 empty2 = {0};
        manifold_impl_create_halfedges(impl, &tv, &empty2);
        manifold_impl_initialize_original(impl);
        manifold_impl_calculate_bbox(impl);
        manifold_impl_set_epsilon(impl, -1.0, false);
        // Don't call sort_geometry — it may collapse degenerate faces
        vec_ivec3_free(&tv);
        vec_ivec3_free(&empty2);
        free(vertUsed);
        qhs_free(&state);
        free(paddedPoints);
        return;
      }
    }
    free(vertUsed);
    qhs_free(&state);
    free(paddedPoints);
    manifold_impl_make_empty(impl, MANIFOLD_ERROR_INVALID_CONSTRUCTION);
    return;
  }

  // Build vertex map: old index -> new index
  int vertCount = 0;
  int *vertMap = (int*)malloc(numPoints * sizeof(int));
  for (size_t i = 0; i < numPoints; i++) {
    if (vertUsed[i]) vertMap[i] = vertCount++;
    else vertMap[i] = -1;
  }

  // Build triangles
  ManifoldVecIVec3 triVerts = {0};
  for (size_t fi = 0; fi < state.mesh.faceLen; fi++) {
    if (qhface_is_disabled(&state.mesh.faces[fi])) continue;
    if (state.mesh.he[state.mesh.faces[fi].he].pairedHalfedge < 0) continue;
    int v[3]; qhm_get_face_verts(&state.mesh, &state.mesh.faces[fi], v);
    if (v[0] < 0 || v[1] < 0 || v[2] < 0) continue;
    if ((size_t)v[0] >= numPoints || (size_t)v[1] >= numPoints || (size_t)v[2] >= numPoints) continue;
    // Compute face startVert from halfedge structure
    int hes[3]; qhm_get_face_he(&state.mesh, &state.mesh.faces[fi], hes);
    int sv0 = state.mesh.he[state.mesh.heNext[state.mesh.heNext[hes[0]]]].endVert;
    int sv1 = state.mesh.he[hes[0]].endVert;
    int sv2 = state.mesh.he[state.mesh.heNext[hes[0]]].endVert;
    vec_ivec3_push(&triVerts, manifold_ivec3(vertMap[sv0], vertMap[sv1], vertMap[sv2]));
  }

  // Copy vertices
  impl->vertPos = vec_vec3_create_n((size_t)vertCount);
  for (size_t i = 0; i < numPoints; i++) {
    if (vertMap[i] >= 0) impl->vertPos.data[vertMap[i]] = points[i];
  }

  ManifoldVecIVec3 emptyTriVert = {0};
  manifold_impl_create_halfedges(impl, &triVerts, &emptyTriVert);
  manifold_impl_initialize_original(impl);
  manifold_impl_calculate_bbox(impl);
  manifold_impl_set_epsilon(impl, -1.0, false);
  manifold_impl_sort_geometry(impl);
  manifold_impl_set_normals_and_coplanar(impl);

  vec_ivec3_free(&triVerts);
  vec_ivec3_free(&emptyTriVert);
  free(vertMap);
  free(vertUsed);
  qhs_free(&state);
  free(paddedPoints);
}
