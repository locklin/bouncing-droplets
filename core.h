/*
 * core.h — Shared types, constants, and function declarations
 * for the bouncing droplet simulator suite.
 */
#ifndef CORE_H
#define CORE_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "raylib.h"

/* ----- Types ----- */
typedef struct { double x, y; } V2;

/* ----- Layout ----- */
#define GRID   200
#define PANEL  400
#define GAP     20
#define WIN_W  (PANEL*2 + GAP*3)
#define WIN_H  (PANEL + GAP*2 + 60)

/* ----- Geometry ----- */
enum { GEO_STADIUM, GEO_DSHAPE, GEO_CIRCLE, GEO_RECT, GEO_NONE, GEO_COUNT };

/* All geometry functions take an explicit geo parameter so they
 * don't depend on a global.  Callers pass their own `geometry` var. */
double      geo_sdf(int geo, double x, double y);
int         geo_inside(int geo, double x, double y);
void        geo_normal(int geo, double x, double y, double *nx, double *ny);
void        geo_draw_outline(int geo, int ox, int oy, int sz);
double      geo_arclength(int geo, double x, double y);
const char *geo_name(int geo);

/* ----- Colormaps ----- */
Color cmap_bwr(double t);
Color cmap_hot(double t);

/* ----- Histogram ----- */
/* Gaussian-splat a position into a GRID x GRID histogram.
 * coords in [-scale, scale] mapped to [0, GRID). */
void hist_splat(int hist[][GRID], int *hpeak, double x, double y, double scale);

/* Render histogram to pixel buffer. Geo clipping if geo >= 0, scale for coord mapping. */
void render_hist(Color *pix, int hist[][GRID], int hpeak, int geo, double scale);

/* ----- Moments (Born rule / Khrennikov) ----- */
void render_born(Color *pix, double h2[][GRID], int n, int geo, double scale);
void render_kurtosis(Color *pix, double h2[][GRID], double h4[][GRID], int n, int geo, double scale);

/* ----- Wave field ----- */
/* Render a GRID x GRID wave field (blue-white-red diverging) with geometry clipping */
void render_wave(Color *pix, double wave[][GRID], int geo, double scale);

/* ----- Poincare section ----- */
/* Render auto-ranged density plot from (x, y) point buffer.
 * Ring buffer: pts array of max_pts, n_pts total accumulated, head = write position. */
void render_poincare(Color *pix, V2 *pts, int n_pts, int head, int max_pts);

/* ----- Trail ----- */
/* Draw fading polyline trail. Ring buffer of V2 positions.
 * Coords mapped from [-scale, scale] to panel. */
void draw_trail(V2 *buf, int head, int count, int vis,
                int ox, int oy, int sz, double scale, Color col);

/* ----- Droplet ----- */
void draw_droplet(double x, double y, int ox, int oy, int sz, double scale);

/* ----- Turbo display ----- */
void draw_turbo(int ox, int oy, int sz, double rate, const char *unit);

/* ----- Angular momentum strip ----- */
/* Draw scrolling L time series. Ring buffer of doubles.
 * ox, oy, w, h define the strip rectangle. */
void draw_L_strip(double *L_hist, int L_head, int L_count, int L_max,
                  double L_scale, int ox, int oy, int w, int h);

/* ----- Coordinate conversion ----- */
Vector2 phys2screen(double px, double py, int ox, int oy, int sz, double scale);

/* ----- Layout with L strip ----- */
#define LSTRIP  80
#define WIN_W_L (PANEL*2 + GAP*3)
#define WIN_H_L (PANEL + LSTRIP + GAP*3 + 40)

#endif /* CORE_H */
