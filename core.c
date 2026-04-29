/*
 * core.c — Shared geometry, rendering, and utility functions
 * for the bouncing droplet simulator suite.
 */
#include "core.h"

/* ================================================================
   GEOMETRY
   ================================================================ */

#define DSHAPE_R   0.9
#define DSHAPE_CUT (-0.3)

double geo_sdf(int geo, double x, double y)
{
    switch (geo) {
    case GEO_CIRCLE:  return 1.0 - sqrt(x*x + y*y);
    case GEO_STADIUM: {
        double a = 0.5, r = 0.5, ax = fabs(x);
        return ax <= a ? r - fabs(y) : r - sqrt((ax-a)*(ax-a) + y*y);
    }
    case GEO_RECT:    return fmin(0.8 - fabs(x), 0.5 - fabs(y));
    case GEO_DSHAPE:  return fmin(DSHAPE_R - sqrt(x*x+y*y), x - DSHAPE_CUT);
    case GEO_ELLSTADIUM: {
        double ax = fabs(x);
        double ymax;
        if (ax <= ESTD_D) ymax = ESTD_G;
        else if (ax <= 1.0) {
            double u = (ax - ESTD_D) / (1.0 - ESTD_D);
            ymax = ESTD_G * sqrt(fmax(0, 1.0 - u*u));
        } else return -fabs(y);
        return ymax - fabs(y);
    }
    default:          return 1.0;  /* GEO_NONE: always inside */
    }
}

int geo_inside(int geo, double x, double y)
{
    return geo == GEO_NONE || geo_sdf(geo, x, y) > 0;
}

void geo_normal(int geo, double x, double y, double *nx, double *ny)
{
    double eps = 0.002;
    *nx = -(geo_sdf(geo, x+eps, y) - geo_sdf(geo, x-eps, y)) / (2*eps);
    *ny = -(geo_sdf(geo, x, y+eps) - geo_sdf(geo, x, y-eps)) / (2*eps);
}

double geo_arclength(int geo, double x, double y)
{
    switch (geo) {
    case GEO_CIRCLE: {
        double a = atan2(y, x);
        return (a + M_PI) / (2*M_PI);
    }
    case GEO_STADIUM: {
        double a = 0.5, r = 0.5;
        double perim = 2*a*2 + 2*M_PI*r;
        double s;
        if (fabs(y+r) < 0.02 && fabs(x) <= a) s = (x+a);
        else if (x > a-0.02) { double ang = atan2(y, x-a); s = 2*a + (ang+M_PI/2)*r; }
        else if (fabs(y-r) < 0.02 && fabs(x) <= a) s = 2*a + M_PI*r + (a-x);
        else { double ang = atan2(y, x+a); s = 2*a + M_PI*r + 2*a + (ang-M_PI/2+M_PI)*r; }
        return fmod(s/perim + 1.0, 1.0);
    }
    case GEO_RECT: {
        double hw = 0.8, hh = 0.5, perim = 2*(2*hw+2*hh);
        double s;
        if (fabs(y+hh) < 0.02) s = (x+hw);
        else if (fabs(x-hw) < 0.02) s = 2*hw + (y+hh);
        else if (fabs(y-hh) < 0.02) s = 2*hw + 2*hh + (hw-x);
        else s = 2*hw + 2*hh + 2*hw + (hh-y);
        return fmod(s/perim + 1.0, 1.0);
    }
    case GEO_DSHAPE: {
        double R = DSHAPE_R, yint = sqrt(R*R - DSHAPE_CUT*DSHAPE_CUT);
        double ang_top = atan2(yint, DSHAPE_CUT), ang_bot = atan2(-yint, DSHAPE_CUT);
        double arc_len = R*(ang_top-ang_bot), chord_len = 2*yint, perim = arc_len+chord_len;
        if (fabs(x-DSHAPE_CUT) < 0.02 && fabs(y) < yint)
            return (y+yint) / perim;
        else {
            double ang = atan2(y, x);
            return fmod((chord_len + R*(ang_top-ang)) / perim + 1.0, 1.0);
        }
    }
    case GEO_ELLSTADIUM: {
        /* Approximate arclength using angle on the boundary */
        double d = ESTD_D, g = ESTD_G;
        double ax = fabs(x);
        if (ax <= d) {
            /* On flat top or bottom */
            if (y > 0) return 0.5 * (x + d) / (2*d + 2*(1-d)); /* top flat, rough */
            else return 0.5 + 0.5 * (d - x) / (2*d + 2*(1-d)); /* bottom flat, rough */
        } else {
            /* On ellipse — use angle */
            double u = (ax - d) / (1 - d);
            double v = fabs(y) / g;
            double ang = atan2(v, u);
            return fmod(ang / (2*M_PI) + 0.5, 1.0);
        }
    }
    default: return 0;
    }
}

void geo_draw_outline(int geo, int ox, int oy, int sz)
{
    int N = 200;
    Color col = LIGHTGRAY;
    #define S2P(px, py) (Vector2){ \
        ox + ((float)(px)*0.5f+0.5f)*sz, oy + ((float)(py)*0.5f+0.5f)*sz }

    switch (geo) {
    case GEO_NONE: break;
    case GEO_CIRCLE:
        DrawCircleLinesV((Vector2){ox+sz/2.0f,oy+sz/2.0f}, sz/2.0f, col); break;
    case GEO_STADIUM: {
        float a=0.5f, r=0.5f;
        DrawLineV(S2P(-a,r), S2P(a,r), col); DrawLineV(S2P(-a,-r), S2P(a,-r), col);
        for (int i=0; i<N; i++) {
            float t0=PI/2+PI*i/(float)N, t1=PI/2+PI*(i+1)/(float)N;
            DrawLineV(S2P(-a+r*cosf(t0),r*sinf(t0)),S2P(-a+r*cosf(t1),r*sinf(t1)),col); }
        for (int i=0; i<N; i++) {
            float t0=-PI/2+PI*i/(float)N, t1=-PI/2+PI*(i+1)/(float)N;
            DrawLineV(S2P(a+r*cosf(t0),r*sinf(t0)),S2P(a+r*cosf(t1),r*sinf(t1)),col); }
        break; }
    case GEO_RECT: {
        float hw=0.8f, hh=0.5f;
        DrawLineV(S2P(-hw,-hh),S2P(hw,-hh),col); DrawLineV(S2P(hw,-hh),S2P(hw,hh),col);
        DrawLineV(S2P(hw,hh),S2P(-hw,hh),col); DrawLineV(S2P(-hw,hh),S2P(-hw,-hh),col);
        break; }
    case GEO_DSHAPE: {
        float R=DSHAPE_R, yint=sqrtf(R*R-DSHAPE_CUT*DSHAPE_CUT);
        float a0=atan2f(yint,DSHAPE_CUT), a1=atan2f(-yint,DSHAPE_CUT), span=a0-a1;
        DrawLineV(S2P(DSHAPE_CUT,-yint), S2P(DSHAPE_CUT,yint), col);
        for (int i=0; i<N; i++) {
            float t0=a0-span*i/(float)N, t1=a0-span*(i+1)/(float)N;
            DrawLineV(S2P(R*cosf(t0),R*sinf(t0)),S2P(R*cosf(t1),R*sinf(t1)),col); }
        break; }
    case GEO_ELLSTADIUM: {
        float d=ESTD_D, g=ESTD_G;
        /* Flat top/bottom */
        DrawLineV(S2P(-d,g), S2P(d,g), col);
        DrawLineV(S2P(-d,-g), S2P(d,-g), col);
        /* Right semiellipse: x = d + (1-d)*cos(t), y = g*sin(t), t from pi/2 to -pi/2 */
        for (int i=0; i<N; i++) {
            float t0=PI/2-PI*i/(float)N, t1=PI/2-PI*(i+1)/(float)N;
            DrawLineV(S2P(d+(1-d)*cosf(t0), g*sinf(t0)),
                      S2P(d+(1-d)*cosf(t1), g*sinf(t1)), col); }
        /* Left semiellipse: mirror of right */
        for (int i=0; i<N; i++) {
            float t0=PI/2-PI*i/(float)N, t1=PI/2-PI*(i+1)/(float)N;
            DrawLineV(S2P(-(d+(1-d)*cosf(t0)), g*sinf(t0)),
                      S2P(-(d+(1-d)*cosf(t1)), g*sinf(t1)), col); }
        break; }
    }
    #undef S2P
}

const char *geo_name(int geo)
{
    const char *names[] = {"Stadium","D-shape","Circle","Ell.Stadium","Rectangle","None (harmonic)"};
    return (geo >= 0 && geo < GEO_COUNT) ? names[geo] : "?";
}

/* ================================================================
   COLORMAPS
   ================================================================ */

Color cmap_bwr(double t)
{
    if (t < 0) t = 0;
    if (t > 1) t = 1;
    if (t < 0.5) { unsigned char v=(unsigned char)(t*2*255); return (Color){v,v,255,255}; }
    else { unsigned char v=(unsigned char)((1-t)*2*255); return (Color){255,v,v,255}; }
}

Color cmap_hot(double t)
{
    if (t < 0) t = 0;
    if (t > 1) t = 1;
    unsigned char r, g, b;
    if (t < 0.33) { r=(unsigned char)(t/0.33*255); g=0; b=0; }
    else if (t < 0.66) { r=255; g=(unsigned char)((t-0.33)/0.33*255); b=0; }
    else { r=255; g=255; b=(unsigned char)((t-0.66)/0.34*255); }
    return (Color){r, g, b, 255};
}

/* ================================================================
   COORDINATE CONVERSION
   ================================================================ */

Vector2 phys2screen(double px, double py, int ox, int oy, int sz, double scale)
{
    return (Vector2){
        ox + (float)((px/scale*0.5+0.5) * sz),
        oy + (float)((py/scale*0.5+0.5) * sz)
    };
}

/* ================================================================
   HISTOGRAM
   ================================================================ */

void hist_splat(int hist[][GRID], int *hpeak, double x, double y, double scale)
{
    double hxf = (x/scale*0.5+0.5)*GRID;
    double hyf = (y/scale*0.5+0.5)*GRID;
    int cx = (int)hxf, cy = (int)hyf;
    for (int dy = -2; dy <= 2; dy++)
        for (int dx = -2; dx <= 2; dx++) {
            int px = cx+dx, py = cy+dy;
            if (px<0 || px>=GRID || py<0 || py>=GRID) continue;
            double gg = (px-hxf)*(px-hxf) + (py-hyf)*(py-hyf);
            hist[py][px] += (int)(100 * exp(-gg));
            if (hist[py][px] > *hpeak) *hpeak = hist[py][px];
        }
}

/* ================================================================
   RENDERING FUNCTIONS
   ================================================================ */

void render_hist(Color *pix, int hist[][GRID], int hpeak, int geo, double scale)
{
    double lmin = 1e30, lmax = 0;
    for (int i = 0; i < GRID*GRID; i++)
        if (hist[i/GRID][i%GRID] > 0) {
            double lv = log(1.0 + hist[i/GRID][i%GRID]);
            if (lv < lmin) lmin = lv;
            if (lv > lmax) lmax = lv;
        }
    double lr = lmax - lmin; if (lr < 1e-10) lr = 1;
    for (int iy=0; iy<GRID; iy++) {
        double y = (2.0*iy/(GRID-1)-1.0)*scale;
        for (int ix=0; ix<GRID; ix++) {
            double x = (2.0*ix/(GRID-1)-1.0)*scale;
            int idx = iy*GRID+ix;
            if (!geo_inside(geo, x, y)) pix[idx] = (Color){30,30,30,255};
            else if (hist[iy][ix] == 0) pix[idx] = BLACK;
            else pix[idx] = cmap_hot((log(1.0+hist[iy][ix]) - lmin) / lr);
        }
    }
}

void render_born(Color *pix, double h2[][GRID], int n, int geo, double scale)
{
    if (n < 1) return;
    double peak = 0;
    for (int iy=0; iy<GRID; iy++)
        for (int ix=0; ix<GRID; ix++) {
            double v = h2[iy][ix]/n; if (v > peak) peak = v; }
    if (peak < 1e-20) peak = 1;
    for (int iy=0; iy<GRID; iy++) {
        double y = (2.0*iy/(GRID-1)-1.0)*scale;
        for (int ix=0; ix<GRID; ix++) {
            double x = (2.0*ix/(GRID-1)-1.0)*scale;
            int idx = iy*GRID+ix;
            pix[idx] = !geo_inside(geo, x, y) ? (Color){30,30,30,255}
                      : cmap_hot(sqrt(h2[iy][ix]/n/peak));
        }
    }
}

void render_kurtosis(Color *pix, double h2[][GRID], double h4[][GRID], int n,
                     int geo, double scale)
{
    if (n < 10) return;
    double vmin=0, vmax=0;
    for (int iy=0; iy<GRID; iy++)
        for (int ix=0; ix<GRID; ix++) {
            double m2=h2[iy][ix]/n, m4=h4[iy][ix]/n;
            double ek = m4 - 3.0*m2*m2;
            if (ek<vmin) vmin=ek; if (ek>vmax) vmax=ek;
        }
    double amax = fmax(fabs(vmin),fabs(vmax)); if (amax<1e-20) amax=1;
    for (int iy=0; iy<GRID; iy++) {
        double y = (2.0*iy/(GRID-1)-1.0)*scale;
        for (int ix=0; ix<GRID; ix++) {
            double x = (2.0*ix/(GRID-1)-1.0)*scale;
            int idx = iy*GRID+ix;
            if (!geo_inside(geo, x, y)) { pix[idx]=(Color){30,30,30,255}; continue; }
            double m2=h2[iy][ix]/n, m4=h4[iy][ix]/n;
            pix[idx] = cmap_bwr(0.5 + 0.5*(m4-3.0*m2*m2)/amax);
        }
    }
}

void render_wave(Color *pix, double wave[][GRID], int geo, double scale)
{
    double vmin=0, vmax=0;
    for (int iy=0; iy<GRID; iy++)
        for (int ix=0; ix<GRID; ix++) {
            if (wave[iy][ix] < vmin) vmin = wave[iy][ix];
            if (wave[iy][ix] > vmax) vmax = wave[iy][ix];
        }
    double range = vmax-vmin; if (range < 1e-10) range = 1;
    for (int iy=0; iy<GRID; iy++) {
        double y = (2.0*iy/(GRID-1)-1.0)*scale;
        for (int ix=0; ix<GRID; ix++) {
            double x = (2.0*ix/(GRID-1)-1.0)*scale;
            int idx = iy*GRID+ix;
            pix[idx] = !geo_inside(geo, x, y) ? (Color){30,30,30,255}
                      : cmap_bwr((wave[iy][ix]-vmin)/range);
        }
    }
}

void render_poincare(Color *pix, V2 *pts, int n_pts, int head, int max_pts)
{
    for (int i = 0; i < GRID*GRID; i++) pix[i] = (Color){15,15,20,255};
    int np = n_pts < max_pts ? n_pts : max_pts;
    if (np == 0) return;

    /* Find data range */
    double x_lo=1e30, x_hi=-1e30, v_lo=1e30, v_hi=-1e30;
    for (int i = 0; i < np; i++) {
        int bi = (n_pts <= max_pts) ? i : (head+i) % max_pts;
        if (pts[bi].x < x_lo) x_lo = pts[bi].x;
        if (pts[bi].x > x_hi) x_hi = pts[bi].x;
        if (pts[bi].y < v_lo) v_lo = pts[bi].y;
        if (pts[bi].y > v_hi) v_hi = pts[bi].y;
    }
    double xm = (x_hi-x_lo)*0.1 + 0.01; x_lo -= xm; x_hi += xm;
    double va = fmax(fabs(v_lo), fabs(v_hi)) * 1.1 + 0.001;
    v_lo = -va; v_hi = va;

    /* Zero axes */
    if (x_lo < 0 && x_hi > 0) {
        int zx = (int)((0-x_lo)/(x_hi-x_lo)*(GRID-1));
        if (zx>=0 && zx<GRID)
            for (int iy=0; iy<GRID; iy++) pix[iy*GRID+zx] = (Color){40,40,50,255};
    }
    for (int ix=0; ix<GRID; ix++) pix[(GRID/2)*GRID+ix] = (Color){40,40,50,255};

    /* Bin into density grid */
    static int pd[GRID][GRID];
    memset(pd, 0, sizeof(pd));
    int pm = 1;
    for (int i = 0; i < np; i++) {
        int bi = (n_pts <= max_pts) ? i : (head+i) % max_pts;
        int gx = (int)((pts[bi].x - x_lo) / (x_hi - x_lo) * (GRID-1));
        int gy = (int)((1.0 - (pts[bi].y - v_lo) / (v_hi - v_lo)) * (GRID-1));
        if (gx>=0 && gx<GRID && gy>=0 && gy<GRID) {
            pd[gy][gx]++;
            if (pd[gy][gx] > pm) pm = pd[gy][gx];
        }
    }
    double lm = log(1.0 + pm);
    for (int iy=0; iy<GRID; iy++)
        for (int ix=0; ix<GRID; ix++)
            if (pd[iy][ix] > 0)
                pix[iy*GRID+ix] = cmap_hot(log(1.0 + pd[iy][ix]) / lm);
}

/* ================================================================
   DRAWING HELPERS
   ================================================================ */

void draw_trail(V2 *buf, int head, int count, int vis,
                int ox, int oy, int sz, double scale, Color col)
{
    int tlen = count < vis ? count : vis;
    int max = count;  /* ring buffer size inferred from usage — caller must ensure buf is big enough */
    for (int k = 1; k < tlen; k++) {
        int i0 = (head - k     + max) % max;
        int i1 = (head - k - 1 + max) % max;
        Vector2 a = phys2screen(buf[i0].x, buf[i0].y, ox, oy, sz, scale);
        Vector2 b = phys2screen(buf[i1].x, buf[i1].y, ox, oy, sz, scale);
        unsigned char al = (unsigned char)(200 * (1.0 - (double)k / tlen));
        DrawLineEx(a, b, 1.5f, (Color){col.r, col.g, col.b, al});
    }
}

void draw_droplet(double x, double y, int ox, int oy, int sz, double scale)
{
    Vector2 dp = phys2screen(x, y, ox, oy, sz, scale);
    DrawCircleV(dp, 5, (Color){255,255,200,120});
    DrawCircleV(dp, 3, WHITE);
}

void draw_turbo(int ox, int oy, int sz, double rate, const char *unit)
{
    DrawRectangle(ox, oy, sz, sz, (Color){30,30,30,255});
    DrawText("TURBO MODE", ox+sz/2-70, oy+sz/2-20, 24, (Color){255,100,100,255});
    char tb[128];
    snprintf(tb, sizeof(tb), "%.0f %s/sec", rate, unit);
    DrawText(tb, ox+sz/2-60, oy+sz/2+15, 18, LIGHTGRAY);
    DrawText("[T] to exit turbo", ox+sz/2-75, oy+sz/2+40, 14, GRAY);
}

void draw_L_strip(double *L_hist, int L_head, int L_count, int L_max,
                  double L_scale, int ox, int oy, int w, int h)
{
    DrawRectangle(ox, oy, w, h, (Color){15,15,20,255});
    int mid = oy + h/2;
    DrawLineV((Vector2){ox,mid}, (Vector2){ox+w,mid}, (Color){40,40,50,255});

    if (L_count > 1) {
        int npts = L_count < w ? L_count : w;
        for (int k = 1; k < npts; k++) {
            int i0 = (L_head - k + L_max) % L_max;
            int i1 = (L_head - k - 1 + L_max) % L_max;
            float x0 = ox + w - k, x1 = ox + w - k - 1;
            float y0 = mid - (float)(L_hist[i0] / L_scale * h * 0.45);
            float y1 = mid - (float)(L_hist[i1] / L_scale * h * 0.45);
            Color lc = L_hist[i0] > 0 ? (Color){100,150,255,200} : (Color){255,100,100,200};
            DrawLineV((Vector2){x0,y0}, (Vector2){x1,y1}, lc);
        }
    }

    char lb[64];
    double Lcur = L_count > 0 ? L_hist[(L_head-1+L_max)%L_max] : 0;
    snprintf(lb, sizeof(lb), "L = %.4f", Lcur);
    DrawText(lb, ox+4, oy+2, 12, LIGHTGRAY);
    DrawText("Angular Momentum", ox+w-150, oy+2, 12, GRAY);
}
