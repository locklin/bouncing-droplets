/*
 * Bouncing Droplet — Fourier-Bessel Mode Simulator
 *
 * Oza-Rosales-Bush trajectory equation with Fourier-Bessel mode
 * decomposition. Supports harmonic confinement and/or corral walls.
 *
 *   x_ddot = -eta*x_dot - chi*x - mu*grad_x(h) + wall_force
 *   C_n_dot = -C_n/Me + J_n(r)*cos(n*theta)
 *   S_n_dot = -S_n/Me + J_n(r)*sin(n*theta)
 *
 * Reference: Budanur & Fleury, Chaos 29, 013122 (2019)
 *
 * Controls:
 *   Space       pause / resume
 *   Up/Down     memory Me
 *   Left/Right  coupling mu (*1.1)
 *   H           toggle harmonic potential on/off
 *   G           cycle geometry (none / circle / stadium / rect / D-shape)
 *   V           cycle right panel view
 *   T           toggle turbo
 *   +/-         simulation speed
 *   R           reset
 *   Click       place droplet
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "raylib.h"

/* ----- Configuration ----- */

#define NMODES    25
#define NDIM      (2*NMODES + 5)
#define GRID      200
#define MAX_TRAIL 4000
#define TRAIL_VIS 1000
#define DT        0.01

typedef struct { double x, y; } V2;

/* ----- Geometry ----- */

enum { GEO_STADIUM, GEO_DSHAPE, GEO_CIRCLE, GEO_RECT, GEO_NONE, GEO_COUNT };
static int geometry = GEO_STADIUM;

#define DSHAPE_R   0.9
#define DSHAPE_CUT (-0.3)

/* SDF returns positive inside, negative outside. Scale factor maps
 * geometry coords [-1,1] to simulation coords [-RMAX,RMAX]. */
static double RMAX = 12.0;

static double geo_sdf_unit(double x, double y)
{
    switch (geometry) {
    case GEO_CIRCLE:  return 1.0 - sqrt(x*x + y*y);
    case GEO_STADIUM: {
        double a=0.5, r=0.5, ax=fabs(x);
        return ax <= a ? r-fabs(y) : r-sqrt((ax-a)*(ax-a)+y*y);
    }
    case GEO_RECT:    return fmin(0.8-fabs(x), 0.5-fabs(y));
    case GEO_DSHAPE:  return fmin(DSHAPE_R-sqrt(x*x+y*y), x-DSHAPE_CUT);
    default:          return 1.0;  /* GEO_NONE: always inside */
    }
}

/* SDF in simulation coordinates */
static double geo_sdf(double sx, double sy)
{
    return geo_sdf_unit(sx / RMAX, sy / RMAX) * RMAX;
}

static int geo_inside(double sx, double sy)
{
    return geometry == GEO_NONE || geo_sdf(sx, sy) > 0;
}

static void geo_normal(double sx, double sy, double *nx, double *ny)
{
    double eps = RMAX * 0.002;
    *nx = -(geo_sdf(sx+eps,sy) - geo_sdf(sx-eps,sy)) / (2*eps);
    *ny = -(geo_sdf(sx,sy+eps) - geo_sdf(sx,sy-eps)) / (2*eps);
}

#define S2P(px, py) (Vector2){ \
    ox + ((float)((px)/RMAX)*0.5f+0.5f)*sz, \
    oy + ((float)((py)/RMAX)*0.5f+0.5f)*sz }

static void geo_draw_outline(int ox, int oy, int sz)
{
    int N = 200; Color col = LIGHTGRAY;
    /* Draw in unit coords scaled to screen */
    #define U2P(ux, uy) (Vector2){ \
        ox + ((float)(ux)*0.5f+0.5f)*sz, oy + ((float)(uy)*0.5f+0.5f)*sz }

    switch (geometry) {
    case GEO_NONE: break;
    case GEO_CIRCLE:
        DrawCircleLinesV((Vector2){ox+sz/2.0f,oy+sz/2.0f}, sz/2.0f, col); break;
    case GEO_STADIUM: {
        float a=0.5f, r=0.5f;
        DrawLineV(U2P(-a,r), U2P(a,r), col); DrawLineV(U2P(-a,-r), U2P(a,-r), col);
        for (int i=0;i<N;i++) {
            float t0=PI/2+PI*i/(float)N, t1=PI/2+PI*(i+1)/(float)N;
            DrawLineV(U2P(-a+r*cosf(t0),r*sinf(t0)),U2P(-a+r*cosf(t1),r*sinf(t1)),col);}
        for (int i=0;i<N;i++) {
            float t0=-PI/2+PI*i/(float)N, t1=-PI/2+PI*(i+1)/(float)N;
            DrawLineV(U2P(a+r*cosf(t0),r*sinf(t0)),U2P(a+r*cosf(t1),r*sinf(t1)),col);}
        break; }
    case GEO_RECT: {
        float hw=0.8f, hh=0.5f;
        DrawLineV(U2P(-hw,-hh),U2P(hw,-hh),col); DrawLineV(U2P(hw,-hh),U2P(hw,hh),col);
        DrawLineV(U2P(hw,hh),U2P(-hw,hh),col); DrawLineV(U2P(-hw,hh),U2P(-hw,-hh),col);
        break; }
    case GEO_DSHAPE: {
        float R=DSHAPE_R, yint=sqrtf(R*R-DSHAPE_CUT*DSHAPE_CUT);
        float a0=atan2f(yint,DSHAPE_CUT),a1=atan2f(-yint,DSHAPE_CUT),span=a0-a1;
        DrawLineV(U2P(DSHAPE_CUT,-yint),U2P(DSHAPE_CUT,yint),col);
        for (int i=0;i<N;i++) {
            float t0=a0-span*i/(float)N, t1=a0-span*(i+1)/(float)N;
            DrawLineV(U2P(R*cosf(t0),R*sinf(t0)),U2P(R*cosf(t1),R*sinf(t1)),col);}
        break; }
    }
    #undef U2P
}
#undef S2P

static const char *geo_name(void) {
    const char *n[] = {"Stadium","D-shape","Circle","Rectangle","None (harmonic)"};
    return geometry < GEO_COUNT ? n[geometry] : "?";
}

/* ----- Parameters ----- */

static struct {
    double eta, chi, mu, Me;
    int harmonic;   /* harmonic potential on/off */
} P = { 0.2, 0.008, 0.0375482401231, 15.0, 1 };

/* ----- State ----- */

static double state[NDIM];
static V2     trail[MAX_TRAIL];
static int    n_trail = 0, trail_head = 0;
static double sim_time = 0;
static int    paused = 0, speed = 50, turbo = 0;
static double turbo_sps = 0;

static inline int iC(int n) { return n == 0 ? 4 : 3 + 2*n; }
static inline int iS(int n) { return 4 + 2*n; }

/* ----- Precomputed Grid ----- */

static double grid_r[GRID * GRID];
static double grid_ct[GRID * GRID];     /* cos(theta) */
static double grid_st[GRID * GRID];     /* sin(theta) */
static double grid_jn[NMODES+1][GRID * GRID];

static void precompute_grid(void)
{
    for (int iy = 0; iy < GRID; iy++) {
        double y = (2.0*iy/(GRID-1) - 1.0) * RMAX;
        for (int ix = 0; ix < GRID; ix++) {
            double x = (2.0*ix/(GRID-1) - 1.0) * RMAX;
            int idx = iy * GRID + ix;
            double r = sqrt(x*x + y*y);
            grid_r[idx] = r;
            grid_ct[idx] = r > 1e-12 ? x/r : 1.0;
            grid_st[idx] = r > 1e-12 ? y/r : 0.0;
            for (int n = 0; n <= NMODES; n++)
                grid_jn[n][idx] = jn(n, r);
        }
    }
    printf("Grid precomputed (%dx%d, %d modes)\n", GRID, GRID, NMODES+1);
}

/* ----- Display ----- */

static double wave_grid[GRID][GRID];
static int    hist[GRID][GRID];
static int    hpeak = 1;
static double h2_sum[GRID][GRID], h4_sum[GRID][GRID];
static int    moment_n = 0;

/* Poincare section: y=0 upward crossings, record (x, vx) */
#define MAX_POINCARE 200000
static V2     poincare_pts[MAX_POINCARE];
static int    n_poincare = 0, poincare_head = 0;
static double poincare_vmax = 1.0;

/* Angular momentum time series */
#define L_HISTORY 4000
static double L_hist[L_HISTORY];
static double L_max = 0.01;       /* auto-scaling */
static int    L_head = 0, L_count = 0;

/* Orbit diagram: (Me, y_tilde_P) from Poincare crossings */
#define MAX_ORBIT 200000
static V2     orbit_pts[MAX_ORBIT];
static int    n_orbit = 0;
static double orbit_Me_min = 100, orbit_Me_max = 0; /* auto range */

/* Auto-sweep */
static int    sweeping = 0;
static double sweep_from = 10.0, sweep_to = 20.0, sweep_step = 0.2;
static int    sweep_settle = 2000;   /* steps to settle per Me */
static int    sweep_record = 50;     /* Poincare crossings to record per Me */
static int    sweep_phase = 0;       /* 0=settling, 1=recording */
static int    sweep_counter = 0;

enum { VIEW_HIST, VIEW_BORN, VIEW_KURTOSIS, VIEW_POINCARE, VIEW_ORBIT, VIEW_COUNT };
static int view_mode = VIEW_HIST;

static Color wpix[GRID*GRID], rpix[GRID*GRID];

/* ----- Colormaps ----- */

static Color cmap_bwr(double t) {
    if (t < 0) t = 0; if (t > 1) t = 1;
    if (t < 0.5) { unsigned char v=(unsigned char)(t*2*255); return (Color){v,v,255,255}; }
    else { unsigned char v=(unsigned char)((1-t)*2*255); return (Color){255,v,v,255}; }
}

static Color cmap_hot(double t) {
    if (t < 0) t = 0; if (t > 1) t = 1;
    unsigned char r,g,b;
    if (t<0.33) { r=(unsigned char)(t/0.33*255); g=0; b=0; }
    else if (t<0.66) { r=255; g=(unsigned char)((t-0.33)/0.33*255); b=0; }
    else { r=255; g=255; b=(unsigned char)((t-0.66)/0.34*255); }
    return (Color){r,g,b,255};
}

/* ----- ODE RHS ----- */

static void rhs(const double *s, double *ds)
{
    double x = s[0], y = s[1], vx = s[2], vy = s[3];
    double r = sqrt(x*x + y*y);
    double ct, st;
    if (r > 1e-12) { ct = x/r; st = y/r; }
    else           { ct = 1.0;  st = 0.0; }

    double Jn[NMODES + 2];
    Jn[0] = j0(r); Jn[1] = j1(r);
    for (int n = 2; n <= NMODES; n++) Jn[n] = jn(n, r);

    double cosnt[NMODES+2], sinnt[NMODES+2];
    cosnt[0] = 1; sinnt[0] = 0;
    cosnt[1] = ct; sinnt[1] = st;
    for (int n = 2; n <= NMODES+1; n++) {
        cosnt[n] = 2*ct*cosnt[n-1] - cosnt[n-2];
        sinnt[n] = 2*ct*sinnt[n-1] - sinnt[n-2];
    }

    double fx = 0, fy = 0;
    for (int n = 0; n <= NMODES; n++) {
        double Cn = s[iC(n)];
        double Sn = (n >= 1) ? s[iS(n)] : 0.0;
        double w = (n == 0) ? 1.0 : 2.0;
        double Kn = Cn*cosnt[n] + Sn*sinnt[n];
        double Ln = Cn*cosnt[n+1] + Sn*sinnt[n+1];
        double Mn = Cn*cosnt[n+1] - Sn*sinnt[n+1];
        double Jnm1 = (n == 0) ? -Jn[1] : Jn[n-1];
        double nJnr;
        if (r > 1e-12) nJnr = (double)n * Jn[n] / r;
        else if (n==1) nJnr = 0.5;
        else           nJnr = 0.0;
        fx += w * (Jnm1*Kn*ct - nJnr*Ln);
        fy += w * (Jnm1*Kn*st + nJnr*Mn);
    }

    double chi_eff = P.harmonic ? P.chi : 0.0;
    ds[0] = vx;
    ds[1] = vy;
    ds[2] = -P.eta*vx - chi_eff*x - P.mu*fx;
    ds[3] = -P.eta*vy - chi_eff*y - P.mu*fy;

    double Me_inv = 1.0 / P.Me;
    ds[iC(0)] = -Me_inv * s[iC(0)] + Jn[0];
    for (int n = 1; n <= NMODES; n++) {
        ds[iC(n)] = -Me_inv * s[iC(n)] + Jn[n]*cosnt[n];
        ds[iS(n)] = -Me_inv * s[iS(n)] + Jn[n]*sinnt[n];
    }
}

/* ----- RK4 ----- */

static void rk4_step(double *y, double dt)
{
    double k1[NDIM], k2[NDIM], k3[NDIM], k4[NDIM], tmp[NDIM];
    rhs(y, k1);
    for (int i=0;i<NDIM;i++) tmp[i] = y[i] + 0.5*dt*k1[i];
    rhs(tmp, k2);
    for (int i=0;i<NDIM;i++) tmp[i] = y[i] + 0.5*dt*k2[i];
    rhs(tmp, k3);
    for (int i=0;i<NDIM;i++) tmp[i] = y[i] + dt*k3[i];
    rhs(tmp, k4);
    for (int i=0;i<NDIM;i++)
        y[i] += dt/6.0 * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]);
}

/* ----- Step with wall confinement ----- */

static void sim_step(void)
{
    double old_x = state[0], old_y = state[1], old_vx = state[2];
    rk4_step(state, DT);
    sim_time += DT;

    /* Wall confinement (when geometry != NONE) */
    if (geometry != GEO_NONE) {
        double d = geo_sdf(state[0], state[1]);
        double wall_zone = RMAX * 0.15;
        if (d < wall_zone) {
            double nx, ny;
            geo_normal(state[0], state[1], &nx, &ny);
            double nm = sqrt(nx*nx + ny*ny);
            if (nm > 1e-10) {
                nx /= nm; ny /= nm;
                double s = (wall_zone - d) / wall_zone;
                double push = 0.5 * s * s;
                state[2] -= push * nx;  /* apply as velocity change */
                state[3] -= push * ny;
                /* Damp velocity into wall */
                double vn = state[2]*nx + state[3]*ny;
                if (vn > 0) { state[2] -= 1.2*vn*nx; state[3] -= 1.2*vn*ny; }
            }
        }
        /* Hard clamp */
        d = geo_sdf(state[0], state[1]);
        if (d < RMAX * 0.02) {
            double nx, ny;
            geo_normal(state[0], state[1], &nx, &ny);
            double nm = sqrt(nx*nx + ny*ny);
            if (nm > 1e-10) {
                state[0] -= (RMAX*0.02 - d) * nx/nm;
                state[1] -= (RMAX*0.02 - d) * ny/nm;
            }
        }
    }

    /* Angular momentum: L = x*vy - y*vx */
    {
        double L = state[0]*state[3] - state[1]*state[2];
        L_hist[L_head] = L;
        L_head = (L_head + 1) % L_HISTORY;
        if (L_count < L_HISTORY) L_count++;
        double aL = fabs(L);
        if (aL > L_max) L_max = aL * 1.1;
        /* Slow decay of L_max for auto-scaling */
        L_max *= 0.99999;
        if (L_max < 0.01) L_max = 0.01;
    }

    /* Poincare section: y crosses zero upward */
    if (old_y <= 0 && state[1] > 0) {
        double dy = state[1] - old_y;
        double frac = (dy > 1e-15) ? -old_y / dy : 0.5;
        double x_cross = old_x + frac * (state[0] - old_x);
        double vx_cross = old_vx + frac * (state[2] - old_vx);
        double av = fabs(vx_cross);
        if (av > poincare_vmax) poincare_vmax = av * 1.1;
        poincare_pts[poincare_head] = (V2){ x_cross, vx_cross };
        poincare_head = (poincare_head + 1) % MAX_POINCARE;
        if (n_poincare < MAX_POINCARE) n_poincare++;

        /* Orbit diagram: record (Me, y_P) where y_P is symmetry-reduced y
         * For the symmetry-reduced section, y_P = y at Poincare crossing
         * (which is ~0 by definition), so use x instead as in Budanur Fig 3.
         * Actually Budanur uses y_tilde_P from the reduced coordinates.
         * Simple version: record x at crossing. */
        double y_cross = 0; /* by definition at y=0 crossing */
        (void)y_cross;
        if (n_orbit < MAX_ORBIT) {
            orbit_pts[n_orbit++] = (V2){ P.Me, x_cross };
            if (P.Me < orbit_Me_min) orbit_Me_min = P.Me;
            if (P.Me > orbit_Me_max) orbit_Me_max = P.Me;
        }

        /* Auto-sweep: count crossings in record phase */
        if (sweeping && sweep_phase == 1) {
            sweep_counter++;
            if (sweep_counter >= sweep_record) {
                /* Advance Me */
                P.Me += sweep_step;
                if (P.Me > sweep_to + 0.01) {
                    sweeping = 0;  /* sweep done */
                } else {
                    sweep_phase = 0;
                    sweep_counter = 0;
                }
            }
        }
    }

    /* Auto-sweep: count settle steps */
    if (sweeping && sweep_phase == 0) {
        sweep_counter++;
        if (sweep_counter >= sweep_settle) {
            sweep_phase = 1;
            sweep_counter = 0;
        }
    }

    /* Trail */
    trail[trail_head] = (V2){ state[0], state[1] };
    trail_head = (trail_head + 1) % MAX_TRAIL;
    if (n_trail < MAX_TRAIL) n_trail++;

    /* Histogram (subsample) */
    static int hcnt = 0;
    if (++hcnt >= 10) {
        hcnt = 0;
        double hxf = (state[0]/RMAX*0.5+0.5)*GRID;
        double hyf = (state[1]/RMAX*0.5+0.5)*GRID;
        int cx=(int)hxf, cy=(int)hyf;
        for (int dy=-2; dy<=2; dy++)
            for (int dx=-2; dx<=2; dx++) {
                int px=cx+dx, py=cy+dy;
                if (px<0||px>=GRID||py<0||py>=GRID) continue;
                double gg = (px-hxf)*(px-hxf)+(py-hyf)*(py-hyf);
                hist[py][px] += (int)(100*exp(-gg));
                if (hist[py][px]>hpeak) hpeak=hist[py][px];
            }
    }
}

/* ----- Wave Reconstruction (uses precomputed grid) ----- */

static void reconstruct_wave(void)
{
    for (int idx = 0; idx < GRID*GRID; idx++) {
        double ct = grid_ct[idx], st = grid_st[idx];
        /* cos(n*theta), sin(n*theta) via recurrence */
        double cn0=1, cn1=ct, sn0=0, sn1=st;
        double h = 0;
        for (int n = 0; n <= NMODES; n++) {
            double cosn, sinn;
            if      (n == 0) { cosn = 1;   sinn = 0;   }
            else if (n == 1) { cosn = ct;  sinn = st;  }
            else {
                cosn = 2*ct*cn1 - cn0; sinn = 2*ct*sn1 - sn0;
                cn0 = cn1; cn1 = cosn; sn0 = sn1; sn1 = sinn;
            }
            double Cn = state[iC(n)];
            double Sn = (n >= 1) ? state[iS(n)] : 0;
            double w = (n == 0) ? 1.0 : 2.0;
            h += w * grid_jn[n][idx] * (Cn*cosn + Sn*sinn);
        }
        int iy = idx / GRID, ix = idx % GRID;
        wave_grid[iy][ix] = h;
        double h2 = h*h;
        h2_sum[iy][ix] += h2;
        h4_sum[iy][ix] += h2*h2;
    }
    moment_n++;
}

/* ----- Rendering ----- */

static void render_wave_tex(Texture2D tex)
{
    double vmin=0, vmax=0;
    for (int i=0; i<GRID*GRID; i++) {
        double v = wave_grid[i/GRID][i%GRID];
        if (v<vmin) vmin=v; if (v>vmax) vmax=v;
    }
    double range = vmax-vmin; if (range<1e-10) range=1;
    for (int iy=0; iy<GRID; iy++) {
        double y = (2.0*iy/(GRID-1)-1.0)*RMAX;
        for (int ix=0; ix<GRID; ix++) {
            double x = (2.0*ix/(GRID-1)-1.0)*RMAX;
            int idx = iy*GRID+ix;
            if (!geo_inside(x,y))
                wpix[idx] = (Color){30,30,30,255};
            else
                wpix[idx] = cmap_bwr((wave_grid[iy][ix]-vmin)/range);
        }
    }
    UpdateTexture(tex, wpix);
}

static void render_right_tex(Texture2D tex)
{
    if (view_mode == VIEW_HIST) {
        double lmin = 1e30, lmax = 0;
        for (int i = 0; i < GRID*GRID; i++)
            if (hist[i/GRID][i%GRID] > 0) {
                double lv = log(1.0 + hist[i/GRID][i%GRID]);
                if (lv < lmin) lmin = lv;
                if (lv > lmax) lmax = lv;
            }
        double lrange = lmax - lmin; if (lrange < 1e-10) lrange = 1;
        for (int iy=0;iy<GRID;iy++) {
            double y=(2.0*iy/(GRID-1)-1.0)*RMAX;
            for (int ix=0;ix<GRID;ix++) {
                double x=(2.0*ix/(GRID-1)-1.0)*RMAX;
                int idx=iy*GRID+ix;
                if (!geo_inside(x,y)) rpix[idx] = (Color){30,30,30,255};
                else if (hist[iy][ix] == 0) rpix[idx] = BLACK;
                else rpix[idx] = cmap_hot((log(1.0+(double)hist[iy][ix]) - lmin) / lrange);
            }
        }
    } else if (view_mode == VIEW_BORN) {
        if (moment_n<1) { UpdateTexture(tex,rpix); return; }
        double peak=0;
        for (int i=0;i<GRID*GRID;i++) {
            double v=h2_sum[i/GRID][i%GRID]/moment_n; if (v>peak) peak=v; }
        if (peak<1e-20) peak=1;
        for (int iy=0;iy<GRID;iy++) {
            double y=(2.0*iy/(GRID-1)-1.0)*RMAX;
            for (int ix=0;ix<GRID;ix++) {
                double x=(2.0*ix/(GRID-1)-1.0)*RMAX;
                int idx=iy*GRID+ix;
                rpix[idx] = !geo_inside(x,y) ? (Color){30,30,30,255}
                           : cmap_hot(sqrt(h2_sum[iy][ix]/moment_n/peak));
            }
        }
    } else if (view_mode == VIEW_KURTOSIS) {
        if (moment_n<10) { UpdateTexture(tex,rpix); return; }
        double vmin=0,vmax=0;
        for (int iy=0;iy<GRID;iy++) for (int ix=0;ix<GRID;ix++) {
            double m2=h2_sum[iy][ix]/moment_n, m4=h4_sum[iy][ix]/moment_n;
            double ek=m4-3.0*m2*m2; if (ek<vmin)vmin=ek; if(ek>vmax)vmax=ek; }
        double amax=fmax(fabs(vmin),fabs(vmax)); if (amax<1e-20) amax=1;
        for (int iy=0;iy<GRID;iy++) {
            double y=(2.0*iy/(GRID-1)-1.0)*RMAX;
            for (int ix=0;ix<GRID;ix++) {
                double x=(2.0*ix/(GRID-1)-1.0)*RMAX;
                int idx=iy*GRID+ix;
                if (!geo_inside(x,y)) { rpix[idx]=(Color){30,30,30,255}; continue; }
                double m2=h2_sum[iy][ix]/moment_n, m4=h4_sum[iy][ix]/moment_n;
                rpix[idx]=cmap_bwr(0.5+0.5*(m4-3.0*m2*m2)/amax);
            }
        }
    } else if (view_mode == VIEW_POINCARE) {
        /* Poincare section: x on horizontal, vx on vertical (y=0 crossings)
         * Auto-range from actual data so nothing falls off-screen. */
        for (int i = 0; i < GRID*GRID; i++)
            rpix[i] = (Color){15, 15, 20, 255};

        if (n_poincare > 0) {
            /* Find data range */
            double x_lo=1e30, x_hi=-1e30, v_lo=1e30, v_hi=-1e30;
            for (int i = 0; i < n_poincare; i++) {
                int bi = (poincare_head - n_poincare + i + MAX_POINCARE) % MAX_POINCARE;
                double px = poincare_pts[bi].x, pv = poincare_pts[bi].y;
                if (px < x_lo) x_lo = px; if (px > x_hi) x_hi = px;
                if (pv < v_lo) v_lo = pv; if (pv > v_hi) v_hi = pv;
            }
            /* Add 10% margin, ensure symmetric around zero for vx */
            double xm = (x_hi-x_lo)*0.1+0.1; x_lo-=xm; x_hi+=xm;
            double va = fmax(fabs(v_lo), fabs(v_hi)) * 1.1 + 0.01;
            v_lo = -va; v_hi = va;

            /* Axes at zero */
            if (x_lo < 0 && x_hi > 0) {
                int zx = (int)((0-x_lo)/(x_hi-x_lo)*(GRID-1));
                if (zx>=0 && zx<GRID)
                    for (int iy=0; iy<GRID; iy++)
                        rpix[iy*GRID+zx] = (Color){40,40,50,255};
            }
            { int zy = GRID/2; /* vx=0 is always at center since symmetric */
              for (int ix=0; ix<GRID; ix++)
                  rpix[zy*GRID+ix] = (Color){40,40,50,255};
            }

            /* Bin into density grid */
            static int pden[GRID][GRID];
            memset(pden, 0, sizeof(pden));
            int pmax = 1;
            for (int i = 0; i < n_poincare; i++) {
                int bi = (poincare_head - n_poincare + i + MAX_POINCARE) % MAX_POINCARE;
                int gx = (int)((poincare_pts[bi].x - x_lo) / (x_hi - x_lo) * (GRID-1));
                int gy = (int)((1.0 - (poincare_pts[bi].y - v_lo) / (v_hi - v_lo)) * (GRID-1));
                if (gx>=0 && gx<GRID && gy>=0 && gy<GRID) {
                    pden[gy][gx]++;
                    if (pden[gy][gx] > pmax) pmax = pden[gy][gx];
                }
            }
            double llmax = log(1.0 + pmax);
            for (int iy=0; iy<GRID; iy++)
                for (int ix=0; ix<GRID; ix++)
                    if (pden[iy][ix] > 0)
                        rpix[iy*GRID+ix] = cmap_hot(log(1.0+pden[iy][ix]) / llmax);
        }
    } else if (view_mode == VIEW_ORBIT) {
        /* Orbit diagram: Me on x-axis, x_P on y-axis */
        for (int i = 0; i < GRID*GRID; i++)
            rpix[i] = (Color){15, 15, 20, 255};
        if (n_orbit > 0) {
            double me_lo = orbit_Me_min - 0.5, me_hi = orbit_Me_max + 0.5;
            if (me_hi - me_lo < 2) { me_lo -= 1; me_hi += 1; }
            /* Find x range */
            double x_lo = 1e30, x_hi = -1e30;
            for (int i = 0; i < n_orbit; i++) {
                if (orbit_pts[i].y < x_lo) x_lo = orbit_pts[i].y;
                if (orbit_pts[i].y > x_hi) x_hi = orbit_pts[i].y;
            }
            double x_margin = (x_hi - x_lo) * 0.1 + 0.1;
            x_lo -= x_margin; x_hi += x_margin;

            /* Bin into density grid */
            static int oden[GRID][GRID];
            memset(oden, 0, sizeof(oden));
            int omax = 1;
            for (int i = 0; i < n_orbit; i++) {
                int px = (int)((orbit_pts[i].x - me_lo) / (me_hi - me_lo) * (GRID-1));
                int py = (int)((1.0 - (orbit_pts[i].y - x_lo) / (x_hi - x_lo)) * (GRID-1));
                if (px >= 0 && px < GRID && py >= 0 && py < GRID) {
                    oden[py][px]++;
                    if (oden[py][px] > omax) omax = oden[py][px];
                }
            }
            double lmax = log(1.0 + omax);
            for (int iy = 0; iy < GRID; iy++)
                for (int ix = 0; ix < GRID; ix++)
                    if (oden[iy][ix] > 0)
                        rpix[iy*GRID+ix] = cmap_hot(log(1.0 + oden[iy][ix]) / lmax);
            /* Zero line */
            if (x_lo < 0 && x_hi > 0) {
                int zy = (int)((1.0 - (0 - x_lo) / (x_hi - x_lo)) * (GRID-1));
                if (zy >= 0 && zy < GRID)
                    for (int ix = 0; ix < GRID; ix++)
                        if (oden[zy][ix] == 0)
                            rpix[zy*GRID+ix] = (Color){40, 40, 50, 255};
            }
        }
    }
    UpdateTexture(tex, rpix);
}

static const char *view_name(void) {
    const char *n[]={"Histogram","Born <h^2>","Excess Kurtosis","Poincare (x,vx)","Orbit Diagram"};
    return view_mode<VIEW_COUNT ? n[view_mode] : "?";
}

static Vector2 p2s(double px, double py, int ox, int oy, int sz) {
    return (Vector2){ ox+(float)((px/RMAX*0.5+0.5)*sz),
                      oy+(float)((py/RMAX*0.5+0.5)*sz) };
}

static void reset_sim(void)
{
    memset(state, 0, sizeof(state));
    if (geometry == GEO_NONE) {
        state[0] = 5.0; state[3] = 0.5;  /* circular-ish orbit */
    } else {
        state[0] = RMAX * 0.25; state[3] = 0.3;
    }
    n_trail=0; trail_head=0; sim_time=0;
    hpeak=1; moment_n=0;
    n_poincare=0; poincare_head=0; poincare_vmax=1.0;
    L_head=0; L_count=0; L_max=0.01;
    n_orbit=0; orbit_Me_min=100; orbit_Me_max=0;
    sweeping=0;
    memset(hist,0,sizeof(hist));
    memset(h2_sum,0,sizeof(h2_sum));
    memset(h4_sum,0,sizeof(h4_sum));
    memset(wave_grid,0,sizeof(wave_grid));
}

/* ----- Main ----- */

#define PANEL  400
#define GAP     20
#define LSTRIP  80      /* angular momentum strip height */
#define WIN_W  (PANEL*2+GAP*3)
#define WIN_H  (PANEL+LSTRIP+GAP*3+40)

int main(void)
{
    InitWindow(WIN_W, WIN_H, "Bouncing Droplet - Fourier-Bessel Modes");
    SetTargetFPS(30);

    precompute_grid();

    Image img = GenImageColor(GRID, GRID, BLACK);
    Texture2D wtex = LoadTextureFromImage(img);
    Texture2D rtex = LoadTextureFromImage(img);
    UnloadImage(img);

    reset_sim();

    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_SPACE)) paused=!paused;
        if (IsKeyPressed(KEY_R))     reset_sim();
        if (IsKeyPressed(KEY_G))     { geometry=(geometry+1)%GEO_COUNT; reset_sim(); }
        if (IsKeyPressed(KEY_H))     { P.harmonic=!P.harmonic; }
        if (IsKeyPressed(KEY_V))     view_mode=(view_mode+1)%VIEW_COUNT;
        if (IsKeyPressed(KEY_T))     { turbo=!turbo; SetTargetFPS(turbo?10:30); }
        if (IsKeyPressed(KEY_UP))  {
            P.Me = fmin(P.Me+0.5, 50.0);
            n_poincare=0; poincare_head=0; poincare_vmax=1.0;
        }
        if (IsKeyPressed(KEY_DOWN)) {
            P.Me = fmax(P.Me-0.5, 1.0);
            n_poincare=0; poincare_head=0; poincare_vmax=1.0;
        }
        if (IsKeyPressed(KEY_RIGHT)) P.mu *= 1.1;
        if (IsKeyPressed(KEY_LEFT))  P.mu /= 1.1;
        if (IsKeyPressed(KEY_EQUAL)) speed = speed<500 ? speed+50 : speed;
        if (IsKeyPressed(KEY_MINUS)) speed = speed>50  ? speed-50 : speed;
        if (IsKeyPressed(KEY_O)) {
            /* Start orbit diagram sweep from current Me to Me+10 */
            sweeping = 1;
            sweep_from = P.Me;
            sweep_to = P.Me + 10.0;
            sweep_step = 0.2;
            sweep_phase = 0;
            sweep_counter = 0;
            view_mode = VIEW_ORBIT;
            speed = 500;  /* crank speed for sweep */
            printf("Orbit sweep: Me %.1f -> %.1f (step %.1f)\n", sweep_from, sweep_to, sweep_step);
        }

        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) && !turbo) {
            Vector2 m = GetMousePosition();
            if (m.x>=GAP && m.x<GAP+PANEL && m.y>=GAP && m.y<GAP+PANEL) {
                double px = ((m.x-GAP)/PANEL*2-1)*RMAX;
                double py = ((m.y-GAP)/PANEL*2-1)*RMAX;
                if (geo_inside(px,py)) {
                    memset(state,0,sizeof(state));
                    state[0]=px; state[1]=py;
                    double r=sqrt(px*px+py*py);
                    if (r>0.1) { state[2]=-0.5*py/r; state[3]=0.5*px/r; }
                    n_trail=0; trail_head=0; sim_time=0;
                }
            }
        }

        if (!paused) {
            if (turbo) {
                struct timespec t0,t1;
                clock_gettime(CLOCK_MONOTONIC,&t0);
                int count=0;
                for (;;) { sim_step(); count++;
                    if (count%500==0) {
                        clock_gettime(CLOCK_MONOTONIC,&t1);
                        double el=(t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)*1e-9;
                        if (el>0.2) { turbo_sps=count/el; break; }
                    }
                }
                reconstruct_wave();
            } else {
                for (int i=0; i<speed; i++) sim_step();
                reconstruct_wave();
            }
        }

        if (!turbo) render_wave_tex(wtex);
        render_right_tex(rtex);

        BeginDrawing();
        ClearBackground((Color){20,20,25,255});
        int lx=GAP, ly=GAP, rx=GAP*2+PANEL, ry=GAP;
        float scale=(float)PANEL/GRID;

        if (turbo) {
            DrawRectangle(lx,ly,PANEL,PANEL,(Color){30,30,30,255});
            geo_draw_outline(lx,ly,PANEL);
            DrawText("TURBO MODE",lx+130,ly+160,24,(Color){255,100,100,255});
            char tb[128]; snprintf(tb,sizeof(tb),"%.0f steps/sec",turbo_sps);
            DrawText(tb,lx+140,ly+195,18,LIGHTGRAY);
            DrawText("[T] exit turbo",lx+140,ly+225,14,GRAY);
        } else {
            DrawTextureEx(wtex,(Vector2){lx,ly},0,scale,WHITE);
            geo_draw_outline(lx,ly,PANEL);

            int tlen = n_trail<TRAIL_VIS ? n_trail : TRAIL_VIS;
            for (int k=1;k<tlen;k++) {
                int i0=(trail_head-k+MAX_TRAIL)%MAX_TRAIL;
                int i1=(trail_head-k-1+MAX_TRAIL)%MAX_TRAIL;
                Vector2 a=p2s(trail[i0].x,trail[i0].y,lx,ly,PANEL);
                Vector2 b=p2s(trail[i1].x,trail[i1].y,lx,ly,PANEL);
                unsigned char al=(unsigned char)(200*(1.0-(double)k/tlen));
                DrawLineEx(a,b,1.5f,(Color){255,255,100,al});
            }
            Vector2 dp=p2s(state[0],state[1],lx,ly,PANEL);
            DrawCircleV(dp,5,(Color){255,255,200,120});
            DrawCircleV(dp,3,WHITE);

            /* Origin marker */
            Vector2 org=p2s(0,0,lx,ly,PANEL);
            DrawLineV((Vector2){org.x-4,org.y},(Vector2){org.x+4,org.y},GRAY);
            DrawLineV((Vector2){org.x,org.y-4},(Vector2){org.x,org.y+4},GRAY);
        }
        DrawText(turbo?"TURBO [T]":"Wave Field",lx,ly+PANEL+4,16,LIGHTGRAY);

        DrawTextureEx(rtex,(Vector2){rx,ry},0,scale,WHITE);
        if (view_mode == VIEW_ORBIT) {
            DrawRectangleLines(rx,ry,PANEL,PANEL,GRAY);
            if (sweeping) {
                char sb[64]; snprintf(sb,sizeof(sb),"Sweeping Me=%.1f",P.Me);
                DrawText(sb,rx+4,ry+4,14,(Color){255,200,100,255});
            }
        } else if (view_mode == VIEW_POINCARE) {
            DrawRectangleLines(rx,ry,PANEL,PANEL,GRAY);
            DrawText("x",rx+PANEL-12,ry+PANEL+4,12,GRAY);
            DrawText("vx",rx+2,ry+2,12,GRAY);
        } else {
            /* Histogram, Born, Kurtosis — show geometry outline */
            geo_draw_outline(rx,ry,PANEL);
        }
        { char l[96]; snprintf(l,sizeof(l),"%s [V]",view_name());
          DrawText(l,rx,ry+PANEL+4,16,LIGHTGRAY); }

        /* Angular momentum strip */
        int Ly = GAP + PANEL + 20;
        int Lw = PANEL*2 + GAP;  /* full width */
        DrawRectangle(lx, Ly, Lw, LSTRIP, (Color){15,15,20,255});
        /* Zero line */
        int Lmid = Ly + LSTRIP/2;
        DrawLineV((Vector2){lx,Lmid},(Vector2){lx+Lw,Lmid},(Color){40,40,50,255});
        /* Plot L(t) */
        if (L_count > 1) {
            for (int k = 1; k < L_count && k < Lw; k++) {
                int i0 = (L_head - k + L_HISTORY) % L_HISTORY;
                int i1 = (L_head - k - 1 + L_HISTORY) % L_HISTORY;
                float x0 = lx + Lw - k;
                float x1 = lx + Lw - k - 1;
                float y0 = Lmid - (float)(L_hist[i0] / L_max * LSTRIP * 0.45);
                float y1 = Lmid - (float)(L_hist[i1] / L_max * LSTRIP * 0.45);
                Color lc = L_hist[i0] > 0 ? (Color){100,150,255,200} : (Color){255,100,100,200};
                DrawLineV((Vector2){x0,y0},(Vector2){x1,y1},lc);
            }
        }
        /* Labels */
        { char lb[64];
          double Lcur = L_count > 0 ? L_hist[(L_head-1+L_HISTORY)%L_HISTORY] : 0;
          snprintf(lb,sizeof(lb),"L = %.4f", Lcur);
          DrawText(lb, lx+4, Ly+2, 12, LIGHTGRAY);
          DrawText("Angular Momentum", lx+Lw-150, Ly+2, 12, GRAY);
        }

        /* Info bar */
        int iy = Ly + LSTRIP + 4;
        char buf[256];
        double r=sqrt(state[0]*state[0]+state[1]*state[1]);
        double v=sqrt(state[2]*state[2]+state[3]*state[3]);
        snprintf(buf,sizeof(buf),
            "t: %.1f  |  Me: %.1f [Up/Dn]  |  r: %.2f  v: %.3f  |  "
            "speed: %d [+/-]  |  %s%s",
            sim_time, P.Me, r, v, speed,
            paused?"PAUSED":"[Space]",
            sweeping?" | SWEEPING [O]":"");
        DrawText(buf,GAP,iy,14,LIGHTGRAY);
        snprintf(buf,sizeof(buf),
            "%s [G]  |  harm: %s [H]  |  Poincare: %d  |  Orbit: %d  |  [O] sweep  |  [T] turbo",
            geo_name(), P.harmonic?"ON":"OFF", n_poincare, n_orbit);
        DrawText(buf,GAP,iy+16,14,GRAY);

        EndDrawing();
    }

    UnloadTexture(wtex); UnloadTexture(rtex);
    CloseWindow();
    return 0;
}
