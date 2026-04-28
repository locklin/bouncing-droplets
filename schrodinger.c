/*
 * Quantum Billiard — Schrodinger Eigenstates
 *
 * Solves -nabla^2 psi = k^2 psi with Dirichlet BCs on the corral
 * boundary. Displays eigenstates and |psi|^2 for comparison with
 * the classical billiard histogram and pilot-wave walker histogram.
 *
 * Uses Lanczos algorithm on a coarse grid to find eigenmodes, then
 * interpolates to fine grid for display.
 *
 * Controls:
 *   Up/Down     cycle through eigenstates
 *   G           cycle geometry
 *   V           cycle view (psi / |psi|^2 / sum |psi_n|^2)
 *   R           recompute eigenmodes
 *   N/M         number of modes in thermal sum
 */

#include "core.h"

#define EGRID    120
#define MAX_MODES 200
#define LANCZOS_M 600

static int geometry = GEO_STADIUM;

/* ----- Eigenmodes ----- */

static double *modes[MAX_MODES];     /* each GRID*GRID on fine grid */
static double  mode_k[MAX_MODES];   /* wavenumber k = sqrt(eigenvalue) */
static int     n_modes = 0;
static int     cur_mode = 0;         /* currently displayed mode */
static int     thermal_n = 10;       /* modes in thermal sum */

/* ----- Linear Algebra ----- */

static int  einterior[EGRID*EGRID];

static double vdot(const double *a, const double *b, int n)
{
    double s = 0;
    for (int i = 0; i < n; i++) s += a[i]*b[i];
    return s;
}

static void neg_lap_coarse(double *out, const double *in)
{
    double dx = 2.0/(EGRID-1), idx2 = 1.0/(dx*dx);
    memset(out, 0, EGRID*EGRID*sizeof(double));
    for (int iy=1; iy<EGRID-1; iy++)
        for (int ix=1; ix<EGRID-1; ix++) {
            int i = iy*EGRID+ix;
            if (!einterior[i]) continue;
            out[i] = (4.0*in[i]-in[i-1]-in[i+1]-in[i-EGRID]-in[i+EGRID])*idx2;
        }
}

/* Sturm count: eigenvalues of tridiag < x */
static int sturm_count(const double *d, const double *e, int n, double x)
{
    int count = 0;
    double q = d[0]-x;
    if (q < 0) count++;
    for (int i=1; i<n; i++) {
        if (q == 0) q = 1e-30;
        q = (d[i]-x) - e[i]*e[i]/q;
        if (q < 0) count++;
    }
    return count;
}

static double tridiag_eigval(const double *d, const double *e, int n, int k)
{
    double lo, hi;
    lo = hi = d[0];
    for (int i=0; i<n; i++) {
        double r = (i>0?fabs(e[i]):0) + (i<n-1?fabs(e[i+1]):0);
        if (d[i]-r<lo) lo=d[i]-r; if (d[i]+r>hi) hi=d[i]+r;
    }
    lo -= 1; hi += 1;
    for (int iter=0; iter<100; iter++) {
        double mid=(lo+hi)/2;
        if (sturm_count(d,e,n,mid) <= k) lo=mid; else hi=mid;
    }
    return (lo+hi)/2;
}

static void tridiag_solve(const double *d, const double *e, int n,
                           double shift, const double *b, double *x)
{
    double *dd=malloc(n*sizeof(double));
    double *ll=malloc(n*sizeof(double));
    double *yy=malloc(n*sizeof(double));
    dd[0]=d[0]-shift; yy[0]=b[0];
    for (int i=1;i<n;i++) {
        if (fabs(dd[i-1])<1e-30) dd[i-1]=(dd[i-1]>=0?1e-30:-1e-30);
        ll[i]=e[i]/dd[i-1];
        dd[i]=(d[i]-shift)-ll[i]*e[i];
        yy[i]=b[i]-ll[i]*yy[i-1];
    }
    if (fabs(dd[n-1])<1e-30) dd[n-1]=1e-30;
    x[n-1]=yy[n-1]/dd[n-1];
    for (int i=n-2;i>=0;i--) x[i]=(yy[i]-e[i+1]*x[i+1])/dd[i];
    free(dd); free(ll); free(yy);
}

/* Interpolate coarse EGRID to fine GRID */
static void interp_c2f(const double *coarse, double *fine)
{
    for (int iy=0;iy<GRID;iy++) {
        double fy = (2.0*iy/(GRID-1)-1.0)*0.5+0.5;
        fy *= (EGRID-1);
        int jy=(int)fy; if (jy<0) jy=0; if (jy>=EGRID-1) jy=EGRID-2;
        double sy=fy-jy;
        for (int ix=0;ix<GRID;ix++) {
            double fx = (2.0*ix/(GRID-1)-1.0)*0.5+0.5;
            fx *= (EGRID-1);
            int jx=(int)fx; if (jx<0) jx=0; if (jx>=EGRID-1) jx=EGRID-2;
            double sx=fx-jx;
            fine[iy*GRID+ix] =
                (1-sx)*(1-sy)*coarse[jy*EGRID+jx] + sx*(1-sy)*coarse[jy*EGRID+jx+1]
              + (1-sx)*sy*coarse[(jy+1)*EGRID+jx] + sx*sy*coarse[(jy+1)*EGRID+jx+1];
        }
    }
}

/* ----- Compute Eigenmodes ----- */

static void compute_modes(void)
{
    int EN = EGRID*EGRID;
    printf("Computing eigenmodes for %s on %dx%d grid...\n", geo_name(geometry), EGRID, EGRID);

    /* Build interior mask */
    int n_int = 0;
    for (int iy=0; iy<EGRID; iy++) {
        double y = 2.0*iy/(EGRID-1)-1.0;
        for (int ix=0; ix<EGRID; ix++) {
            double x = 2.0*ix/(EGRID-1)-1.0;
            int i = iy*EGRID+ix;
            einterior[i] = (iy>0 && iy<EGRID-1 && ix>0 && ix<EGRID-1
                            && geo_inside(geometry,x,y)) ? 1 : 0;
            if (einterior[i]) n_int++;
        }
    }
    printf("  Interior points: %d\n", n_int);

    /* Free old modes */
    for (int i=0; i<n_modes; i++) { free(modes[i]); modes[i]=NULL; }
    n_modes = 0; cur_mode = 0;

    int M = LANCZOS_M;
    if (M > n_int) M = n_int;
    double *diag = calloc(M, sizeof(double));
    double *offd = calloc(M+1, sizeof(double));
    double **V = malloc(M * sizeof(double*));
    for (int i=0; i<M; i++) V[i] = calloc(EN, sizeof(double));
    double *w = malloc(EN * sizeof(double));

    srand(42);
    for (int i=0; i<EN; i++)
        V[0][i] = einterior[i] ? ((double)rand()/RAND_MAX-0.5) : 0;
    double nrm = sqrt(vdot(V[0],V[0],EN));
    for (int i=0; i<EN; i++) V[0][i] /= nrm;

    int m_actual = M;
    for (int j=0; j<M; j++) {
        neg_lap_coarse(w, V[j]);
        diag[j] = vdot(V[j], w, EN);
        for (int i=0; i<EN; i++) {
            w[i] -= diag[j]*V[j][i];
            if (j>0) w[i] -= offd[j]*V[j-1][i];
        }
        for (int k=0; k<=j; k++) {
            double c = vdot(V[k],w,EN);
            for (int i=0; i<EN; i++) w[i] -= c*V[k][i];
        }
        double beta = sqrt(vdot(w,w,EN));
        if (beta<1e-12 || j+1>=M) { m_actual=j+1; break; }
        offd[j+1] = beta;
        for (int i=0; i<EN; i++) V[j+1][i] = w[i]/beta;
    }
    printf("  Lanczos: %d iterations\n", m_actual);

    /* Find the first MAX_MODES eigenvalues from the bottom */
    int n_want = MAX_MODES;
    int n_avail = sturm_count(diag, offd, m_actual, 1e15);
    if (n_want > n_avail) n_want = n_avail;

    double *tvec = malloc(m_actual * sizeof(double));
    double *rhs  = malloc(m_actual * sizeof(double));
    double *coarse = malloc(EN * sizeof(double));

    for (int idx=0; idx<n_want; idx++) {
        double lam = tridiag_eigval(diag, offd, m_actual, idx);
        if (lam < 0) continue;  /* skip negative eigenvalues (numerical artifact) */

        /* Inverse iteration for eigenvector */
        for (int i=0; i<m_actual; i++) rhs[i] = (double)rand()/RAND_MAX-0.5;
        tridiag_solve(diag, offd, m_actual, lam, rhs, tvec);
        tridiag_solve(diag, offd, m_actual, lam, tvec, rhs);
        double tn = sqrt(vdot(rhs,rhs,m_actual));
        for (int i=0; i<m_actual; i++) rhs[i] /= tn;

        /* Transform to coarse grid */
        memset(coarse, 0, EN*sizeof(double));
        for (int j=0; j<m_actual; j++) {
            double c = rhs[j];
            for (int i=0; i<EN; i++) coarse[i] += c*V[j][i];
        }
        double pn = sqrt(vdot(coarse,coarse,EN));
        if (pn < 1e-15) continue;
        for (int i=0; i<EN; i++) coarse[i] /= pn;

        /* Interpolate to fine grid */
        double *fine = calloc(GRID*GRID, sizeof(double));
        interp_c2f(coarse, fine);

        modes[n_modes] = fine;
        mode_k[n_modes] = sqrt(fabs(lam));
        n_modes++;
        if (n_modes >= MAX_MODES) break;
    }

    printf("  Found %d eigenmodes (k = %.2f to %.2f)\n",
           n_modes, n_modes>0?mode_k[0]:0, n_modes>0?mode_k[n_modes-1]:0);

    free(tvec); free(rhs); free(coarse); free(w);
    for (int i=0; i<M; i++) free(V[i]);
    free(V); free(diag); free(offd);
}

/* ----- Rendering ----- */

static Color wpix[GRID*GRID], rpix[GRID*GRID];

enum { VIEW_PSI, VIEW_PSI2, VIEW_THERMAL, VIEW_COUNT };
static int view_mode = VIEW_PSI2;

static void render_panels(Texture2D ltex, Texture2D rtex)
{
    if (n_modes == 0) return;
    int cm = cur_mode < n_modes ? cur_mode : n_modes-1;

    /* Left: psi (blue-white-red) */
    double vmin=0, vmax=0;
    double *phi = modes[cm];
    for (int i=0; i<GRID*GRID; i++) {
        if (phi[i]<vmin) vmin=phi[i];
        if (phi[i]>vmax) vmax=phi[i];
    }
    double amax = fmax(fabs(vmin),fabs(vmax));
    if (amax < 1e-20) amax = 1;
    for (int iy=0; iy<GRID; iy++) {
        double y=2.0*iy/(GRID-1)-1.0;
        for (int ix=0; ix<GRID; ix++) {
            double x=2.0*ix/(GRID-1)-1.0;
            int idx=iy*GRID+ix;
            wpix[idx] = !geo_inside(geometry,x,y) ? (Color){30,30,30,255}
                       : cmap_bwr(0.5 + 0.5*phi[idx]/amax);
        }
    }
    UpdateTexture(ltex, wpix);

    /* Right: depends on view mode */
    if (view_mode == VIEW_PSI || view_mode == VIEW_PSI2) {
        /* |psi|^2 for current mode */
        double peak = 0;
        for (int i=0; i<GRID*GRID; i++) {
            double v = phi[i]*phi[i];
            if (v > peak) peak = v;
        }
        if (peak < 1e-20) peak = 1;
        for (int iy=0; iy<GRID; iy++) {
            double y=2.0*iy/(GRID-1)-1.0;
            for (int ix=0; ix<GRID; ix++) {
                double x=2.0*ix/(GRID-1)-1.0;
                int idx=iy*GRID+ix;
                rpix[idx] = !geo_inside(geometry,x,y) ? (Color){30,30,30,255}
                           : cmap_hot(sqrt(phi[idx]*phi[idx]/peak));
            }
        }
    } else if (view_mode == VIEW_THERMAL) {
        /* Sum of |psi_n|^2 for first thermal_n modes */
        int nm = thermal_n < n_modes ? thermal_n : n_modes;
        double peak = 0;
        /* First pass: find peak */
        for (int i=0; i<GRID*GRID; i++) {
            double s = 0;
            for (int m=0; m<nm; m++) s += modes[m][i]*modes[m][i];
            if (s > peak) peak = s;
        }
        if (peak < 1e-20) peak = 1;
        for (int iy=0; iy<GRID; iy++) {
            double y=2.0*iy/(GRID-1)-1.0;
            for (int ix=0; ix<GRID; ix++) {
                double x=2.0*ix/(GRID-1)-1.0;
                int idx=iy*GRID+ix;
                if (!geo_inside(geometry,x,y)) { rpix[idx]=(Color){30,30,30,255}; continue; }
                double s = 0;
                for (int m=0; m<nm; m++) s += modes[m][idx]*modes[m][idx];
                rpix[idx] = cmap_hot(sqrt(s/peak));
            }
        }
    }
    UpdateTexture(rtex, rpix);
}


/* ----- Main ----- */

#define PANEL 400
#define GAP    20
#define WIN_W (PANEL*2+GAP*3)
#define WIN_H (PANEL+GAP*2+60)

int main(void)
{
    InitWindow(WIN_W, WIN_H, "Quantum Billiard - Schrodinger Eigenstates");
    SetTargetFPS(30);

    Image img = GenImageColor(GRID, GRID, BLACK);
    Texture2D ltex = LoadTextureFromImage(img);
    Texture2D rtex = LoadTextureFromImage(img);
    UnloadImage(img);

    compute_modes();

    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_UP))    cur_mode = (cur_mode+1) % (n_modes>0?n_modes:1);
        if (IsKeyPressed(KEY_DOWN))  cur_mode = (cur_mode-1+n_modes) % (n_modes>0?n_modes:1);
        if (IsKeyPressed(KEY_G))     { geometry=(geometry+1)%GEO_COUNT; compute_modes(); }
        if (IsKeyPressed(KEY_V))     view_mode=(view_mode+1)%VIEW_COUNT;
        if (IsKeyPressed(KEY_R))     compute_modes();
        if (IsKeyPressed(KEY_N))     thermal_n = thermal_n<n_modes ? thermal_n+1 : thermal_n;
        if (IsKeyPressed(KEY_M))     thermal_n = thermal_n>1 ? thermal_n-1 : thermal_n;

        render_panels(ltex, rtex);

        BeginDrawing();
        ClearBackground((Color){20,20,25,255});
        int lx=GAP, ly=GAP, rx=GAP*2+PANEL, ry=GAP;
        float scale=(float)PANEL/GRID;

        DrawTextureEx(ltex,(Vector2){lx,ly},0,scale,WHITE);
        geo_draw_outline(geometry,lx,ly,PANEL);
        { char l[64];
          snprintf(l,sizeof(l),"psi_%d  (k=%.2f) [Up/Dn]",
                   cur_mode, n_modes>0?mode_k[cur_mode]:0);
          DrawText(l,lx,ly+PANEL+4,16,LIGHTGRAY); }

        DrawTextureEx(rtex,(Vector2){rx,ry},0,scale,WHITE);
        geo_draw_outline(geometry,rx,ry,PANEL);
        { char l[64];
          if (view_mode == VIEW_THERMAL)
              snprintf(l,sizeof(l),"sum |psi_n|^2 (n=0..%d) [N/M]",thermal_n-1);
          else
              snprintf(l,sizeof(l),"|psi_%d|^2 [V]",cur_mode);
          DrawText(l,rx,ry+PANEL+4,16,LIGHTGRAY); }

        int iy=GAP+PANEL+24; char buf[256];
        snprintf(buf,sizeof(buf),
            "Mode %d/%d  |  k = %.3f  |  k^2 = %.2f  |  %s [G]  |  [V] view  |  [R] recompute",
            cur_mode, n_modes, n_modes>0?mode_k[cur_mode]:0,
            n_modes>0?mode_k[cur_mode]*mode_k[cur_mode]:0, geo_name(geometry));
        DrawText(buf,GAP,iy,14,LIGHTGRAY);
        snprintf(buf,sizeof(buf),
            "Thermal modes: %d [N/M]  |  Grid: %dx%d  |  Lanczos: %d iter",
            thermal_n, EGRID, EGRID, LANCZOS_M);
        DrawText(buf,GAP,iy+18,14,GRAY);

        EndDrawing();
    }

    for (int i=0; i<n_modes; i++) free(modes[i]);
    UnloadTexture(ltex); UnloadTexture(rtex);
    CloseWindow();
    return 0;
}
