/*
 * Walker — Level 0 (Direct Bessel, first-order discrete map)
 */
#include "core.h"

#define MAX_MEM   512
#define TRAIL_VIS 200
#define MAX_POINCARE 100000

static struct { double kF, mu, beta; } P = { 12.0, 0.95, 0.003 };
static int geometry = GEO_STADIUM;

static V2     mem[MAX_MEM];
static int    nmem = 0, mhead = 0;
static V2     drop = {0.25, 0.0}, prev_drop = {0.25, 0.0};
static int    total = 0, paused = 0, speed = 10, turbo = 0;
static double turbo_bps = 0;

static double wave[GRID][GRID];
static int    hist[GRID][GRID];
static int    hpeak = 1;
static double h2_sum[GRID][GRID], h4_sum[GRID][GRID];
static int    moment_n = 0;

static V2     poincare[MAX_POINCARE];
static int    n_poincare = 0, poincare_head = 0;

#define L_HISTORY 4000
static double L_hist[L_HISTORY];
static double L_max_val = 0.001;
static int    L_head = 0, L_count = 0;

enum { VIEW_HIST, VIEW_BORN, VIEW_KURTOSIS, VIEW_POINCARE, VIEW_COUNT };
static int view_mode = VIEW_HIST;
static Color wpix[GRID*GRID], rpix[GRID*GRID];

/* ----- Physics ----- */

/* Pending wave sources: bounces since last wave update */
#define MAX_PENDING 256
static V2  pending[MAX_PENDING];
static int n_pending = 0;

/* Flush pending bounces into the wave field (batched incremental update) */
static void wave_flush(void)
{
    if (n_pending == 0) return;
    for (int iy = 0; iy < GRID; iy++) {
        double y = 2.0*iy/(GRID-1)-1.0;
        for (int ix = 0; ix < GRID; ix++) {
            /* Decay for all pending bounces at once: mu^n_pending */
            wave[iy][ix] *= pow(P.mu, n_pending);
            if (geo_inside(geometry, 2.0*ix/(GRID-1)-1.0, y)) {
                /* Add each pending source with correct relative decay */
                double decay = 1.0;
                for (int p = n_pending-1; p >= 0; p--) {
                    double dx = 2.0*ix/(GRID-1)-1.0 - pending[p].x;
                    double dy = y - pending[p].y;
                    wave[iy][ix] += decay * j0(P.kF * sqrt(dx*dx+dy*dy));
                    decay *= P.mu;
                }
            }
            double h2 = wave[iy][ix]*wave[iy][ix];
            h2_sum[iy][ix] += h2;
            h4_sum[iy][ix] += h2*h2;
        }
    }
    moment_n++;
    n_pending = 0;
}

/* Record a bounce source for later flushing */
static void wave_record(double bx, double by)
{
    if (n_pending >= MAX_PENDING) wave_flush();  /* auto-flush if full */
    pending[n_pending++] = (V2){bx, by};
}

static void bounce(void)
{
    double gx = 0, gy = 0;
    int n = nmem < MAX_MEM ? nmem : MAX_MEM;
    double decay = 1.0;
    for (int k = 0; k < n; k++) {
        int i = (mhead-1-k+MAX_MEM)%MAX_MEM;
        double dx = drop.x-mem[i].x, dy = drop.y-mem[i].y;
        double d = sqrt(dx*dx+dy*dy);
        if (d < 1e-15) { decay *= P.mu; continue; }
        double f = -P.kF * j1(P.kF*d) / d;
        gx += decay*f*dx; gy += decay*f*dy;
        decay *= P.mu;
    }

    drop.x -= P.beta*gx;
    drop.y -= P.beta*gy;

    /* Wall repulsion */
    double d = geo_sdf(geometry, drop.x, drop.y);
    if (d < 0.15) {
        double nx, ny;
        geo_normal(geometry, drop.x, drop.y, &nx, &ny);
        double nm = sqrt(nx*nx+ny*ny);
        if (nm > 1e-10) { nx/=nm; ny/=nm; double s=(0.15-d)/0.15;
            drop.x -= 0.03*s*s*nx; drop.y -= 0.03*s*s*ny; }
    }
    d = geo_sdf(geometry, drop.x, drop.y);
    if (d < 0.02) {
        double nx, ny;
        geo_normal(geometry, drop.x, drop.y, &nx, &ny);
        double nm = sqrt(nx*nx+ny*ny);
        if (nm > 1e-10) { drop.x -= (0.02-d)*nx/nm; drop.y -= (0.02-d)*ny/nm; }
    }

    /* Poincare: y=0 upward crossing, record (x, dx) */
    if (prev_drop.y <= 0 && drop.y > 0) {
        double dy = drop.y-prev_drop.y;
        double frac = (dy>1e-15) ? -prev_drop.y/dy : 0.5;
        double xc = prev_drop.x + frac*(drop.x-prev_drop.x);
        double dxc = drop.x - prev_drop.x;
        if (n_poincare < MAX_POINCARE) poincare[n_poincare++] = (V2){xc,dxc};
        else { poincare[poincare_head] = (V2){xc,dxc}; poincare_head = (poincare_head+1)%MAX_POINCARE; }
    }
    /* Angular momentum: L = x*dy - y*dx (displacement as velocity proxy) */
    {
        double dx = drop.x - prev_drop.x, dy = drop.y - prev_drop.y;
        double L = drop.x*dy - drop.y*dx;
        L_hist[L_head] = L;
        L_head = (L_head+1) % L_HISTORY;
        if (L_count < L_HISTORY) L_count++;
        double aL = fabs(L);
        if (aL > L_max_val) L_max_val = aL*1.1;
        L_max_val *= 0.99999;
        if (L_max_val < 0.001) L_max_val = 0.001;
    }
    prev_drop = drop;

    mem[mhead] = drop;
    mhead = (mhead+1)%MAX_MEM;
    if (nmem < MAX_MEM) nmem++;
    total++;

    wave_record(drop.x, drop.y);
    hist_splat(hist, &hpeak, drop.x, drop.y, 1.0);
}

static void reset(void)
{
    nmem=0; mhead=0; total=0; hpeak=1; moment_n=0;
    n_poincare=0; poincare_head=0; n_pending=0;
    L_head=0; L_count=0; L_max_val=0.001;
    drop = (V2){0.25,0.0}; prev_drop = drop;
    memset(wave,0,sizeof(wave)); memset(hist,0,sizeof(hist));
    memset(h2_sum,0,sizeof(h2_sum)); memset(h4_sum,0,sizeof(h4_sum));
}

/* ----- Main ----- */

int main(void)
{
    InitWindow(WIN_W_L, WIN_H_L, "Walker - Level 0 (Direct Bessel)");
    SetTargetFPS(30);
    Image img = GenImageColor(GRID,GRID,BLACK);
    Texture2D wtex = LoadTextureFromImage(img);
    Texture2D rtex = LoadTextureFromImage(img);
    UnloadImage(img);

    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_SPACE)) paused=!paused;
        if (IsKeyPressed(KEY_R))     reset();
        if (IsKeyPressed(KEY_G))     { geometry=(geometry+1)%GEO_COUNT; reset(); }
        if (IsKeyPressed(KEY_V))     view_mode=(view_mode+1)%VIEW_COUNT;
        if (IsKeyPressed(KEY_T))     { turbo=!turbo; SetTargetFPS(turbo?10:30); }
        if (IsKeyPressed(KEY_UP))    P.mu = fmin(P.mu+0.005, 0.999);
        if (IsKeyPressed(KEY_DOWN))  P.mu = fmax(P.mu-0.005, 0.5);
        if (IsKeyPressed(KEY_RIGHT)) P.beta *= 1.15;
        if (IsKeyPressed(KEY_LEFT))  P.beta /= 1.15;
        if (IsKeyPressed(KEY_L))     { P.kF+=1.0; reset(); }
        if (IsKeyPressed(KEY_K))     { P.kF=fmax(P.kF-1.0,2.0); reset(); }
        if (IsKeyPressed(KEY_EQUAL)) speed = speed<50 ? speed+1 : speed;
        if (IsKeyPressed(KEY_MINUS)) speed = speed>1  ? speed-1 : speed;

        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) && !turbo) {
            Vector2 m = GetMousePosition();
            if (m.x>=GAP && m.x<GAP+PANEL && m.y>=GAP && m.y<GAP+PANEL) {
                double px=((m.x-GAP)/PANEL*2-1), py=((m.y-GAP)/PANEL*2-1);
                if (geo_inside(geometry,px,py)) { drop.x=px; drop.y=py; }
            }
        }

        if (!paused) {
            if (turbo) {
                struct timespec t0,t1; clock_gettime(CLOCK_MONOTONIC,&t0);
                int count=0;
                for (;;) { bounce(); count++;
                    if (count%100==0) { clock_gettime(CLOCK_MONOTONIC,&t1);
                        double el=(t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)*1e-9;
                        if (el>0.2) { turbo_bps=count/el; break; } } }
            } else { for (int i=0;i<speed;i++) bounce(); }
        }

        /* Flush wave sources and render */
        wave_flush();
        if (!turbo) { render_wave(wpix, wave, geometry, 1.0); UpdateTexture(wtex, wpix); }
        switch (view_mode) {
        case VIEW_HIST:     render_hist(rpix, hist, hpeak, geometry, 1.0); break;
        case VIEW_BORN:     render_born(rpix, h2_sum, moment_n, geometry, 1.0); break;
        case VIEW_KURTOSIS: render_kurtosis(rpix, h2_sum, h4_sum, moment_n, geometry, 1.0); break;
        case VIEW_POINCARE: render_poincare(rpix, poincare, n_poincare, poincare_head, MAX_POINCARE); break;
        }
        UpdateTexture(rtex, rpix);

        /* Draw */
        BeginDrawing();
        ClearBackground((Color){20,20,25,255});
        int lx=GAP, ly=GAP, rx=GAP*2+PANEL, ry=GAP;
        float scale=(float)PANEL/GRID;

        if (turbo) { draw_turbo(lx,ly,PANEL,turbo_bps,"bounces"); geo_draw_outline(geometry,lx,ly,PANEL); }
        else {
            DrawTextureEx(wtex,(Vector2){lx,ly},0,scale,WHITE);
            geo_draw_outline(geometry,lx,ly,PANEL);
            draw_trail(mem, mhead, nmem<MAX_MEM?nmem:MAX_MEM, TRAIL_VIS, lx,ly,PANEL,1.0,(Color){255,255,100,255});
            draw_droplet(drop.x,drop.y,lx,ly,PANEL,1.0);
        }
        DrawText(turbo?"TURBO [T]":"Wave Field",lx,ly+PANEL+4,16,LIGHTGRAY);

        DrawTextureEx(rtex,(Vector2){rx,ry},0,scale,WHITE);
        if (view_mode == VIEW_POINCARE) DrawRectangleLines(rx,ry,PANEL,PANEL,GRAY);
        else geo_draw_outline(geometry,rx,ry,PANEL);
        { const char *vn[]={"Histogram","Born <h^2>","Kurtosis","Poincare (x,dx)"};
          DrawText(vn[view_mode],rx,ry+PANEL+4,16,LIGHTGRAY); }

        /* Angular momentum strip */
        int Ly = GAP+PANEL+20, Lw = PANEL*2+GAP;
        draw_L_strip(L_hist, L_head, L_count, L_HISTORY, L_max_val, lx, Ly, Lw, LSTRIP);

        int iy = Ly+LSTRIP+4; char buf[256];
        snprintf(buf,sizeof(buf),"Bounces: %d | mu: %.3f [Up/Dn] | beta: %.4f [Lt/Rt] | speed: %d [+/-] | %s",
            total,P.mu,P.beta,speed,paused?"PAUSED":"[Space]");
        DrawText(buf,GAP,iy,14,LIGHTGRAY);
        snprintf(buf,sizeof(buf),"kF: %.1f [K/L] | %s [G] | [V] view | [T] turbo | [R] reset",
            P.kF,geo_name(geometry));
        DrawText(buf,GAP,iy+18,14,GRAY);

        EndDrawing();
    }
    UnloadTexture(wtex); UnloadTexture(rtex); CloseWindow();
    return 0;
}
