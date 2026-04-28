/*
 * Bouncing Droplet Simulator — Level 1 (Oza-Rosales-Bush)
 *
 * Stroboscopic trajectory equation for a walking droplet:
 *
 *   m * x_ddot + D * x_dot = -F * grad h(x, t)
 *   h(x, t) = sum_{n<N} J_0(kF * |x - x_n|) * mu^{N-n}
 *
 * The droplet has inertia (velocity) and drag. The wave field is the
 * same Bessel kernel as Level 0 but the dynamics are second-order.
 * This is the discrete stroboscopic form of the Oza et al. (2013)
 * integro-differential trajectory equation.
 *
 * Reference: Oza, Rosales & Bush, JFM 737, 552-570 (2013)
 *
 * Controls:
 *   Space       pause / resume
 *   Up/Down     memory parameter mu
 *   Left/Right  kick strength beta
 *   D/F         drag (D=more inertia, F=more friction)
 *   +/-         simulation speed
 *   K/L         Faraday wavenumber kF (resets)
 *   G           cycle geometry
 *   V           cycle right panel view
 *   T           toggle turbo mode
 *   R           reset simulation
 *   Click       place droplet (left panel)
 */

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "raylib.h"

#define MAX_MEM   512
#define GRID      200
#define TRAIL_VIS 200

typedef struct { double x, y; } V2;

/* ----- Geometry (same as Level 0) ----- */

enum { GEO_STADIUM, GEO_DSHAPE, GEO_CIRCLE, GEO_RECT, GEO_COUNT };
static int geometry = GEO_STADIUM;

#define DSHAPE_R    0.9
#define DSHAPE_CUT (-0.3)

static double geo_sdf(double x, double y)
{
    switch (geometry) {
    case GEO_CIRCLE:  return 1.0 - sqrt(x*x + y*y);
    case GEO_STADIUM: {
        double a=0.5, r=0.5, ax=fabs(x);
        return ax<=a ? r-fabs(y) : r-sqrt((ax-a)*(ax-a)+y*y);
    }
    case GEO_RECT:    return fmin(0.8-fabs(x), 0.5-fabs(y));
    case GEO_DSHAPE:  return fmin(DSHAPE_R-sqrt(x*x+y*y), x-DSHAPE_CUT);
    default: return -1;
    }
}

static int geo_inside(double x, double y) { return geo_sdf(x,y) > 0; }

static void geo_normal(double x, double y, double *nx, double *ny)
{
    double eps = 0.002;
    *nx = -(geo_sdf(x+eps,y)-geo_sdf(x-eps,y))/(2*eps);
    *ny = -(geo_sdf(x,y+eps)-geo_sdf(x,y-eps))/(2*eps);
}

#define S2P(px, py) (Vector2){ \
    ox+((float)(px)*0.5f+0.5f)*sz, oy+((float)(py)*0.5f+0.5f)*sz }

static void geo_draw_outline(int ox, int oy, int sz)
{
    int N=200; Color col=LIGHTGRAY;
    switch (geometry) {
    case GEO_CIRCLE:
        DrawCircleLinesV((Vector2){ox+sz/2.0f,oy+sz/2.0f},sz/2.0f,col); break;
    case GEO_STADIUM: {
        float a=0.5f,r=0.5f;
        DrawLineV(S2P(-a,r),S2P(a,r),col); DrawLineV(S2P(-a,-r),S2P(a,-r),col);
        for (int i=0;i<N;i++) {
            float t0=PI/2+PI*i/(float)N,t1=PI/2+PI*(i+1)/(float)N;
            DrawLineV(S2P(-a+r*cosf(t0),r*sinf(t0)),S2P(-a+r*cosf(t1),r*sinf(t1)),col); }
        for (int i=0;i<N;i++) {
            float t0=-PI/2+PI*i/(float)N,t1=-PI/2+PI*(i+1)/(float)N;
            DrawLineV(S2P(a+r*cosf(t0),r*sinf(t0)),S2P(a+r*cosf(t1),r*sinf(t1)),col); }
        break;
    }
    case GEO_RECT: {
        float hw=0.8f,hh=0.5f;
        DrawLineV(S2P(-hw,-hh),S2P(hw,-hh),col); DrawLineV(S2P(hw,-hh),S2P(hw,hh),col);
        DrawLineV(S2P(hw,hh),S2P(-hw,hh),col); DrawLineV(S2P(-hw,hh),S2P(-hw,-hh),col);
        break;
    }
    case GEO_DSHAPE: {
        float R=DSHAPE_R, yint=sqrtf(R*R-DSHAPE_CUT*DSHAPE_CUT);
        float a0=atan2f(yint,DSHAPE_CUT),a1=atan2f(-yint,DSHAPE_CUT),span=a0-a1;
        DrawLineV(S2P(DSHAPE_CUT,-yint),S2P(DSHAPE_CUT,yint),col);
        for (int i=0;i<N;i++) {
            float t0=a0-span*i/(float)N,t1=a0-span*(i+1)/(float)N;
            DrawLineV(S2P(R*cosf(t0),R*sinf(t0)),S2P(R*cosf(t1),R*sinf(t1)),col); }
        break;
    }
    }
}
#undef S2P

static const char *geo_name(void) {
    const char *n[]={"Stadium","D-shape","Circle","Rectangle"};
    return geometry<GEO_COUNT ? n[geometry] : "?";
}

/* ----- Parameters & State ----- */

static struct { double kF, mu, beta; } P = { 12.0, 0.97, 0.005 };
static double drag = 0.95;  /* velocity retention per bounce (1 = no drag) */

static V2     mem[MAX_MEM];
static int    nmem = 0, mhead = 0;
static V2     drop = { 0.25, 0.0 };
static V2     vel  = { 0.0, 0.0 };
static int    total = 0, paused = 0, speed = 3, turbo = 0;
static double turbo_bps = 0;

static double wave[GRID][GRID];
static int    hist[GRID][GRID];
static int    hpeak = 1;
static double h2_sum[GRID][GRID];
static double h4_sum[GRID][GRID];
static int    moment_n = 0;

enum { VIEW_HIST, VIEW_BORN, VIEW_KURTOSIS, VIEW_COUNT };
static int view_mode = VIEW_HIST;

static Color wpix[GRID * GRID];
static Color rpix[GRID * GRID];

/* ----- Colormaps ----- */

static Color cmap_bwr(double t) {
    if (t<0) t=0; if (t>1) t=1;
    if (t<0.5) { unsigned char v=(unsigned char)(t*2*255); return (Color){v,v,255,255}; }
    else { unsigned char v=(unsigned char)((1-t)*2*255); return (Color){255,v,v,255}; }
}

static Color cmap_hot(double t) {
    if (t<0) t=0; if (t>1) t=1;
    unsigned char r,g,b;
    if (t<0.33) { r=(unsigned char)(t/0.33*255); g=0; b=0; }
    else if (t<0.66) { r=255; g=(unsigned char)((t-0.33)/0.33*255); b=0; }
    else { r=255; g=255; b=(unsigned char)((t-0.66)/0.34*255); }
    return (Color){r,g,b,255};
}

/* ----- Physics ----- */

static void wave_update(double bx, double by)
{
    for (int iy = 0; iy < GRID; iy++) {
        double y = 2.0*iy/(GRID-1)-1.0;
        for (int ix = 0; ix < GRID; ix++) {
            double x = 2.0*ix/(GRID-1)-1.0;
            wave[iy][ix] *= P.mu;
            if (geo_inside(x,y)) {
                double dx=x-bx, dy=y-by;
                wave[iy][ix] += j0(P.kF * sqrt(dx*dx+dy*dy));
            }
            double h2 = wave[iy][ix]*wave[iy][ix];
            h2_sum[iy][ix] += h2;
            h4_sum[iy][ix] += h2*h2;
        }
    }
    moment_n++;
}

static void bounce(void)
{
    /* Save bounce position (wave records where droplet WAS) */
    double bx = drop.x, by = drop.y;

    /* Wave gradient at current position (analytic from memory) */
    double gx = 0, gy = 0;
    int n = nmem < MAX_MEM ? nmem : MAX_MEM;
    double decay = 1.0;
    for (int k = 0; k < n; k++) {
        int i = (mhead-1-k+MAX_MEM) % MAX_MEM;
        double dx = drop.x - mem[i].x;
        double dy = drop.y - mem[i].y;
        double d = sqrt(dx*dx + dy*dy);
        if (d < 1e-15) { decay *= P.mu; continue; }
        double f = -P.kF * j1(P.kF*d) / d;
        gx += decay * f * dx;
        gy += decay * f * dy;
        decay *= P.mu;
    }

    /* Second-order dynamics: velocity with drag */
    vel.x = drag * vel.x - P.beta * gx;
    vel.y = drag * vel.y - P.beta * gy;
    drop.x += vel.x;
    drop.y += vel.y;

    /* Wall: repulsion + velocity damping */
    double d = geo_sdf(drop.x, drop.y);
    if (d < 0.15) {
        double nx, ny;
        geo_normal(drop.x, drop.y, &nx, &ny);
        double nm = sqrt(nx*nx+ny*ny);
        if (nm > 1e-10) {
            nx /= nm; ny /= nm;
            double s = (0.15-d)/0.15;
            double push = 0.03*s*s;
            drop.x -= push*nx; drop.y -= push*ny;
            /* Reflect velocity off wall */
            double vn = vel.x*nx + vel.y*ny;
            if (vn > 0) { vel.x -= 1.5*vn*nx; vel.y -= 1.5*vn*ny; }
        }
    }
    d = geo_sdf(drop.x, drop.y);
    if (d < 0.02) {
        double nx, ny;
        geo_normal(drop.x, drop.y, &nx, &ny);
        double nm = sqrt(nx*nx+ny*ny);
        if (nm > 1e-10) {
            drop.x -= (0.02-d)*nx/nm;
            drop.y -= (0.02-d)*ny/nm;
        }
    }

    /* Record at OLD position (where bounce actually occurred) */
    mem[mhead] = (V2){bx, by};
    mhead = (mhead+1) % MAX_MEM;
    if (nmem < MAX_MEM) nmem++;
    total++;

    /* Wave update at old position */
    wave_update(bx, by);

    /* Histogram */
    double hxf = (drop.x*0.5+0.5)*GRID, hyf = (drop.y*0.5+0.5)*GRID;
    int cx = (int)hxf, cy = (int)hyf;
    for (int dy = -2; dy <= 2; dy++)
        for (int dx = -2; dx <= 2; dx++) {
            int px=cx+dx, py=cy+dy;
            if (px<0||px>=GRID||py<0||py>=GRID) continue;
            double gg = (px-hxf)*(px-hxf)+(py-hyf)*(py-hyf);
            hist[py][px] += (int)(100*exp(-gg));
            if (hist[py][px] > hpeak) hpeak = hist[py][px];
        }
}

/* ----- Rendering (identical to Level 0) ----- */

static void render_wave_tex(Texture2D tex)
{
    double vmin=0, vmax=0;
    for (int iy=0;iy<GRID;iy++) for (int ix=0;ix<GRID;ix++) {
        if (wave[iy][ix]<vmin) vmin=wave[iy][ix];
        if (wave[iy][ix]>vmax) vmax=wave[iy][ix];
    }
    double range=vmax-vmin; if (range<1e-10) range=1;
    for (int iy=0;iy<GRID;iy++) {
        double y=2.0*iy/(GRID-1)-1.0;
        for (int ix=0;ix<GRID;ix++) {
            double x=2.0*ix/(GRID-1)-1.0;
            int idx=iy*GRID+ix;
            wpix[idx] = !geo_inside(x,y) ? (Color){30,30,30,255}
                       : cmap_bwr((wave[iy][ix]-vmin)/range);
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
        for (int iy=0;iy<GRID;iy++) { double y=2.0*iy/(GRID-1)-1.0;
            for (int ix=0;ix<GRID;ix++) { double x=2.0*ix/(GRID-1)-1.0;
                int idx=iy*GRID+ix;
                if (!geo_inside(x,y)) rpix[idx] = (Color){30,30,30,255};
                else if (hist[iy][ix] == 0) rpix[idx] = BLACK;
                else rpix[idx] = cmap_hot((log(1.0+(double)hist[iy][ix]) - lmin) / lrange);
            }
        }
    } else if (view_mode == VIEW_BORN) {
        if (moment_n<1) { UpdateTexture(tex,rpix); return; }
        double peak=0;
        for (int iy=0;iy<GRID;iy++) for (int ix=0;ix<GRID;ix++) {
            double v=h2_sum[iy][ix]/moment_n; if (v>peak) peak=v; }
        if (peak<1e-20) peak=1;
        for (int iy=0;iy<GRID;iy++) { double y=2.0*iy/(GRID-1)-1.0;
            for (int ix=0;ix<GRID;ix++) { double x=2.0*ix/(GRID-1)-1.0;
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
            double ek=m4-3.0*m2*m2; if (ek<vmin) vmin=ek; if (ek>vmax) vmax=ek; }
        double amax=fmax(fabs(vmin),fabs(vmax)); if (amax<1e-20) amax=1;
        for (int iy=0;iy<GRID;iy++) { double y=2.0*iy/(GRID-1)-1.0;
            for (int ix=0;ix<GRID;ix++) { double x=2.0*ix/(GRID-1)-1.0;
                int idx=iy*GRID+ix;
                if (!geo_inside(x,y)) { rpix[idx]=(Color){30,30,30,255}; continue; }
                double m2=h2_sum[iy][ix]/moment_n, m4=h4_sum[iy][ix]/moment_n;
                rpix[idx]=cmap_bwr(0.5+0.5*(m4-3.0*m2*m2)/amax);
            }
        }
    }
    UpdateTexture(tex, rpix);
}

static const char *view_name(void) {
    const char *n[]={"Histogram","Born <h^2>","Excess Kurtosis"};
    return view_mode<VIEW_COUNT ? n[view_mode] : "?";
}

static Vector2 p2s(double px, double py, int ox, int oy, int sz) {
    return (Vector2){ ox+(float)((px*0.5+0.5)*sz), oy+(float)((py*0.5+0.5)*sz) };
}

static void reset_sim(void)
{
    nmem=0; mhead=0; total=0; hpeak=1; moment_n=0;
    drop=(V2){0.25,0.0}; vel=(V2){0.0,0.0};
    memset(wave,0,sizeof(wave)); memset(hist,0,sizeof(hist));
    memset(h2_sum,0,sizeof(h2_sum)); memset(h4_sum,0,sizeof(h4_sum));
}

/* ----- Main ----- */

#define PANEL 400
#define GAP    20
#define WIN_W (PANEL*2+GAP*3)
#define WIN_H (PANEL+GAP*2+60)

int main(void)
{
    InitWindow(WIN_W, WIN_H, "Bouncing Droplet - Level 1 (Oza-Rosales-Bush)");
    SetTargetFPS(30);

    Image img = GenImageColor(GRID, GRID, BLACK);
    Texture2D wtex = LoadTextureFromImage(img);
    Texture2D rtex = LoadTextureFromImage(img);
    UnloadImage(img);

    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_SPACE)) paused=!paused;
        if (IsKeyPressed(KEY_R))     reset_sim();
        if (IsKeyPressed(KEY_G))     { geometry=(geometry+1)%GEO_COUNT; reset_sim(); }
        if (IsKeyPressed(KEY_V))     view_mode=(view_mode+1)%VIEW_COUNT;
        if (IsKeyPressed(KEY_T))     { turbo=!turbo; SetTargetFPS(turbo?10:30); }
        if (IsKeyPressed(KEY_UP))    P.mu=fmin(P.mu+0.005,0.999);
        if (IsKeyPressed(KEY_DOWN))  P.mu=fmax(P.mu-0.005,0.5);
        if (IsKeyPressed(KEY_RIGHT)) P.beta*=1.15;
        if (IsKeyPressed(KEY_LEFT))  P.beta/=1.15;
        if (IsKeyPressed(KEY_D))     drag=fmin(drag+0.01,0.999);
        if (IsKeyPressed(KEY_F))     drag=fmax(drag-0.01,0.0);
        if (IsKeyPressed(KEY_L))     { P.kF+=1.0; reset_sim(); }
        if (IsKeyPressed(KEY_K))     { P.kF=fmax(P.kF-1.0,2.0); reset_sim(); }
        if (IsKeyPressed(KEY_EQUAL)) speed=speed<50?speed+1:speed;
        if (IsKeyPressed(KEY_MINUS)) speed=speed>1?speed-1:speed;

        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT) && !turbo) {
            Vector2 m=GetMousePosition();
            if (m.x>=GAP && m.x<GAP+PANEL && m.y>=GAP && m.y<GAP+PANEL) {
                double px=((m.x-GAP)/PANEL*2-1), py=((m.y-GAP)/PANEL*2-1);
                if (geo_inside(px,py)) { drop.x=px; drop.y=py; vel.x=0; vel.y=0; }
            }
        }

        if (!paused) {
            if (turbo) {
                struct timespec t0,t1;
                clock_gettime(CLOCK_MONOTONIC,&t0);
                int count=0;
                for (;;) { bounce(); count++;
                    if (count%100==0) {
                        clock_gettime(CLOCK_MONOTONIC,&t1);
                        double el=(t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)*1e-9;
                        if (el>0.2) { turbo_bps=count/el; break; }
                    }
                }
            } else {
                for (int i=0; i<speed; i++) bounce();
            }
        }

        if (!turbo) render_wave_tex(wtex);
        render_right_tex(rtex);

        BeginDrawing();
        ClearBackground((Color){20,20,25,255});
        int lx=GAP,ly=GAP,rx=GAP*2+PANEL,ry=GAP;
        float scale=(float)PANEL/GRID;

        if (turbo) {
            DrawRectangle(lx,ly,PANEL,PANEL,(Color){30,30,30,255});
            geo_draw_outline(lx,ly,PANEL);
            DrawText("TURBO MODE",lx+130,ly+160,24,(Color){255,100,100,255});
            char tb[128]; snprintf(tb,sizeof(tb),"%.0f bounces/sec",turbo_bps);
            DrawText(tb,lx+140,ly+195,18,LIGHTGRAY);
            DrawText("[T] to exit turbo",lx+135,ly+225,14,GRAY);
        } else {
            DrawTextureEx(wtex,(Vector2){lx,ly},0,scale,WHITE);
            geo_draw_outline(lx,ly,PANEL);
            int tlen=nmem<TRAIL_VIS?nmem:TRAIL_VIS;
            for (int k=1;k<tlen;k++) {
                int i0=(mhead-k+MAX_MEM)%MAX_MEM, i1=(mhead-k-1+MAX_MEM)%MAX_MEM;
                Vector2 a=p2s(mem[i0].x,mem[i0].y,lx,ly,PANEL);
                Vector2 b=p2s(mem[i1].x,mem[i1].y,lx,ly,PANEL);
                unsigned char al=(unsigned char)(180*(1.0-(double)k/tlen));
                DrawLineEx(a,b,1.5f,(Color){255,255,100,al});
            }
            Vector2 dp=p2s(drop.x,drop.y,lx,ly,PANEL);
            DrawCircleV(dp,5,(Color){255,255,200,120});
            DrawCircleV(dp,3,WHITE);
        }
        DrawText(turbo?"TURBO [T]":"Wave Field",lx,ly+PANEL+4,16,LIGHTGRAY);

        DrawTextureEx(rtex,(Vector2){rx,ry},0,scale,WHITE);
        geo_draw_outline(rx,ry,PANEL);
        { char l[64]; snprintf(l,sizeof(l),"%s [V]",view_name()); DrawText(l,rx,ry+PANEL+4,16,LIGHTGRAY); }

        int iy=GAP+PANEL+24; char buf[256];
        snprintf(buf,sizeof(buf),
            "Bounces: %d  |  mu: %.3f [Up/Dn]  |  beta: %.4f [Lt/Rt]  |  "
            "drag: %.3f [D/F]  |  speed: %d [+/-]  |  %s",
            total,P.mu,P.beta,drag,speed,paused?"PAUSED":"[Space]");
        DrawText(buf,GAP,iy,14,LIGHTGRAY);
        snprintf(buf,sizeof(buf),
            "kF: %.1f [K/L]  |  %s [G]  |  vel: %.4f  |  Samples: %d  |  [T] turbo",
            P.kF,geo_name(),sqrt(vel.x*vel.x+vel.y*vel.y),moment_n);
        DrawText(buf,GAP,iy+18,14,GRAY);

        EndDrawing();
    }

    UnloadTexture(wtex); UnloadTexture(rtex);
    CloseWindow();
    return 0;
}
