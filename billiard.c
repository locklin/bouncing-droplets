/*
 * Classical Billiard — Point particle bouncing elastically off walls.
 *
 * Displays trajectory, position histogram, and Poincare section
 * (boundary map in Birkhoff coordinates: arclength s vs sin(alpha)
 * where alpha is the reflection angle).
 *
 * Controls:
 *   Space       pause / resume
 *   G           cycle geometry
 *   V           cycle right panel (histogram / Poincare)
 *   T           turbo mode
 *   +/-         speed
 *   R           reset (random initial conditions)
 *   Click       set position (velocity tangential)
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "raylib.h"

#define GRID      200
#define MAX_TRAIL 8000
#define TRAIL_VIS 2000
#define DT        0.002
typedef struct { double x, y; } V2;

/* ----- Geometry ----- */

enum { GEO_STADIUM, GEO_DSHAPE, GEO_CIRCLE, GEO_RECT, GEO_COUNT };
static int geometry = GEO_STADIUM;

#define DSHAPE_R   0.9
#define DSHAPE_CUT (-0.3)

static double geo_sdf(double x, double y)
{
    switch (geometry) {
    case GEO_CIRCLE:  return 1.0 - sqrt(x*x+y*y);
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
    double eps = 0.001;
    *nx = -(geo_sdf(x+eps,y)-geo_sdf(x-eps,y))/(2*eps);
    *ny = -(geo_sdf(x,y+eps)-geo_sdf(x,y-eps))/(2*eps);
}

/* Compute normalized arclength s in [0,1] for a point on the boundary.
 * Uses the angle/parameter that naturally parameterizes each geometry. */
static double geo_arclength(double x, double y)
{
    switch (geometry) {
    case GEO_CIRCLE: {
        double a = atan2(y, x);
        return (a + M_PI) / (2*M_PI);  /* [0,1] */
    }
    case GEO_STADIUM: {
        /* Boundary: bottom flat, right semicircle, top flat, left semicircle.
         * Total perimeter = 2*1.0 + 2*pi*0.5 = 2 + pi ≈ 5.14 */
        double a = 0.5, r = 0.5;
        double perim = 2*a*2 + 2*M_PI*r;  /* 2 + pi */
        double s;
        if (fabs(y + r) < 0.02 && fabs(x) <= a) {
            /* Bottom flat: x from -a to a */
            s = (x + a) / (2*a) * (2*a);
        } else if (x > a - 0.02) {
            /* Right semicircle */
            double ang = atan2(y, x - a);
            s = 2*a + (ang + M_PI/2) * r;
        } else if (fabs(y - r) < 0.02 && fabs(x) <= a) {
            /* Top flat: x from a to -a */
            s = 2*a + M_PI*r + (a - x);
        } else {
            /* Left semicircle */
            double ang = atan2(y, x + a);
            s = 2*a + M_PI*r + 2*a + (ang - M_PI/2 + M_PI) * r;
        }
        return fmod(s / perim + 1.0, 1.0);
    }
    case GEO_RECT: {
        double hw = 0.8, hh = 0.5;
        double perim = 2*(2*hw + 2*hh);
        double s;
        if (fabs(y + hh) < 0.02) s = (x + hw);               /* bottom */
        else if (fabs(x - hw) < 0.02) s = 2*hw + (y + hh);   /* right */
        else if (fabs(y - hh) < 0.02) s = 2*hw + 2*hh + (hw - x); /* top */
        else s = 2*hw + 2*hh + 2*hw + (hh - y);               /* left */
        return fmod(s / perim + 1.0, 1.0);
    }
    case GEO_DSHAPE: {
        /* Chord + arc. Approximate with angle-based parameterization. */
        double R = DSHAPE_R;
        double yint = sqrt(R*R - DSHAPE_CUT*DSHAPE_CUT);
        double ang_top = atan2(yint, DSHAPE_CUT);
        double ang_bot = atan2(-yint, DSHAPE_CUT);
        double arc_len = R * (ang_top - ang_bot);
        double chord_len = 2 * yint;
        double perim = arc_len + chord_len;

        if (fabs(x - DSHAPE_CUT) < 0.02 && fabs(y) < yint) {
            /* On the chord */
            double s = (y + yint) / chord_len * chord_len;
            return s / perim;
        } else {
            /* On the arc */
            double ang = atan2(y, x);
            double s = chord_len + R * (ang_top - ang);
            return fmod(s / perim + 1.0, 1.0);
        }
    }
    default: return 0;
    }
}

#define S2P(px,py) (Vector2){ ox+((float)(px)*0.5f+0.5f)*sz, \
                               oy+((float)(py)*0.5f+0.5f)*sz }
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
            float t0=PI/2+PI*i/(float)N, t1=PI/2+PI*(i+1)/(float)N;
            DrawLineV(S2P(-a+r*cosf(t0),r*sinf(t0)),S2P(-a+r*cosf(t1),r*sinf(t1)),col);}
        for (int i=0;i<N;i++) {
            float t0=-PI/2+PI*i/(float)N, t1=-PI/2+PI*(i+1)/(float)N;
            DrawLineV(S2P(a+r*cosf(t0),r*sinf(t0)),S2P(a+r*cosf(t1),r*sinf(t1)),col);}
        break; }
    case GEO_RECT: {
        float hw=0.8f,hh=0.5f;
        DrawLineV(S2P(-hw,-hh),S2P(hw,-hh),col); DrawLineV(S2P(hw,-hh),S2P(hw,hh),col);
        DrawLineV(S2P(hw,hh),S2P(-hw,hh),col); DrawLineV(S2P(-hw,hh),S2P(-hw,-hh),col);
        break; }
    case GEO_DSHAPE: {
        float R=DSHAPE_R, yint=sqrtf(R*R-DSHAPE_CUT*DSHAPE_CUT);
        float a0=atan2f(yint,DSHAPE_CUT),a1=atan2f(-yint,DSHAPE_CUT),span=a0-a1;
        DrawLineV(S2P(DSHAPE_CUT,-yint),S2P(DSHAPE_CUT,yint),col);
        for (int i=0;i<N;i++) {
            float t0=a0-span*i/(float)N, t1=a0-span*(i+1)/(float)N;
            DrawLineV(S2P(R*cosf(t0),R*sinf(t0)),S2P(R*cosf(t1),R*sinf(t1)),col);}
        break; }
    }
}
#undef S2P

static const char *geo_name(void) {
    const char *n[]={"Stadium","D-shape","Circle","Rectangle"};
    return geometry<GEO_COUNT ? n[geometry] : "?";
}

/* ----- State ----- */

static V2 pos, vel;
static V2 trail[MAX_TRAIL];
static int n_trail = 0, trail_head = 0;
static int total_bounces = 0, total_steps = 0;
static int paused = 0, speed = 200, turbo = 0;
static double turbo_sps = 0;

static int    hist[GRID][GRID];
static int    hpeak = 1;

/* Poincare section: density grid in (s, sin_alpha) */
static int    poincare_hist[GRID][GRID];
static int    ppeak = 1;
static int    n_poincare = 0;

enum { VIEW_HIST, VIEW_POINCARE, VIEW_COUNT };
static int view_mode = VIEW_HIST;

static Color rpix[GRID * GRID];

static Color cmap_hot(double t) {
    if (t < 0) t = 0;
    if (t > 1) t = 1;
    unsigned char r,g,b;
    if (t<0.33) { r=(unsigned char)(t/0.33*255); g=0; b=0; }
    else if (t<0.66) { r=255; g=(unsigned char)((t-0.33)/0.33*255); b=0; }
    else { r=255; g=255; b=(unsigned char)((t-0.66)/0.34*255); }
    return (Color){r,g,b,255};
}

/* ----- Billiard Step ----- */

static void billiard_step(void)
{
    double dt_remain = DT;
    for (int bounce = 0; bounce < 10 && dt_remain > 1e-10; bounce++) {
        double nx = pos.x + vel.x * dt_remain;
        double ny = pos.y + vel.y * dt_remain;

        if (geo_sdf(nx, ny) > 0) {
            pos.x = nx; pos.y = ny;
            break;
        }

        /* Binary search for crossing time */
        double t_lo = 0, t_hi = dt_remain;
        for (int i = 0; i < 40; i++) {
            double t_mid = (t_lo + t_hi) / 2;
            double mx = pos.x + vel.x * t_mid;
            double my = pos.y + vel.y * t_mid;
            if (geo_sdf(mx, my) > 0) t_lo = t_mid;
            else t_hi = t_mid;
        }

        pos.x += vel.x * t_lo;
        pos.y += vel.y * t_lo;
        dt_remain -= t_lo;

        /* Wall normal */
        double wnx, wny;
        geo_normal(pos.x, pos.y, &wnx, &wny);
        double nm = sqrt(wnx*wnx + wny*wny);
        if (nm > 1e-10) {
            wnx /= nm; wny /= nm;

            /* Record Poincare section BEFORE reflecting */
            double vn = vel.x*wnx + vel.y*wny; /* component into wall */
            double vt_x = vel.x - vn*wnx;      /* tangential component */
            double vt_y = vel.y - vn*wny;
            double vt = sqrt(vt_x*vt_x + vt_y*vt_y);
            double v_total = sqrt(vel.x*vel.x + vel.y*vel.y);

            /* sin(alpha) = tangential / total (alpha = angle from normal) */
            double sin_alpha = (v_total > 1e-12) ? vt / v_total : 0;
            /* Sign: positive if tangential is "counterclockwise" */
            double cross = wnx * vt_y - wny * vt_x;
            if (cross < 0) sin_alpha = -sin_alpha;

            double s = geo_arclength(pos.x, pos.y);

            {
                int bx = (int)(s * (GRID-1));
                int by = (int)((1.0 - (sin_alpha*0.5+0.5)) * (GRID-1));
                if (bx >= 0 && bx < GRID && by >= 0 && by < GRID) {
                    poincare_hist[by][bx]++;
                    if (poincare_hist[by][bx] > ppeak) ppeak = poincare_hist[by][bx];
                }
                n_poincare++;
            }

            /* Reflect */
            if (vn > 0) {
                vel.x -= 2*vn*wnx;
                vel.y -= 2*vn*wny;
            }
        }

        pos.x -= 0.001*wnx;
        pos.y -= 0.001*wny;
        total_bounces++;
    }

    total_steps++;

    /* Trail */
    trail[trail_head] = pos;
    trail_head = (trail_head+1) % MAX_TRAIL;
    if (n_trail < MAX_TRAIL) n_trail++;

    /* Histogram */
    static int hcnt = 0;
    if (++hcnt >= 5) {
        hcnt = 0;
        double hxf=(pos.x*0.5+0.5)*GRID, hyf=(pos.y*0.5+0.5)*GRID;
        int cx=(int)hxf, cy=(int)hyf;
        for (int dy=-2; dy<=2; dy++)
            for (int dx=-2; dx<=2; dx++) {
                int px=cx+dx, py=cy+dy;
                if (px<0||px>=GRID||py<0||py>=GRID) continue;
                double gg=(px-hxf)*(px-hxf)+(py-hyf)*(py-hyf);
                hist[py][px] += (int)(100*exp(-gg));
                if (hist[py][px]>hpeak) hpeak=hist[py][px];
            }
    }
}

/* ----- Rendering ----- */

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
        for (int iy=0; iy<GRID; iy++) {
            double y=2.0*iy/(GRID-1)-1.0;
            for (int ix=0; ix<GRID; ix++) {
                double x=2.0*ix/(GRID-1)-1.0;
                int idx=iy*GRID+ix;
                if (!geo_inside(x,y)) rpix[idx] = (Color){30,30,30,255};
                else if (hist[iy][ix] == 0) rpix[idx] = BLACK;
                else rpix[idx] = cmap_hot((log(1.0+(double)hist[iy][ix]) - lmin) / lrange);
            }
        }
    } else if (view_mode == VIEW_POINCARE) {
        /* Poincare section: s on x-axis [0,1], sin(alpha) on y-axis [-1,1] */
        double lmin = 1e30, lmax = 0;
        for (int i = 0; i < GRID*GRID; i++)
            if (poincare_hist[i/GRID][i%GRID] > 0) {
                double lv = log(1.0 + poincare_hist[i/GRID][i%GRID]);
                if (lv < lmin) lmin = lv;
                if (lv > lmax) lmax = lv;
            }
        double lrange = lmax - lmin; if (lrange < 1e-10) lrange = 1;
        for (int iy=0; iy<GRID; iy++)
            for (int ix=0; ix<GRID; ix++) {
                int idx = iy*GRID+ix;
                rpix[idx] = poincare_hist[iy][ix] > 0
                    ? cmap_hot((log(1.0 + poincare_hist[iy][ix]) - lmin) / lrange)
                    : (Color){15, 15, 20, 255};
            }
        /* Axis line at sin(alpha)=0 */
        int mid_y = GRID/2;
        for (int ix=0; ix<GRID; ix++)
            if (poincare_hist[mid_y][ix] == 0)
                rpix[mid_y*GRID+ix] = (Color){40, 40, 50, 255};
    }
    UpdateTexture(tex, rpix);
}

static Vector2 p2s(double px, double py, int ox, int oy, int sz) {
    return (Vector2){ ox+(float)((px*0.5+0.5)*sz), oy+(float)((py*0.5+0.5)*sz) };
}

static void reset_sim(void)
{
    srand(time(NULL));
    do {
        pos.x = ((double)rand()/RAND_MAX-0.5)*1.4;
        pos.y = ((double)rand()/RAND_MAX-0.5)*1.4;
    } while (!geo_inside(pos.x, pos.y));

    double angle = (double)rand()/RAND_MAX * 2*M_PI;
    vel.x = 0.5*cos(angle);
    vel.y = 0.5*sin(angle);

    n_trail=0; trail_head=0; total_bounces=0; total_steps=0;
    hpeak=1; n_poincare=0; ppeak=1;
    memset(hist, 0, sizeof(hist));
    memset(poincare_hist, 0, sizeof(poincare_hist));
}

/* ----- Main ----- */

#define PANEL 400
#define GAP    20
#define WIN_W (PANEL*2+GAP*3)
#define WIN_H (PANEL+GAP*2+60)

int main(void)
{
    InitWindow(WIN_W, WIN_H, "Classical Billiard");
    SetTargetFPS(30);

    Image img = GenImageColor(GRID, GRID, BLACK);
    Texture2D rtex = LoadTextureFromImage(img);
    UnloadImage(img);

    reset_sim();

    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_SPACE)) paused=!paused;
        if (IsKeyPressed(KEY_R))     reset_sim();
        if (IsKeyPressed(KEY_G))     { geometry=(geometry+1)%GEO_COUNT; reset_sim(); }
        if (IsKeyPressed(KEY_V))     view_mode=(view_mode+1)%VIEW_COUNT;
        if (IsKeyPressed(KEY_T))     { turbo=!turbo; SetTargetFPS(turbo?10:30); }
        if (IsKeyPressed(KEY_EQUAL)) speed = speed<5000 ? speed+200 : speed;
        if (IsKeyPressed(KEY_MINUS)) speed = speed>200  ? speed-200 : speed;

        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
            Vector2 m=GetMousePosition();
            if (m.x>=GAP && m.x<GAP+PANEL && m.y>=GAP && m.y<GAP+PANEL) {
                double px=((m.x-GAP)/PANEL*2-1), py=((m.y-GAP)/PANEL*2-1);
                if (geo_inside(px,py)) {
                    pos.x=px; pos.y=py;
                    double r=sqrt(px*px+py*py);
                    if (r>0.01) { vel.x=-0.5*py/r; vel.y=0.5*px/r; }
                    n_trail=0; trail_head=0;
                }
            }
        }

        if (!paused) {
            if (turbo) {
                struct timespec t0,t1;
                clock_gettime(CLOCK_MONOTONIC,&t0);
                int count=0;
                for (;;) { billiard_step(); count++;
                    if (count%2000==0) {
                        clock_gettime(CLOCK_MONOTONIC,&t1);
                        double el=(t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)*1e-9;
                        if (el>0.2) { turbo_sps=count/el; break; }
                    }
                }
            } else {
                for (int i=0; i<speed; i++) billiard_step();
            }
        }

        render_right_tex(rtex);

        BeginDrawing();
        ClearBackground((Color){20,20,25,255});
        int lx=GAP, ly=GAP, rx=GAP*2+PANEL, ry=GAP;
        float scale=(float)PANEL/GRID;

        /* Left: trajectory */
        if (turbo) {
            DrawRectangle(lx,ly,PANEL,PANEL,(Color){30,30,30,255});
            geo_draw_outline(lx,ly,PANEL);
            DrawText("TURBO",lx+160,ly+180,24,(Color){255,100,100,255});
            char tb[64]; snprintf(tb,sizeof(tb),"%.0f steps/sec",turbo_sps);
            DrawText(tb,lx+145,ly+215,16,LIGHTGRAY);
        } else {
            DrawRectangle(lx,ly,PANEL,PANEL,(Color){15,15,20,255});
            geo_draw_outline(lx,ly,PANEL);
            int tlen = n_trail<TRAIL_VIS ? n_trail : TRAIL_VIS;
            for (int k=1;k<tlen;k++) {
                int i0=(trail_head-k+MAX_TRAIL)%MAX_TRAIL;
                int i1=(trail_head-k-1+MAX_TRAIL)%MAX_TRAIL;
                Vector2 a=p2s(trail[i0].x,trail[i0].y,lx,ly,PANEL);
                Vector2 b=p2s(trail[i1].x,trail[i1].y,lx,ly,PANEL);
                unsigned char al=(unsigned char)(200*(1.0-(double)k/tlen));
                DrawLineEx(a,b,1.0f,(Color){100,200,255,al});
            }
            Vector2 dp=p2s(pos.x,pos.y,lx,ly,PANEL);
            DrawCircleV(dp,4,(Color){255,255,200,200});
            DrawCircleV(dp,2,WHITE);
        }
        DrawText("Trajectory",lx,ly+PANEL+4,16,LIGHTGRAY);

        /* Right panel */
        DrawTextureEx(rtex,(Vector2){rx,ry},0,scale,WHITE);
        if (view_mode == VIEW_POINCARE) {
            /* Axis labels */
            DrawText("s",rx+PANEL-12,ry+PANEL+4,14,GRAY);
            DrawText("sin a",rx+2,ry+2,12,GRAY);
            /* Draw thin border */
            DrawRectangleLines(rx,ry,PANEL,PANEL,GRAY);
        } else {
            geo_draw_outline(rx,ry,PANEL);
        }
        { const char *vn[] = {"Histogram [V]","Poincare (s, sin a) [V]"};
          DrawText(vn[view_mode],rx,ry+PANEL+4,16,LIGHTGRAY); }

        /* Info */
        int iy=GAP+PANEL+24; char buf[256];
        snprintf(buf,sizeof(buf),
            "Steps: %d  |  Bounces: %d  |  Poincare pts: %d  |  %s [G]  |  speed: %d [+/-]  |  %s",
            total_steps, total_bounces, n_poincare, geo_name(), speed,
            paused?"PAUSED":"[Space]");
        DrawText(buf,GAP,iy,14,LIGHTGRAY);
        snprintf(buf,sizeof(buf),
            "pos: (%.3f, %.3f)  |  [V] view  |  [T] turbo  |  [R] reset  |  Click to place",
            pos.x, pos.y);
        DrawText(buf,GAP,iy+18,14,GRAY);

        EndDrawing();
    }

    UnloadTexture(rtex);
    CloseWindow();
    return 0;
}
