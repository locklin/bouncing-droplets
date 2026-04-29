/*
 * Classical Billiard — elastic point particle, no waves.
 * Poincare section in Birkhoff coordinates (s, sin alpha).
 */
#include "core.h"

#define MAX_TRAIL 8000
#define TRAIL_VIS 2000
#define DT        0.002
#define MAX_POINCARE 100000

static int geometry = GEO_STADIUM;

static V2     pos, vel;
static V2     trail[MAX_TRAIL];
static int    n_trail = 0, trail_head = 0;
static int    total_bounces = 0, total_steps = 0;
static int    paused = 0, speed = 200, turbo = 0;
static double turbo_sps = 0;

static int    hist[GRID][GRID];
static int    hpeak = 1;

/* Birkhoff Poincare: (s, sin_alpha) binned into density grid */
static int    poincare_hist[GRID][GRID];
static int    ppeak = 1;
static int    n_poincare = 0;

/* (X, Vx) Poincare: y=0 upward crossings, like walker/oza */
#define MAX_POINCARE_XV 100000
static V2     poincare_xv[MAX_POINCARE_XV];
static int    n_poincare_xv = 0, poincare_xv_head = 0;

enum { VIEW_HIST, VIEW_POINCARE, VIEW_POINCARE_XV, VIEW_COUNT };
static int view_mode = VIEW_HIST;
static Color rpix[GRID*GRID];

/* ----- Billiard Step ----- */

static void billiard_step(void)
{
    double old_y = pos.y;
    double dt_remain = DT;
    for (int bounce = 0; bounce < 10 && dt_remain > 1e-10; bounce++) {
        double nx = pos.x + vel.x*dt_remain, ny = pos.y + vel.y*dt_remain;
        if (geo_sdf(geometry, nx, ny) > 0) { pos.x=nx; pos.y=ny; break; }

        double t_lo=0, t_hi=dt_remain;
        for (int i=0; i<40; i++) {
            double t_mid=(t_lo+t_hi)/2;
            double mx=pos.x+vel.x*t_mid, my=pos.y+vel.y*t_mid;
            if (geo_sdf(geometry,mx,my)>0) t_lo=t_mid; else t_hi=t_mid;
        }
        pos.x += vel.x*t_lo; pos.y += vel.y*t_lo;
        dt_remain -= t_lo;

        double wnx, wny;
        geo_normal(geometry, pos.x, pos.y, &wnx, &wny);
        double nm = sqrt(wnx*wnx+wny*wny);
        if (nm > 1e-10) {
            wnx/=nm; wny/=nm;
            double vn = vel.x*wnx+vel.y*wny;
            double vt_x=vel.x-vn*wnx, vt_y=vel.y-vn*wny;
            double vt=sqrt(vt_x*vt_x+vt_y*vt_y);
            double vt_total=sqrt(vel.x*vel.x+vel.y*vel.y);
            double sin_alpha = (vt_total>1e-12) ? vt/vt_total : 0;
            if (wnx*vt_y - wny*vt_x < 0) sin_alpha = -sin_alpha;

            double s = geo_arclength(geometry, pos.x, pos.y);
            int bx=(int)(s*(GRID-1)), by=(int)((1.0-(sin_alpha*0.5+0.5))*(GRID-1));
            if (bx>=0&&bx<GRID&&by>=0&&by<GRID) {
                poincare_hist[by][bx]++;
                if (poincare_hist[by][bx]>ppeak) ppeak=poincare_hist[by][bx];
            }
            n_poincare++;

            if (vn > 0) { vel.x -= 2*vn*wnx; vel.y -= 2*vn*wny; }
        }
        pos.x -= 0.001*wnx; pos.y -= 0.001*wny;
        total_bounces++;
    }
    total_steps++;

    /* (X, Vx) Poincare: y=0 upward crossing */
    if (old_y <= 0 && pos.y > 0) {
        double dy = pos.y - old_y;
        double frac = (dy > 1e-15) ? -old_y / dy : 0.5;
        double xc = pos.x - vel.x*dt_remain*(1-frac); /* approximate crossing x */
        (void)xc; xc = pos.x; /* simpler: just use current x */
        if (n_poincare_xv < MAX_POINCARE_XV)
            poincare_xv[n_poincare_xv++] = (V2){pos.x, vel.x};
        else {
            poincare_xv[poincare_xv_head] = (V2){pos.x, vel.x};
            poincare_xv_head = (poincare_xv_head+1) % MAX_POINCARE_XV;
        }
    }

    trail[trail_head] = pos;
    trail_head = (trail_head+1)%MAX_TRAIL;
    if (n_trail < MAX_TRAIL) n_trail++;

    static int hcnt = 0;
    if (++hcnt >= 5) { hcnt=0; hist_splat(hist, &hpeak, pos.x, pos.y, 1.0); }
}

static void render_right(Texture2D tex)
{
    if (view_mode == VIEW_HIST) {
        render_hist(rpix, hist, hpeak, geometry, 1.0);
    } else if (view_mode == VIEW_POINCARE) {
        /* Birkhoff Poincare: pre-binned density */
        double lmin=1e30, lmax=0;
        for (int i=0; i<GRID*GRID; i++)
            if (poincare_hist[i/GRID][i%GRID]>0) {
                double lv=log(1.0+poincare_hist[i/GRID][i%GRID]);
                if (lv<lmin) lmin=lv; if (lv>lmax) lmax=lv; }
        double lr=lmax-lmin; if (lr<1e-10) lr=1;
        for (int iy=0;iy<GRID;iy++) for (int ix=0;ix<GRID;ix++) {
            int idx=iy*GRID+ix;
            rpix[idx] = poincare_hist[iy][ix]>0
                ? cmap_hot((log(1.0+poincare_hist[iy][ix])-lmin)/lr)
                : (Color){15,15,20,255};
        }
        int mid=GRID/2;
        for (int ix=0;ix<GRID;ix++)
            if (poincare_hist[mid][ix]==0) rpix[mid*GRID+ix]=(Color){40,40,50,255};
    } else if (view_mode == VIEW_POINCARE_XV) {
        /* (X, Vx) Poincare from y=0 crossings — same as walker/oza */
        render_poincare(rpix, poincare_xv, n_poincare_xv, poincare_xv_head, MAX_POINCARE_XV);
    }
    UpdateTexture(tex, rpix);
}

static void reset_sim(void)
{
    srand(time(NULL));
    do { pos.x=((double)rand()/RAND_MAX-0.5)*1.4;
         pos.y=((double)rand()/RAND_MAX-0.5)*1.4;
    } while (!geo_inside(geometry,pos.x,pos.y));
    double angle=(double)rand()/RAND_MAX*2*M_PI;
    vel.x=0.5*cos(angle); vel.y=0.5*sin(angle);
    n_trail=0; trail_head=0; total_bounces=0; total_steps=0;
    hpeak=1; n_poincare=0; ppeak=1;
    n_poincare_xv=0; poincare_xv_head=0;
    memset(hist,0,sizeof(hist)); memset(poincare_hist,0,sizeof(poincare_hist));
}

int main(void)
{
    InitWindow(WIN_W, WIN_H, "Classical Billiard");
    SetTargetFPS(30);
    Image img=GenImageColor(GRID,GRID,BLACK);
    Texture2D rtex=LoadTextureFromImage(img);
    UnloadImage(img);
    reset_sim();

    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_SPACE)) paused=!paused;
        if (IsKeyPressed(KEY_R))     reset_sim();
        if (IsKeyPressed(KEY_G))     { geometry=(geometry+1)%GEO_COUNT; reset_sim(); }
        if (IsKeyPressed(KEY_V))     view_mode=(view_mode+1)%VIEW_COUNT;
        if (IsKeyPressed(KEY_T))     { turbo=!turbo; SetTargetFPS(turbo?10:30); }
        if (IsKeyPressed(KEY_EQUAL)) speed=speed<5000?speed+200:speed;
        if (IsKeyPressed(KEY_MINUS)) speed=speed>200?speed-200:speed;

        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)) {
            Vector2 m=GetMousePosition();
            if (m.x>=GAP&&m.x<GAP+PANEL&&m.y>=GAP&&m.y<GAP+PANEL) {
                double px=((m.x-GAP)/PANEL*2-1), py=((m.y-GAP)/PANEL*2-1);
                if (geo_inside(geometry,px,py)) {
                    pos.x=px; pos.y=py;
                    double r=sqrt(px*px+py*py);
                    if (r>0.01) { vel.x=-0.5*py/r; vel.y=0.5*px/r; }
                    n_trail=0; trail_head=0;
                }
            }
        }

        if (!paused) {
            if (turbo) {
                struct timespec t0,t1; clock_gettime(CLOCK_MONOTONIC,&t0);
                int count=0;
                for (;;) { billiard_step(); count++;
                    if (count%2000==0) { clock_gettime(CLOCK_MONOTONIC,&t1);
                        double el=(t1.tv_sec-t0.tv_sec)+(t1.tv_nsec-t0.tv_nsec)*1e-9;
                        if (el>0.2) { turbo_sps=count/el; break; } } }
            } else { for (int i=0;i<speed;i++) billiard_step(); }
        }

        render_right(rtex);
        BeginDrawing();
        ClearBackground((Color){20,20,25,255});
        int lx=GAP,ly=GAP,rx=GAP*2+PANEL,ry=GAP;
        float sc=(float)PANEL/GRID;

        if (turbo) { draw_turbo(lx,ly,PANEL,turbo_sps,"steps"); geo_draw_outline(geometry,lx,ly,PANEL); }
        else {
            DrawRectangle(lx,ly,PANEL,PANEL,(Color){15,15,20,255});
            geo_draw_outline(geometry,lx,ly,PANEL);
            draw_trail(trail,trail_head,n_trail,TRAIL_VIS,lx,ly,PANEL,1.0,(Color){100,200,255,255});
            draw_droplet(pos.x,pos.y,lx,ly,PANEL,1.0);
        }
        DrawText("Trajectory",lx,ly+PANEL+4,16,LIGHTGRAY);

        DrawTextureEx(rtex,(Vector2){rx,ry},0,sc,WHITE);
        if (view_mode==VIEW_POINCARE) {
            DrawRectangleLines(rx,ry,PANEL,PANEL,GRAY);
            DrawText("s",rx+PANEL-12,ry+PANEL+4,12,GRAY);
            DrawText("sin a",rx+2,ry+2,12,GRAY);
        } else if (view_mode==VIEW_POINCARE_XV) {
            DrawRectangleLines(rx,ry,PANEL,PANEL,GRAY);
            DrawText("x",rx+PANEL-12,ry+PANEL+4,12,GRAY);
            DrawText("vx",rx+2,ry+2,12,GRAY);
        } else geo_draw_outline(geometry,rx,ry,PANEL);
        { const char *vn[]={"Histogram","Birkhoff (s,sin a)","Poincare (x,vx)"};
          DrawText(vn[view_mode],rx,ry+PANEL+4,16,LIGHTGRAY); }

        int iy=GAP+PANEL+24; char buf[256];
        snprintf(buf,sizeof(buf),"Steps: %d | Bounces: %d | Poincare: %d | %s [G] | speed: %d | %s",
            total_steps,total_bounces,n_poincare,geo_name(geometry),speed,paused?"PAUSED":"[Space]");
        DrawText(buf,GAP,iy,14,LIGHTGRAY);
        snprintf(buf,sizeof(buf),"[V] view | [T] turbo | [R] reset | Click to place");
        DrawText(buf,GAP,iy+18,14,GRAY);
        EndDrawing();
    }
    UnloadTexture(rtex); CloseWindow();
    return 0;
}
