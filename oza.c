/*
 * Oza — Level 1 (Oza-Rosales-Bush, second-order with velocity/drag)
 */
#include "core.h"

#define MAX_MEM   512
#define TRAIL_VIS 200
#define MAX_POINCARE 100000

static struct { double kF, mu, beta; } P = { 12.0, 0.97, 0.005 };
static double drag = 0.95;
static int geometry = GEO_STADIUM;

static V2     mem[MAX_MEM];
static int    nmem = 0, mhead = 0;
static V2     drop = {0.25, 0.0}, vel = {0,0};
static int    total = 0, paused = 0, speed = 3, turbo = 0;
static double turbo_bps = 0;

static double wave[GRID][GRID];
static int    hist[GRID][GRID];
static int    hpeak = 1;
static double h2_sum[GRID][GRID], h4_sum[GRID][GRID];
static int    moment_n = 0;

static V2     poincare[MAX_POINCARE];
static int    n_poincare = 0, poincare_head = 0;

/* Angular momentum */
#define L_HISTORY 4000
static double L_hist[L_HISTORY];
static double L_max_val = 0.001;
static int    L_head = 0, L_count = 0;

/* Orbit diagram: (mu, x_P) from Poincare crossings */
#define MAX_ORBIT 200000
static V2     orbit_pts[MAX_ORBIT];
static int    n_orbit = 0;
static double orbit_mu_min = 1, orbit_mu_max = 0;

/* Auto-sweep */
static int    sweeping = 0;
static double sweep_from, sweep_to, sweep_step = 0.005;
static int    sweep_settle = 2000, sweep_record = 50;
static int    sweep_phase = 0, sweep_counter = 0;

enum { VIEW_HIST, VIEW_BORN, VIEW_KURTOSIS, VIEW_POINCARE, VIEW_ORBIT, VIEW_COUNT };
static int view_mode = VIEW_HIST;
static Color wpix[GRID*GRID], rpix[GRID*GRID];

static void wave_update(double bx, double by)
{
    for (int iy=0; iy<GRID; iy++) {
        double y = 2.0*iy/(GRID-1)-1.0;
        for (int ix=0; ix<GRID; ix++) {
            double x = 2.0*ix/(GRID-1)-1.0;
            wave[iy][ix] *= P.mu;
            if (geo_inside(geometry,x,y)) {
                double dx=x-bx, dy=y-by;
                wave[iy][ix] += j0(P.kF*sqrt(dx*dx+dy*dy));
            }
            double h2=wave[iy][ix]*wave[iy][ix];
            h2_sum[iy][ix]+=h2; h4_sum[iy][ix]+=h2*h2;
        }
    }
    moment_n++;
}

static void bounce(void)
{
    double bx=drop.x, by=drop.y;
    double gx=0, gy=0;
    int n = nmem<MAX_MEM ? nmem : MAX_MEM;
    double decay = 1.0;
    for (int k=0; k<n; k++) {
        int i=(mhead-1-k+MAX_MEM)%MAX_MEM;
        double dx=drop.x-mem[i].x, dy=drop.y-mem[i].y;
        double d=sqrt(dx*dx+dy*dy);
        if (d<1e-15) { decay*=P.mu; continue; }
        double f=-P.kF*j1(P.kF*d)/d;
        gx+=decay*f*dx; gy+=decay*f*dy; decay*=P.mu;
    }

    vel.x = drag*vel.x - P.beta*gx;
    vel.y = drag*vel.y - P.beta*gy;
    drop.x += vel.x; drop.y += vel.y;

    /* Angular momentum: L = x*vy - y*vx */
    {
        double L = drop.x*vel.y - drop.y*vel.x;
        L_hist[L_head] = L;
        L_head = (L_head+1) % L_HISTORY;
        if (L_count < L_HISTORY) L_count++;
        double aL = fabs(L);
        if (aL > L_max_val) L_max_val = aL*1.1;
        L_max_val *= 0.99999;
        if (L_max_val < 0.001) L_max_val = 0.001;
    }

    /* Wall */
    double d = geo_sdf(geometry, drop.x, drop.y);
    if (d < 0.15) {
        double nx, ny; geo_normal(geometry, drop.x, drop.y, &nx, &ny);
        double nm=sqrt(nx*nx+ny*ny);
        if (nm>1e-10) { nx/=nm; ny/=nm; double s=(0.15-d)/0.15;
            drop.x-=0.03*s*s*nx; drop.y-=0.03*s*s*ny;
            double vn=vel.x*nx+vel.y*ny;
            if (vn>0) { vel.x-=1.5*vn*nx; vel.y-=1.5*vn*ny; } }
    }
    d = geo_sdf(geometry, drop.x, drop.y);
    if (d<0.02) { double nx,ny; geo_normal(geometry,drop.x,drop.y,&nx,&ny);
        double nm=sqrt(nx*nx+ny*ny);
        if (nm>1e-10) { drop.x-=(0.02-d)*nx/nm; drop.y-=(0.02-d)*ny/nm; } }

    /* Poincare */
    if (by<=0 && drop.y>0) {
        double dy=drop.y-by, frac=(dy>1e-15)?-by/dy:0.5;
        double xc=bx+frac*(drop.x-bx);
        if (n_poincare<MAX_POINCARE) poincare[n_poincare++]=(V2){xc,vel.x};
        else { poincare[poincare_head]=(V2){xc,vel.x}; poincare_head=(poincare_head+1)%MAX_POINCARE; }

        /* Orbit diagram */
        if (n_orbit < MAX_ORBIT) {
            orbit_pts[n_orbit++] = (V2){P.mu, xc};
            if (P.mu < orbit_mu_min) orbit_mu_min = P.mu;
            if (P.mu > orbit_mu_max) orbit_mu_max = P.mu;
        }
        /* Auto-sweep: count crossings in record phase */
        if (sweeping && sweep_phase == 1) {
            if (++sweep_counter >= sweep_record) {
                P.mu += sweep_step;
                if (P.mu > sweep_to + 0.0001) sweeping = 0;
                else { sweep_phase = 0; sweep_counter = 0; }
            }
        }
    }
    /* Auto-sweep: count settle steps */
    if (sweeping && sweep_phase == 0) {
        if (++sweep_counter >= sweep_settle) { sweep_phase = 1; sweep_counter = 0; }
    }

    mem[mhead]=(V2){bx,by}; mhead=(mhead+1)%MAX_MEM;
    if (nmem<MAX_MEM) nmem++; total++;
    wave_update(bx, by);
    hist_splat(hist, &hpeak, drop.x, drop.y, 1.0);
}

static void reset_sim(void)
{
    nmem=0; mhead=0; total=0; hpeak=1; moment_n=0;
    n_poincare=0; poincare_head=0;
    n_orbit=0; orbit_mu_min=1; orbit_mu_max=0; sweeping=0;
    L_head=0; L_count=0; L_max_val=0.001;
    drop=(V2){0.25,0.0}; vel=(V2){0,0};
    memset(wave,0,sizeof(wave)); memset(hist,0,sizeof(hist));
    memset(h2_sum,0,sizeof(h2_sum)); memset(h4_sum,0,sizeof(h4_sum));
}

int main(void)
{
    InitWindow(WIN_W_L, WIN_H_L, "Oza - Level 1 (velocity + drag)");
    SetTargetFPS(30);
    Image img=GenImageColor(GRID,GRID,BLACK);
    Texture2D wtex=LoadTextureFromImage(img), rtex=LoadTextureFromImage(img);
    UnloadImage(img);

    while (!WindowShouldClose()) {
        if (IsKeyPressed(KEY_SPACE)) paused=!paused;
        if (IsKeyPressed(KEY_R))     reset_sim();
        if (IsKeyPressed(KEY_G))     { geometry=(geometry+1)%GEO_COUNT; reset_sim(); }
        if (IsKeyPressed(KEY_V))     view_mode=(view_mode+1)%VIEW_COUNT;
        if (IsKeyPressed(KEY_T))     { turbo=!turbo; SetTargetFPS(turbo?10:30); }
        if (IsKeyPressed(KEY_UP))   { P.mu=fmin(P.mu+0.005,0.999); n_poincare=0; poincare_head=0; }
        if (IsKeyPressed(KEY_DOWN)) { P.mu=fmax(P.mu-0.005,0.5); n_poincare=0; poincare_head=0; }
        if (IsKeyPressed(KEY_RIGHT)) P.beta*=1.15;
        if (IsKeyPressed(KEY_LEFT))  P.beta/=1.15;
        if (IsKeyPressed(KEY_D))     drag=fmin(drag+0.01,0.999);
        if (IsKeyPressed(KEY_F))     drag=fmax(drag-0.01,0.0);
        if (IsKeyPressed(KEY_L))     { P.kF+=1.0; reset_sim(); }
        if (IsKeyPressed(KEY_K))     { P.kF=fmax(P.kF-1.0,2.0); reset_sim(); }
        if (IsKeyPressed(KEY_EQUAL)) speed=speed<50?speed+1:speed;
        if (IsKeyPressed(KEY_MINUS)) speed=speed>1?speed-1:speed;
        if (IsKeyPressed(KEY_O)) {
            sweeping=1; sweep_from=P.mu; sweep_to=fmin(P.mu+0.1,0.999);
            sweep_step=0.005; sweep_phase=0; sweep_counter=0;
            view_mode=VIEW_ORBIT; speed=50;
            printf("Orbit sweep: mu %.3f -> %.3f\n",sweep_from,sweep_to);
        }

        if (IsMouseButtonPressed(MOUSE_BUTTON_LEFT)&&!turbo) {
            Vector2 m=GetMousePosition();
            if (m.x>=GAP&&m.x<GAP+PANEL&&m.y>=GAP&&m.y<GAP+PANEL) {
                double px=((m.x-GAP)/PANEL*2-1), py=((m.y-GAP)/PANEL*2-1);
                if (geo_inside(geometry,px,py)) { drop.x=px; drop.y=py; vel.x=0; vel.y=0; }
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

        if (!turbo) { render_wave(wpix,wave,geometry,1.0); UpdateTexture(wtex,wpix); }
        switch (view_mode) {
        case VIEW_HIST:     render_hist(rpix,hist,hpeak,geometry,1.0); break;
        case VIEW_BORN:     render_born(rpix,h2_sum,moment_n,geometry,1.0); break;
        case VIEW_KURTOSIS: render_kurtosis(rpix,h2_sum,h4_sum,moment_n,geometry,1.0); break;
        case VIEW_POINCARE: render_poincare(rpix,poincare,n_poincare,poincare_head,MAX_POINCARE); break;
        case VIEW_ORBIT: {
            /* Orbit diagram: mu on x, x_P on y */
            for (int i=0;i<GRID*GRID;i++) rpix[i]=(Color){15,15,20,255};
            if (n_orbit > 0) {
                double mu_lo=orbit_mu_min-0.01, mu_hi=orbit_mu_max+0.01;
                if (mu_hi-mu_lo<0.02) { mu_lo-=0.01; mu_hi+=0.01; }
                double x_lo=1e30,x_hi=-1e30;
                for (int i=0;i<n_orbit;i++) {
                    if (orbit_pts[i].y<x_lo) x_lo=orbit_pts[i].y;
                    if (orbit_pts[i].y>x_hi) x_hi=orbit_pts[i].y;
                }
                double xm=(x_hi-x_lo)*0.1+0.01; x_lo-=xm; x_hi+=xm;
                static int od[GRID][GRID]; memset(od,0,sizeof(od)); int om=1;
                for (int i=0;i<n_orbit;i++) {
                    int px=(int)((orbit_pts[i].x-mu_lo)/(mu_hi-mu_lo)*(GRID-1));
                    int py=(int)((1.0-(orbit_pts[i].y-x_lo)/(x_hi-x_lo))*(GRID-1));
                    if (px>=0&&px<GRID&&py>=0&&py<GRID) {
                        od[py][px]++; if (od[py][px]>om) om=od[py][px]; }
                }
                double lm=log(1.0+om);
                for (int iy=0;iy<GRID;iy++) for (int ix=0;ix<GRID;ix++)
                    if (od[iy][ix]>0) rpix[iy*GRID+ix]=cmap_hot(log(1.0+od[iy][ix])/lm);
            }
            break;
        }
        }
        UpdateTexture(rtex,rpix);

        BeginDrawing();
        ClearBackground((Color){20,20,25,255});
        int lx=GAP,ly=GAP,rx=GAP*2+PANEL,ry=GAP;
        float sc=(float)PANEL/GRID;

        if (turbo) { draw_turbo(lx,ly,PANEL,turbo_bps,"bounces"); geo_draw_outline(geometry,lx,ly,PANEL); }
        else {
            DrawTextureEx(wtex,(Vector2){lx,ly},0,sc,WHITE);
            geo_draw_outline(geometry,lx,ly,PANEL);
            draw_trail(mem,mhead,nmem<MAX_MEM?nmem:MAX_MEM,TRAIL_VIS,lx,ly,PANEL,1.0,(Color){255,255,100,255});
            draw_droplet(drop.x,drop.y,lx,ly,PANEL,1.0);
        }
        DrawText(turbo?"TURBO [T]":"Wave Field",lx,ly+PANEL+4,16,LIGHTGRAY);

        DrawTextureEx(rtex,(Vector2){rx,ry},0,sc,WHITE);
        if (view_mode==VIEW_POINCARE || view_mode==VIEW_ORBIT) {
            DrawRectangleLines(rx,ry,PANEL,PANEL,GRAY);
            if (sweeping && view_mode==VIEW_ORBIT) {
                char sb[64]; snprintf(sb,sizeof(sb),"Sweeping mu=%.3f",P.mu);
                DrawText(sb,rx+4,ry+4,14,(Color){255,200,100,255});
            }
        } else geo_draw_outline(geometry,rx,ry,PANEL);
        { const char *vn[]={"Histogram","Born <h^2>","Kurtosis","Poincare (x,vx)","Orbit Diagram"};
          DrawText(vn[view_mode],rx,ry+PANEL+4,16,LIGHTGRAY); }

        /* Angular momentum strip */
        int Ly = GAP+PANEL+20, Lw = PANEL*2+GAP;
        draw_L_strip(L_hist, L_head, L_count, L_HISTORY, L_max_val, lx, Ly, Lw, LSTRIP);

        int iy = Ly+LSTRIP+4; char buf[256];
        snprintf(buf,sizeof(buf),"Bounces: %d | mu: %.3f | beta: %.4f | drag: %.3f [D/F] | speed: %d | %s",
            total,P.mu,P.beta,drag,speed,paused?"PAUSED":"[Space]");
        DrawText(buf,GAP,iy,14,LIGHTGRAY);
        snprintf(buf,sizeof(buf),"kF: %.1f [K/L] | %s [G] | vel: %.4f | Orbit: %d | [O] sweep | [V] view | [T] turbo",
            P.kF,geo_name(geometry),sqrt(vel.x*vel.x+vel.y*vel.y),n_orbit);
        DrawText(buf,GAP,iy+18,14,GRAY);

        EndDrawing();
    }
    UnloadTexture(wtex); UnloadTexture(rtex); CloseWindow();
    return 0;
}
