#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "datadef.h"
#include "init.h"

#define max(x,y) ((x)>(y)?(x):(y))
#define min(x,y) ((x)<(y)?(x):(y))

/* Modified slightly by D. Orchard (2010) from the classic code from: 

    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
    Numerical Simulation in Fluid Dynamics,
    SIAM, 1998.

    http://wissrech.ins.uni-bonn.de/research/projects/NaSt2D/index.html

*/

/* Computation of tentative velocity field (f, g) */
void compute_tentative_velocity(
    int jMAX, 
    DoubleMatrix& u, 
    DoubleMatrix& v, 
    DoubleMatrix& f, 
    DoubleMatrix& g,
    CharMatrix&   flag, 
    int      iMAX,  
    double   del_t, 
    double   delx, 
    double   dely,
    double   gamma, 
    double   Re)
{
    //f array written with [i][j] pattern
    //u array read with [i][j] [i+1][j] [i-1][j] pattern
    //v array read with [i][j] [i+1][j] [i+1][j-1] pattern
    for (int i=1; i<=iMAX-1; i++) {
        for (int j=1; j<=jMAX; j++) {
            double du2dx, duvdy, laplu;
            /* only if both adjacent cells are fluid cells */
            if ((flag[i][j] & C_F) && (flag[i+1][j] & C_F)) {
                du2dx = ((u[i][j]+u[i+1][j])*(u[i][j]+u[i+1][j])+
                    gamma*fabs(u[i][j]+u[i+1][j])*(u[i][j]-u[i+1][j])-
                    (u[i-1][j]+u[i][j])*(u[i-1][j]+u[i][j])-
                    gamma*fabs(u[i-1][j]+u[i][j])*(u[i-1][j]-u[i][j]))
                    /(4.0*delx);
                duvdy = ((v[i][j]+v[i+1][j])*(u[i][j]+u[i][j+1])+
                    gamma*fabs(v[i][j]+v[i+1][j])*(u[i][j]-u[i][j+1])-
                    (v[i][j-1]+v[i+1][j-1])*(u[i][j-1]+u[i][j])-
                    gamma*fabs(v[i][j-1]+v[i+1][j-1])*(u[i][j-1]-u[i][j]))
                    /(4.0*dely);
                laplu = (u[i+1][j]-2.0*u[i][j]+u[i-1][j])/delx/delx+
                    (u[i][j+1]-2.0*u[i][j]+u[i][j-1])/dely/dely;

                f[i][j] = u[i][j]+del_t*(laplu/Re-du2dx-duvdy);
            } else {
                f[i][j] = u[i][j];
            }
        }
    }

    for (int i=1; i<=iMAX; i++) {
        for (int j=1; j<=jMAX-1; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[i][j] & C_F) && (flag[i][j+1] & C_F)) {
                double duvdx, dv2dy, laplv;
                duvdx = ((u[i][j]+u[i][j+1])*(v[i][j]+v[i+1][j])+
                    gamma*fabs(u[i][j]+u[i][j+1])*(v[i][j]-v[i+1][j])-
                    (u[i-1][j]+u[i-1][j+1])*(v[i-1][j]+v[i][j])-
                    gamma*fabs(u[i-1][j]+u[i-1][j+1])*(v[i-1][j]-v[i][j]))
                    /(4.0*delx);
                dv2dy = ((v[i][j]+v[i][j+1])*(v[i][j]+v[i][j+1])+
                    gamma*fabs(v[i][j]+v[i][j+1])*(v[i][j]-v[i][j+1])-
                    (v[i][j-1]+v[i][j])*(v[i][j-1]+v[i][j])-
                    gamma*fabs(v[i][j-1]+v[i][j])*(v[i][j-1]-v[i][j]))
                    /(4.0*dely);

                laplv = (v[i+1][j]-2.0*v[i][j]+v[i-1][j])/delx/delx+
                    (v[i][j+1]-2.0*v[i][j]+v[i][j-1])/dely/dely;

                g[i][j] = v[i][j]+del_t*(laplv/Re-duvdx-dv2dy);
            } else {
                g[i][j] = v[i][j];
            }
        }
    }

    /* f & g at external boundaries */
    for (int j=1; j<=jMAX; j++) {
        f[0][j]    = u[0][j];
        f[iMAX][j] = u[iMAX][j];
    }
    for (int i=1; i<=iMAX; i++) {
        g[i][0]    = v[i][0];
        g[i][jMAX] = v[i][jMAX];
    }
}


/* Calculate the right hand side of the pressure equation */
void compute_rhs(
    int jMAX, 
    DoubleMatrix& f, 
    DoubleMatrix& g, 
    DoubleMatrix& rhs, 
    CharMatrix&   flag, 
    int      iMAX,
    double   del_t, 
    double   delx, 
    double   dely)
{
    for (int i=1;i<=iMAX;i++) {
        for (int j=1;j<=jMAX;j++) {
            if (flag[i][j] & C_F) {
                /* only for fluid and non-surface cells */
                rhs[i][j] = (
                   (f[i][j]-f[i-1][j])/delx +
                   (g[i][j]-g[i][j-1])/dely) / del_t;
            }
        }
    }
}


/* Red/Black SOR to solve the poisson equation */
int poisson(
    int jMAX, 
    DoubleMatrix& p, 
    DoubleMatrix& rhs, 
    CharMatrix&   flag, 
    double&  res, 
    int      ifull)
{
    int iter;
    double beta_2;
    double p0 = 0.0;

    double rdx2 = 1.0/(delx*delx);
    double rdy2 = 1.0/(dely*dely);
    beta_2 = -omega/(2.0*(rdx2+rdy2));

    /* Calculate sum of squares */
    for (int i = 1; i <= iMAX; i++) {
        for (int j=1; j<=jMAX; j++) {
            if (flag[i][j] & C_F) { p0 += p[i][j]*p[i][j]; }
        }
    }

    p0 = sqrt(p0/ifull);
    if (p0 < 0.0001) { p0 = 1.0; }


    /* Red/Black SOR-iteration */
    for (iter = 0; iter < itermax; iter++) {
        for (int rb = 0; rb <= 1; rb++) {
            for (int i = 1; i <= iMAX; i++) {
                for (int j = 1; j <= jMAX; j++) {
                    /* This is tricky, try and think about what happens here */
                    if ((i+j) % 2 != rb) { continue; }
                    if (flag[i][j] == (C_F | B_NSEW)) {
                        /* five point star for interior fluid cells */
                        p[i][j] = (1.-omega)*p[i][j] - 
                            beta_2*(
                                (p[i+1][j]+p[i-1][j])*rdx2
                                + (p[i][j+1]+p[i][j-1])*rdy2
                                -  rhs[i][j]);
                    } else if (flag[i][j] & C_F) { 
                        /* modified star near boundary */
                        double beta_mod = -omega/((eps_E+eps_W)*rdx2+(eps_N+eps_S)*rdy2);
                        p[i][j] = (1.-omega)*p[i][j] -
                            beta_mod*(
                                (eps_E*p[i+1][j]+eps_W*p[i-1][j])*rdx2
                                + (eps_N*p[i][j+1]+eps_S*p[i][j-1])*rdy2
                                - rhs[i][j]);
                    }
                } /* end of j */
            } /* end of i */
        } /* end of rb */

        /* Partial computation of residual */
        double resl = 0.0;
        for (int i = 1; i <= iMAX; i++) {
            for (int j = 1; j <= jMAX; j++) {
                if (flag[i][j] & C_F) {
        /* only fluid cells */
                    double add = (eps_E*(p[i+1][j]-p[i][j]) - 
                        eps_W*(p[i][j]-p[i-1][j])) * rdx2  +
                        (eps_N*(p[i][j+1]-p[i][j]) -
                        eps_S*(p[i][j]-p[i][j-1])) * rdy2  -  rhs[i][j];
                    resl += add*add;
                }
            }
        }
        resl = sqrt((resl)/ifull)/p0;
        res = resl;
        /* convergence? */
        if (res<eps) break;
    } /* end of iter */

    return iter;
}


/* Update the velocity values based on the tentative
 * velocity values and the new pressure matrix
 */
void update_velocity(int jMAX, 
                     DoubleMatrix& u; 
                     DoubleMatrix& v;
                     DoubleMatrix& f 
                     DoubleMatrix& g, 
                     DoubleMatrix& p,
                     CharMatrix&   flag 
)
{

    for (int i=1; i<=iMAX-1; i++) {
        for (int j=1; j<=jMAX; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[i][j] & C_F) && (flag[i+1][j] & C_F)) {
                u[i][j] = f[i][j]-(p[i+1][j]-p[i][j])*del_t/delx;
            }
        }
    }

    for (int i=1; i<=iMAX; i++) {
        for (int j=1; j<=jMAX-1; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag[i][j] & C_F) && (flag[i][j+1] & C_F)) {
               v[i][j] = g[i][j]-(p[i][j+1]-p[i][j])*del_t/dely;
           }
       }
   }
}


/* Set the timestep size so that we satisfy the Courant-Friedrichs-Lewy
 * conditions (ie no particle moves more than one cell width in one
 * timestep). Otherwise the simulation becomes unstable.
 */
void set_timestep_interval(
    DoubleMatrix& u, 
    DoubleMatrix& v 
    )
{
    /* del_t satisfying CFL conditions */
    if (tau >= 1.0e-10) { /* else no time stepsize control */
        double umax = 1.0e-10;
        double vmax = 1.0e-10; 
        for (int i=0; i<=iMAX+1; i++) {
            for (int j=1; j<=jMAX+1; j++) {
                umax = max(fabs(u[i][j]), umax);
            }
        }
        for (int i=1; i<=iMAX+1; i++) {
            for (int j=0; j<=jMAX+1; j++) {
                vmax = max(fabs(v[i][j]), vmax);
            }
        }

        double deltu = delx/umax;
        double deltv = dely/vmax; 
        double deltRe = 1/(1/(delx*delx)+1/(dely*dely))*Re/2.0;

        if (deltu<deltv) {
            *del_t = min(deltu, deltRe);
        } else {
            *del_t = min(deltv, deltRe);
        }
            *del_t = tau * (*del_t); /* multiply by safety factor */
    }
}
