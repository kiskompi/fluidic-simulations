#include "simulation.h"
#define max(x,y) ((x)>(y)?(x):(y))
#define min(x,y) ((x)<(y)?(x):(y))

/* Modified slightly by D. Orchard (2010) from the classic code from:

    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
    Numerical Simulation in Fluid Dynamics,
    SIAM, 1998.

    http://wissrech.ins.uni-bonn.de/research/projects/NaSt2D/index.html

*/

Simulation::Simulation()
{
}

/* Computation of tentative velocity field (f, g) */
void Simulation::compute_tentative_velocity()
{
    //f array written with [i][j] pattern
    //u array read with [i][j] [i+1][j] [i-1][j] pattern
    //v array read with [i][j] [i+2][j] [i+1][j-1] pattern
    Kokkos::parallel_for (
        iMAX*jMAX,
	KOKKOS_LAMBDA(size_t idx) {
                const int i = idx/iMAX + 1;
                const int j = idx%iMAX + 1;
                double du2dx = 0, duvdy = 0, laplu = 0;
                /* only if both adjacent cells are fluid cells */
                if ((flag(i,j) & C_F) && (flag(i+1,j) & C_F)) {
                    du2dx = (   (u(  i  , j) + u(i + 1, j))  * (u(  i  , j)   +  u(i + 1, j)) +
                                gamma * fabs(u(i,  j)   +  u(i + 1, j))  * (u(  i  , j) - u(i + 1,j)) -
                                (u(i - 1, j) + u(  i  , j))  * (u(i - 1, j)   +  u(  i  , j)) -
                                gamma * fabs(u(i - 1,j) +  u(  i  , j) ) * (u(i - 1, j) - u(i    ,j))
                            ) / (4.0 * delx);
                    duvdy = ((v(i,j)+v(i+1,j))*(u(i,j)+u(i,j+1))+
                             gamma*fabs(v(i,j)+v(i+1,j))*(u(i,j)-u(i,j+1))-
                             (v(i,j-1)+v(i+1,j-1))*(u(i,j-1)+u(i,j))-
                             gamma*fabs(v(i,j-1)+v(i+1,j-1))*(u(i,j-1)-u(i,j)))
                            /(4.0*dely);
                    laplu = (u(i+1,j)-2.0*u(i,j)+u(i-1,j))/delx/delx+
                            (u(i,j+1)-2.0*u(i,j)+u(i,j-1))/dely/dely;

                    f(i,j) = u(i,j)+step_delta*(laplu/reynolds-du2dx-duvdy);
                } else {
                    f(i,j) = u(i,j);
                }
    });

    Kokkos::parallel_for (
        iMAX*jMAX,
	KOKKOS_LAMBDA(size_t idx) {
                const int i = idx/iMAX+1;
                const int j = idx%iMAX+1;
            /* only if both adjacent cells are fluid cells */
            if ((flag(i,j) & C_F) && (flag(i,j+1) & C_F)) {
                double duvdx = 0, dv2dy = 0, laplv = 0;
                duvdx = ((u(i,j)+u(i,j+1))*(v(i,j)+v(i+1,j))+
                         gamma*fabs(u(i,j)+u(i,j+1))*(v(i,j)-v(i+1,j))-
                         (u(i-1,j)+u(i-1,j+1))*(v(i-1,j)+v(i,j))-
                         gamma*fabs(u(i-1,j)+u(i-1,j+1))*(v(i-1,j)-v(i,j)))
                        /(4.0*delx);
                dv2dy = ((v(i,j)+v(i,j+1))*(v(i,j)+v(i,j+1))+
                         gamma*fabs(v(i,j)+v(i,j+1))*(v(i,j)-v(i,j+1))-
                         (v(i,j-1)+v(i,j))*(v(i,j-1)+v(i,j))-
                         gamma*fabs(v(i,j-1)+v(i,j))*(v(i,j-1)-v(i,j)))
                        /(4.0*dely);

                laplv = (v(i+1,j)-2.0*v(i,j)+v(i-1,j))/delx/delx+
                        (v(i,j+1)-2.0*v(i,j)+v(i,j-1))/dely/dely;

                g(i,j) = v(i,j)+step_delta*(laplv/reynolds-duvdx-dv2dy);
            } else {
                g(i,j) = v(i,j);
            }
    });

    /* f & g at external boundaries */
    Kokkos::parallel_for (
        jMAX,
	KOKKOS_LAMBDA(size_t j) {
        f(0,j)    = u(0,j);
        f(iMAX,j) = u(iMAX,j);
    });
    Kokkos::parallel_for (
        iMAX,
	KOKKOS_LAMBDA(size_t i) {
        g(i,0)    = v(i,0);
        g(i,jMAX) = v(i,jMAX);
    });
}


/* Calculate the right hand side of the pressure equation */
void Simulation::compute_rhs()
{
    Kokkos::parallel_for (
        iMAX*jMAX,
        KOKKOS_LAMBDA(size_t idx) {
        const int i = idx/iMAX+1;
        const int j = idx%iMAX+1;
        if (flag(i,j) & C_F) {
            /* only for fluid and non-surface cells */
            rhs(i,j) = ((f(i,j)-f(i-1,j))/delx +
                         (g(i,j)-g(i,j-1))/dely) / step_delta;
        }
    });
}


/* Red/Black SOR to solve the poisson equation */
int Simulation::poisson()
{
    int iter = 0;
    double p0 = 0.0;

    const double rdx2 = 1.0/(delx*delx);
    const double rdy2 = 1.0/(dely*dely);
    const double beta_2 = -omega/(2.0*(rdx2+rdy2));

    /* Calculate sum of squares */
    Kokkos::parallel_reduce (
        iMAX*jMAX,
	KOKKOS_LAMBDA(size_t idx, double& inner_p0) {
        const int i = idx/iMAX + 1;
        const int j = idx%iMAX + 1;
            if (flag(i,j) & C_F) {
                inner_p0 += p(i,j)*p(i,j);
            }
    }, p0);

    p0 = sqrt(p0/ifluid);
    if (p0 < 0.0001) {
        p0 = 1.0;
    }


    /* Red/Black SOR-iteration */
    for (iter = 0; iter < itermax; iter++) {
        for (int rb = 0; rb <= 1; rb++) {
        Kokkos::parallel_for (
            iMAX*jMAX,
            KOKKOS_LAMBDA(size_t idx) {
            const int i = idx/iMAX+1;
            const int j = idx%iMAX+1;
                    /* This is tricky, try and think about what happens here */
                    if ((i+j) % 2 != rb) {
                        return;
                    }
                    if (flag(i,j) == (C_F | B_NSEW)) {
                        /* five point star for interior fluid cells */
                        p(i,j) = (1.-omega)*p(i,j) -
                                  beta_2*(
                                      (p(i+1,j)+p(i-1,j))*rdx2
                                      + (p(i,j+1)+p(i,j-1))*rdy2
                                      -  rhs(i,j));
                    } else if (flag(i,j) & C_F) {
                        /* modified star near boundary */
                        double beta_mod = -omega/((eps_E+eps_W)*rdx2+(eps_N+eps_S)*rdy2);
                        p(i,j) = (1.-omega)*p(i,j) -
                                  beta_mod*(
                                      (eps_E*p(i+1,j)+eps_W*p(i-1,j))*rdx2
                                      + (eps_N*p(i,j+1)+eps_S*p(i,j-1))*rdy2
                                     - rhs(i,j));
                    }
            }); /* end of i */
        } /* end of rb */

        /* Partial computation of residual */
        double resl = 0.0;
        Kokkos::parallel_reduce (
            iMAX,
            KOKKOS_LAMBDA(size_t i, double& resl_lmb) {
            for (int j = 1; j <= jMAX; j++) {
                if (flag(i,j) & C_F) {
                    /* only fluid cells */
                    double add = (eps_E*(p(i+1,j)-p(i,j)) -
                                  eps_W*(p(i,j)-p(i-1,j))) * rdx2  +
                                 (eps_N*(p(i,j+1)-p(i,j)) -
                                  eps_S*(p(i,j)-p(i,j-1))) * rdy2  -  rhs(i,j);
                    resl_lmb += add*add;
                }
            }
        }, resl);
        resl = sqrt((resl)/ifluid)/p0;
        res = resl;
        /* convergence? */
        if (res<eps) break;
    } /* end of iter */

    return iter;
}


/* Update the velocity values based on the tentative
 * velocity values and the new pressure matrix
 */
void Simulation::update_velocity()
{
    Kokkos::parallel_for (
        iMAX,
	KOKKOS_LAMBDA(size_t i) {
        for (int j=1; j<=jMAX; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag(i,j) & C_F) && (flag(i+1,j) & C_F)) {
                u(i,j) = f(i,j)-(p(i+1,j)-p(i,j))*step_delta/delx;
            }
        }
    });

    Kokkos::parallel_for (
        iMAX,
	KOKKOS_LAMBDA(size_t i) {
        for (int j=1; j<=jMAX-1; j++) {
            /* only if both adjacent cells are fluid cells */
            if ((flag(i,j) & C_F) && (flag(i,j+1) & C_F)) {
                v(i,j) = g(i,j)-(p(i,j+1)-p(i,j))*step_delta/dely;
            }
        }
    });
}


/* 
 * Set the timestep size so that we satisfy the Courant-Friedrichs-Lewy
 * conditions (ie no particle moves more than one cell width in one
 * timestep). Otherwise the simulation becomes unstable.
 */
void Simulation::set_timestep_interval()
{
    /* step_delta satisfying CFL conditions */
    if (tau >= 1.0e-10) { /* else no time stepsize control */
        double umax = 1.0e-10;
        double vmax = 1.0e-10;
        for (int i=1; i<=iMAX+1; i++) {
            for (int j=1; j<=jMAX+1; j++) {
                umax= max(fabs(u(i,j)), umax);
            }
        }
        for (int i=1; i<=iMAX+1; i++) {
            for (int j=0; j<=jMAX+1; j++) {
                vmax= max(fabs(v(i,j)), vmax);
            }
        }

        const double deltu = delx/umax;
        const double deltv = dely/vmax;
        const double deltRe = 1/(1/(delx*delx)+1/(dely*dely))*reynolds/2.0;

        if (deltu<deltv) {
            step_delta = min(deltu, deltRe);
        } else {
            step_delta = min(deltv, deltRe);
        }
        step_delta = tau * (step_delta); /* multiply by safety factor */
    }
}

/* Initialize the flag array, marking any obstacle cells and the edge cells
 * as boundaries. The cells adjacent to boundary cells have their relevant
 * flags set too.
 */
void Simulation::init_flag()
{
    /* Mask of a circular obstacle */
    const double mx   = 20.0/41.0 * jMAX * dely;
    const double my   = mx;
    const double rad1 = 5.00/41.0 * jMAX * dely;

    printf("%g, %g, %g, %g, %g\n", delx, dely, mx, my, rad1);
    host_view_type h_flag = Kokkos::create_mirror_view (flag);
    for (int i = 0; i <= iMAX + 1; ++i) {
        for (int j = 0; j <= jMAX + 1; ++j) {
           h_flag(i,j) = 0; 
        }
    }

    int checksum = 0;
    for (int i = 0; i <= iMAX + 1; ++i) {
        for (int j = 0; j <= jMAX + 1; ++j) {
           checksum += h_flag(i,j)*i; 
        }
    }

    for (int i = 1; i <= iMAX; ++i) {
        for (int j = 1; j <= jMAX; ++j) {
            const int x = (i - 0.5) * delx - mx;
            const int y = (j - 0.5) * dely - my;
            h_flag(i,j) = (x * x + y * y <= rad1 * rad1) ? C_B : C_F;
        }
    }
    
    printf("Checksum: %d\n", checksum);

    /* Mark the north & south boundary cells */
    for (int i = 0; i <= iMAX + 1; ++i) {
        h_flag(i,0)      = C_B;
        h_flag(i,jMAX+1) = C_B;
    }
    /* Mark the east and west boundary cells */
    for (int j = 1; j <= jMAX; ++j) {
        h_flag(0,j)      = C_B;
        h_flag(iMAX+1,j) = C_B;
    }

    /* flags for boundary cells */
    // nem kell kokkosozni
    ibound = 0;
    for (int i = 1; i <= iMAX; ++i) {
        for (int j = 1; j <= jMAX; ++j) {
            if (!(h_flag(i,j) & C_F)) {
                ++ibound;
                if (h_flag(i-1,j) & C_F) h_flag(i,j) = h_flag(i,j) | B_W;
                if (h_flag(i+1,j) & C_F) h_flag(i,j) = h_flag(j,j) | B_E;
                if (h_flag(i,j-1) & C_F) h_flag(i,j) = h_flag(j,j) | B_S;
                if (h_flag(i,j+1) & C_F) h_flag(i,j) = h_flag(j,j) | B_N;
            }
        }
    }
    ifluid = (iMAX * jMAX) - ibound;
    printf("%d\n", ibound);
    Kokkos::deep_copy (flag, h_flag); // Copy from host to device.
}

void Simulation::apply_boundary_conditions()
{
    for (int j=0; j<=jMAX+1; j++) {
        /* Fluid freely flows in from the west */
        u(0,j) = u(1,j);
        v(0,j) = v(1,j);

        /* Fluid freely flows out to the east */
        u(iMAX,j)   = u(iMAX-1,j);
        v(iMAX+1,j) = v(iMAX,j);
    }

    for (int i=0; i<=iMAX+1; i++) {
        /* The vertical velocity approaches 0 at the north and south
        * boundaries, but fluid flows freely in the horizontal direction */
        v(i,jMAX)   = 0.0;
        u(i,jMAX+1) = u(i,jMAX);

        v(i,0) = 0.0;
        u(i,0) = u(i,1);
    }

    /* Apply no-slip boundary conditions to cells that are adjacent to
    * internal obstacle cells. This forces the u and v velocity to
    * tend towards zero in these cells.
    */
    Kokkos::parallel_for (
        iMAX*jMAX,                     
        KOKKOS_LAMBDA(size_t idx) {
        const int i = idx/iMAX + 1;
        const int j = idx%iMAX + 1;    
            if (flag(i,j) & Cell::B_NSEW) {
                switch (flag(i,j)) {
                case Cell::B_N:
                    u(i,j)   = -u(i,j+1);
                    break;
                case Cell::B_E:
                    u(i,j)   = 0.0;
                    break;
                case Cell::B_NE:
                    u(i,j)   = 0.0;
                    break;
                case Cell::B_SE:
                    u(i,j)   = 0.0;
                    break;
                case Cell::B_NW:
                    u(i,j)   = -u(i,j+1);
                    break;
                case Cell::B_S:
                    u(i,j)   = -u(i,j-1);
                    break;
                case Cell::B_SW:
                    u(i,j)   = -u(i,j-1);
                    break;
                }
            }
    });
    Kokkos::parallel_for (
        iMAX*jMAX,                     
        KOKKOS_LAMBDA(size_t idx) {
        const int i = idx/iMAX + 1;
        const int j = idx%iMAX + 1;         
        if (flag(i+1,j) & Cell::B_NSEW) {
            switch (flag(i+1,j)) {
            case Cell::B_N:
                u(i,j) = -u(i,j+1);
                break;
            case Cell::B_W:
                u(i,j) = 0.0;
                break;
            case Cell::B_NE:
                u(i,j) = -u(i,j+1);
                break;
            case Cell::B_SW:
                u(i,j) = 0.0;
                break;
            case Cell::B_NW:
                u(i,j) = 0.0;
                break;
            case Cell::B_S:
                u(i,j) = -u(i,j-1);
                break;
            case Cell::B_SE:
                u(i,j) = -u(i,j-1);
                break;
            }
        }
    });

    Kokkos::parallel_for (
        iMAX*jMAX,                     
        KOKKOS_LAMBDA(size_t idx) {
        const int i = idx/iMAX + 1;
        const int j = idx%iMAX + 1;         

        if (flag(i,j) & Cell::B_NSEW) {
            switch (flag(i,j)) {
            case Cell::B_N:
                v(i,j)   = 0.0;
                break;
            case Cell::B_E:
                v(i,j)   = -v(i+1,j);
                break;
            case Cell::B_NE:
                v(i,j)   = 0.0;
                break;
            case Cell::B_SE:
                v(i,j)   = -v(i+1,j);
                break;
            case Cell::B_NW:
                v(i,j)   = 0.0;
                break;
            case Cell::B_W:
                v(i,j)   = -v(i-1,j);
                break;
            case Cell::B_SW:
                v(i,j)   = -v(i-1,j);
                break;
            }
        }
    });

     Kokkos::parallel_for (
        iMAX*jMAX,                     
        KOKKOS_LAMBDA(size_t idx) {
        const int i = idx/iMAX + 1;
        const int j = idx%iMAX + 1;         

        if (flag(i,j+1) & Cell::B_NSEW) {
            switch (flag(i,j+1)) {
            case Cell::B_E:
                v(i,j) = -v(i+1,j);
                break;
            case Cell::B_S:
                v(i,j) = 0.0;
                break;
            case Cell::B_NE:
                v(i,j) = -v(i+1,j);
                break;
            case Cell::B_SE:
                v(i,j) = 0.0;
                break;
            case Cell::B_SW:
                v(i,j) = 0.0;
                break;
            case Cell::B_W:
                v(i,j) = -v(i-1,j);
                break;
            case Cell::B_NW:
                v(i,j) = -v(i-1,j);
                break;
            }
        }
    });

    /* Finally, fix the horizontal velocity at the  western edge to have
    * a continual flow of fluid into the simulation.
    */
    v(0,0) = 2*vi-v(1,0);
    for (int j=1; j<=jMAX; j++) {
        u(0,j) = ui;
        v(0,j) = 2*vi-v(1,j);
    }
}


/* Computation of stream function and vorticity */
void  Simulation::calc_psi_zeta(DoubleMatrix zeta) const
{
    /* Computation of the vorticity zeta at the upper right corner     */
    /* of cell (i,j) (only if the corner is surrounded by fluid cells) */
     Kokkos::parallel_for (
        iMAX*jMAX,                     
        KOKKOS_LAMBDA(size_t idx) {
        const int i = idx/iMAX + 1;
        const int j = idx%iMAX + 1;         
        if (is_surrounded(i, j)) {
            zeta(i,j) = (u(i,j+1)-u(i,j))/dely-
                         (v(i+1,j)-v(i,j))/delx;
        } else {
            zeta(i,j) = 0.0;
        }
    });
}

// ====================== FRIEND FUNCTIONS ====================== //

void write_ppm(const Simulation& sim,
               char*         outname,
               int           iters,
               int           freq)
{
    double    zmax = -1e10, zmin = 1e10;
    double    pmax = -1e10, pmin = 1e10;
    double    vmax = -1e10, vmin = 1e10;
    double    umax = -1e10, umin = 1e10;
    char      outfile[64];
    char      outmode  = 0;
    Simulation::DoubleMatrix psi  = Simulation::DoubleMatrix();
    Simulation::DoubleMatrix zeta = Simulation::DoubleMatrix();

    sprintf(outfile, "%s/%06d.ppm", outname, (iters/freq));

    __mode_t process_mask = umask(0);
    mkdir(outname, S_IRWXU | S_IRWXG | S_IRWXO);
    umask(process_mask);
    FILE *fout = fopen(outfile, "wb");

    if (!fout) {
        fprintf (stderr, "Could not open '%s'\n", outfile);
        return;
    }



    sim.calc_psi_zeta(zeta);

    fprintf(fout, "P6 %d %d 255\n", iMAX, jMAX);


    for (size_t j = 1; j < jMAX+1 ; j++) {
        for (size_t i = 1; i < iMAX+1 ; i++) {
            int r = 0, g = 0, b = 0;
            if (!(sim.get_flag()(i,j) & C_F)) {
                r = 0;
                b = 0;
                g = 255;
            } else {
                if (outmode == 0) {
                    double z = (i < iMAX && j < jMAX) ? zeta(i,j) : 0.0;
                    r = g = b = pow(fabs(z/12.6),.4) * 255;
                } else if (outmode == 1) {
                    double p = (i < iMAX && j < jMAX) ? psi(i,j) : 0.0;
                    r = g = b = (p + 3.0) / 7.5 * 255;
                } else if (outmode == 2) {
                    r = g = b = (sim.get_p()(i,j)-pmin) / (pmax-pmin) * 255;
                }
            }
            fprintf(fout, "%c%c%c", r, g, b);
        }
    }

    fclose(fout);
}
