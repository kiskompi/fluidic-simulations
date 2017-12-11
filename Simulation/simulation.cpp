#include "simulation.h"
#include <Kokkos_Macros.hpp>
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
    double_host_view h_p= Kokkos::create_mirror_view (p);
    for (int i = 0; i <= iMAX + 1; ++i) {
        for (int j = 0; j <= jMAX + 1; ++j) {
           h_p(i,j) = 0.0; 
        }
    }
    Kokkos::deep_copy (p, h_p); // Copy from host to device.
 
    double_host_view h_g= Kokkos::create_mirror_view (g);
    for (int i = 0; i <= iMAX + 1; ++i) {
        for (int j = 0; j <= jMAX + 1; ++j) {
           h_g(i,j) = 0.0; 
        }
    }
    Kokkos::deep_copy (g, h_g); // Copy from host to device.
   
    double_host_view h_rhs= Kokkos::create_mirror_view (rhs);
    for (int i = 0; i <= iMAX + 1; ++i) {
        for (int j = 0; j <= jMAX + 1; ++j) {
           h_rhs(i,j) = 0.0; 
        }
    }
    Kokkos::deep_copy (rhs, h_rhs); // Copy from host to device.

}

/* Computation of tentative velocity field (f, g) */
void Simulation::compute_tentative_velocity()
{
    //f array written with [i][j] pattern
    //u array read with [i][j] [i+1][j] [i-1][j] pattern
    //v array read with [i][j] [i+2][j] [i+1][j-1] pattern
    auto u_t = u;
    auto v_t = v;
    auto f_t = f;
    auto g_t = g;
    auto p_t = p;
    double step_delta_t = step_delta;
    double delx_t = delx;
    double dely_t = dely;
    auto flag_t = flag;
    double gamma_t = gamma;
    double reynolds_t = reynolds;

    Kokkos::parallel_for (
        iMAX*jMAX,
	KOKKOS_LAMBDA(size_t idx) {
                const int i = idx/iMAX + 1;
                const int j = idx%iMAX + 1;
                /* only if both adjacent cells are fluid cells */
                if ((flag_t(i,j) & C_F) && (flag_t(i+1,j) & C_F)) {
                    const double du2dx = 
                            ((u_t(i,j)+u_t(i+1,j))*(u_t(i,j)+ u_t(i+1,j))+
                              gamma_t*fabs(u_t(i,j)+u_t(i+1, j))*(u_t(i, j)-u_t(i+1,j)) -
                              (u_t(i-1,j)+u_t(i,j))*(u_t(i-1,j)+u_t(i,j))-
                              gamma_t*fabs(u_t(i-1,j)+u_t(i,j))*(u_t(i-1,j)-u_t(i,j))
                            )/(4.0*delx_t);
                    const double duvdy = 
                            ((v_t(i,j)+v_t(i+1,j))*(u_t(i,j)+u_t(i,j+1))+
                             gamma_t*fabs(v_t(i,j)+v_t(i+1,j))*(u_t(i,j)-u_t(i,j+1))-
                             (v_t(i,j-1)+v_t(i+1,j-1))*(u_t(i,j-1)+u_t(i,j))-
                             gamma_t*fabs(v_t(i,j-1)+v_t(i+1,j-1))*(u_t(i,j-1)-u_t(i,j)))
                            /(4.0*dely_t);
                    const double laplu = 
                            (u_t(i+2,j)-2.0*u_t(i,j)+u_t(i-1,j))
                            /delx_t/delx_t+
                            (u_t(i,j+1)-2.0*u_t(i,j)+u_t(i,j-1))
                            /dely_t/dely_t;

                    f_t(i,j) = u_t(i,j)+step_delta_t*(laplu/reynolds_t-du2dx-duvdy);
                } else {
                    f_t(i,j) = u_t(i,j);
                }
    });

    Kokkos::parallel_for (
        iMAX*jMAX,
	KOKKOS_LAMBDA(size_t idx) {
            //  const int i = idx/iMAX+1;
            //  const int j = idx%iMAX+1;
            /* only if both adjacent cells are fluid cells */
            /*
            if ((flag_t(i,j) & C_F) && (flag_t(i,j+1) & C_F)) {
                const double duvdx = 
                        ((u_t(i,j)+u_t(i,j+1))     * (v_t(i,j)+v_t(i+1,j))+
                         gamma_t*fabs(u_t(i,j)+u_t(i,j+1))*(v_t(i,j)-v_t(i+1,j))-
                         (u_t(i-1,j)+u_t(i-1,j+1)) * (v_t(i-1,j)+v_t(i,j))-
                         gamma_t*fabs(u_t(i-1,j)   +  u_t(i-1,j+1))*(v_t(i-1,j)-v_t(i,j)))
                        /(4.0*delx_t);
                const double dv2dy = 
                        ((v_t(i,j)+v_t(i,j+1)) * (v_t(i,j)+v_t(i,j+1))+
                         gamma_t*fabs(v_t(i,j)+v_t(i,j+1)) * (v_t(i,j)-v_t(i,j+1))-
                         (v_t(i,j-1)+v_t(i,j)) * (v_t(i,j-1)+v_t(i,j))-
                         gamma_t*fabs(v_t(i,j-1)+v_t(i,j)) * (v_t(i,j-1)-v_t(i,j)))
                        /(4.0*dely_t);

                const double laplv = (v_t(i+1,j)-2.0*v_t(i,j)+v_t(i-1,j))/delx_t/delx_t+
                        (v_t(i,j+1)-2.0*v_t(i,j)+v_t(i,j-1))/dely_t/dely_t;

                g_t(i,j) = v_t(i,j)+step_delta_t*(laplv/reynolds_t-duvdx-dv2dy);
            } else {
                g_t(i,j) = v_t(i,j);
            }*/
    });

    /* f & g at external boundaries */
    Kokkos::parallel_for (
        jMAX,
	KOKKOS_LAMBDA(size_t j) {
        f_t(0,j)    = u_t(0,j);
        f_t(iMAX,j) = u_t(iMAX,j);
    });
    Kokkos::parallel_for (
        iMAX,
	KOKKOS_LAMBDA(size_t i) {
        g_t(i,0)    = v_t(i,0);
        g_t(i,jMAX) = v_t(i,jMAX);
    });
}

/* Calculate the right hand side of the pressure equation */
void Simulation::compute_rhs()
{
    auto g_t = g;
    auto f_t = f;
    auto rhs_t = rhs;
    auto flag_t = flag;
    auto dely_t = delx;
    auto delx_t = dely;
    auto step_delta_t = step_delta;
    Kokkos::parallel_for (
        iMAX*jMAX,
        KOKKOS_LAMBDA(size_t idx) {
        const int i = idx/iMAX+1;
        const int j = idx%iMAX+1;
        if (flag_t(i,j) & C_F) {
            /* only for fluid and non-surface cells */
            rhs_t(i,j) = ((f_t(i,j)-f_t(i-1,j))/delx_t +
                         (g_t(i,j)-g_t(i,j-1))/dely_t) / step_delta_t;
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

    /* Calculate sum of squares */
    /*
    Kokkos::parallel_reduce (
        iMAX*jMAX,
	KOKKOS_LAMBDA(size_t idx, double& inner_p0) {
            const int i = idx/iMAX + 1;
            const int j = idx%iMAX + 1;
            if (flag_t(i,j) & C_F) {
                //printf("p(%d,%d): %f\ninner_p0: %f\n", i, j, p(i,j), inner_p0);
                inner_p0 += p(i,j)*p(i,j);
            }
        }, p0
    );
    */
    double_host_view h_p = Kokkos::create_mirror_view (p); 
    char_host_view   h_f = Kokkos::create_mirror_view (flag); 
    double_host_view   h_rhs_t = Kokkos::create_mirror_view (rhs); 
    for (int i = 1; i <=iMAX; ++i) {
        for (int j = 1; j <= jMAX; ++j) {
            if (h_f(i,j) & C_F) {
                p0 += h_p(i,j)*h_p(i,j);
            }
        }
    }
    p0 = sqrt(p0/ifluid);
    printf("p0: %f\n", p0);
    if (p0 < 0.0001) {
        p0 = 1.0;
    }

    auto rhs_t = rhs;
    auto flag_t= flag;
    /* Red/Black SOR-iteration */
    for (iter = 0; iter < itermax; ++iter) {
        for (int rb = 0; rb <= 1; ++rb) {
        Kokkos::parallel_for (
            iMAX*jMAX,
            KOKKOS_LAMBDA(size_t idx) {
            const int i = idx/iMAX + 1;
            const int j = idx%jMAX + 1;
            /* This is tricky, try and think about what happens here */
            if ((i+j) % 2 != rb) {
                return;
            }
            if (flag(i,j) == (C_F | B_NSEW)) {
                /* five point star for interior fluid cells */
                const double beta_2 = -omega/(2.0*(rdx2+rdy2));
                double pijprev = p(i,j);
                p(i,j) = (1.0-omega)*p(i,j) -
                          beta_2*(
                                (p(i+1,j)+p(i-1,j))*rdx2
                              + (p(i,j+1)+p(i,j-1))*rdy2
                              -  rhs_t(i,j));
                double pij = p(i,j);
                double rhs_tij = rhs_t(i,j);
            } else if (flag_t(i,j) & C_F) {
                /* modified star near boundary */
                double pijprev = p(i,j);
                const double beta_mod = -omega/((eps_E+eps_W)*rdx2+(eps_N+eps_S)*rdy2);
                p(i,j) = (1.-omega)*p(i,j) -
                          beta_mod*(
                                (eps_E*p(i+1,j)+eps_W*p(i-1,j))*rdx2
                              + (eps_N*p(i,j+1)+eps_S*p(i,j-1))*rdy2
                              - rhs_t(i,j));
                double pij = p(i,j);
                double rhs_tij = rhs_t(i,j);

            }
            }); /* end of i */
        } /* end of rb */

        /* Partial computation of residual */
        double resl = 0.0;
        /*
        Kokkos::parallel_reduce (
            iMAX*jMAX,
            KOKKOS_LAMBDA(size_t idx, double& resl_lmb) {
            const int i = idx/iMAX+1;
            const int j = idx%iMAX+1;

            if (flag_t(i,j) & C_F) {
                // only fluid cells 
                const  double add = (eps_E*(p(i+1,j)-p(i,j)) -
                              eps_W*(p(i, j )-p(i-1,j))) * rdx2  +
                             (eps_N*(p(i,j+1)-p( i ,j)) -
                              eps_S*(p(i, j )-p(i,j-1))) * rdy2  -  rhs_t(i,j);
                resl_lmb += add*add;
            }
        }, resl);
        */
        for (int i = 0; i <= iMAX; ++i) {
            for (int j = 0; j <= jMAX; ++j) {
            if (flag(i,j) & C_F) {
                // only fluid cells 
                const  double add = (eps_E*(h_p(i+1,j)-h_p(i,j)) -
                              eps_W*(h_p(i, j )-h_p(i-1,j))) * rdx2  +
                             (eps_N*(h_p(i,j+1)-h_p( i ,j)) -
                              eps_S*(h_p(i, j )-h_p(i,j-1))) * rdy2  -  h_rhs_t(i,j);
                resl += add*add;
            }
 
            }
        }
        resl = sqrt((resl)/ifluid)/p0;
        // printf("%f, %f, %f\n", resl, ifluid, p0);
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
    auto u_t = u;
    auto v_t = v;
    auto f_t = f;
    auto g_t = g;
    auto p_t = p;
    auto step_delta_t = step_delta;
    auto delx_t = delx;
    auto dely_t = dely;
    Kokkos::parallel_for (
            iMAX*jMAX,
            KOKKOS_LAMBDA(size_t idx) {
            const int i = idx/iMAX+1;
            const int j = idx%jMAX+1;
            u_t(i,j) = f_t(i,j)-(p_t(i+1,j)-p_t(i,j))*step_delta_t/delx_t;
            v_t(i,j) = g_t(i,j)-(p_t(i,j+1)-p_t(i,j))*step_delta_t/dely_t;
    });
}


/* 
 * Set the timestep size so that we satisfy the Courant-Friedrichs-Lewy
 * conditions (ie no particle moves more than one cell width in one
 * timestep). Otherwise the simulation becomes unstable.
 */
void Simulation::set_timestep_interval()
{
    auto u_t = u;
    auto v_t = v;
    /* step_delta satisfying CFL conditions */
    if (tau >= 1.0e-10) { /* else no time stepsize control */
        double umax = 1.0e-10;
        double vmax = 1.0e-10;
        Kokkos::parallel_reduce (
            iMAX,
            KOKKOS_LAMBDA(size_t idx, double& max_ph) {
                const int i = idx/iMAX+1;
                const int j = idx%iMAX+1;
                /* only if both adjacent cells are fluid cells */
                max_ph = max(fabs(u_t(i,j)), max_ph);
        }, umax);

        Kokkos::parallel_reduce (
            iMAX*jMAX,
            KOKKOS_LAMBDA(size_t idx, double& max_ph) {
                const int i = idx/iMAX+1;
                const int j = idx%iMAX+1;
                /* only if both adjacent cells are fluid cells */
                max_ph = max(fabs(v_t(i,j)), max_ph);
        }, vmax);

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

/* Initialize the flag_t array, marking any obstacle cells and the edge cells
 * as boundaries. The cells adjacent to boundary cells have their relevant
 * flag_ts set too.
 */
void Simulation::init_flag()
{
    /* Mask of a circular obstacle */
    const double mx   = 20.0/41.0 * jMAX * dely;
    const double my   = mx;
    const double rad1 = 5.00/41.0 * jMAX * dely;

    printf("%g, %g, %g, %g, %g\n", delx, dely, mx, my, rad1);
    char_host_view h_flag = Kokkos::create_mirror_view (flag);
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
            const double x = (i - 0.5) * delx - mx;
            const double y = (j - 0.5) * dely - my;
            h_flag(i,j) = (x * x + y * y <= rad1 * rad1) ? C_B : C_F;
        }
    }
    
    printf("Checksum: %d\n", checksum);

    /* Mark the north & south boundary cells */
    for (int i = 0; i <= iMAX+1; ++i) {
        h_flag(i,0)      = C_B;
        h_flag(i,jMAX+1) = C_B;
    }
    /* Mark the east and west boundary cells */
    for (int j = 1; j <= jMAX; ++j) {
        h_flag(0,j)      = C_B;
        h_flag(iMAX+1,j) = C_B;
    }

    /* flag_ts for boundary cells */
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
    Kokkos::deep_copy (flag,h_flag); // Copy from host to device.
}

void Simulation::apply_boundary_conditions()
{
     auto u_t = u;
     auto v_t = v;
     Kokkos::parallel_for (
        jMAX+1,                     
        KOKKOS_LAMBDA(size_t j) {
        /* Fluid freely flows in from the west */
        u_t(0,j) = u_t(1,j);
        v_t(0,j) = v_t(1,j);

        /* Fluid freely flows out to the east */
        u_t(iMAX,j)   = u_t(iMAX-1,j);
        v_t(iMAX+1,j) = v_t(iMAX,j);
    });
    Kokkos::parallel_for (
        iMAX+1,                     
        KOKKOS_LAMBDA(size_t i) {
        /* The vertical velocity approaches 0 at the north and south
        * boundaries, but fluid flows freely in the horizontal direction */
        v_t(i,jMAX)   = 0.0;
        u_t(i,jMAX+1) = u_t(i,jMAX);

        v_t(i,0) = 0.0;
        u_t(i,0) = u_t(i,1);
    });

    /* Apply no-slip boundary conditions to cells that are adjacent to
    * internal obstacle cells. This forces the u and v velocity to
    * tend towards zero in these cells.
    */
    auto flag_t = flag;
    Kokkos::parallel_for (
        iMAX*jMAX,                     
        KOKKOS_LAMBDA(size_t idx) {
        const int i = idx/iMAX + 1;
        const int j = idx%iMAX + 1;    
            if (flag_t(i,j) & B_NSEW) {
                switch (flag_t(i,j)) {
                case B_N:
                    u_t(i,j)   = -u_t(i,j+1);
                    break;
                case B_E:
                    u_t(i,j)   = 0.0;
                    break;
                case B_NE:
                    u_t(i,j)   = 0.0;
                    break;
                case B_SE:
                    u_t(i,j)   = 0.0;
                    break;
                case B_NW:
                    u_t(i,j)   = -u_t(i,j+1);
                    break;
                case B_S:
                    u_t(i,j)   = -u_t(i,j-1);
                    break;
                case B_SW:
                    u_t(i,j)   = -u_t(i,j-1);
                    break;
                }
            }
    });
    Kokkos::parallel_for (
        iMAX*jMAX,                     
        KOKKOS_LAMBDA(size_t idx) {
        const int i = idx/iMAX + 1;
        const int j = idx%iMAX + 1;         
        if (flag_t(i+1,j) & B_NSEW) {
            switch (flag_t(i+1,j)) {
            case B_N:
                u_t(i,j) = -u_t(i,j+1);
                break;
            case B_W:
                u_t(i,j) = 0.0;
                break;
            case B_NE:
                u_t(i,j) = -u_t(i,j+1);
                break;
            case B_SW:
                u_t(i,j) = 0.0;
                break;
            case B_NW:
                u_t(i,j) = 0.0;
                break;
            case B_S:
                u_t(i,j) = -u_t(i,j-1);
                break;
            case B_SE:
                u_t(i,j) = -u_t(i,j-1);
                break;
            }
        }
    });

    Kokkos::parallel_for (
        iMAX*jMAX,                     
        KOKKOS_LAMBDA(size_t idx) {
        const int i = idx/iMAX + 1;
        const int j = idx%iMAX + 1;         

        if (flag_t(i,j) & B_NSEW) {
            switch (flag_t(i,j)) {
            case B_N:
                v_t(i,j)   = 0.0;
                break;
            case B_E:
                v_t(i,j)   = -v_t(i+1,j);
                break;
            case B_NE:
                v_t(i,j)   = 0.0;
                break;
            case B_SE:
                v_t(i,j)   = -v_t(i+1,j);
                break;
            case B_NW:
                v_t(i,j)   = 0.0;
                break;
            case B_W:
                v_t(i,j)   = -v_t(i-1,j);
                break;
            case B_SW:
                v_t(i,j)   = -v_t(i-1,j);
                break;
            }
        }
    });

     Kokkos::parallel_for (
        iMAX*jMAX,                     
        KOKKOS_LAMBDA(size_t idx) {
        const int i = idx/iMAX + 1;
        const int j = idx%iMAX + 1;         

        if (flag_t(i,j+1) & B_NSEW) {
            switch (flag_t(i,j+1)) {
            case B_E:
                v_t(i,j) = -v_t(i+1,j);
                break;
            case B_S:
                v_t(i,j) = 0.0;
                break;
            case B_NE:
                v_t(i,j) = -v_t(i+1,j);
                break;
            case B_SE:
                v_t(i,j) = 0.0;
                break;
            case B_SW:
                v_t(i,j) = 0.0;
                break;
            case B_W:
                v_t(i,j) = -v_t(i-1,j);
                break;
            case B_NW:
                v_t(i,j) = -v_t(i-1,j);
                break;
            }
        }
    });

    /* Finally, fix the horizontal velocity at the  western edge to have
    * a continual flow of fluid into the simulation.
    */
    auto vi_t = vi;
    auto ui_t = ui;
    Kokkos::parallel_for (
        jMAX,                     
        KOKKOS_LAMBDA(size_t idx) {
        const int j = idx + 1;
        v_t(0,0) = 2*vi_t-v_t(1,0);
        u_t(0,j) = ui_t;
        v_t(0,j) = 2*vi_t-v_t(1,j);
    });
}


/* Computation of stream function and vorticity */
void  Simulation::calc_psi_zeta(DoubleMatrix zeta) const
{
    /* Computation of the vorticity zeta at the upper right corner     */
    /* of cell (i,j) (only if the corner is surrounded by fluid cells) */
     auto v_t = v;
     auto u_t = u;
     auto delx_t = delx;
     auto dely_t = dely;
     auto flag_t = flag;
     Kokkos::parallel_for (
        iMAX*jMAX,                     
        KOKKOS_LAMBDA(size_t idx) {
        const int i = idx/iMAX + 1;
        const int j = idx%iMAX + 1;         
        if ((flag_t(i-1,j-1) & C_F) &&
               (flag_t(i-1, j ) & C_F) &&
               (flag_t(i-1,j+1) & C_F) &&
               (flag_t( i ,j-1) & C_F) &&
               (flag_t( i ,j+1) & C_F) &&
               (flag_t(i+1,j-1) & C_F) &&
               (flag_t(i+1, j ) & C_F) &&
               (flag_t(i+1,j+1) & C_F)) {
            zeta(i,j) = (u_t(i,j+1)-u_t(i,j))/dely_t-
                         (v_t(i+1,j)-v_t(i,j))/delx_t;
        } else {
            zeta(i,j) = 0.0;
        }
    });
}

// ====================== FRIEND FUNCTIONS ====================== //
