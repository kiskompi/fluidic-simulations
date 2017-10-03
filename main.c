#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>

#include "boundary.h"
#include "datadef.h"
#include "init.h"
#include "simulation.h"

/* Modified slightly by D. Orchard (2010) from the classic code from: 

    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
    Numerical Simulation in Fluid Dynamics,
    SIAM, 1998.

    http://people.sc.fsu.edu/~jburkardt/cpp_src/nast2d/nast2d.html

*/

int main(int argc, char *argv[])
{
    int verbose = 1;          /* Verbosity level */
    double xlength = 22.0;     /* Width of simulated domain */
    double ylength = 4.1;      /* Height of simulated domain */
    int imax = 660;           /* Number of cells horizontally */
    int jmax = 120;           /* Number of cells vertically */

    char *outname;
    int output = 0;
    int output_frequency = 0;

    double t_end = 40; //2.1       /* Simulation runtime */
    double del_t = 0.003;      /* Duration of each timestep */
    double tau = 0.5;          /* Safety factor for timestep control */

    int itermax = 100;        /* Maximum number of iterations in SOR */
    double eps = 0.001;        /* Stopping error threshold for SOR */
    double omega = 1.7;        /* Relaxation parameter for SOR */
    double gamma = 0.9;        /* Upwind differencing factor in PDE
                                 discretisation */

    double Re = 150.0;         /* Reynolds number */
    double ui = 1.0;           /* Initial X velocity */
    double vi = 0.0;           /* Initial Y velocity */

    double t, delx, dely;
    int  i, j, itersor = 0, ifluid = 0, ibound = 0;
    double res;
    double (*u)[jmax+2];
    double (*v)[jmax+2];
    double (*p)[jmax+2];
    double (*rhs)[jmax+2];
    double (*f)[jmax+2];
    double (*g)[jmax+2];
    char  (*flag)[jmax+2];
    int init_case, iters = 0;
    int show_help = 0, show_usage = 0, show_version = 0;

    if (argc > 1) {
      output = 1;
      outname = argv[1];
      output_frequency = 1;
    }

    if (argc > 2) {
      output_frequency = atoi(argv[2]);
    }

    delx = xlength/imax;
    dely = ylength/jmax;

    /* Allocate arrays */
    u    = (double (*)[jmax+2])malloc((imax+2) * (jmax+2) * sizeof(double));
    v    = (double (*)[jmax+2])malloc((imax+2) * (jmax+2) * sizeof(double));
    f    = (double (*)[jmax+2])malloc((imax+2) * (jmax+2) * sizeof(double));
    g    = (double (*)[jmax+2])malloc((imax+2) * (jmax+2) * sizeof(double));
    p    = (double (*)[jmax+2])malloc((imax+2) * (jmax+2) * sizeof(double));
    rhs  = (double (*)[jmax+2])malloc((imax+2) * (jmax+2) * sizeof(double)); 
    flag = (char (*)[jmax+2])malloc((imax+2) * (jmax+2) * sizeof(char));

    if (!u || !v || !f || !g || !p || !rhs || !flag) {
        fprintf(stderr, "Couldn't allocate memory for matrices.\n");
        return 1;
    }

    // Set up initial values
    for (i=0;i<=imax+1;i++) {
         for (j=0;j<=jmax+1;j++) {
             u[i][j] = ui;
             v[i][j] = vi;
             p[i][j] = 0.0;
         }
     }

    init_flag(jmax, flag, imax, delx, dely, &ibound);
    apply_boundary_conditions(jmax, u, v, flag, imax, ui, vi);
    
    // Main loop

    for (t = 0.0; t < t_end; t += del_t, iters++) {
        set_timestep_interval(jmax, &del_t, imax, delx, dely, u, v, Re, tau);

        ifluid = (imax * jmax) - ibound;

        compute_tentative_velocity(jmax, u, v, f, g, flag, imax,
            del_t, delx, dely, gamma, Re);

        compute_rhs(jmax, f, g, rhs, flag, imax, del_t, delx, dely);

        if (ifluid > 0) {
            itersor = poisson(jmax, p, rhs, flag, imax, delx, dely,
                        eps, itermax, omega, &res, ifluid);
        } else {
            itersor = 0;
        }

         printf("%d t:%g, del_t:%g, SOR iters:%3d, res:%e, bcells:%d\n",
                iters, t+del_t, del_t, itersor, res, ibound);

	
        update_velocity(jmax, u, v, f, g, p, flag, imax, del_t, delx, dely);

        apply_boundary_conditions(jmax, u, v, flag, imax, ui, vi);

    	if (output && (iters % output_frequency == 0)) {
    	  write_ppm(jmax, u, v, p, flag, imax, xlength, ylength, outname,
    		    iters, output_frequency);
    	}
    }

    free(u);
    free(v);
    free(f);
    free(g);
    free(p);
    free(rhs);
    free(flag);

    return 0;
}

// Used for comparing computations when debugging other implementations

unsigned int simplest_checksum_char(char** in, int imax, int jmax) {
  unsigned int checksum = 0;
  int i;
  int j;
  for (i=0; i<(imax+2); i++){
    for (j=0; j<(jmax+2); j++){
      checksum+=in[i][j]*(i);
    }
  }
  return checksum;
}

double simplest_checksum(double** in, int imax, int jmax) {
  double checksum = 0.0;
  int i;
  int j;
  for (i=0; i<(imax+2); i++){
    for (j=0; j<(jmax+2); j++){
      checksum+=in[i][j]*((double)(i*jmax)+j);
    }
  }
  return checksum;
}
