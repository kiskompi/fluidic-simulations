#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <errno.h>

#include "boundary.h"
#include "datadef.h"
#include "init.h"
#include "simulation.h"
#include "typedefs.hpp"
/* Modified slightly by D. Orchard (2010) from the classic code from: 

    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
    Numerical Simulation in Fluid Dynamics,
    SIAM, 1998.

    http://people.sc.fsu.edu/~jburkardt/cpp_src/nast2d/nast2d.html

*/


int main(int argc, char *argv[])
{
   
    int    output           = 0;
    int    output_frequency = 0;
    char*  outname          = nullptr;
    double res              = 0;

    DoubleMatrix u    = DoubleMatrix(); 
    DoubleMatrix v    = DoubleMatrix();
    DoubleMatrix p    = DoubleMatrix();
    DoubleMatrix rhs  = DoubleMatrix(); 
    DoubleMatrix f    = DoubleMatrix();
    DoubleMatrix g    = DoubleMatrix();
    CharMatrix   flag = CharMatrix();

    if (argc > 1) {
      output = 1;
      outname = argv[1];
      output_frequency = 1;
    }

    if (argc > 2) {
      output_frequency = atoi(argv[2]);
    }
    int ibound = init_flag(flag);
    apply_boundary_conditions(u, v, flag);
    
    // Main loop

    int    iters  = 0;
    int    itersor= 0;
    double time_v = 0.0;
    for (; time_v < t_end; time_v += step_delta, iters++) {
        set_timestep_interval(u, v);

        const double ifluid = (iMAX * jMAX) - ibound;

        compute_tentative_velocity(u, v, f, g, flag);

        compute_rhs(f, g, rhs, flag);

        itersor = (ifluid > 0) ? poisson(p, rhs, flag, res, ifluid) : 0;

        printf("%d time_v:%g, step_delta:%g, SOR iters:%3d, res:%e, bcells:%d\n",
               iters, time_v+step_delta, step_delta, itersor, res, ibound);

	
        update_velocity(u, v, f, g, p, flag);

        apply_boundary_conditions(u, v, flag);

    	if (output && (iters % output_frequency == 0)) {
    	  write_ppm(u, v, p, flag, outname, iters, output_frequency);
    	}
    }


    return 0;
}

// Used for comparing computations when debugging other implementations

unsigned int simplest_checksum_char(CharMatrix& in) {
  unsigned int checksum = 0;
  int i;
  int j;
  for (i=0; i<(iMAX+2); i++){
    for (j=0; j<(jMAX+2); j++){
      checksum+=in[i][j]*(i);
    }
  }
  return checksum;
}

double simplest_checksum(DoubleMatrix& in) {
  double checksum = 0.0;
  int i;
  int j;
  for (i=0; i<(iMAX+2); i++){
    for (j=0; j<(jMAX+2); j++){
      checksum+=in[i][j]*((double)(i*jMAX)+j);
    }
  }
  return checksum;
}
