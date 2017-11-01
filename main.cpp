#include <errno.h>
#include "output.c"
#include <iostream>
#include "Simulation/simulation.h"
/* Modified slightly by D. Orchard (2010) from the classic code from: 

    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
    Numerical Simulation in Fluid Dynamics,
    SIAM, 1998.

    http://people.sc.fsu.edu/~jburkardt/cpp_src/nast2d/nast2d.html

*/


int main(int argc, char *argv[])
{   
    Kokkos::initialize();
    int    output           = 0;
    int    output_frequency = 0;
    char*  outname          = nullptr;

    Simulation sim = Simulation();
    sim.init_flag();
    if (argc > 1) {
      output = 1;
      outname = argv[1];
      output_frequency = 1;
    }

    if (argc > 2) {
      output_frequency = atoi(argv[2]);
    }

    sim.apply_boundary_conditions();
    
    // Main loop
    int    iters  = 0;
    int    itersor= 0;
    
    for (double time_v = 0.0; time_v < sim.time_end; time_v += sim.get_step_delta(), iters++) {
        sim.set_timestep_interval();
        sim.compute_tentative_velocity();
        sim.compute_rhs();

        itersor = (sim.ifluid > 0) ? sim.poisson() : 0;

        printf("%d time_v:%g, step_delta:%g, SOR iters:%3d, res:%e, bcells:%d\n",
               iters, time_v+sim.get_step_delta(), sim.get_step_delta(), itersor, sim.get_res(), sim.get_ibound());

	
        sim.update_velocity();
        sim.apply_boundary_conditions();

    	// if (output && (iters % output_frequency == 0)) {
    	//   write_ppm(sim, outname, iters, output_frequency);
    	// }
    }
    Kokkos::finalize();
    return 0;
}

// Used for comparing computations when debugging other implementations

unsigned int simplest_checksum_char(Simulation::CharMatrix& in) {
  unsigned int checksum = 0;
  int i;
  int j;
  for (i=0; i<(Simulation::iMAX+2); i++){
    for (j=0; j<(Simulation::jMAX+2); j++){
      checksum+=in[i][j]*(i);
    }
  }
  return checksum;
}

double simplest_checksum(Simulation::DoubleMatrix& in) {
  double checksum = 0.0;
  int i;
  int j;
  for (i=0; i<(Simulation::iMAX+2); i++){
    for (j=0; j<(Simulation::jMAX+2); j++){
      checksum+=in[i][j]*((double)(i*Simulation::jMAX)+j);
    }
  }
  return checksum;
}
