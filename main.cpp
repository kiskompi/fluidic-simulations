#include <errno.h>
#include "output.c"
#include <iostream>
#include "Simulation/simulation.h"
#include <time.h>

/* Modified slightly by D. Orchard (2010) from the classic code from: 

    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
    Numerical Simulation in Fluid Dynamics,
    SIAM, 1998.

    http://people.sc.fsu.edu/~jburkardt/cpp_src/nast2d/nast2d.html

*/


int main(int argc, char *argv[])
{
    clock_t start, end;
    double cpu_time_used;
     
    start = clock();
    Kokkos::initialize();

    Simulation sim = Simulation();
    sim.init_flag();

    sim.apply_boundary_conditions();
    
    // Main loop
    int    iters  = 0;
    int    itersor= 0;
    
    for (double time_v = 0.0; time_v < sim.time_end; time_v += sim.get_step_delta(), ++iters) {
        sim.set_timestep_interval();
        sim.compute_tentative_velocity();
        sim.compute_rhs();

        itersor = (sim.ifluid > 0) ? sim.poisson() : 0;

        printf("%d time_v:%g, step_delta:%g, SOR iters:%3d, res:%e, bcells:%d\n",
               iters, 
               time_v+sim.get_step_delta(), 
               sim.get_step_delta(), 
               itersor, 
               sim.get_res(), 
               sim.get_ibound());

	
        sim.update_velocity();
        sim.apply_boundary_conditions();
        if (iters == 100) break;
    	// if (output && (iters % output_frequency == 0)) {
    	//   write_ppm(sim, outname, iters, output_frequency);
    	// }
    }
    Kokkos::finalize();
    printf("%.2f\n", (double)(time(NULL) - start));
    end = clock();
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    std::cout<<"\nChrono Runtime"<<cpu_time_used<<"\n";
    return 0;
}

