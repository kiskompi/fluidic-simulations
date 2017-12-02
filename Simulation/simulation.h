#ifndef SIMULATION_H
#define SIMULATION_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <fcntl.h>
#include <getopt.h>
#include <sys/stat.h>

#include <Kokkos_Core.hpp>
#include "datadef.h"

#define iMAX 660
#define jMAX 120

class Simulation {
public:
//    constexpr static const int    iMAX       = 660;     /* Number of cells horizontally */
//    constexpr static const int    jMAX       = 120;     /* Number of cells vertically   */

    typedef Kokkos::View<double**>    DoubleMatrix;
    typedef Kokkos::View<char**  >    CharMatrix;
    typedef CharMatrix::HostMirror    char_host_view;
    typedef DoubleMatrix::HostMirror  double_host_view;

    const double        xlength    = 22.0;    /* Width of simulated domain    */
    const double        ylength    = 4.10;     /* Height of simulated domain   */

    const double        time_end   = 40;      /* Simulation runtime 2.1             */
    const double        tau        = 0.5;     /* Safety factor for timestep control */

    const int           itermax    = 100;     /* Maximum number of iterations in SOR         */
    const double        eps        = 0.001;   /* Stopping error threshold for SOR            */
    const double        omega      = 1.7;     /* Relaxation parameter for SOR discretisation */
    const double        gamma      = 0.9;
    
    const double        reynolds   = 150.0;   /* Reynolds number    */
    const double        ui         = 1.0;     /* Initial X velocity */
    const double        vi         = 0.0;     /* Initial Y velocity */
    
    const double        delx       = xlength/iMAX;
    const double        dely       = ylength/jMAX;
    double              ifluid     = 0.0;
    
private:    
    double res        = 0.0;
    int    ibound     = 0;
    double step_delta = 0.003;   /* Duration of each timestep */
 
    DoubleMatrix u    = DoubleMatrix  ("u",    iMAX+2, jMAX+2); 
    DoubleMatrix v    = DoubleMatrix  ("v",    iMAX+2, jMAX+2);
    DoubleMatrix p    = DoubleMatrix  ("p",    iMAX+2, jMAX+2);
    DoubleMatrix f    = DoubleMatrix  ("f",    iMAX+2, jMAX+2);
    DoubleMatrix g    = DoubleMatrix  ("g",    iMAX+2, jMAX+2);
    DoubleMatrix rhs  = DoubleMatrix  ("rhs",  iMAX+2, jMAX+2); 
    CharMatrix   flag = CharMatrix    ("flag", iMAX+2, jMAX+2);

    bool is_surrounded(const size_t i, const size_t j) const{
        return (flag(i-1,j-1) & C_F) &&
               (flag(i-1, j ) & C_F) &&
               (flag(i-1,j+1) & C_F) &&
               (flag( i ,j-1) & C_F) &&
               (flag( i ,j+1) & C_F) &&
               (flag(i+1,j-1) & C_F) &&
               (flag(i+1, j ) & C_F) &&
               (flag(i+1,j+1) & C_F);
    }

public:
    Simulation();
    void compute_tentative_velocity();
    void compute_rhs();
    void update_velocity();
    void set_timestep_interval();
    void apply_boundary_conditions();
    int  poisson();
    void calc_psi_zeta(DoubleMatrix zeta) const;
    void init_flag();       // Initialize the flag array

    inline const CharMatrix get_flag() const {
        return flag;
    }
    
    inline const DoubleMatrix get_p() const {
        return p;
    }

    inline const double& get_res() const {
        return res;
    }
    
    inline const int& get_ibound() const {
        return ibound;
    }
    
    inline const double& get_step_delta() const {
        return step_delta;
    }

    friend void write_ppm(
        const Simulation&   sim,
        char*         outname, 
        int           iters, 
        int           freq);
    friend unsigned int simplest_checksum_char(Simulation::CharMatrix in);
    friend double simplest_checksum(Simulation::DoubleMatrix in);
};


#endif
