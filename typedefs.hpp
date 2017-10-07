#ifndef TYPEDEFS_H
#define TYPEDEFS_H
#include "Matrix.hpp"
const int    iMAX       = 660;           /* Number of cells horizontally */
const int    jMAX       = 120;           /* Number of cells vertically */
const double xlength    = 22.0;    /* Width of simulated domain  */
const double ylength    = 4.1;     /* Height of simulated domain */

const double t_end      = 40;      /* Simulation runtime 2.1             */
const double step_delta = 0.003;   /* Duration of each timestep          */
const double tau        = 0.5;     /* Safety factor for timestep control */

const int    itermax    = 100;     /* Maximum number of iterations in SOR         */
const double eps        = 0.001;   /* Stopping error threshold for SOR            */
const double omega      = 1.7;     /* Relaxation parameter for SOR discretisation */
const double gamma      = 0.9;

const double Re         = 150.0;   /* Reynolds number    */
const double ui         = 1.0;     /* Initial X velocity */
const double vi         = 0.0;     /* Initial Y velocity */

const double delx = xlength/iMAX;
const double dely = ylength/jMAX;


typedef Matrix<double, iMAX + 2, jMAX + 2>  DoubleMatrix;
typedef Matrix<char,   iMAX + 2, jMAX + 2>  CharMatrix;
#endif