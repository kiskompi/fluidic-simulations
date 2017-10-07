#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <algorithm>
#include <sys/stat.h>
#include "datadef.h"

/* Modified slightly by D. Orchard (2010) from the classic code from: 

    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
    Numerical Simulation in Fluid Dynamics,
    SIAM, 1998.

    http://people.sc.fsu.edu/~jburkardt/cpp_src/nast2d/nast2d.html

*/

/* Computation of stream function and vorticity */
void inline calc_psi_zeta(DoubleMatrix& u,
                          DoubleMatrix& v,
                          DoubleMatrix& psi,
                          DoubleMatrix& zeta,
                          CharMatrix&   flag)
{
    /* Computation of the vorticity zeta at the upper right corner     */
    /* of cell (i,j) (only if the corner is surrounded by fluid cells) */
    for (size_t i=1;i<=iMAX-1;i++) {
        for (size_t j=1;j<=jMAX-1;j++) {
            if (is_surrounded(flag[i, j])) {
                zeta[i][j] = (u[i][j+1]-u[i][j])/dely-
                             (v[i+1][j]-v[i][j])/delx;
            } else {
                zeta[i][j] = 0.0;
            }
        }
    }
}

void write_ppm(DoubleMatrix& u,
               DoubleMatrix& v,
               DoubleMatrix& p,
               CharMatrix&   flag,
               char*  outname, 
               int    iters, 
               int    freq) 
{
    double    zmax = -1e10, zmin = 1e10;
    double    pmax = -1e10, pmin = 1e10;
    double    vmax = -1e10, vmin = 1e10;
    double    umax = -1e10, umin = 1e10;
    char      outfile[64];
    char      outmode  = 0;    
    DoubleMatrix&   psi  = new double[jMAX+2];  
    DoubleMatrix&   zeta = new double[jMAX+2];     

    sprintf(outfile, "%s/%06d.ppm", outname, (iters/freq));

    __mode_t process_mask = umask(0);
    mkdir(outname, S_IRWXU | S_IRWXG | S_IRWXO);
    umask(process_mask);
    FILE *fout = fopen(outfile, "wb");
    
    if (!fout) {
      fprintf (stderr, "Could not open '%s'\n", outfile);
      return;
    }

    DoubleMatrix zeta = DoubleMatrix();

    calc_psi_zeta(u, v, psi, zeta, flag);

    fprintf(fout, "P6 %d %d 255\n", iMAX, jMAX);    

    for (size_t j = 1; j < jMAX+1 ; j++) {
        for (size_t i = 1; i < iMAX+1 ; i++) {
            int r, g, b;
            if (!(flag[i][j] & C_F)) {
                r = 0; b = 0; g = 255;
            } else {
	            if (outmode == 0) {
                    double z = (i < iMAX && j < jMAX) ? zeta[i][j] : 0.0;
                    r = g = b = pow(fabs(z/12.6),.4) * 255;
                } else if (outmode == 1) {
                    double p = (i < iMAX && j < jMAX) ? psi[i][j] : 0.0;
                    r = g = b = (p + 3.0) / 7.5 * 255; 
                } else if (outmode == 2) {
	    	        r = g = b = (p[i][j]-pmin) / (pmax-pmin) * 255;
	            }
            }
            fprintf(fout, "%c%c%c", r, g, b);
        }
    } 

    fclose(fout);    
    if (outmode==1) {
        free(psi);
    }
    
    if (outmode==0) {
        free(zeta);
    } 
}

static inline bool is_surrounded(char* flag){
    return (flag[0][0] & C_F) && 
           (flag[1][0] & C_F) &&
           (flag[0][1] & C_F) && 
           (flag[1][1] & C_F);
}


