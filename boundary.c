#include <stdio.h>
#include <string.h>
#include "datadef.h"
/*
 * SSH: akomporday@haswell
 */
/* Modified slightly by D. Orchard (2010) from the classic code from: 

    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
    Numerical Simulation in Fluid Dynamics,
    SIAM, 1998.

    http://people.sc.fsu.edu/~jburkardt/cpp_src/nast2d/nast2d.html

*/

/* Given the boundary conditions defined by the flag matrix, update
 * the u and v velocities. Also enforce the boundary conditions at the
 * edges of the matrix.
 */
void apply_boundary_conditions(
                               DoubleMatrix& u,
                               DoubleMatrix& v,
                               CharMatrix&   flag)
{
    for (int j=0; j<=jMAX+1; j++) {
        /* Fluid freely flows in from the west */
        u[0][j] = u[1][j];
        v[0][j] = v[1][j];

        /* Fluid freely flows out to the east */
        u[iMAX][j]   = u[iMAX-1][j];
        v[iMAX+1][j] = v[iMAX][j];
    }

    for (int i=0; i<=iMAX+1; i++) {
        /* The vertical velocity approaches 0 at the north and south
         * boundaries, but fluid flows freely in the horizontal direction */
        v[i][jMAX]   = 0.0;
        u[i][jMAX+1] = u[i][jMAX];

        v[i][0] = 0.0;
        u[i][0] = u[i][1];
    }

    /* Apply no-slip boundary conditions to cells that are adjacent to
     * internal obstacle cells. This forces the u and v velocity to
     * tend towards zero in these cells.
     */

   for (int i = 1; i <= iMAX; i++) {
        for (int j = 1; j <= jMAX; j++) {
            if (flag[i][j] & Cell::B_NSEW) {
                switch (flag[i][j]) {
                    case Cell::B_N: 
                        u[i][j]   = -u[i][j+1];
                        break;
                    case Cell::B_E: 
                        u[i][j]   = 0.0;
                        break;
                    case Cell::B_NE:
                        u[i][j]   = 0.0;
                        break;
                    case Cell::B_SE:
                        u[i][j]   = 0.0;                      
                        break;
                    case Cell::B_NW:
                        u[i][j]   = -u[i][j+1];
                        break;
                    case Cell::B_S:
                        u[i][j]   = -u[i][j-1];
                        break;
                    case Cell::B_SW:
                        u[i][j]   = -u[i][j-1];
                        break;
                }
            }
	}
    } 
    for (int i=0; i<=(iMAX-1); i++) {
        for (int j=1; j<=jMAX; j++) {
            if (flag[i+1][j] & Cell::B_NSEW) {
                switch (flag[i+1][j]) {
                    case Cell::B_N: 
                        u[i][j] = -u[i][j+1];
                        break;
                    case Cell::B_W: 
                        u[i][j] = 0.0;
                        break;
                    case Cell::B_NE:
                        u[i][j] = -u[i][j+1];
                        break;
                    case Cell::B_SW:
                        u[i][j] = 0.0;
                        break;
                    case Cell::B_NW:
                        u[i][j] = 0.0;
                        break;
                    case Cell::B_S:
                        u[i][j] = -u[i][j-1];
                        break;
                    case Cell::B_SE:
                        u[i][j] = -u[i][j-1];
                        break;
                }
            }
	}
    } 


    for (int i=1; i<=iMAX; i++) {
        for (int j=1; j<=jMAX; j++) {
            if (flag[i][j] & Cell::B_NSEW) {
                switch (flag[i][j]) {
                    case Cell::B_N: 
                        v[i][j]   = 0.0;
                        break;
                    case Cell::B_E: 
                        v[i][j]   = -v[i+1][j];
                        break;
                    case Cell::B_NE:
                        v[i][j]   = 0.0;
                        break;
                    case Cell::B_SE:
                        v[i][j]   = -v[i+1][j];
                        break;
                    case Cell::B_NW:
                        v[i][j]   = 0.0;
                        break;
                    case Cell::B_W: 
                        v[i][j]   = -v[i-1][j];
                        break;
                    case Cell::B_SW:
                        v[i][j]   = -v[i-1][j];
                        break;
                }
            }
	 }
      } 

    for (int i=1; i<=iMAX; i++) {
      for (int j=0; j<=(jMAX-1); j++) {
            if (flag[i][j+1] & Cell::B_NSEW) {
                switch (flag[i][j+1]) {
                    case Cell::B_E: 
                        v[i][j] = -v[i+1][j];
                        break;
                    case Cell::B_S:
                        v[i][j] = 0.0;
                        break;
                    case Cell::B_NE:
                        v[i][j] = -v[i+1][j];
                        break;
                    case Cell::B_SE:
                        v[i][j] = 0.0;
                        break;
                    case Cell::B_SW:
                        v[i][j] = 0.0;
			break;
                    case Cell::B_W: 
                        v[i][j] = -v[i-1][j];
                        break;
                    case Cell::B_NW:
                        v[i][j] = -v[i-1][j];
                        break;
                }
            }
	}
     } 

    /* Finally, fix the horizontal velocity at the  western edge to have
     * a continual flow of fluid into the simulation.
     */
    v[0][0] = 2*vi-v[1][0];
    for (int j=1;j<=jMAX;j++) {
        u[0][j] = ui;
        v[0][j] = 2*vi-v[1][j];
    }
}
