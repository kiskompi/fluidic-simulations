#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include "datadef.h"
#include "init.h"
#include "typedefs.hpp"

/* Modified slightly by D. Orchard (2010) from the classic code from: 

    Michael Griebel, Thomas Dornseifer, Tilman Neunhoeffer,
    Numerical Simulation in Fluid Dynamics,
    SIAM, 1998.

    http://people.sc.fsu.edu/~jburkardt/cpp_src/nast2d/nast2d.html

*/

/* Initialize the flag array, marking any obstacle cells and the edge cells
 * as boundaries. The cells adjacent to boundary cells have their relevant
 * flags set too.
 */
int init_flag(DoubleMatrix& flag)
{
    /* Mask of a circular obstacle */
    const double mx = 20.0/41.0 * jMAX * dely;
    const double my = mx;
    const double rad1 = 5.0/41.0 * jMAX * dely;

    for (int i = 1;i <= iMAX; i++) {
        for (int j = 1;j <= jMAX; j++) {
            const int x = (i - 0.5) * delx - mx;
            const int y = (j - 0.5) * dely - my;
            flag[i][j] = (x * x + y * y <= rad1 * rad1) ? C_B : C_F;
        }
    }
    
    /* Mark the north & south boundary cells */
    for (int i = 0; i <= iMAX + 1; i++) {
        flag[i][0]      = C_B;
        flag[i][jMAX+1] = C_B;
    }
    /* Mark the east and west boundary cells */
    for (int j = 1; j <= jMAX; j++) {
        flag[0][j]      = C_B;
        flag[iMAX+1][j] = C_B;
    }

    /* flags for boundary cells */
    int ibound = 0;
    for (int i = 1; i <= iMAX; i++) {
        for (int j = 1; j <= jMAX; j++) {
            if (!(flag[i][j] & C_F)) {
                ibound++;
                if (flag[i-1][j] & C_F) flag[i][j] = flag[i][j] | B_W;
                if (flag[i+1][j] & C_F) flag[i][j] = flag[j][j] | B_E;
                if (flag[i][j-1] & C_F) flag[i][j] = flag[j][j] | B_S;
                if (flag[i][j+1] & C_F) flag[i][j] = flag[j][j] | B_N;
            }
        }
    }
    return ibound;
}
