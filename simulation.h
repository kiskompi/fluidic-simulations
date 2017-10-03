void compute_tentative_velocity(int jmax, double (*u)[jmax+2], double (*v)[jmax+2], double (*f)[jmax+2], double (*g)[jmax+2],
    char (*flag)[jmax+2], int imax, double del_t, double delx, double dely,
    double gamma, double Re);

void compute_rhs(int jmax, double (*f)[jmax+2], double (*g)[jmax+2], double (*rhs)[jmax+2], char (*flag)[jmax+2], int imax,
    double del_t, double delx, double dely);

int poisson(int jmax, double (*p)[jmax+2], double (*rhs)[jmax+2], char (*flag)[jmax+2], int imax, 
    double delx, double dely, double eps, int itermax, double omega,
    double *res, int ifull);

void update_velocity(int jmax, double (*u)[jmax+2], double (*v)[jmax+2], double (*f)[jmax+2], double (*g)[jmax+2], double (*p)[jmax+2],
    char (*flag)[jmax+2], int imax, double del_t, double delx, double dely);

void set_timestep_interval(int jmax, double *del_t, int imax, double delx,
    double dely, double (*u)[jmax+2], double (*v)[jmax+2], double Re, double tau);
