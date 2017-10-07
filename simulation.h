void compute_tentative_velocity(int jMAX, 
    DoubleMatrix& u, 
    DoubleMatrix& v, 
    DoubleMatrix& f, 
    DoubleMatrix& g,
    CharMatrix&   flag, 
    int      iMAX,  
    double   del_t, 
    double   delx, 
    double   dely,
    double   gamma, 
    double   Re);

void compute_rhs(int jMAX, 
    DoubleMatrix& f, 
    DoubleMatrix& g, 
    DoubleMatrix& rhs, 
    CharMatrix&   flag, 
    int      iMAX,
    double   del_t, 
    double   delx, 
    double   dely);

int poisson(int jMAX, 
    DoubleMatrix& p, 
    DoubleMatrix& rhs, 
    CharMatrix&   flag, 
    int      iMAX, 
    double   delx, 
    double   dely, 
    double   eps, 
    int      itermax, 
    double   omega,
    double&  res, 
    int      ifull);

void update_velocity(
    int jMAX, 
    DoubleMatrix& u, 
    DoubleMatrix& v,
    DoubleMatrix& f, 
    DoubleMatrix& g, 
    DoubleMatrix& p,
    CharMatrix&   flag, 
    int      iMAX,  
    double   del_t, 
    double   delx, 
    double   dely);


void set_timestep_interval(
    int      jMAX, 
    double*  del_t, 
    int      iMAX,  
    double   delx,
    double   dely, 
    DoubleMatrix& u, 
    DoubleMatrix& v, 
    double   Re, 
    double   tau);
