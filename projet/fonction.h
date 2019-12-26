

#define u_sys 1
#define l_sys 0


void init_mat(double **mat, double n, double v);
void init_mat_test(double **mat, int n, double a, double b);
void init_mat_test2(double **mat, int n);
void init_v(double *q, int n, double v);
void affiche_mat(double **mat,int n);
void affiche_v(double *v,int n);
double norm(double *q, int n);
double* dot(double *q, double **mat, int n);
double* prod_v_scal(double *q, int n,double scal );
double* dev_v_scal(double *q, int n,double scal );
double prod_scal(double* q1, double* q2, int n);
double* somme_vec(double* q1, double* q2, int n );
double* sub_vec(double* q1, double* q2, int n );
double* extract_v(double** mat, int i, int n);
void insert(double** mat, double* q, int i, int n);
void arnoldi(double **mat, double *v,int k, int n, double **V, double **H);
double absolute_value(double value);
void seq_LU_fact(double **A , double ** LU , size_t ** permutation , int size );
void res_sys_tri (double** A , double* X , double* b  ,int size , int stat);
void invers_power(double **A , int size , double E );











