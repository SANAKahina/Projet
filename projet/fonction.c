#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fonction.h"



void init_mat(double **mat, double n, double v)
{ 
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            mat[i][j]=v;
}

void init_mat_test(double **mat, int n, double a, double b)
{ 
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
        	if (i == j)
        	{
            	mat[i][j]=a;
            	if (i<n-1)
            	{
            		mat[i][j+1]=b;
            	}
            	if (j>0)
            	{
            		mat[i][j-1]=b;
            	}
        	}
}

void init_mat_test2(double **mat, int n)
{ 

    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
          if (i == j)
          {
              mat[i][j]=i+1;
              if (i<n-1)
              {
                mat[i][j+1]=-0.1;
              }
              if (j>0)
              {
                mat[i][j-1]=0.1;
              }
          }
}

void init_v(double *q, int n, double v)
{ 
        for(int j=0;j<n;j++)
            q[j]=v;
}

void affiche_mat(double **mat,int n)
{ 
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++)
            printf("%f ",mat[i][j]);
        printf("\n");
    }
}
void affiche_v(double *v,int n)
{ 
    for(int i=0;i<n;i++){
            printf("%f ",v[i]);
        printf("\n");
    }
}

double norm(double *q, int n)
{
  double tmp=0;
 
    for (int i = 0 ; i < n ; i++)
    {
      tmp = q[i] * q[i] + tmp;
    }
  return sqrt(tmp);
  
}

double* dot(double *q, double **mat, int n)
{
  double* res;
  res =malloc(sizeof(double*)*n);
  init_v(res, n, 0);

  for(int i=0;i<n;i++)
	  {
	  for( int j=0; j<n; j++)
	    {
	    res[i] = res[i] +mat[i][j]*q[j];

	    }
	  }
  return res;
}

double* prod_v_scal(double *q, int n,double scal )
{ 
      for(int j=0;j<n;j++)
      q[j]=q[j]*scal;
      return q;
	 
}

double* dev_v_scal(double *q, int n,double scal )
{ 
  prod_v_scal(q, n, 1/scal);	 
}


double prod_scal(double* q1, double* q2, int n)
{
  double tmp = 0;
  
  for( int j=0; j<n; j++)
	    {
	    tmp = tmp+ q1[j]*q2[j];
	    }
          return tmp;
}

double* somme_vec(double* q1, double* q2, int n )
{
  double* vt= NULL;
  vt=malloc(sizeof(double*)*n);
  init_v(vt, n, 0);

      for(int j=0;j<n;j++)
      vt[j]=q1[j]+q2[j];

  return vt;
	 
}

double* sub_vec(double* q1, double* q2, int n )
{
  double* vt= NULL;
  vt=malloc(sizeof(double*)*n);
  init_v(vt, n, 0);

      for(int j=0;j<n;j++)
      vt[j]=q1[j]-q2[j];

  return vt;
	 
}

double* extract_v(double** mat, int i, int n)
{
  double* vt= NULL;
  vt=malloc(sizeof(double*)*n);
  init_v(vt, n, 0);

      for(int j=0;j<n;j++)
      vt[j]=mat[i][j];
      return vt;
  
}

void insert(double** mat, double* q, int i, int n)
{
        for(int j=0;j<n;j++)
        mat[i][j]= q[j];
}


void arnoldi(double **mat, double *v, 
       int k, int n, double **V, double **H){
    double* vt= NULL;
    vt=malloc(sizeof(double*)*n);
    init_v(vt, n, 0);

    insert(V, dev_v_scal(v, n, norm(v,n)), 0,n);
    // affiche_v( extract_v(V, 0, n),n);
     
    for (unsigned int m=0; m<k; m++) {
      vt = dot(extract_v(V,m,n) ,mat ,n);
        for (unsigned int j=0; j<m+1; j++) {
    H[j][m] = prod_scal(vt,extract_v(V,j,n) , n);
    vt = sub_vec(vt,prod_v_scal(extract_v(V,j,n), n, H[j][m]), n );
        }
        H[m+1][m] = norm(vt, n);
        if (m != k-1)
    insert(V, dev_v_scal(vt, n,H[m+1][m] )  ,m+1, n);
    }
    
}



double absolute_value(double value) {
  return (value < 0) ? -value : value;
}

void seq_LU_fact(double **A , double ** LU ,size_t** permutation , int size ) {

  //A matrice a factorisé
  //LU la factorisation LU de A
  //b la taille AOO
  //size taille de la matrice carrée
  *permutation = (size_t*)calloc(size, sizeof(size_t));

  //*permutation = (double*)calloc(size, sizeof(double));
  for (size_t i = 0; i < size; ++i) {
    (*permutation)[i] = i;
  }

  double ** B ;
  double * temp;
  double pivot;
  size_t highest_pivot;
  size_t temp_size;

  for (size_t j = 0; j < size-1 ; j++) {

    highest_pivot = j;

    if (j == 0 ) {
      B = A;
    } else {
      B = LU;
    }

    for (size_t i = j + 1; i < size; ++i) {
      if (absolute_value(B[highest_pivot][j]) < absolute_value(B[i][j]))   highest_pivot = i;
    }

    if (highest_pivot != j) {
      temp_size = (*permutation)[highest_pivot];
      (*permutation)[highest_pivot] = (*permutation)[j];
      (*permutation)[j] = temp_size;
      temp = B[highest_pivot];
      B[highest_pivot] = B[j];
      B[j] = temp;
    }

    if (j == 0) {
      //intialisation de la 1er ligne de LU
      for (size_t i = 0; i <  size; i++) {
        LU[0][i] = A[0][i];
      }
    }

    if (B[j][j] == 0 && j != 0)   continue;

    for (size_t i = j + 1; i < size; i++) {
      if (B[j][j] != 0)   pivot = B[i][j] / B[j][j];
      else                pivot = 0;
    //  printf("%f ", pivot);
      LU[i][j] = pivot;
      if (pivot != 0 || j == 0) {
        for (size_t k = j + 1; k < size; k++) {
          LU[i][k] = B[i][k] - pivot*B[j][k];
        }
      }
    }

    if (j == 0 && highest_pivot != j) {
      temp = B[highest_pivot];
      B[highest_pivot] = B[j];
      B[j] = temp;
    }
   // printf("\n");
  }
}

void res_sys_tri (double** A , double* X , double* b  ,int size , int stat) {

  // A : L00 ou U00 AX=b
  // X : sert a stocké les résultats soit par ligne (L10)ou par colone(U01)
  // b : la partie de la matrice A(louwla) qui represente le second memebre (ligne cas U00 ou colone cas du L00)
  // n : size du vecteur resultat
  // state : 0 -> L , 1 -> U
  // index : indice ou stcoker le resultat dans LU et ou prendre le b dans A

  if (stat) {
    size_t index_i;
    //le cas ou un a une forme A = U * X
    if (A[size-1][size-1] != 0.0)  X[size-1] = b[size-1] / A[size-1][size-1];
    else                                     X[size-1] = 0;
    // printf("\n\n\nres[%d] = %f\n", size-1, X[size-1]);
    for (size_t i = 1; i < size; i++) {
      index_i = size -1 - i;
      X[index_i] = b[index_i];
      // printf("res[%ld] = %f", size -1 - i, X[size -1 - i]);
      for (size_t j = size - 1; j > index_i; j--) {
        X[index_i] -= A[index_i][j] * X[j];
        // printf("- %f * %f ", A[size -1 - i][j], X[j]);
      }
      if (A[index_i][index_i] != 0.0)          X[index_i] /= A[index_i][index_i];
      else                       X[index_i] = 0;
      // printf("\nX[%ld] /= %f -> X[%ld] = %f\n", size -1 - i, A[size -1 - i][size -1 - i], size -1 - i, X[size -1 - i]);
    }
  } else {
    //le cas ou un a une forme A = L * X
    X[0] = b[0];
    for (size_t i = 1; i < size; i++) {
      X[i] = b[i];
      for (size_t j = 0; j < i; j++) {
        X[i] -= A[i][j] * X[j];
      }
    }
  }
}

void invers_power(double **A , int size , double E ){

  double  **LU;

  LU=malloc(sizeof(double*)*size);

     for(int i=0;i<size;i++)
        LU[i]=malloc(sizeof(double)*size);

 // mat_init( &LU , size ,size);
  init_mat(LU,size, 0);

  size_t *perm;
  double old_norme ;

  double *X = (double*)malloc(size*sizeof(double));
  double *X_start = (double*)malloc(size*sizeof(double));
  double *X_end = (double*)malloc(size*sizeof(double));
  double *temp = (double*)malloc(size*sizeof(double));


  double y = 0;
  double p = 0;
  double res = 0;

  // calcule de la LU
  seq_LU_fact(A , LU , &perm, size );

  // remplire X_start randomly
  for (size_t i = 0; i < size; i++) {
      X_start[i] = rand()%100;
  }

  old_norme = norm(X_start , size);
  for (int i = 0; i < size; i++) {
     X_start[i] = X_start[i] / old_norme;
  }
  //printf("%lf ", old_norme );

  do {

    y=p;
    p=0;

    //inversement de X_start
    for (int i = 0; i < size; ++i)  {
      temp[i] = X_start[i];
    }

    for (int i = 0; i < size; ++i)  {
      X_start[perm[i]] = temp[i];
    }

  res_sys_tri (LU, X, X_start, size, l_sys);
  res_sys_tri (LU, X_end, X, size, u_sys);

  p = norm(X_end , size);


  for (int i = 0; i < size; i++) {
      X_end[i] = X_end[i] / p;

  }

  for (size_t i = 0; i < size; i++) {
      X_start[i] = X_end[i] ;
  }

} while(absolute_value(p-y) > E);

//Write_data(size , X_end);
printf("vecteur propre : \n");
  affiche_v(X_end,size);


}