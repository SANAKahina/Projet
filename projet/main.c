#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include"fonction.h"

#define u_sys 1
#define l_sys 0

int main(int argc, char const *argv[])
{
    int n; 
    double** p=NULL;
    double* q=NULL;
    double* q2=NULL;

    double** V=NULL;
    double** H=NULL;
    double** test=NULL;
    double** LU=NULL;
    double** test2=NULL;

    n=atoi(argv[1]);
    int k = n/2;
    
    p=malloc(sizeof(double*)*n);
    q=malloc(sizeof(double*)*n);
    q2=malloc(sizeof(double*)*n);
    LU=malloc(sizeof(double*)*n);
    test2=malloc(sizeof(double*)*n);
    test=malloc(sizeof(double*)*n);

    V=(double **) malloc(n*sizeof(double*));
    for (int i=0; i<n; i++)
      V[i]=(double *) malloc(k*sizeof(double));

    
    H=(double **) malloc((k+1)*sizeof(double*));
    for (int i=0; i<n; i++)
      H[i]=(double *) malloc((k+1)*sizeof(double));
    
    for(int i=0;i<n;i++)
        p[i]=malloc(sizeof(double)*n);

    for(int i=0;i<n;i++)
        test[i]=malloc(sizeof(double)*n);

    for(int i=0;i<n;i++)
        test2[i]=malloc(sizeof(double)*n);
         
      for(int i=0;i<n;i++)
        LU[i]=malloc(sizeof(double)*n);

    init_mat(p,n, 0);
    init_mat(LU,n, 0);

    init_v(q,n, 0);
    q[0]=1;
    q2[0]=55;

printf("Vecteur Test :\n");
affiche_v(q,n);

printf("Matrices Test :\n");
    init_mat_test(test, n,2,-1);
    affiche_mat(test,n);
    printf("\n");


    init_mat_test2(test2,n);
    affiche_mat(test2,n);
    printf("\n");

printf("****arnoldi Methode******\n");
//arnoldi(test, q2, k, n, V,H);

printf("Matrice H :\n");
affiche_mat(H,k+1);

printf("Matrice Q:\n");
affiche_mat(V,k+1);

printf("*****Methode puissance inverse****** :\n");

printf("factorisation LU de la Matrice test :\n");


size_t *perm;
seq_LU_fact(test ,LU, &perm,n);
affiche_mat(LU,n);

const double E=0.00001;
printf("\n");

invers_power(test , n , E);  

    //pour le free
    for(int i=0;i<n;i++)free(p[i]);
    free(p);

    return 0;
}


