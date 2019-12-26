#include <stdio.h>
#include <stdlib.h>
#include <math.h>


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

void  gradient_conjugue(double **A ,double *b ,double *d,double *r,int n)
{	
	
	int iter =10;
	double E=0.00001;
	double alpha=0,beta=0;
	double *rr,*dd;
	double*x,*xx;
	double *AX;
	x=malloc(sizeof(double)*n);
	AX=malloc(sizeof(double)*n);
	rr=malloc(sizeof(double)*n);
	xx=malloc(sizeof(double)*n);
	dd=malloc(sizeof(double)*n);
	
	init_v(AX, n, 0);
	init_v(rr, n, 0);
	init_v(dd, n, 0);
	init_v(x, n,0);
  //affiche_v(x,n);
	AX=dot(x,A,n);
	r=sub_vec(b,AX,n);
  //printf("initialisation de x :\n");
  //affiche_v(x,n);

  printf("début \n");
do {
	//L’algorithme démarre avec le choix d0=r0=b−A x0 
	for(int i=0;i<n;i++)
	{
		d[i]=r[i];
	}
	for (int j=0;j<iter;j++)
	{
		alpha=(norm(r,n)/prod_scal(dot(d,A,n),d,n));
		//mise a jour de x 
		xx=somme_vec(x,prod_v_scal(d,n,alpha),n);
		
		//mise a jour du residus
		rr=sub_vec(r,prod_v_scal(dot(d,A,n),n,alpha),n);
		beta=norm(rr,n)/norm(r,n);
		dd=somme_vec(rr,prod_v_scal(d,n,beta),n);

          x[j]=xx[j];
          r[j]=rr[j];
          d[j]=dd[j];
printf("x_iter[%d]\n",j);
affiche_v(xx,n);
  printf("\n");

	}
  

}while (norm(r,n)<E);

free(dd);
free(xx);
free(rr);
printf("xx resultant 'solution': \n");
affiche_v(xx,n);
printf("x resultant 'solution': \n");

affiche_v(x,n);

printf("fin\n");
}
int main(int argc, char const *argv[])
{
    double** p=NULL;int n; 
    double* q=NULL;
    double* q2=NULL;

    double** V=NULL;
    double** H=NULL;
    double** test=NULL;
    double** LU=NULL;
    double** test2=NULL;


    //scanf("%d",&n);
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
    
    
    double *d=malloc(sizeof(double*)*n);
		double *b=malloc(sizeof(double*)*n);
		double *r=malloc(sizeof(double*)*n);
		
    double *d2=malloc(sizeof(double*)*n);
    double *b2=malloc(sizeof(double*)*n);
    double *r2=malloc(sizeof(double*)*n);
    
    
    

	
  init_mat_test2(test2,n);
  printf("\n");


	 init_v(q,n,0);
    init_v(d,n,0);
    init_v(r,n,0);


   init_v(q,n,0);
    init_v(d2,n,0);
    init_v(r2,n,0);

    init_v(b,n,0);
    init_v(b2,n,0);

    for (int i=0;i<n;i++)
    {    
		b[i]=i;
    b2[i]=i;


	}

  //affiche_v(b,n);

    q[0]=1;
    q2[0]=1;

  printf("\n");

    init_mat_test(test, n,2,-1);
    affiche_mat(test2,n);


printf("affiche gradient sur test2\n");
gradient_conjugue(test2 ,b,d,r,n);
printf("affiche gradient sur test\n");
gradient_conjugue(test ,b2,d2,r2,n);


    for(int i=0;i<n;i++)free(p[i]);
    free(p);

    return 0;
}


