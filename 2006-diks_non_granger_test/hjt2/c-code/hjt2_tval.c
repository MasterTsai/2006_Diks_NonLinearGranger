/*
 * This program calculates p-values for the Hiemstra-Jones statistic
 * as well as for the Diks-Panchenko statistic. For details see the 
 * README file and our paper:
 *
 * Diks, C. and Panchenko, V. (2006)
 * A new statistic and practical guidelines for nonparametric Granger causality
 * testing, Journal of Economic Dynamics and Control 30 (9-10), 1647-1669
 *
 * Cees Diks -- June 2008
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define max(a,b)	a>b?a:b

#define EPS             1.5

int Ndat=15000;
int Mmax=10;     /* maximum emb. dim. */
int n, K;
double T2,Q,Qad,**A,C[4],*h,*ohm,*cov;
double p_HJ, p_T2, p_HJad, HJ_TVAL, T2_TVAL;

/* determine the log of ratios of correlation integrals */

void redun(double *x, double *y, int N, int m, int mmax, double epsilon)
{

  int i, j, s;
  int IYij, IXYij, IYZij, IXYZij;
  double disx, disy, disz, *Cy, *Cxy, *Cyz, *Cxyz, tCy=0, tCxy=0, tCyz=0, tCxyz=0;

  Cy = (double *) malloc(N*sizeof(double));
  Cxy = (double *) malloc(N*sizeof(double));
  Cyz = (double *) malloc(N*sizeof(double));
  Cxyz = (double *) malloc(N*sizeof(double));
  
  for (i=0;i!=N;i++)
    h[i] = Cy[i] = Cxy[i] = Cyz[i] = Cxyz[i] = 0.0;
  
  T2=Q=0.0;

  n = N - mmax;
  
  for (i=mmax;i!=N;i++)
  {

    Cy[i]=Cxy[i]=Cyz[i]=Cxyz[i]=0.0;

    for (j=mmax;j!=N;j++)
    if (j!=i)
    {
	    
      disx = disy = 0.0;
      for (s=1;s!=m+1;s++)
        disx = max(fabs(x[i-s]-x[j-s]),disx);
      
      for (s=1;s!=mmax+1;s++)
        disy = max(fabs(y[i-s]-y[j-s]),disy);

      if (disy <= epsilon)
      {
	Cy[i]++;
        A[3][i]++;

        if (disx <= epsilon)
        {
          Cxy[i]++;
          A[1][i]++;
        }

        disz = max(fabs(y[i]-y[j]),disy);

        if (disz <= epsilon)
        {
          Cyz[i]++;
          A[2][i]++;
          if (disx <= epsilon)
          {
            Cxyz[i]++;
            A[0][i]++;
          }
        }
      }   /* end condition |Yi - Yj| < epsilon */
    }   /* end loop over j */


    Cy[i] /= (double)(n);
    Cxy[i] /= (double)(n);
    Cyz[i] /= (double)(n);
    Cxyz[i] /= (double)(n);

    h[i] += 2.0*(Cxyz[i]*Cy[i] - Cxy[i]*Cyz[i])/6.0;

    for (j=mmax;j!=N;j++)
    if (j!=i)
    {
	    
      IYij = IXYij = IYZij = IXYZij = 0;
      
      disx = disy = 0.0;

      for (s=1;s!=m+1;s++)
        disx = max(fabs(x[i-s]-x[j-s]),disx);
     
      for (s=1;s!=mmax+1;s++)
        disy = max(fabs(y[i-s]-y[j-s]),disy);

      if (disy <= epsilon)
      {
        IYij=1;
        if (disx <= epsilon)
        {
	  IXYij = 1;
        }

        disz = max(fabs(y[i]-y[j]),disy);

        if (disz <= epsilon)
        {
	  IYZij = 1;	
        if (disx <= epsilon)
            IXYZij = 1;
        }
      }   /* end condition |Yi - Yj| < epsilon */

      h[j] += 2.0*(Cy[i]*IXYZij - Cyz[i]*IXYij)/(double)(6*n);
      h[j] += 2.0*(Cxyz[i]*IYij - Cxy[i]*IYZij)/(double)(6*n);
       
    }   /* end second loop over j */

    tCy+=n*Cy[i]; tCxy+=n*Cxy[i], tCyz+=n*Cyz[i], tCxyz+=n*Cxyz[i];

    T2 += Cxyz[i]*Cy[i] - Cyz[i]*Cxy[i];
     
  } /* end loop over i */

/*
 Alternatively one might determine T2 as the sum of the h[i]

  for (i=mmax;i!=N;i++)
    T2 += h[i];
*/

  T2 /= (double)(n);

  for (i=mmax;i!=N;i++)
  {
    h[i] -= T2;
  }

  Q = (double) tCxyz/tCxy - (double) tCyz/tCy; 

  C[0] = tCxyz/(double)(n*(n-1));
  C[1] = tCxy/(double)(n*(n-1));
  C[2] = tCyz/(double)(n*(n-1));
  C[3] = tCy/(double)(n*(n-1));

  Qad = C[0]*C[3] - C[1]*C[2];

  for (i=0;i!=4;i++)
   for (j=mmax;j!=N;j++)
   {
     A[i][j] /= (double)(n-1);
     A[i][j] -= C[i];
   } 

  free(Cy); free(Cxy); free(Cyz); free(Cxyz);
}

void InsertionSort(double *X, int *S, int M)
{
    int i, *I;
    int j;
    int r;
    double R;

    I= (int*) malloc (M*sizeof(int));

    for (i=0;i<M;i++)
      I[i]=i;

    for (i=1; i<M; i++)
      {
        R = X[i];
        r = i;
	for (j=i-1; (j>=0) && (X[j]>R); j--)
        {
	  X[j+1] = X[j];
          I[j+1] = I[j];
        }
	X[j+1] = R;
        I[j+1] = r;
      }
    for (i=0; i<M; i++)
      S[I[i]]=i;

}


void  uniform (double *X, int M)
{
  int *I, i;

  I = (int*) malloc (M*sizeof(int));
  InsertionSort(X, I, M);

  for (i=0;i<M;i++)
    X[i] = (double) I[i]/M*3.464101615;        // to make unit variance

}


/* normalize the time series to unit std. dev. */

void normalise(double *x, int N)
{
  int i;
  double mean=0.0, var=0.0;

  for (i=0;i!=N;i++)
  {
    mean += x[i];
    var += x[i]*x[i];
  }

  mean /= (double)(N);
  var /= (double)(N);
  var -= mean*mean;

  for (i=0;i!=N;i++)
    x[i] = (x[i]-mean)/sqrt(var);

  return;
}

//double rand01(void)
//{
//  return ((double)(rand())/(double)(RAND_MAX));
//}

void testcaus(double *x, double *y, int N, double epsilon, int m, int U)
{

  int i, j, k, l, mmax;
  double HJad_TVAL, S2, S2ad, VT2, d[4], sigma[4][4];
  
  for (j=0;j!=4;j++)
  {
    C[j] = 0.0;
    for (i=0;i!=N;i++)
    {
      A[j][i] = 0.0;
    }
  }

  if (U<=0)
  {
    normalise(x, N);
    normalise(y, N);
  }
  else
  {
    uniform(x, N);
    uniform(y, N);
  }

  redun(x,y,N,m,m,epsilon);

  mmax = m;
    
  for (i=0;i!=4;i++)
  {
    for (j=0;j!=4;j++)
    {
      sigma[i][j] = 0.0;
      for (k=0;k!=K;k++)
      {
        for (l=mmax+k;l!=N;l++)
          sigma[i][j] += 4.0*ohm[k]*(A[i][l]*A[j][l-k]+A[i][l-k]*A[j][l])/(double)(2*(n-k));
      }
    }
  }

  S2=S2ad=0.0;
    
  d[0] = 1.0/C[1];
  d[1] = -C[0]/(C[1]*C[1]);
  d[2] = -1.0/C[3];
  d[3] = C[2]/(C[3]*C[3]);

  for (i=0;i!=4;i++)
    for (j=0;j!=4;j++)
      S2 += d[i]*sigma[i][j]*d[j];

   HJ_TVAL = Q*sqrt(n)/sqrt(S2);

  d[0] = C[3];
  d[1] = -C[2];
  d[2] = -C[1];
  d[3] = C[0];

  for (i=0;i!=4;i++)
     for (j=0;j!=4;j++)
        S2ad += d[i]*sigma[i][j]*d[j];

  HJad_TVAL = Qad*sqrt(n)/sqrt(S2ad);
    
    
  /* determine autocovariance of h[i] */
    
  for (k=0;k!=K;k++)
  {
    cov[k] = 0.0;
    
    for (i=mmax+k;i!=N;i++)
      cov[k] += h[i]*h[i-k];

    cov[k] /= (double)(n-k);
  } 
    
  VT2=0.0;
    
/* variance of T2 */
	
  for (k=0;k!=K;k++)
    VT2 += 9.0*ohm[k]*cov[k];
  
  T2_TVAL = T2*sqrt(n)/sqrt(VT2);

  if (HJ_TVAL>0)
     p_HJ = 0.5 - 0.5*erf(HJ_TVAL/sqrt(2.0));
  else
     p_HJ = 0.5 + 0.5*erf(-HJ_TVAL/sqrt(2.0));

  if (T2_TVAL>0)
     p_T2 = 0.5 - 0.5*erf(T2_TVAL/sqrt(2.0));
  else
     p_T2 = 0.5 + 0.5*erf(-T2_TVAL/sqrt(2.0));

  if (HJad_TVAL>0)
     p_HJad = 0.5 - 0.5*erf(HJad_TVAL/sqrt(2.0));
  else
     p_HJad = 0.5 + 0.5*erf(-HJad_TVAL/sqrt(2.0));

}


int main()
{
  char filename1[128], filename2[128], filename[128];
  double x[Ndat], y[Ndat], tmp, epsilon=EPS;
  int i, j, k, m, N;
//  long seed;
  FILE *infil, *outfil;

  printf("Input file 1: "); scanf("%s", filename1);

  if ( (infil=fopen(filename1,"r")) == NULL)
  {
    fprintf(stderr,"Error: unable to open file. Exiting...\n");
    exit(1);
  }

  i = 0;
  while (fscanf(infil,"%lf", &tmp) != EOF)
  {
      x[i] = tmp;
      i++;
  }

  N = i;

  printf("%d data read\n", i);
  
  printf("Input file 2: "); scanf("%s", filename2);

  fclose(infil);

  if ( (infil=fopen(filename2,"r")) == NULL)
  {
    fprintf(stderr,"Error: unable to open file. Exiting...\n");
    exit(1);
  }

  i = 0;
  while (fscanf(infil,"%lf", &tmp) != EOF)
  {
      y[i] = tmp;
      i++;
  }

  fclose(infil);

  printf("%d data read\n", i);

  if (i!=N)
  {
    fprintf(stderr,"Error: data number mismatch. Exiting...\n");
    exit(1);
  }

  A = (double **) malloc(4*sizeof(double *));
  for (i=0;i!=4;i++)
    A[i] = (double *) malloc(N*sizeof(double));

  
  h = (double *) malloc(N*sizeof(double));
  
  K = (int)(sqrt(sqrt(N)));
  ohm = (double *) malloc(K*sizeof(double));
  cov = (double *) malloc(K*sizeof(double));

  ohm[0] = 1.0;
  for (k=1;k<K;k++)
    ohm[k] = 2.0*(1.0-k/(double)(K));

  /* Generate output file name */

  //Generate output file name
  i=0;

  while (filename1[i]!='\0' && filename1[i]!='.')
  {
    filename[i]=filename1[i];
    i++;
  }
  filename1[i]='\0';
  filename[i]='_';

  j=0;

  while (filename2[j]!='\0' && filename2[j]!='.')
  {
    filename[i+j+1]=filename2[j];
    j++;
  }

  filename2[j]='\0';
  filename[i+j+1]='.';  filename[i+j+2]='o'; filename[i+j+3]='u'; 
  filename[i+j+4]='t'; filename[i+j+5]='\0';


  if ( (outfil=fopen(filename,"w")) == NULL)
  {
    fprintf(stderr,"Error: unable to open file for writing.\
	Exiting...\n");
    exit(1);
  }

  fprintf(outfil,"file 1: %s\n", filename1);
  fprintf(outfil,"file 2: %s\n", filename2);
  fprintf(outfil,"series length: %d\n",N);
  fprintf(outfil,"epsilon=%f\n",epsilon);


  A = (double **) malloc(4*sizeof(double *));
  for (i=0;i!=4;i++)
   A[i] = (double *) malloc(N*sizeof(double));

  h = (double *) malloc(N*sizeof(double));
  
 // printf("seed: "); scanf("%ld",&seed);
 // srand(seed);


  fprintf(stderr,"epsilon=%f\n",epsilon);
  fprintf(outfil,"epsilon=%f\n",epsilon);
  /* fprintf(stderr,"       T2        TVAL       p           Q       TVAL      p          Qad      TVAL      p \n");
  */

  fprintf(stderr,"\n%s does not granger cause %s\n", filename1, filename2);
  fprintf(outfil,"\n%s does not granger cause %s\n", filename1, filename2);
  
  for (m=1;m<=Mmax;m++)
  { testcaus(x,y,N,epsilon,m,0);

    fprintf(stderr,"lX=lY=%d, p_HJ=%1.6f T_HJ=%1.6f p_T2=%1.6f T_T2=%1.6f\n",m,\
		    p_HJ, HJ_TVAL, p_T2, T2_TVAL);
    fprintf(outfil,"lX=lY=%d, p_HJ=%1.6f T_HJ=%1.6f p_T2=%1.6f  T_T2=%1.6f\n",m,\
		    p_HJ, HJ_TVAL, p_T2, T2_TVAL);
  }
  fprintf(stderr,"\n%s does not granger cause %s\n", filename2, filename1);
  fprintf(outfil,"\n%s does not granger cause %s\n", filename2, filename1);

  for (m=1;m<=Mmax;m++)
  {
    testcaus(y,x,N,epsilon,m,0);

    fprintf(stderr,"lX=lY=%d, p_HJ=%1.6f T_HJ=%1.6f p_T2=%1.6f T_T2=%1.6f\n",m,\
		    p_HJ, HJ_TVAL, p_T2, T2_TVAL);
    fprintf(outfil,"lX=lY=%d, p_HJ=%1.6f T_HJ=%1.6f p_T2=%1.6f  T_T2=%1.6f\n",m,\
		    p_HJ, HJ_TVAL, p_T2, T2_TVAL);
  }

  fprintf(stderr,"\n%s does not granger cause %s (UNIF)\n", filename1, 
filename2);
  fprintf(outfil,"\n%s does not granger cause %s (UNIF)\n", filename1, 
filename2);

  for (m=1;m<=Mmax;m++)
  {
    testcaus(x,y,N,epsilon,m,1);

    fprintf(stderr,"lX=lY=%d, p_HJ=%1.6f T_HJ=%1.6f p_T2=%1.6f T_T2=%1.6f\n",m,\
		    p_HJ, HJ_TVAL, p_T2, T2_TVAL);
    fprintf(outfil,"lX=lY=%d, p_HJ=%1.6f T_HJ=%1.6f p_T2=%1.6f  T_T2=%1.6f\n",m,\
		    p_HJ, HJ_TVAL, p_T2, T2_TVAL);
  }

  fprintf(stderr,"\n%s does not granger cause %s (UNIF)\n", filename2, 
filename1);
  fprintf(outfil,"\n%s does not granger cause %s (UNIF)\n", filename2, 
filename1);

  for (m=1;m<=Mmax;m++)
  {
    testcaus(y,x,N,epsilon,m,1);

    fprintf(stderr,"lX=lY=%d, p_HJ=%1.6f T_HJ=%1.6f p_T2=%1.6f T_T2=%1.6f\n",m,\
		    p_HJ, HJ_TVAL, p_T2, T2_TVAL);
    fprintf(outfil,"lX=lY=%d, p_HJ=%1.6f T_HJ=%1.6f p_T2=%1.6f  T_T2=%1.6f\n",m,\
		    p_HJ, HJ_TVAL, p_T2, T2_TVAL);
  }

  fclose(outfil);

  return(0);
}
