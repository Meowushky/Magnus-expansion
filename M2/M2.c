#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <complex.h>

enum
{
  FLAVS=3
};

double ex(double);

double ex(double t)
{
  double g=65956.0, n=10.54;
  return g*exp(-n*t);
}

double lambda(double p, double q, int k)
{
  double L;
  int m; //знак
  
  if(q<=0) m=1;
  else m=-1;
	L=m*2*sqrt(p/3)*cos(acos((3*q*sqrt(3))/(2*p*sqrt(p)))/3-2*M_PI*k/3);
  return L; 	
}

int main()
{ 
  double k=435.196,b0=0.030554;
  double s12=sqrt(0.308),s13=sqrt(0.0234),c12=sqrt(1-s12*s12),c13=sqrt(1-s13*s13);  
  int One[FLAVS][FLAVS]= //Единичная матрица
  {
    {1, 0, 0},
    {0, 1, 0},
    {0, 0, 1}
  };
  double H0[FLAVS][FLAVS]=
  {
    {0, 0, 0},
    {0, k*b0, 0},
    {0, 0, k}
  };
  double W[FLAVS][FLAVS]=
  {
    {c13*c13*c12*c12, c12*s12*c13*c13, c12*c13*s13},
    {c12*s12*c13*c13, s12*s12*c13*c13, s12*c13*s13},
    {c12*s13*c13, s12*c13*s13, s13*s13}
  };
  double L[FLAVS],L_r[FLAVS]; //L Собственные значения обесшпуренной матрицы, L_r -- те же значения, но переобозначенные (r-relabeled)
  double Psi_0[FLAVS]=   //Начальное значение
  {
    c12*c13,
	s12*c13,
	s13
  };
  double complex Psi[FLAVS]; //все промежуточные значения  	
  double f,Tr,Tr2,Det; // Tr, Det - след, определитель матрицы, Tr2 - след квадрата матрицы
  double d0=0.1,d1=1.0,x,h,P; //d0, d1 - интервал интегрирования, P - сумма квадратов модулей пси
  double p,q,z,A[FLAVS][FLAVS],A2[FLAVS][FLAVS],a,b; //Для манипуляций с матричной экспонентой, A2 - квадрат матрицы А
  double complex Eom2[FLAVS][FLAVS]; //Матрица эволюции(?) для 2 порядка (e^omega)
  double (*v)(double);
    
  v=ex;
  
  int N=10,i;
  double complex r0,r1;
  
  h=(d1-d0)/(N-1);
  
  for(int j1=0;j1<FLAVS;j1++) 
    Psi[j1]=Psi_0[j1]; //начальные условия
  
  for(i=0,x=d0;i<N-1;i++)
  {
   f=v(x+h/2);
   x+=h;
   Tr=0;
   Tr2=0;
   
    for(int j1=0;j1<FLAVS;j1++)
    {
	  for(int j2=0;j2<FLAVS;j2++)
	  { 
		  A[j1][j2]=-H0[j1][j2]-f*W[j1][j2];
		  //Omega[j1][j2]=A[j1][j2]*h;  
      }
      
      Tr+=A[j1][j1]; 
    }
    z=Tr/3.;
    
    //обесшпуривание
	 A[0][0]=A[0][0]-z;
	 A[1][1]=A[1][1]-z;
     A[2][2]=A[2][2]-z;  
 
    //По Крамеру находит переменную q=Det
    Det = A[0][0]*A[1][1]*A[2][2];
    Det+=A[0][1]*A[1][2]*A[2][0];
    Det+=A[0][2]*A[1][0]*A[2][1];
    Det-=A[0][0]*A[1][2]*A[2][1];
    Det-=A[0][1]*A[1][0]*A[2][2];
    Det-=A[0][2]*A[1][1]*A[2][0];
    
    //Для переменной p
    for(int j1=0;j1<FLAVS;j1++)
	  for(int j2=0;j2<FLAVS;j2++)
        Tr2+=A[j1][j2]*A[j2][j1]; 
     
    p=Tr2/2; q=Det;
       
	for(int t=0;t<FLAVS;t++)
	  L[t]=lambda(p,q,t);
	  
    //переобозначаем, чтобы собственное значение возрастало взависимости от индекса 0<1<2
    if((L[0]<L[1])&&(L[1]<L[2])) //012 уже нужное распределение
      for(int t=0;t<FLAVS;t++)
	    L_r[t]=L[t]; 
	if((L[0]<L[2])&&(L[2]<L[1]))   //021, нужно поменять местами 1 и 2
	{ 
	   L_r[0]=L[0];
	   L_r[1]=L[2];
	   L_r[2]=L[1];
    }
    if((L[1]<L[0])&&(L[0]<L[2]))  //102, нужно поменять местами 0 и 1
    { 
	   L_r[0]=L[1];
	   L_r[1]=L[0];
	   L_r[2]=L[2];
    }
    if((L[1]<L[2])&&(L[2]<L[0]))  //120,  0->1, 1->2, 2->0
    { 
	   L_r[0]=L[1];
	   L_r[1]=L[2];
	   L_r[2]=L[0];
    }
    if((L[2]<L[0])&&(L[0]<L[1]))  //201, 0->2, 1->0, 2->1
    { 
	   L_r[0]=L[2];
	   L_r[1]=L[0];
	   L_r[2]=L[1];
    }
    if((L[2]<L[1])&&(L[1]<L[0]))  //210, нужно поменять местами 0 и 2
    { 
	   L_r[0]=L[2];
	   L_r[1]=L[1];
	   L_r[2]=L[0];
    }
    
  a=L_r[1]-L_r[0];
  b=L_r[2]-L_r[0];
  //printf("a=%lf	b=%lf	",a,b);
  //printf("Exp=%lf+i%lf\n",creal(cexp(I*a*h)),cimag(cexp(I*a*h)));
  r0=-(1-cexp(I*a*h))/a;
  r1=-(-r0-(1-cexp(I*b*h)/b))/(a-b);
 // printf("r0=%lf+i%lf	r1=%lf+i%lf\n",creal(r0),cimag(r0),creal(r1),cimag(r1));
 // continue;
 
  //обнуление значений матрицы А2
  for(int j1=0;j1<FLAVS;j1++)
	for(int j2=0;j2<FLAVS;j2++)
		A2[j1][j2]=0.;
  
  for(int j1=0;j1<FLAVS;j1++)
	for(int j2=0;j2<FLAVS;j2++)
	{   
	    for(int j3=0;j3<FLAVS;j3++)
	      A2[j1][j2]+=A[j1][j3]*A[j3][j2];
        Eom2[j1][j2]=cexp(I*h*z)*cexp(I*L_r[0]*h)*((1-L_r[0]*(r0-L_r[1]*r1))*One[j1][j2]+(r0+L_r[2]*r1)*A[j1][j2]+r1*A2[j1][j2]);
    }
  
  P=0.;  
  for(int j1=0;j1<FLAVS;j1++)
  {
      Psi[j1]=Eom2[j1][0]*Psi[0];
      Psi[j1]+=Eom2[j1][1]*Psi[1];
      Psi[j1]+=Eom2[j1][2]*Psi[2];
      P+=creal(Psi[j1])*creal(Psi[j1])+cimag(Psi[j1])*cimag(Psi[j1]);
   }     
  printf("[%d] %lf \n",i,P);
}
  
  return 0;
}
