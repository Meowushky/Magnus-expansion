#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <complex.h>
#include <stdint.h>
#include <stdbool.h>
#include <string.h>
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

enum
{
  FLAVS=3,
  MAX_LEN=256,
  MAX_PAR_LIST=5,
  MAX_READ_LEN=1024,
  MODELS=2
};

enum
{
 SUN=0,
 SN
};

enum
{
  P_MODEL=0,
  P_A,
  P_B,
  P_E,
  P_TOL,
  P_S12,
  P_S13,
  P_S23,
  P_OUT,
  P_PSI0,
  P_EOPL,
  PARAMS=10
};

enum
{
  CPD_NUMBER=0,
  CPD_CX_NUMBER,
  CPD_STRING,
  CPD_VECTOR,
  CPD_CX_VECTOR,
  CPD_MATRIX,
  CPD_CX_MATRIX
};

static const double norm_accuracy=1e-10;
static double model_par[MODELS][MAX_PAR_LIST]=
{{65956.,10.54,0,0,0},{52.934,0,0,0,0}};
uint8_t chosen_model=0;

typedef struct
{
  double a,b;
  char name[MAX_LEN];
  double (*dens)(double);
  double par[MAX_PAR_LIST];
} model_info;

typedef struct
{
  double d0,d1,tol;
  double (*dens)(double);
  double (*H0)[FLAVS],(*W)[FLAVS],(*cH0W)[FLAVS],(*cH0H0W)[FLAVS],(*cWH0W)[FLAVS];  
} wf_ctx;

typedef struct
{
  double complex Psi[FLAVS];
  uint64_t calls;
  double last_step;
  double prev_step;
} rwf_ctx;

typedef struct
{
  uint8_t tag;
  union
  {
    double v; //value
    char n[MAX_LEN]; //name
    double complex z[FLAVS]; //psi0
  } par;
  char name[MAX_LEN];
} conf_data;

void EigenV(double*, double, double);
void aWF_calc(wf_ctx*, rwf_ctx*);
void init_model_data(model_info*);
double sun_model(double);
double sn_model(double);
void init_conf(conf_data*, model_info*);
void print_conf(FILE*, conf_data*);
void get_conf(char*,char*, conf_data*, model_info*);

int main(int argc,char **argv)
{ 
  if(P_EOPL!=PARAMS)
  {
    fprintf(stderr, "[ERROR] Непредвиденная внутренняя ошибка программы: несогласованное состояние настроечных параметров: ожидаемое число параметров %d, запрограммированное: %d\n",
    P_EOPL, PARAMS);
    return 1;
  }
  
  int m_len;
  double mod2;
  char tmp[MAX_READ_LEN];
  
  FILE *stream;
  
  model_info mode[MODELS];
  init_model_data(mode);
  
  //параметры по умолчанию
  
  conf_data cfg[PARAMS];
  init_conf(cfg,mode);
  
  //чтение параметров из lua
  
  if(argc==1)
  {
    m_len=snprintf(tmp,MAX_LEN,"%s.lua",argv[0]);
    
    if(m_len>=MAX_LEN)
    {
      fprintf(stderr,"Ошибка! Файл %s имеет слишком длинное имя: количество символов превышает допустимое значение %d.\n",
      tmp,MAX_LEN);
      return 2;
    }
    
    FILE *f=fopen(tmp,"r");
    
    if(f==NULL)
    {
     fprintf(stderr,"Ошибка! Указанный файл c параметрами %s не найден.\n",
     tmp);
    }
    else
    {
      fclose(f);
    } 
  }

  if(argc>1)
  {
    m_len=snprintf(tmp,MAX_LEN,"%s",argv[1]);
    
    if(m_len>=MAX_LEN)
    {
      fprintf(stderr,"Ошибка! Файл %s имеет слишком длинное имя: количество символов превышает допустимое значение %d.\n",
      tmp,MAX_LEN);
      return 2;
    }
    FILE *f=fopen(tmp,"r");
    
    if(f==NULL)
    {
     fprintf(stderr,"Ошибка! Указанный файл с параметрами %s не найден.\n",
     tmp);
     return 1;
    }
    else
    {
       fclose(f);
    }
    
    get_conf(argv[0],tmp,cfg,mode);
  }
  
  //проверка параметров
  
  //выбор модели
  
  //a0<a
  if(mode[chosen_model].a>cfg[P_A].par.v)
  {
    fprintf(stderr,"Ошибка! Заданная начальная точка интервала слишком мала %s=%lf: заходит за границы применимости используемой модели с %s=%lf.\n",
    cfg[P_A].name, cfg[P_A].par.v, cfg[P_A].name, mode[chosen_model].a);
    return 1;
  }
  
  //b<=b0
  if(mode[chosen_model].b<cfg[P_B].par.v)
  {
    fprintf(stderr,"Ошибка! Заданная конечная точка интервала слишком большая %s=%lf: заходит за границы применимости используемой модели с %s=%lf.\n",
    cfg[P_B].name, cfg[P_B].par.v, cfg[P_B].name, mode[chosen_model].b);
    return 1;
  }
  
  //a<b
  if(cfg[P_A].par.v>cfg[P_B].par.v)
  {
    fprintf(stderr,"Ошибка! Неверно заданы параметры: конечная точка интервала %s=%lf не может быть меньше начальной %s=%lf.\n",
    cfg[P_B].name, cfg[P_B].par.v, cfg[P_A].name, cfg[P_A].par.v);
    return 1;
  }

  //E>0
  if(cfg[P_E].par.v<=0)
  {
    fprintf(stderr,"Ошибка! Задана отрицательная или нулевая энергия %s=%lf.\n",
    cfg[P_E].name, cfg[P_E].par.v);
    return 1;
  }
  
  //TOL>0
  if(cfg[P_TOL].par.v<=0)
  {
    fprintf(stderr,"Ошибка! Задан отрицательный параметр чувствительности %s=%lf.\n",
    cfg[P_TOL].name, cfg[P_TOL].par.v);
    return 1;
  }
  
  //|S12|<=1
  if(fabs(cfg[P_S12].par.v)>1)
  {
    fprintf(stderr,"Ошибка! Неправильно задан угол: абсолютное значение синуса_12 |%s|=%lf не может принимать значения, больше 1.\n",
    cfg[P_S12].name, cfg[P_S12].par.v);
    return 1;
  }
  
  //|S13|<=1
  if(fabs(cfg[P_S13].par.v)>1)
  {
    fprintf(stderr,"Ошибка! Неправильно задан угол: абсолютное значение синуса_13 |%s|=%lf не может принимать значения, больше 1.\n",
    cfg[P_S13].name, cfg[P_S13].par.v);
    return 1;
  }
  
  //|S23|<=1
  if(fabs(cfg[P_S23].par.v)>1)
  {
    fprintf(stderr,"Ошибка! Неправильно задан угол: абсолютное значение синуса_23 |%s|=%lf не может принимать значения, больше 1.\n",
    cfg[P_S23].name, cfg[P_S23].par.v);
    return 1;
  }
  
  //|PSI0|^2-1>10^{-10}
  mod2=cfg[P_PSI0].par.z[0]*conj(cfg[P_PSI0].par.z[0]);
  mod2+=cfg[P_PSI0].par.z[1]*conj(cfg[P_PSI0].par.z[1]);
  mod2+=cfg[P_PSI0].par.z[2]*conj(cfg[P_PSI0].par.z[2]);
  
  if(fabs(mod2-1.)>norm_accuracy)
  {
    fprintf(stderr,"Ошибка! Квадрат модуля заданного вектора |%s|^2-1=%4.3e больше запрограммированного параметра %4.3e.\n",
    cfg[P_PSI0].name, mod2-1., norm_accuracy);
    return 1;
  }
  
  if(strcmp(cfg[P_OUT].par.n,"SCREEN")==0)
  {
   stream=stderr; 
  }
  else
  {
    strncpy(tmp,cfg[P_OUT].par.n,MAX_LEN-1);
    tmp[MAX_LEN]='\0';
/*  snprintf(tmp,MAX_LEN,"%s-%s_E%4.3lf_a%4.3lf_b%4.3lf_t%4.3lf.dat",
    argv[0],
    mode[chosen_model].name,
    cfg[P_E].par.v,
    cfg[P_A].par.v,
    cfg[P_B].par.v,
    cfg[P_TOL].par.v);
*/
    stream=fopen(tmp,"w");
  }

  double a=4351960.,b=0.030554,E=cfg[P_E].par.v,tol=cfg[P_TOL].par.v;
  double s12=cfg[P_S12].par.v,
    s13=cfg[P_S13].par.v,
    c12=sqrt(1.-s12*s12),
    c13=sqrt(1.-s13*s13);
  double W[FLAVS][FLAVS]=
  {
    {c13*c13*c12*c12, c12*s12*c13*c13, c12*c13*s13},
    {c12*s12*c13*c13, s12*s12*c13*c13, s12*c13*s13},
    {c12*s13*c13, s12*c13*s13, s13*s13}
  };

  double d0=cfg[P_A].par.v,
    d1=cfg[P_B].par.v;
  double Pee;
  
  double H0[FLAVS][FLAVS]=
  {
    {0., 0., 0.},
    {0., a*b/E, 0.},
    {0., 0., a/E}
  };
  double cH0W[FLAVS][FLAVS],cH0H0W[FLAVS][FLAVS],cWH0W[FLAVS][FLAVS];
  
  //вычисление коммутатора H0 с W
  for(int j1=0;j1<FLAVS;j1++)
    for(int j2=0;j2<FLAVS;j2++)
      for(int j3=0;j3<FLAVS;j3++)
        cH0W[j1][j2]=H0[j1][j3]*W[j3][j2]-W[j1][j3]*H0[j3][j2];
      
  //вычисление коммутатора [H0,[H0,W]]
  for(int j1=0;j1<FLAVS;j1++)
    for(int j2=0;j2<FLAVS;j2++)
      for(int j3=0;j3<FLAVS;j3++)
        cH0H0W[j1][j2]=H0[j1][j3]*cH0W[j3][j2]-cH0W[j1][j3]*H0[j3][j2];
      
  //вычисление коммутатора [W,[H0,W]]
  for(int j1=0;j1<FLAVS;j1++)
    for(int j2=0;j2<FLAVS;j2++)
      for(int j3=0;j3<FLAVS;j3++)
        cWH0W[j1][j2]=W[j1][j3]*cH0W[j3][j2]-cH0W[j1][j3]*W[j3][j2];
  
  wf_ctx ctx;
  ctx.d0=d0;
  ctx.d1=d1;
  ctx.dens=mode[chosen_model].dens;
  ctx.H0=H0;
  ctx.W=W;
  ctx.cH0W=cH0W;
  ctx.cH0H0W=cH0H0W;
  ctx.cWH0W=cWH0W;
  ctx.tol=tol;
  
  rwf_ctx res;
  res.Psi[0]=cfg[P_PSI0].par.z[0];
  res.Psi[1]=cfg[P_PSI0].par.z[1];
  res.Psi[2]=cfg[P_PSI0].par.z[2];
  res.calls=0;
  res.last_step=-1.;
  res.prev_step=-1.;

  print_conf(stream,cfg);
  
  aWF_calc(&ctx,&res);
  
  Pee=c12*c12*c13*c13*res.Psi[0]*conj(res.Psi[0]);
  Pee+=s12*s12*c13*c13*res.Psi[1]*conj(res.Psi[1]);
  Pee+=s13*s13*res.Psi[2]*conj(res.Psi[2]);
  
  mod2=res.Psi[0]*conj(res.Psi[0]);
  mod2+=res.Psi[1]*conj(res.Psi[1]);
  mod2+=res.Psi[2]*conj(res.Psi[2]);;
  
  fprintf(stream,"## |psi1|=%9.8lf\n", sqrt(creal(res.Psi[0])*creal(res.Psi[0])+cimag(res.Psi[0])*cimag(res.Psi[0])));
  fprintf(stream,"## |psi2|=%9.8lf\n", sqrt(creal(res.Psi[1])*creal(res.Psi[1])+cimag(res.Psi[1])*cimag(res.Psi[1])));
  fprintf(stream,"## |psi3|=%9.8lf\n", sqrt(creal(res.Psi[2])*creal(res.Psi[2])+cimag(res.Psi[2])*cimag(res.Psi[2])));
  
  fprintf(stream,"# calls=%ld\n", res.calls);
  fprintf(stream,"# prev. step=%4.3e\n", res.prev_step);
  fprintf(stream,"# last step=%4.3e\n", res.last_step);
  fprintf(stream,"# psi={{%11.10lf,%11.10lf},{%11.10lf,%11.10lf},{%11.10lf,%11.10lf}}\n",
    creal(res.Psi[0]),cimag(res.Psi[0]),
    creal(res.Psi[1]),cimag(res.Psi[1]),
    creal(res.Psi[2]),cimag(res.Psi[2]));
  fprintf(stream,"# |psi|^2-1=%11.10lf\n", mod2-1.);
  
  fprintf(stream,"# a b E Pee\n");  
  fprintf(stream,"%lf\t%lf\t%lf\t%lf\n", d0, d1, E, Pee);
  
  if(stream!=stderr)
  {
    fclose(stream);
  }
  
return 0;
}

void aWF_calc(wf_ctx *ctx, rwf_ctx *res)
{
  bool flag=false;
  double h,x,Er,S1[FLAVS][FLAVS];
  x=ctx->d0;
  h=ctx->tol/2.;
  double xi_m,xi_p,f_m,f_p;
  double z,p,q,F;
  double L[2*FLAVS-1];
  double complex A[FLAVS][FLAVS],Eom4[FLAVS][FLAVS],S2[FLAVS][FLAVS];
  double complex r0,r1,psi_0[FLAVS],psi_1[FLAVS];
  double One[FLAVS][FLAVS]=
  {
    {1., 0., 0.},
    {0., 1., 0.},
    {0., 0., 1.}
  };
  
  while(x<ctx->d1) 
  {
    if(true==flag)
    {
      h=ctx->d1-x;
    }
    
    if((x+h<ctx->d1)&&(x+2.*h>ctx->d1)&&(false==flag)) 
      {
        res->prev_step=h;
        h=(ctx->d1-x)/2.;
          
        flag=true;
      }
    
    res->calls+=1;
    
    xi_m=x+(1.-1./sqrt(3.))*h/2.;
    xi_p=x+(1.+1./sqrt(3.))*h/2.;
    f_m=ctx->dens(xi_m);
    f_p=ctx->dens(xi_p);
    x+=h;
  
    z=0.;
    
    for(int j1=0;j1<FLAVS;j1++) 
      for(int j2=0;j2<FLAVS;j2++)
      {
        S1[j1][j2]=-sqrt(3.)*(f_p-f_m)*ctx->cH0W[j1][j2]/12.;
        S2[j1][j2]=I*sqrt(3.)*(f_p-f_m)*(ctx->cH0H0W[j1][j2]+(f_p+f_m)*ctx->cWH0W[j1][j2]/2.)/24.;
      }
    for(int j1=0;j1<FLAVS;j1++)
    {
      for(int j2=0;j2<FLAVS;j2++)
      {
        A[j1][j2]=-ctx->H0[j1][j2]-(f_p+f_m)*ctx->W[j1][j2]/2.+
        I*S1[j1][j2]*h;
      }
      z+= (double) A[j1][j1]/3.;
    }
    
    //обесшпуривание
    A[0][0]=A[0][0]-(double complex) z;
    A[1][1]=A[1][1]-(double complex) z;
    A[2][2]=A[2][2]-(double complex) z;  
 
    //По Крамеру находит переменную q=Det
    q = (double) A[0][0]*A[1][1]*A[2][2];
    q+= (double) A[0][1]*A[1][2]*A[2][0];
    q+= (double) A[0][2]*A[1][0]*A[2][1];
    q-= (double) A[0][0]*A[1][2]*A[2][1];
    q-= (double) A[0][1]*A[1][0]*A[2][2];
    q-= (double) A[0][2]*A[1][1]*A[2][0];
    
    p=0.;

    for(int j1=0;j1<FLAVS;j1++)
      p+= (double) (A[j1][0]*A[0][j1]+A[j1][1]*A[1][j1]+A[j1][2]*A[2][j1])/2.;

    EigenV(L,p,q);
    
    F=2.*sqrt(p/3.);

    r0=-(2.*sin(L[3]*F*h/2.)*sin(L[3]*F*h/2.)-I*sin(L[3]*F*h))/L[3]; // r0/F
    r1=-(-r0-(2.*sin(L[4]*F*h/2.)*sin(L[4]*F*h/2.)-I*sin(L[4]*F*h))/L[4])/(L[3]-L[4]);// r1/F^2

    for(int j1=0;j1<FLAVS;j1++)
      for(int j2=0;j2<FLAVS;j2++)
      {
        Eom4[j1][j2]=cexp(I*h*z)*cexp(I*L[0]*F*h)*
        ((1.-L[0]*(r0-L[1]*r1))*One[j1][j2]+
        (r0+L[2]*r1)*A[j1][j2]/F+
        r1*(A[j1][0]*A[0][j2]+A[j1][1]*A[1][j2]+A[j1][2]*A[2][j2])/(F*F));
      }
    
    psi_1[0]=res->Psi[0]; 
    psi_1[1]=res->Psi[1]; 
    psi_1[2]=res->Psi[2]; //сохранение предыдущих значений
    
    for(int j1=0;j1<FLAVS;j1++)
    {
      psi_0[j1]=Eom4[j1][0]*res->Psi[0];
      psi_0[j1]+=Eom4[j1][1]*res->Psi[1];
      psi_0[j1]+=Eom4[j1][2]*res->Psi[2];
    }

    for(int j1=0;j1<FLAVS;j1++)
      res->Psi[j1]=psi_0[j1];

    for(int j1=0;j1<FLAVS;j1++)
      for(int j2=0;j2<FLAVS;j2++)
        for(int j3=0;j3<FLAVS;j3++)
          psi_0[j1]=(S1[j1][j2]+h*S2[j1][j2]+h*h*S1[j1][j3]*S1[j3][j2])*res->Psi[j2];
      
    Er=h*h*sqrt( (double) psi_0[0]*conj(psi_0[0])+psi_0[1]*conj(psi_0[1])+psi_0[2]*conj(psi_0[2]));
    
    if((Er>ctx->tol)&&(false==flag))
    { 
      x=x-h;
      h=h*0.8*pow(ctx->tol/Er,1./3.);
      res->Psi[0]=psi_1[0]; 
      res->Psi[1]=psi_1[1]; 
      res->Psi[2]=psi_1[2]; 
    }
  }
  res->last_step=h;
}

void init_model_data(model_info *mode)
{
  strncpy(mode[SUN].name,"sun",MAX_LEN);
  mode[SUN].a=0.1;
  mode[SUN].b=1.;
  mode[SUN].par[0]=model_par[SUN][0];
  mode[SUN].par[1]=model_par[SUN][1];
  mode[SUN].dens=sun_model;
  
  strncpy(mode[SN].name,"sn",MAX_LEN);
  mode[SN].a=0.02;
  mode[SN].b=20.;
  mode[SN].par[0]=model_par[SN][0];
  mode[SN].dens=sn_model;
}

double sun_model(double t)
{
  return model_par[SUN][0]*exp(-model_par[SUN][1]*t);
}

double sn_model(double t)
{
  return model_par[SN][0]/(t*t*t);
}

void EigenV(double *L,double p, double q)
{
  double arg,u;
  arg=acos((3.*q*sqrt(3.))/(2.*p*sqrt(p)))/3.;

  L[0]=cos(arg);
  L[1]=cos(arg-2.*M_PI/3.);

  if(L[1]<L[0])
  {
    u=L[0];
    L[0]=L[1];
    L[1]=u;
  }
  
  L[2]=cos(arg-4.*M_PI/3.);

  if(L[2]<L[0])
  {
    u=L[0];
    L[0]=L[2];
    L[2]=u;

    u=L[1];
    L[1]=L[2];
    L[2]=u;
  }
  else
    if(L[2]<L[1])
      {
        u=L[1];
        L[1]=L[2];
        L[2]=u;
      }
  L[3]=L[1]-L[0]; //a
  L[4]=L[2]-L[0]; //b
}

void init_conf(conf_data *cfg, model_info *mode)
{
  cfg[P_A].par.v=mode[chosen_model].a;
  strncpy(cfg[P_A].name,"a",MAX_LEN);
  cfg[P_A].tag=CPD_NUMBER;
  
  cfg[P_B].par.v=mode[chosen_model].b;
  strncpy(cfg[P_B].name,"b",MAX_LEN);
  cfg[P_B].tag=CPD_NUMBER;
  
  cfg[P_E].par.v=1.0;
  strncpy(cfg[P_E].name,"E",MAX_LEN);
  cfg[P_E].tag=CPD_NUMBER;
  
  cfg[P_TOL].par.v=1e-5;
  strncpy(cfg[P_TOL].name,"tol",MAX_LEN);
  cfg[P_TOL].tag=CPD_NUMBER;
  
  cfg[P_S12].par.v=sqrt(0.308);
  strncpy(cfg[P_S12].name,"s12",MAX_LEN);
  cfg[P_S12].tag=CPD_NUMBER;
  
  cfg[P_S13].par.v=sqrt(0.0234);
  strncpy(cfg[P_S13].name,"s13",MAX_LEN);
  cfg[P_S13].tag=CPD_NUMBER;
  
  cfg[P_S23].par.v=sqrt(0.437);
  strncpy(cfg[P_S23].name,"s23",MAX_LEN);
  cfg[P_S23].tag=CPD_NUMBER;
  
  strncpy(cfg[P_OUT].par.n,"SCREEN",MAX_LEN);
  strncpy(cfg[P_OUT].name,"out",MAX_LEN);
  cfg[P_OUT].tag=CPD_STRING;
  
  strncpy(cfg[P_MODEL].par.n,"sun",MAX_LEN);
  strncpy(cfg[P_MODEL].name,"model",MAX_LEN);
  cfg[P_MODEL].tag=CPD_STRING;
  
  cfg[P_PSI0].par.z[0]=sqrt((1-cfg[P_S12].par.v*cfg[P_S12].par.v)*(1-cfg[P_S13].par.v*cfg[P_S13].par.v));
  cfg[P_PSI0].par.z[1]=cfg[P_S12].par.v*sqrt(1-cfg[P_S13].par.v*cfg[P_S13].par.v);
  cfg[P_PSI0].par.z[2]=cfg[P_S13].par.v;
  strncpy(cfg[P_PSI0].name,"psi0",MAX_LEN);
  cfg[P_PSI0].tag=CPD_CX_NUMBER;
}

void print_conf(FILE *stream, conf_data *cfg)
{  
  fprintf(stream,"# %s=%lf\n", cfg[P_A].name, cfg[P_A].par.v);
  fprintf(stream,"# %s=%lf\n", cfg[P_B].name, cfg[P_B].par.v);
  fprintf(stream,"# %s=%e\n", cfg[P_E].name, cfg[P_E].par.v);
  fprintf(stream,"# %s=%e\n", cfg[P_TOL].name, cfg[P_TOL].par.v);
  fprintf(stream,"# %s=%lf\n", cfg[P_S12].name, cfg[P_S12].par.v);
  fprintf(stream,"# %s=%lf\n", cfg[P_S13].name, cfg[P_S13].par.v);
  fprintf(stream,"# %s=%lf\n", cfg[P_S23].name, cfg[P_S23].par.v);
  fprintf(stream,"# %s^2=%lf\n", cfg[P_S12].name, cfg[P_S12].par.v*cfg[P_S12].par.v);
  fprintf(stream,"# %s^2=%lf\n", cfg[P_S13].name, cfg[P_S13].par.v*cfg[P_S13].par.v);
  fprintf(stream,"# %s^2=%lf\n", cfg[P_S23].name, cfg[P_S23].par.v*cfg[P_S23].par.v);
  fprintf(stream,"# %s=%s\n", cfg[P_OUT].name, cfg[P_OUT].par.n);
  fprintf(stream,"# %s=%s\n", cfg[P_MODEL].name, cfg[P_MODEL].par.n);
  
  fprintf(stream,"# %s=", cfg[P_PSI0].name);
  fprintf(stream,"%lf\t%lf\t", creal(cfg[P_PSI0].par.z[0]), cimag(cfg[P_PSI0].par.z[0]));
  fprintf(stream,"%lf\t%lf\t", creal(cfg[P_PSI0].par.z[1]), cimag(cfg[P_PSI0].par.z[1]));
  fprintf(stream,"%lf\t%lf\n", creal(cfg[P_PSI0].par.z[2]), cimag(cfg[P_PSI0].par.z[2]));
}

void get_conf(char *prog_name, char *cfgfile, conf_data *cfg, model_info *mode)
{
  double Re,Im;
  int stat;
  
  lua_State *L=luaL_newstate();
  luaL_openlibs(L);
  stat=luaL_loadfile(L,cfgfile);
  
  if(stat!=LUA_OK)
  {
    const char *mess=lua_tostring(L, -1);
    fprintf(stderr, "[%s] SYNTAX ERROR: %s\n ", prog_name, mess);
    lua_pop(L, 1);
    lua_close(L);
    exit(2);
  }
  stat = lua_pcall(L, 0, LUA_MULTRET, 0);
  if( stat != LUA_OK )
  {
    const char *mess=lua_tostring(L, -1);
    fprintf(stderr, "[%s] RUNTIME ERROR: %s\n", prog_name, mess);
    lua_pop(L, 1);
    lua_close(L);
    exit(2);
  }
  
  for(int j1=0;j1<PARAMS;j1++)
  {
    stat=lua_getglobal(L,cfg[j1].name);
    if(stat==LUA_TNONE||stat==LUA_TNIL)
    {
      fprintf(stderr,"Такого параметра нет в файле: (%s).\n",
      cfg[j1].name);
      continue;
    }
    switch(stat)
    {
      case LUA_TNUMBER:
        cfg[j1].par.v=lua_tonumber(L,-1);
        lua_pop(L,1);
        break;
      
      case LUA_TSTRING:
        strncpy(cfg[j1].par.n,lua_tostring(L,-1),MAX_LEN);
        lua_pop(L,1);
          if(strcmp(cfg[j1].name,"model")==0)
          {
              while(chosen_model<MODELS)
              {
                if(strcmp(cfg[P_MODEL].par.n, mode[chosen_model].name)==0)
                {
                  break;
                }
                chosen_model++;
              }
              cfg[P_A].par.v=mode[chosen_model].a;
              cfg[P_B].par.v=mode[chosen_model].b;
          }
        break;
      
      case LUA_TTABLE:
        for(int i=1;i-1<FLAVS;i++)
        {
          lua_pushnumber(L,i);
          stat=lua_gettable(L,-2);

          if(LUA_TTABLE!=stat)
          {
            fprintf(stderr,"Значение для %s в файле имеет неверный формат.\n",
            cfg[j1].name);
            exit(1);
          }
    
          lua_pushnumber(L,1);
          stat=lua_gettable(L,-2);
          Re=lua_tonumber(L,-1);
          lua_pop(L,1);

          lua_pushnumber(L,2);
          stat=lua_gettable(L,-2);
          Im=lua_tonumber(L,-1);
          lua_pop(L,2);
    
          cfg[j1].par.z[i-1]=Re+I*Im;
        }
        break;
      
        default:
          fprintf(stderr,"Неверно задано значение.\n");
          exit(1);
    }
  }

  lua_close(L);
}
