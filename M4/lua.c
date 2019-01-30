#include <stdio.h>
#include <complex.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

enum
{
  MAX_LEN=256,
  MAX_READ_LEN=1024,
  FLAVS=3
};

enum
{
  P_A=0,
  P_B,
  P_E,
  P_TOL,
  P_S12,
  P_S13,
  P_S23,
  P_OUT,
  P_MODEL,
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

typedef union
{
  double v; //value
  char n[MAX_LEN]; //name
  double complex z[FLAVS]; //psi0
} conf_par;

typedef struct
{
  uint8_t tag;
  conf_par par;
  char name[MAX_LEN];
} conf_data;

void init_conf(conf_data*);
void print_conf(FILE*, conf_data*);

int main(int argc, char **argv)
{
  int m_len,stat;
  double Re,Im;
  char tmp[MAX_READ_LEN];
  conf_data cfg[PARAMS];
  init_conf(cfg);
  
  FILE *stream=fopen("params.dat","w");

  if(argc==1)
  {
    m_len=snprintf(tmp,MAX_LEN,"%s.lua",argv[0]);
    
    if(m_len>=MAX_LEN)
    {
      fprintf(stderr,"Ошибка: слишком длинная строка (%s).",tmp);
      return 2;
    }
    
    FILE *f=fopen(tmp,"r");
    
    if(f==NULL)
    {
     fprintf(stderr,"Такого файла нет: %s.",tmp);
     return 1;
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
      fprintf(stderr,"Ошибка: слишком длинная строка (%s).",tmp);
      return 2;
    }
    FILE *f=fopen(tmp,"r");
    
      if(f==NULL)
    {
     fprintf(stderr,"Такого файла нет: %s.",tmp);
     return 1;
    }
    else
    {
       fclose(f);
    } 
  }
  lua_State *L=luaL_newstate();
  luaL_openlibs(L);
  stat=luaL_loadfile(L,tmp);
  
  if(stat!=LUA_OK)
  {
    const char *mess=lua_tostring(L, -1);
    fprintf(stderr, "[%s] SYNTAX ERROR: %s\n ", argv[0], mess);
    lua_pop(L, 1);
    lua_close(L);
    return 2;
  }
  stat = lua_pcall(L, 0, LUA_MULTRET, 0);
  if( stat != LUA_OK )
  {
    const char *mess=lua_tostring(L, -1);
    fprintf(stderr, "[%s] RUNTIME ERROR: %s\n", argv[0], mess);
    lua_pop(L, 1);
    lua_close(L);
    return 2;
  }
  
  for(int j1=0;j1<PARAMS;j1++)
  {
    stat=lua_getglobal(L,cfg[j1].name);
    if(stat==LUA_TNONE||stat==LUA_TNIL)
    {
      fprintf(stderr,"Такого параметра нет в файле: (%s).\n",cfg[j1].name);
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
      break;
      
      case LUA_TTABLE:
      for(int i=1;i-1<FLAVS;i++)
      {
        lua_pushnumber(L,i);
        stat=lua_gettable(L,-2);
    
        if(LUA_TTABLE!=stat)
        {
          fprintf(stderr,"Значение для %s в файле имеет неверный формат.\n",cfg[j1].name);
          return 1;
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
    }
  }

  print_conf(stream,cfg);

  lua_close(L);
  fclose(stream);
  
  return 0;  
}

void init_conf(conf_data *cfg)
{
  cfg[P_A].par.v=0.1;
  strncpy(cfg[P_A].name,"a",MAX_LEN);
  cfg[P_A].tag=CPD_NUMBER;
  
  cfg[P_B].par.v=1.0;
  strncpy(cfg[P_B].name,"b",MAX_LEN);
  cfg[P_B].tag=CPD_NUMBER;
  
  cfg[P_E].par.v=1.0;
  strncpy(cfg[P_E].name,"E",MAX_LEN);
  cfg[P_E].tag=CPD_NUMBER;
  
  cfg[P_TOL].par.v=1e-4;
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
  
  strncpy(cfg[P_OUT].par.n,"screen",MAX_LEN);
  strncpy(cfg[P_OUT].name,"out",MAX_LEN);
  cfg[P_OUT].tag=CPD_STRING;
  
  strncpy(cfg[P_MODEL].par.n,"sun",MAX_LEN);
  strncpy(cfg[P_MODEL].name,"model",MAX_LEN);
  cfg[P_MODEL].tag=CPD_STRING;
  
  cfg[P_PSI0].par.z[0]=-1;
  cfg[P_PSI0].par.z[1]=-1;
  cfg[P_PSI0].par.z[2]=-1;
  strncpy(cfg[P_PSI0].name,"psi0",MAX_LEN);
  cfg[P_PSI0].tag=CPD_CX_NUMBER;
}

void print_conf(FILE *stream, conf_data *cfg)
{  
  fprintf(stream,"%s=%lf\n",cfg[P_A].name,cfg[P_A].par.v);
  fprintf(stream,"%s=%lf\n",cfg[P_B].name,cfg[P_B].par.v);
  fprintf(stream,"%s=%lf\n",cfg[P_E].name,cfg[P_E].par.v);
  fprintf(stream,"%s=%lf\n",cfg[P_TOL].name,cfg[P_TOL].par.v);
  fprintf(stream,"%s=%lf\n",cfg[P_S12].name,cfg[P_S12].par.v);
  fprintf(stream,"%s=%lf\n",cfg[P_S13].name,cfg[P_S13].par.v);
  fprintf(stream,"%s=%lf\n",cfg[P_S23].name,cfg[P_S23].par.v);
  fprintf(stream,"%s=%s\n",cfg[P_OUT].name,cfg[P_OUT].par.n);
  fprintf(stream,"%s=%s\n",cfg[P_MODEL].name,cfg[P_MODEL].par.n);
  
  fprintf(stream,"%s=\n",cfg[P_PSI0].name);
  fprintf(stream,"(%lf,%lf)\n",creal(cfg[P_PSI0].par.z[0]),cimag(cfg[P_PSI0].par.z[0]));
  fprintf(stream,"(%lf,%lf)\n",creal(cfg[P_PSI0].par.z[1]),cimag(cfg[P_PSI0].par.z[1]));
  fprintf(stream,"(%lf,%lf)\n",creal(cfg[P_PSI0].par.z[2]),cimag(cfg[P_PSI0].par.z[2]));
}
