#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>

enum
{
  MAX_LEN=256,
  MAX_READ_LEN=1024,
};

int main(int argc, char **argv)
{
  int m_len,stat;
  double a=0.1,b,E,tol,s12,s13;
  const char *out,*model;
  char tmp[MAX_READ_LEN];
  //char filen[MAX_LEN];

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
  stat=lua_getglobal(L,"a");
  if(stat==LUA_TNONE||stat==LUA_TNIL)
  {
    fprintf(stderr,"Такого параметра нет в файле: (%s).\n","a");
  }
  
  if(stat==LUA_TNUMBER)
  {
    a=lua_tonumber(L,-1);
    lua_pop(L,1);
  }
  fprintf(stderr,"a=%lf\n",a);
  
  stat=lua_getglobal(L,"b");
  if(stat==LUA_TNONE||stat==LUA_TNIL)
  {
    fprintf(stderr,"Такого параметра нет в файле: (%s).\n","b");
  }
  
  if(stat==LUA_TNUMBER)
  {
    b=lua_tonumber(L,-1);
    lua_pop(L,1);
  }
  fprintf(stderr,"b=%lf\n",b);
  
  stat=lua_getglobal(L,"E");
  if(stat==LUA_TNONE||stat==LUA_TNIL)
  {
    fprintf(stderr,"Такого параметра нет в файле: (%s).\n","E");
  }
  
  if(stat==LUA_TNUMBER)
  {
    E=lua_tonumber(L,-1);
    lua_pop(L,1);
  }
  fprintf(stderr,"E=%lf\n",E);
  
  stat=lua_getglobal(L,"tol");
  if(stat==LUA_TNONE||stat==LUA_TNIL)
  {
    fprintf(stderr,"Такого параметра нет в файле: (%s).\n","tol");
  }
  
  if(stat==LUA_TNUMBER)
  {
    tol=lua_tonumber(L,-1);
    lua_pop(L,1);
  }
  fprintf(stderr,"tol=%lf\n",tol);
  
  stat=lua_getglobal(L,"s12");
  if(stat==LUA_TNONE||stat==LUA_TNIL)
  {
    fprintf(stderr,"Такого параметра нет в файле: (%s).\n","s12");
  }
  
  if(stat==LUA_TNUMBER)
  {
    s12=lua_tonumber(L,-1);
    lua_pop(L,1);
  }
  fprintf(stderr,"s12=%lf\n",s12);
  
    stat=lua_getglobal(L,"s13");
  if(stat==LUA_TNONE||stat==LUA_TNIL)
  {
    fprintf(stderr,"Такого параметра нет в файле: (%s).\n","s13");
  }
  
  if(stat==LUA_TNUMBER)
  {
    s13=lua_tonumber(L,-1);
    lua_pop(L,1);
  }
  fprintf(stderr,"s13=%lf\n",s13);
  
    stat=lua_getglobal(L,"out");
  if(stat==LUA_TNONE||stat==LUA_TNIL)
  {
    fprintf(stderr,"Такого параметра нет в файле: (%s).\n","out");
  }
  
  if(stat==LUA_TSTRING)
  {
    out=lua_tostring(L,-1);
    lua_pop(L,1);
  }
  fprintf(stderr,"out=%s\n",out);
  
    stat=lua_getglobal(L,"model");
  if(stat==LUA_TNONE||stat==LUA_TNIL)
  {
    fprintf(stderr,"Такого параметра нет в файле: (%s).\n","model");
  }
  
  if(stat==LUA_TSTRING)
  {
    model=lua_tostring(L,-1);
    lua_pop(L,1);
  }
  fprintf(stderr,"model=%s\n",model);
  
  lua_close(L);

  return 0;  
}
