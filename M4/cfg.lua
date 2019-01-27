a=0.1
b=1.0
E=10.0
tol=1e-3
out='screen'
model='sun'
psi0={{-0.5,-0.1},{0.25,0.1},{-0.25,-0.75}}
s12=math.sqrt(0.308)
s13=math.sqrt(0.0234)
print("a=",a)
io.write(string.format("a=%f\n",a))
FLAVS=3
for i=1,FLAVS do 
io.write(string.format("%i=%f+i%f\n",i,psi0[i][1],psi0[i][2]))
end
--a 0.1
--the_sun_is=
function set_name()
  return string.format("%s_a%1.1e_b%1.1e_E%1.1e_tol%1.1e.dat",model,a,b,E,tol)
end
out=set_name()
io.write(string.format("out='%s'\n",out))
