import sys
import random
import numpy 
bv=float(sys.argv[1])
b=-bv/3 #needed to get original normalization of couplings
mask=numpy.ones(17)
#mask=numpy.zeros(17)

#g is the full 35th order calculation, gl only goes to 34
g=[]
gl=[]

g.append(+0.*b**0+0*b**1+-0.05555555555555525*b**2+-0.0015432098765432445*b**3+-0.0030864197530864387*b**4)
g.append(+0.*b**0+0.16666666666666652*b**1+0*b**2+0.009259259259259231*b**3+-0.0007716049382716084*b**4)
g.append(+0.*b**0+0.16666666666666652*b**1+0*b**2+0.009259259259259231*b**3+-0.0007716049382716084*b**4)
g.append(+0.*b**0+0*b**1+-0.05555555555555558*b**2+-0.0030864197530864057*b**3+-0.007716049382716035*b**4)
g.append(+0.*b**0+0*b**1+0*b**2+0.009259259259259266*b**3+-0.00025720164609054075*b**4)
g.append(+0.*b**0+0*b**1+0*b**2+0.009259259259259266*b**3+-0.00025720164609054075*b**4)
g.append(+0.*b**0+0*b**1+0*b**2+0.004629629629629629*b**3+-0.0010288065843621387*b**4)
g.append(+0.*b**0+0*b**1+0*b**2+0.004629629629629629*b**3+-0.0010288065843621387*b**4)
g.append(+0.*b**0+0*b**1+0*b**2+-2.7755575615628914e-17*b**3+-0.0046296296296296224*b**4)
g.append(+0.*b**0+0*b**1+0*b**2+0*b**3+-0.0015432098765432029*b**4)
g.append(+0.*b**0+0*b**1+0*b**2+0*b**3+-0.0015432098765432029*b**4)
g.append(+0.*b**0+0*b**1+0*b**2+0*b**3+-0.0015432098765432029*b**4)
g.append(+0.*b**0+0*b**1+0*b**2+0*b**3+-0.0015432098765432029*b**4)
g.append(+0.*b**0+0*b**1+0*b**2+0*b**3+-0.000257201646090539*b**4)
g.append(+0.*b**0+0*b**1+0*b**2+0*b**3+-0.000257201646090539*b**4)
g.append(+0.*b**0+0*b**1+0*b**2+0*b**3+0*b**4)
g.append(+0.*b**0+0*b**1+0*b**2+0*b**3+0*b**4)
gl.append(+0.*b**0+0*b**1+-0.05555555555555525*b**2+-0.0015432098765432445*b**3)
gl.append(+0.*b**0+0.16666666666666652*b**1+0*b**2+0.009259259259259231*b**3)
gl.append(+0.*b**0+0.16666666666666652*b**1+0*b**2+0.009259259259259231*b**3)
gl.append(+0.*b**0+0*b**1+-0.05555555555555558*b**2+-0.0030864197530864057*b**3)
gl.append(+0.*b**0+0*b**1+0*b**2+0.009259259259259266*b**3)
gl.append(+0.*b**0+0*b**1+0*b**2+0.009259259259259266*b**3)
gl.append(+0.*b**0+0*b**1+0*b**2+0.004629629629629629*b**3)
gl.append(+0.*b**0+0*b**1+0*b**2+0.004629629629629629*b**3)
gl.append(+0.*b**0+0*b**1+0*b**2+-2.7755575615628914e-17*b**3)
gl.append(+0.*b**0+0*b**1+0*b**2+0*b**3)
gl.append(+0.*b**0+0*b**1+0*b**2+0*b**3)
gl.append(+0.*b**0+0*b**1+0*b**2+0*b**3)
gl.append(+0.*b**0+0*b**1+0*b**2+0*b**3)
gl.append(+0.*b**0+0*b**1+0*b**2+0*b**3)
gl.append(+0.*b**0+0*b**1+0*b**2+0*b**3)
gl.append(+0.*b**0+0*b**1+0*b**2+0*b**3)
gl.append(+0.*b**0+0*b**1+0*b**2+0*b**3)

#print(g)
#print(gl)
#print([a_i - b_i for a_i, b_i in zip(g, gl)])


seed=random.randint(0,99999)
D=int(sys.argv[2])
nt=int(sys.argv[3])
nx=int(sys.argv[4])
group=str(sys.argv[5])
fol=str(sys.argv[6])

outfile=fol+"out_b"+str(bv)+"_g"+group+"_D"+str(D)+"_nt"+str(nt)+"_nx"+str(nx)+"_s"+str(seed)+".log"

command="./dym-mod-metro ./groups/"+group+"ct "+str(D)+" "+str(nt)+" "+str(nx)+" "

for i in range(len(g)):
#    if i == 1:
#        command+=str(bv)
#    else:
#        command+=str(0)
    command+=str(mask[i]*g[i])
    command+=" "

command+=str(seed)+" >> "+outfile

print(command)
#print(outfile)

