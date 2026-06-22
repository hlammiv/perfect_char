import sys
import random
import numpy 
bv=float(sys.argv[1])
b=-bv/3 #needed to get original normalization of couplings
mask=numpy.ones(17)
#mask=numpy.zeros(17)

#g is the full 24th order calculation, gl only goes to 23
g=[]
gl=[]

g.append(+0)
g.append(+0.*b**0+0.33333333333333337*b**1)
g.append(+0.*b**0+0.33333333333333337*b**1)
g.append(+0.*b**0+0*b**1)
g.append(+0.*b**0+0*b**1)
g.append(+0.*b**0+0*b**1)
g.append(+0.*b**0+0*b**1+0*b**2)
g.append(+0.*b**0+0*b**1+0*b**2)
g.append(+0.*b**0+0*b**1+0*b**2)
g.append(+0.*b**0+0*b**1+0*b**2)
g.append(+0.*b**0+0*b**1+0*b**2)
g.append(+0.*b**0+0*b**1+0*b**2)
g.append(+0.*b**0+0*b**1+0*b**2)
g.append(+0.*b**0+0*b**1+0*b**2)
g.append(+0.*b**0+0*b**1+0*b**2)
g.append(+0.*b**0+0*b**1+0*b**2)
g.append(+0.*b**0+0*b**1+0*b**2)

gl.append(+-0.5*b**0+0*b**1)
gl.append(+0.*b**0+0.33333333333333337*b**1)
gl.append(+0.*b**0+0.33333333333333337*b**1)
gl.append(+0.*b**0+0*b**1)
gl.append(+0.*b**0+0*b**1)
gl.append(+0.*b**0+0*b**1)
gl.append(+0.*b**0+0*b**1)
gl.append(+0.*b**0+0*b**1)
gl.append(+0.*b**0+0*b**1)
gl.append(+0.*b**0+0*b**1)
gl.append(+0.*b**0+0*b**1)
gl.append(+0.*b**0+0*b**1)
gl.append(+0.*b**0+0*b**1)
gl.append(+0.*b**0+0*b**1)
gl.append(+0.*b**0+0*b**1)
gl.append(+0.*b**0+0*b**1)
gl.append(+0.*b**0+0*b**1)

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
#        command+=str(b)
#    else:
#        command+=str(0)
    command+=str(mask[i]*g[i])
    command+=" "

command+=str(seed)+" >> "+outfile

print(command)
#print(outfile)

