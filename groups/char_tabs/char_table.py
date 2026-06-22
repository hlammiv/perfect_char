from numpy import *
from numpy.linalg import norm

mu1=(1-sqrt(5))/2
mu2=(1+sqrt(5))/2
w =(1+1j*sqrt(3))/2
ws = (1-1j*sqrt(3))/2
print(mu1,mu2,w,ws)
char_tab = zeros((17,17),dtype=complex)

#character 1
for i in range(17):
    char_tab[0][i]=1

#character 2
char_tab[1][0]=3
char_tab[1][1]=mu2
char_tab[1][2]=1
char_tab[1][3]=w
char_tab[1][4]=ws
char_tab[1][5]=-mu1*w
char_tab[1][6]=-mu1*ws
char_tab[1][7]=0
char_tab[1][8]=0
char_tab[1][9]=-ws
char_tab[1][10]=-w
char_tab[1][11]=mu1
char_tab[1][12]=-mu2*ws
char_tab[1][13]=-mu2*w
char_tab[1][14]=-1
char_tab[1][15]=-3*ws
char_tab[1][16]=-3*w

#character 3
for i in range(17):
        char_tab[2][i]=conjugate(char_tab[1][i])

#character 4
char_tab[3][0]=8
char_tab[3][1]=mu2
char_tab[3][2]=0
char_tab[3][3]=0
char_tab[3][4]=0
char_tab[3][5]=mu1
char_tab[3][6]=mu1
char_tab[3][7]=-1
char_tab[3][8]=-1
char_tab[3][9]=0
char_tab[3][10]=0
char_tab[3][11]=mu1
char_tab[3][12]=mu2
char_tab[3][13]=mu2
char_tab[3][14]=0
char_tab[3][15]=8
char_tab[3][16]=8

#character 5
char_tab[4][0]=6
char_tab[4][1]=1
char_tab[4][2]=0
char_tab[4][3]=-2*ws
char_tab[4][4]=-2*w
char_tab[4][5]=-ws
char_tab[4][6]=-w
char_tab[4][7]=0
char_tab[4][8]=0
char_tab[4][9]=0
char_tab[4][10]=0
char_tab[4][11]=1
char_tab[4][12]=-w
char_tab[4][13]=-ws
char_tab[4][14]=2
char_tab[4][15]=-6*w
char_tab[4][16]=-6*ws

#character 6
for i in range(17):
        char_tab[5][i]=conjugate(char_tab[4][i])

#character 7
char_tab[6][0]=15
char_tab[6][1]=0
char_tab[6][2]=-1
char_tab[6][3]=w
char_tab[6][4]=ws
char_tab[6][5]=0
char_tab[6][6]=0
char_tab[6][7]=0
char_tab[6][8]=0
char_tab[6][9]=ws
char_tab[6][10]=w
char_tab[6][11]=0
char_tab[6][12]=0
char_tab[6][13]=0
char_tab[6][14]=-1
char_tab[6][15]=-15*ws
char_tab[6][16]=-15*w

#character 8
for i in range(17):
        char_tab[7][i]=conjugate(char_tab[6][i])

#character 9
char_tab[8][0]=10
char_tab[8][1]=0
char_tab[8][2]=0
char_tab[8][3]=-2
char_tab[8][4]=-2
char_tab[8][5]=0
char_tab[8][6]=0
char_tab[8][7]=1
char_tab[8][8]=1
char_tab[8][9]=0
char_tab[8][10]=0
char_tab[8][11]=0
char_tab[8][12]=0
char_tab[8][13]=0
char_tab[8][14]=-2
char_tab[8][15]=10
char_tab[8][16]=10

#character 10
char_tab[9][0]=5
char_tab[9][1]=0
char_tab[9][2]=-1
char_tab[9][3]=1
char_tab[9][4]=1
char_tab[9][5]=0
char_tab[9][6]=0
char_tab[9][7]=-1
char_tab[9][8]=2
char_tab[9][9]=-1
char_tab[9][10]=-1
char_tab[9][11]=0
char_tab[9][12]=0
char_tab[9][13]=0
char_tab[9][14]=1
char_tab[9][15]=5
char_tab[9][16]=5

#character 11
char_tab[10][0]=5
char_tab[10][1]=0
char_tab[10][2]=-1
char_tab[10][3]=1
char_tab[10][4]=1
char_tab[10][5]=0
char_tab[10][6]=0
char_tab[10][7]=2
char_tab[10][8]=-1
char_tab[10][9]=-1
char_tab[10][10]=-1
char_tab[10][11]=0
char_tab[10][12]=0
char_tab[10][13]=0
char_tab[10][14]=1
char_tab[10][15]=5
char_tab[10][16]=5


#character 12
char_tab[11][0]=8
char_tab[11][1]=mu1
char_tab[11][2]=0
char_tab[11][3]=0
char_tab[11][4]=0
char_tab[11][5]=mu2
char_tab[11][6]=mu2
char_tab[11][7]=-1
char_tab[11][8]=-1
char_tab[11][9]=0
char_tab[11][10]=0
char_tab[11][11]=mu2
char_tab[11][12]=mu1
char_tab[11][13]=mu1
char_tab[11][14]=0
char_tab[11][15]=8
char_tab[11][16]=8

#character 13
char_tab[12][0]=9
char_tab[12][1]=-1
char_tab[12][2]=1
char_tab[12][3]=1
char_tab[12][4]=1
char_tab[12][5]=-1
char_tab[12][6]=-1
char_tab[12][7]=0
char_tab[12][8]=0
char_tab[12][9]=1
char_tab[12][10]=1
char_tab[12][11]=-1
char_tab[12][12]=-1
char_tab[12][13]=-1
char_tab[12][14]=1
char_tab[12][15]=9
char_tab[12][16]=9

#character 14
char_tab[13][0]=9
char_tab[13][1]=-1
char_tab[13][2]=1
char_tab[13][3]=-ws
char_tab[13][4]=-w
char_tab[13][5]=ws
char_tab[13][6]=w
char_tab[13][7]=0
char_tab[13][8]=0
char_tab[13][9]=-w
char_tab[13][10]=-ws
char_tab[13][11]=-1
char_tab[13][12]=w
char_tab[13][13]=ws
char_tab[13][14]=1
char_tab[13][15]=-9*w
char_tab[13][16]=-9*ws

#character 15
for i in range(17):
        char_tab[14][i]=conjugate(char_tab[13][i])


#character 16
char_tab[15][0]=3
char_tab[15][1]=mu1
char_tab[15][2]=1
char_tab[15][3]=w
char_tab[15][4]=ws
char_tab[15][5]=-mu2*w
char_tab[15][6]=-mu2*ws
char_tab[15][7]=0
char_tab[15][8]=0
char_tab[15][9]=-ws
char_tab[15][10]=-w
char_tab[15][11]=mu2
char_tab[15][12]=-mu1*ws
char_tab[15][13]=-mu1*w
char_tab[15][14]=-1
char_tab[15][15]=-3*ws
char_tab[15][16]=-3*w

#character 17
for i in range(17):
        char_tab[16][i]=conjugate(char_tab[15][i])


set_printoptions(linewidth=260,precision=2)
print(char_tab)

for i in range(17):
    for j in range(17):
        print(real(char_tab[i,j]),end=" ")
    print("")

for i in range(17):
    for j in range(17):
        print(imag(char_tab[i,j]),end=" ")
    print("")

