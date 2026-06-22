import numpy as np
import math as math
from numpy.linalg import norm

def equal(A, B):
        return norm(A-B) < 1e-5

read_dictionary = np.load('result_class_element.npy',allow_pickle='TRUE').item()
def dag(u):
    return np.transpose(u.conjugate())

def tr(u):
    return np.trace(u)

def mult(a,b):
    return np.matmul(a,b)

def chi_2_m1(u):
    return 1/2*(tr(u)**2*tr(dag(u))+tr(mult(u,u))*tr(dag(u)))-tr(u)

def chi_3_1(u):
    return 1/8*(tr(u)**4+2*tr(mult(u,u))*tr(u)**2-tr(mult(u,u))**2-2*tr(mult(u,mult(u,mult(u,u)))))

u = read_dictionary['5'][0]

elems=list(read_dictionary.values())
els = [item for sublist in elems for item in sublist]

print(len(els))

#for A in els:
#    for B in els:
#        M = A @ B
#        for i in range(len(els)):
#            if equal(M, els[i]):
#                print(i, end=' ')
#    print('')

for A in els:
    for clas, Us in read_dictionary.items():  # for name, age in dictionary.iteritems():  (for Python 2.x)
        for i in range(len(Us)):
            if equal(Us[i],A):
                print(int(clas)-1,end=' ')
    print('')

