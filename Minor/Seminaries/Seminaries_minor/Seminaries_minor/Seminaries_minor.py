import numpy as np
import sympy as sp
<<<<<<< HEAD

=======
import scipy as sip
>>>>>>> 11d9c2dffe39db0a3f04a1c838ba68a75f4ae955


#1

L = np.array([[2, 3, -1], 
              [3, -2, 1], 
              [5, 1, 0]])

psi = np.array([2, 5, 0]).reshape((3,1))
G = np.hstack((L, psi))
print('G', G, sep='=')

#2

A = np.arange(0, 20).reshape((4,5))
print('A', A, sep='=')
A = A.T
print('A', A, sep='=')
A = np.vstack((A, A*2))
print('A', A, sep='=')
f = A[: , 0]
l = A[: , -1]
A = np.hstack((A, f.reshape((f.size,1)), l.reshape((l.size,1))))
print('A', A, sep='=')

#3

def print_col(A,  n = 0):
    A = np.array(A)
    if 0 <= n and n  <= A[0,  :].size:
        print(A[: , n].reshape((A[:, 0].size,1)))
    else:
        print("Mistake")

[print_col(A, i) for i in range(A[0,  :].size)]
print_col(A)

#3*

def concate(A, n=0):
    if 0 <= n and n  <= A[0,  :].size:
        h = A[: , n].reshape((A[:, 0].size,1))
        print(A, h, np.hstack((A, h)), sep='\n')
        return np.hstack((A, h.reshape((f.size,1))))
    else:
        print("Mistake")

_ = concate(A, 2)

#4

def n_concate(A, *args):
    ns = args
    for n in ns:
        if 0 <= n and n  <= A[1,  :].size:
            l = A[n , :]
            A = np.vstack((A, l))
        else:
            print("Mistake")
            return None
    return A
print("a) ",  n_concate(L))
print("b) ",  n_concate(L, 1))
print("b) ",  n_concate(L, 0, 2))

#5






