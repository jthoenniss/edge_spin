from math import e
from os import stat_result
from numpy.lib.type_check import real
from scipy.linalg import expm
from scipy.linalg import eigvalsh
from scipy.sparse.linalg import eigsh #diagonalization via Lanczos for initial state
from scipy.linalg import block_diag
import  numpy as np
from scipy import sparse

np.set_printoptions(linewidth=np.nan, precision=5, suppress=True)

L = 7# number of sites of the bath

#Pauli matrices
sigma_x = np.array([[0,1],[1,0]])
sigma_y = np.array([[0,-1j],[1j,0]])
sigma_z = np.array([[1,0],[0,-1]])

#Parameters for Floquet evolution (can handle KIC as well as XY model)
Jx = 0
Jy = 0.31
g = np.pi / 4

ham_XY = Jx * np.kron(sigma_x, sigma_x) + Jy * np.kron(sigma_y, sigma_y)
kick_two_site = g * ( np.kron(sigma_z, np.identity(2)) + np.kron(np.identity(2), sigma_z)  ) 

F_odd = expm(1j * kick_two_site) @ expm(1j * ham_XY)
F_even = expm(1j * ham_XY)

#odd layer
layer_odd = np.identity(2)
for i in range ((L - 1)//2 ):
    layer_odd = np.kron(layer_odd, F_odd)
if L%2 == 0:
    layer_odd = np.kron(layer_odd,np.identity(2))

#even layer
layer_even = F_even
for i in range (L//2 - 1):
    layer_even = np.kron(layer_even, F_even)
if L%2 != 0:
    layer_even = np.kron(layer_even,np.identity(2))


F = layer_even @ layer_odd
#diagonalize Floquet operator

#compute time evolved Pauli-operators

#compute expecation value < Z(t) Z(0) >

