import numpy as np 
import matplotlib as plt 

k_boltz = 1.380649 * 10**-23 # J/K

class atom: 
    def __init__(self, x_coord, y_coord, z_coord,x_p, y_p,z_p):
        self.x_coord = x_coord 
        self.y_coord = y_coord
        self.z_coord = z_coord
        self.x_p = x_p
        self.y_p = y_p
        self.z_p = z_p


def prob_pe_accept(old,new):
     expression = -(old - new) / k_boltz
     return np.exp(expression)

def gen_atom_arr(size): 
    arr = np.empty((size), dtype=object)
    x = 0 
    y = 0 
    for i in range(size):
            x+=1
            y+=1
            arr[i] = atom(x,y,0)
    return arr

def get_pe(atom1,atom2):
    ...
    
arr_atoms = gen_atom_arr(7)

for atom in arr_atoms:
     print(atom)