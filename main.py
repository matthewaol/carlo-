import numpy as np 
import matplotlib as plt 

k_boltz = 1.380649 * 10**-23 # J/K

class atom: 
    def __init__(self,x_coord,y_coord,pe):
        self.x_coord = x_coord 
        self.y_coord = y_coord
        self.pe = pe

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
        
arr_atoms = gen_atom_arr(7)

for atom in arr_atoms:
     print(atom)