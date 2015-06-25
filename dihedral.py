'''Calculates diehdral angles from Jaguar output files
Doug Crandell
9/9/2014'''

import numpy as np
import sys
import os

def main():
    if len(sys.argv) < 6:
        print "Not enough parameters specified"
    else:
        coords = get_coords(sys.argv[1])
        points_array = make_array(sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],coords)
        print round(dihedral(points_array),3)

def dihedral(p):
     #Calculate the vectors between the points
     b = p[:-1] - p[1:]
     b[0] *= -1 #"Flip the fist vector so that eclipsing vectors have dihedral =0
     # Use dot product to find the components of b1 and b3 that are not
     # perpendicular to b2. Subtract those components. The resulting vectors
     # lie in parallel planes.
     v = np.array( [ v - (v.dot(b[1])/b[1].dot(b[1])) * b[1] for v in [b[0], b[2]] ] )
     # Normalize vectors
     v /= np.sqrt(np.einsum('...i,...i', v, v)).reshape(-1,1)
     b1 = b[1] / np.linalg.norm(b[1])
     x = np.dot(v[0], v[1])
     m = np.cross(v[0], b1)
     y = np.dot(m, v[1])
     return np.degrees(np.arctan2( y, x ))

def dihedral_alternate(p):
    b = p[:-1] - p[1:]
    b[0] *= -1
    v = np.array( [np.cross(v,b[1]) for v in [b[0], b[2]] ] )
    # Normalize vectors
    v /= np.sqrt(np.einsum('...i,...i', v, v)).reshape(-1,1)
    return np.degrees(np.arccos( v[0].dot(v[1]) ))

def make_array(atom1,atom2,atom3,atom4,coords):
    '''Make an array of the xyz coordinates''''
    points = np.array([
                [ coords[atom1][0],coords[atom1][1],coords[atom1][2]],
                [ coords[atom2][0],coords[atom2][1],coords[atom2][2]],
                [ coords[atom3][0],coords[atom3][1],coords[atom3][2]],
                [ coords[atom4][0],coords[atom4][1],coords[atom4][2]]
                ])
    return points
    

def get_coords(file):
    '''Return a dictionary of the atoms and their cartesian coordinates'''
    data = open_file(file)
    line_number = 0
    #Find section containing the last geometry
    for line in data.split('\n'):
        if line.strip().startswith("Input geometry:"):
            start_geom = line_number
        if line.strip().startswith("final geometry:"): #If G.O. get the final geometry
            start_geom = line_number
        if line.strip().startswith("principal moments of inertia:"):
            end_geom = line_number
        line_number += 1
    geometry = data.split('\n')[start_geom+3:end_geom-1]
    atom_dict = {}
    for atom in geometry:
        atom = filter(None,atom.strip().split(' '))
        atom_dict[atom[0]] = [float(atom[1]),float(atom[2]),float(atom[3])]
    return atom_dict


def open_file(file):
    '''Read in the text of an output file and return it'''
    if len(file.split('.')) > 1:
        if file.split('.')[1] == 'out':
            with open(file,'r') as f:
                data = f.read()
    else:
        try:
            outfile = os.getcwd() +'/' + file + '.out' #Get output file based on user input of calc_id
            with open (outfile, 'r') as f:
                data = f.read()
        except:
            print "The file does not exist!"
    return data


if __name__ == "__main__":
    main()
