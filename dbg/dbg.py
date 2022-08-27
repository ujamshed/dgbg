import pandas as pd
import re
import mdtraj as md
import nglview
import numpy as np
import math
import os

from IPython.display import display, Image

pd.options.mode.chained_assignment = None  # default='warn'

class DBGException(Exception):
    pass

class binding_box():
    
    def __init__(self, pdb_id_file):
        
        # Sanitize file
        filename, ext = os.path.splitext(pdb_id_file)
        if (ext != ".pdb"):
            raise DBGException("Incorrect file format. Only use PDB files.")
        else:
            self.pdb_id = pdb_id_file
    
    def get_ligand_data(self, ligand_code, chain_id):
        list_of_values = []

        with open(self.pdb_id) as file:
            lines = file.readlines()
            pattern = "(^ATOM.*|^HETATM.*)"    
            for line in lines:
                if (re.search(pattern, line) != None):
                    test = re.search(r'^HETATM[\d]+', line)
                    if (test != None):
                        new_match = re.match(r'(^HETATM)([\d]+)', test.group(0))
                        new_string = new_match.groups()[0] + " " + new_match.groups()[1]
                        new_line = re.sub(r'^HETATM[\d]+', new_string, line)
                        list_of_values.append(new_line)
                    else:
                        list_of_values.append(line)

        df = pd.DataFrame(list_of_values)
        output = df[0].str.split(expand=True).rename({0: 'TYPE', 1: 'SERIAL_NUM', 2: 'ELEMENT', 3: 'AA', 4: 'CHAIN_ID', 5:'RESIDUE_NUM', 6: 'X_COOR', 7: 'Y_COOR', 8: 'Z_COOR', 9: 'OCCUPANCY', 10: 'TEMP_FACTOR', 11: 'ELEMENT'}, axis=1)
        self.output = output[(output['AA']==ligand_code) & (output['CHAIN_ID']==chain_id)]
        
        self.output['X_COOR'] = self.output['X_COOR'].apply(pd.to_numeric, errors='coerce')
        self.output['Y_COOR'] = self.output['Y_COOR'].apply(pd.to_numeric, errors='coerce')
        self.output['Z_COOR'] = self.output['Z_COOR'].apply(pd.to_numeric, errors='coerce')
        
        
        return self.output
    
    def create_box(self, spacing, padding):
        #add spacing and padding into the function parameters
        
        self.spacing = spacing
        
        min_x = self.output['X_COOR'].min()
        max_x = self.output['X_COOR'].max()
            
        min_y = self.output['Y_COOR'].min()
        max_y = self.output['Y_COOR'].max()
            
        min_z = self.output['Z_COOR'].min()
        max_z = self.output['Z_COOR'].max()
        
        x_distance = max_x - min_x
        y_distance = max_y - min_y
        z_distance = max_z - min_z
        
        self.distances = np.array([x_distance, y_distance, z_distance])
        
        self.box_size = self.distances // self.spacing
        
        self.xyz = self.output.loc[:, ['X_COOR', 'Y_COOR', 'Z_COOR']]
        
        self.centroid = np.array([self.output['X_COOR'].mean(), self.output['Y_COOR'].mean(), self.output['Z_COOR'].mean()]) 
            
        return self.centroid, (self.box_size + padding)
    
    def show_result(self, centroid, box_size):
        
        
        protein_mdtraj = md.load_pdb(self.pdb_id)
        v1 = nglview.show_mdtraj(protein_mdtraj)
        
        # Box_size / 2 (for positive and negative directions)
        box_size = ((box_size) / 2) * self.spacing
        
        # Create box points
        pt1 = [centroid[0]-box_size[0], centroid[1]-box_size[1], centroid[2]-box_size[2]]
        pt2 = [centroid[0]-box_size[0], centroid[1]+box_size[1], centroid[2]-box_size[2]]
        pt3 = [centroid[0]+box_size[0], centroid[1]-box_size[1], centroid[2]-box_size[2]]
        pt4 = [centroid[0]+box_size[0], centroid[1]+box_size[1], centroid[2]-box_size[2]]
        pt5 = [centroid[0]-box_size[0], centroid[1]-box_size[1], centroid[2]+box_size[2]]
        pt6 = [centroid[0]-box_size[0], centroid[1]+box_size[1], centroid[2]+box_size[2]]
        pt7 = [centroid[0]+box_size[0], centroid[1]-box_size[1], centroid[2]+box_size[2]]
        pt8 = [centroid[0]+box_size[0], centroid[1]+box_size[1], centroid[2]+box_size[2]]
        
        # Create centroid
        v1.shape.add_sphere(list(self.centroid), [1, 2, 2], 0.5, 'Centroid')

        # Create box
        v1.shape.add_cylinder(pt1, pt2, [ 1, 1, 0 ], 0.1 )
        v1.shape.add_cylinder(pt1, pt3, [ 1, 1, 0 ], 0.1 )
        v1.shape.add_cylinder(pt1, pt5, [ 1, 1, 0 ], 0.1 )
        v1.shape.add_cylinder(pt8, pt4, [ 1, 1, 0 ], 0.1 )
        v1.shape.add_cylinder(pt8, pt7, [ 1, 1, 0 ], 0.1 )
        v1.shape.add_cylinder(pt8, pt6, [ 1, 1, 0 ], 0.1 )

        v1.shape.add_cylinder(pt5, pt6, [ 1, 1, 0 ], 0.1 )
        v1.shape.add_cylinder(pt5, pt7, [ 1, 1, 0 ], 0.1 )
        v1.shape.add_cylinder(pt3, pt7, [ 1, 1, 0 ], 0.1 )
        v1.shape.add_cylinder(pt2, pt6, [ 1, 1, 0 ], 0.1 )
        v1.shape.add_cylinder(pt4, pt2, [ 1, 1, 0 ], 0.1 )
        v1.shape.add_cylinder(pt4, pt3, [ 1, 1, 0 ], 0.1 )
        
        return v1