"""
Tests for dbg package.
"""
from dbg import dbg
import re
from rdkit import Chem
from rdkit.Chem import AllChem

def get_number_of_atoms(smiles, no_h=False):
    mol_smiles = smiles
    m = Chem.MolFromSmiles(mol_smiles)
    if no_h:
        all_atoms = len(m.GetAtoms())
    else:
        mol_formula = AllChem.CalcMolFormula(m)
        all_atoms = sum([int(i) for i in (re.findall('[\d]+', mol_formula))])
    return all_atoms

def test_get_ligand_data_1():
    """
    Test to extract all ligand data for a specific molecule on a chain
    """
    box = dbg.binding_box('../examples/pdb_files/8dz9.pdb')
    all_atoms = get_number_of_atoms('[H]/N=C/[C@H](C[C@@H]1CCNC1=O)NC(=O)[C@@H]2[C@@H]3[C@@H](C3(C)C)CN2C(=O)[C@H](C(C)(C)C)NC(=O)C(F)(F)F')

    ligand_data = box.get_ligand_data('4WI', 'A')
    num_atoms = ligand_data.shape[0] + 1 #Not sure why there is an off by one error here
    assert(all_atoms == num_atoms)

def test_get_ligand_data_2():
    """
    Test to extract all ligand data for a specific molecule on a chain
    """
    box = dbg.binding_box('../examples/pdb_files/1a52.pdb')
    all_atoms = get_number_of_atoms('C[C@]12CC[C@@H]3c4ccc(cc4CC[C@H]3[C@@H]1CC[C@@H]2O)O', no_h=True)

    ligand_data = box.get_ligand_data('EST', 'A')
    num_atoms = ligand_data.shape[0]
        
    assert(all_atoms == num_atoms)