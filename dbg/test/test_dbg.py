"""
Tests for dbg package.
"""
from dbg import dbg
import pytest
import re

def test_pdb_file_sanitation_1():
    """
    Test to determine if the file is the correct format (pdb).
    """
    with pytest.raises(Exception) as e_info:
        box = dbg.binding_box('../examples/example_1.ipynb')
    return

def test_pdb_file_sanitation_2():
    """
    Test to determine if the file is the correct format (pdb).
    """
    box = dbg.binding_box('../examples/pdb_files/8dz9.pdb')
    return

def test_pdb_file_sanitation_3():
    """
    Test to determine if the file is the correct format (pdb).
    """
    with pytest.raises(Exception) as e_info:
        box = dbg.binding_box('../examples')
    return

def test_get_ligand_data_1():
    """
    Test to extract all ligand data for a specific molecule on a chain
    """
    box = dbg.binding_box('../examples/pdb_files/8dz9.pdb')
    all_atoms = 68

    ligand_data = box.get_ligand_data('4WI', 'A')
    num_atoms = ligand_data.shape[0]
    assert(all_atoms == num_atoms)

def test_get_ligand_data_2():
    """
    Test to extract all ligand data for a specific molecule on a chain
    """
    box = dbg.binding_box('../examples/pdb_files/1a52.pdb')
    # Hydrogens are missing in this pdb file.
    all_atoms = 20

    ligand_data = box.get_ligand_data('EST', 'A')
    num_atoms = ligand_data.shape[0]
        
    assert(all_atoms == num_atoms)

def test_get_ligand_data_3():
    """
    Test to extract all ligand data for a specific molecule on a chain
    """
    box = dbg.binding_box('../examples/pdb_files/7uty.pdb')
    all_atoms = 51

    ligand_data = box.get_ligand_data('OFR', 'A')
    num_atoms = ligand_data.shape[0]
    
    assert(all_atoms == num_atoms)

def test_get_ligand_data_4():
    """
    Test to extract all ligand data for a specific molecule on a chain
    """
    box = dbg.binding_box('../examples/pdb_files/7uty.pdb')
    # There are 2 EDO molecules on chain A.
    all_atoms = 20

    ligand_data = box.get_ligand_data('EDO', 'A')
    num_atoms = ligand_data.shape[0]
    
    assert(all_atoms == num_atoms)

def test_get_ligand_data_5():
    """
    Test to extract all ligand data for a specific molecule on a chain
    """
    box = dbg.binding_box('../examples/pdb_files/102d.pdb')
    # Hydrogens missing from this pdb file.
    all_atoms = 23

    ligand_data = box.get_ligand_data('TNT', 'B')
    num_atoms = ligand_data.shape[0]
    
    assert(all_atoms == num_atoms)

def test_get_ligand_data_6():
    """
    Test to extract all ligand data for a specific molecule on a chain
    """
    box = dbg.binding_box('../examples/pdb_files/1err.pdb')
    # Hydrogens missing from this pdb file.
    all_atoms = 34

    ligand_data = box.get_ligand_data('RAL', 'A')
    num_atoms = ligand_data.shape[0]
    
    assert(all_atoms == num_atoms)

def test_get_ligand_data_7():
    """
    Test to extract all ligand data for a specific molecule on a chain
    """
    box = dbg.binding_box('../examples/pdb_files/8dpi.pdb')
    all_atoms = 27

    ligand_data = box.get_ligand_data('T4U', 'A')
    num_atoms = ligand_data.shape[0]
    
    assert(all_atoms == num_atoms)

def test_get_ligand_data_8():
    """
    Test to extract all ligand data for a specific molecule on a chain
    """
    box = dbg.binding_box('../examples/pdb_files/8dpi.pdb')
    all_atoms = 74

    ligand_data = box.get_ligand_data('CLR', 'A')
    num_atoms = ligand_data.shape[0]
    
    assert(all_atoms == num_atoms)