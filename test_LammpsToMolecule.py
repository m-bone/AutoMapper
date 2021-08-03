##############################################################################
# Developed by: Matthew Bone
# Last Updated: 30/07/2021
# Updated by: Matthew Bone
#
# Contact Details:
# Bristol Composites Institute (BCI)
# Department of Aerospace Engineering - University of Bristol
# Queen's Building - University Walk
# Bristol, BS8 1TR
# U.K.
# Email - matthew.bone@bristol.ac.uk
#
# File Description:
# A unit test file designed for PyTest. Tests the LammpsToMolecule function to
# ensure that it produces certain key outputs, where they should be.
##############################################################################

import os
from LammpsToMolecule import lammps_to_molecule
from LammpsTreatmentFuncs import clean_data
from LammpsSearchFuncs import get_data, find_sections

def test_lammps_to_molecule():
    path = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'Test_Cases/Methane_Ethane') # Allows for relative pathing in pytest
    lammps_to_molecule(path, 'cleanedpre_reaction.data', 'pre-', ['1', '6'])
    with open ('pre-molecule.data', 'r') as f:
        mol = f.readlines()
    
    mol = clean_data(mol)
    sectionIndex = find_sections(mol)
    types = get_data('Types', mol, sectionIndex)
    charges = get_data('Charges', mol, sectionIndex)
    coords = get_data('Coords', mol, sectionIndex)
    atomSectionLengths = len(types) == len(charges) and len(types) == len(coords)

    # Number of section keywords (+1 for end of file), Last section header is Angles, type of atom 6 is 2, atom based sections should be equal length, 12 angles in molecule
    checkValues = [len(sectionIndex), mol[sectionIndex[-2]], int(types[5][1]), atomSectionLengths, int(mol[2].split()[0])] 
    expected = [6, 'Angles', 2, True, 12]

    assert checkValues == expected