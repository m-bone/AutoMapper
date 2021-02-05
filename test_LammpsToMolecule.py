from LammpsToMolecule import lammps_to_molecule
from LammpsTreatmentFuncs import get_data, clean_data, find_sections

def test_lammps_to_molecule():
    lammps_to_molecule('/home/matt/Documents/Bond_React_Python/Test_Cases', 'cleanedpre-system.data', 'pre')
    with open ('premolecule.data', 'r') as f:
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