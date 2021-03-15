import os
from natsort import natsorted
from LammpsTreatmentFuncs import clean_data
from LammpsSearchFuncs import get_data, find_sections
from BondDistanceMatrix import breadth_first_search, fill_graph, calc_bond_length, calc_path_distance, get_bond_path

def bond_powers(directory, fileName, bondingAtoms):
    os.chdir(directory)

    # Load molecule file
    with open(fileName, 'r') as f:
        lines = f.readlines()

    # Clean data and get coords and bonds
    data = clean_data(lines)
    sections = find_sections(data)
    coords = get_data('Coords', data, sections)

    coordDict = {row[0]: row[1:] for row in coords}

    atomIDs = [row[0] for row in coords]
    bonds = get_data('Bonds', data, sections)

    # Calculate bond lengths and build in a dictionary
    bondsLengthList = [calc_bond_length(bond, coordDict) for bond in bonds]
    bondLengthDict = {bond[0]: bond[1] for bond in bondsLengthList}

    # Edit bondLengthDict to make bonding atoms bond length zero
    sortBondingAtoms = natsorted(bondingAtoms)
    bondingAtomsBondID = []
    for bond in bonds:
        if bond[2] == sortBondingAtoms[0] and bond[3] == sortBondingAtoms[1]:
            bondingAtomsBondID.append(bond[0])
    if len(bondingAtomsBondID) > 0: # Stops this happening in pre-molecule cases before bond occurs
        bondLengthDict[bondingAtomsBondID[0]] = 0.0

    # Create graph of bonding atom pairs for bond path search
    moleculeGraph = fill_graph(atomIDs, bonds) 