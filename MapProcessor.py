import os
import contextlib
from natsort import natsorted
from PathSearch import map_from_path
from LammpsToMoleculePartial import lammps_to_molecule_partial
from LammpsTreatmentFuncs import save_text_file

# DIR = '/home/matt/Documents/Oct20-Dec20/Bonding_Test/Generic_PU/Reaction'
# EBT = ['H', 'H', 'C', 'C', 'C', 'C', 'N', 'N', 'O', 'O', 'O', 'O']

def map_processor(directory, preDataFileName, postDataFileName, preMoleculeFileName, postMoleculeFileName, preBondingAtoms, postBondingAtoms, deleteAtoms, elementsByType, debug=False):
    # Split delete atoms list, if given
    if deleteAtoms is not None:
        assert len(deleteAtoms) % 2 == 0, 'Error: Different numbers of delete atom IDs for pre- and post-bond supplied.'
        deleteAtomIndex = len(deleteAtoms) // 2
        preDeleteAtoms = deleteAtoms[:deleteAtomIndex]
        postDeleteAtoms = deleteAtoms[deleteAtomIndex:]
    else:
        preDeleteAtoms = None
        postDeleteAtoms = None
    
    # Initial molecule creation
    with restore_dir(): # Allows for relative directory usage
        preRenumberedAtomIDDict = lammps_to_molecule_partial(directory, preDataFileName, preMoleculeFileName, elementsByType, preBondingAtoms, deleteAtoms=preDeleteAtoms)
    
    with restore_dir():
        postRenumberedAtomIDDict = lammps_to_molecule_partial(directory, postDataFileName, postMoleculeFileName, elementsByType, postBondingAtoms, deleteAtoms=postDeleteAtoms)

    # Initial map creation
    with restore_dir():
        mappedIDList, extendPreEdgeAtomDict, commentPreBondingAtoms, commentPreEdgeAtomDict, commentPreDeleteAtoms  = map_from_path(directory, preMoleculeFileName, postMoleculeFileName, elementsByType, debug=debug)

    # If edge atoms don't need to be extended, output map
    if list(extendPreEdgeAtomDict.values()).count(0) == len(extendPreEdgeAtomDict):
        os.chdir(directory)
        save_output(mappedIDList, commentPreBondingAtoms, commentPreEdgeAtomDict, commentPreDeleteAtoms)
    else:
        print('Note: Initial edge atoms were too close to atoms that changed type. Expanding partial molecule and remapping to fix.')

        # Convert edge atom numbers to original number system
        # Get the preAtomID for the edge atoms
        extendPostEdgeAtomDict = {}
        for key, value in extendPreEdgeAtomDict.items():
            for atomPair in mappedIDList:
                if atomPair[0] == key:
                    extendPostEdgeAtomDict[atomPair[1]] = value

        def revert_id_numbering(edgeAtomDict, renumberedAtomDict):
            revertedEdgeAtomDict = {}
            invertedRenumberedAtomDict = {value: key for key, value in renumberedAtomDict.items()}
            for key, value in edgeAtomDict.items():
                revertedEdgeAtomDict[invertedRenumberedAtomDict[key]] = value
            
            return revertedEdgeAtomDict

        revertedExtendPreEdgeAtomDict = revert_id_numbering(extendPreEdgeAtomDict, preRenumberedAtomIDDict)
        revertedExtendPostEdgeAtomDict = revert_id_numbering(extendPostEdgeAtomDict, postRenumberedAtomIDDict)


        with restore_dir():
            preRenumberedAtomIDDict = lammps_to_molecule_partial(directory, preDataFileName, preMoleculeFileName, elementsByType, preBondingAtoms, extendEdgeAtomsDict=revertedExtendPreEdgeAtomDict)
        
        with restore_dir():
            postRenumberedAtomIDDict = lammps_to_molecule_partial(directory, postDataFileName, postMoleculeFileName, elementsByType, postBondingAtoms, extendEdgeAtomsDict=revertedExtendPostEdgeAtomDict)

        with restore_dir():
            mappedIDList, _, commentPreBondingAtoms, commentPreEdgeAtomDict, commentPreDeleteAtoms = map_from_path(directory, preMoleculeFileName, postMoleculeFileName, elementsByType, debug=debug)
        
        # Output map
        os.chdir(directory)
        save_output(mappedIDList, commentPreBondingAtoms, commentPreEdgeAtomDict, commentPreDeleteAtoms)
    
    # Returns mappedIDList for other functions to use e.g. testing
    return mappedIDList

def output_map(mappedIDList, preBondingAtoms, preEdgeAtomDict, preDeleteAtoms):
    # Bonding atoms
    bondingAtoms = [['BondingIDs', '\n']]
    for atom in preBondingAtoms:
        bondingAtoms.extend([[atom]])
    bondingAtoms.extend(['\n'])

    # Delete Atoms
    deleteIDCount = []
    deleteAtoms = []
    if preDeleteAtoms is not None:
        deleteIDCount.extend([[str(len(preDeleteAtoms)) + ' deleteIDs']])
        deleteAtoms.extend([['DeleteIDs', '\n']])
        for atom in preDeleteAtoms:
            deleteAtoms.extend([[atom]])
        deleteAtoms.extend(['\n'])

    # Edge Atoms
    edgeIDCount = []
    edgeAtoms = []
    if preEdgeAtomDict is not None:
        edgeIDCount.extend([[str(len(preEdgeAtomDict)) + ' edgeIDs']])
        edgeAtoms.extend([['EdgeIDs', '\n']])
        for atom in preEdgeAtomDict.keys():
            edgeAtoms.extend([[atom]])
        edgeAtoms.extend(['\n'])
    edgeIDCount.extend('\n')

    # Equivalences
    equivalences = [['#This is an AutoMapper generated map\n'], [str(len(mappedIDList)) + ' equivalences']]
    equivalenceAtoms = [['Equivalences', '\n']]
    for atomPair in mappedIDList:
        equivalenceAtoms.extend([[atomPair[0] + '\t' + atomPair[1]]])

    # Output data
    output = []
    totalOutput = [equivalences, deleteIDCount, edgeIDCount, bondingAtoms, deleteAtoms, edgeAtoms, equivalenceAtoms]

    for section in totalOutput:
        output.extend(section)

    return output

def save_output(mappedIDList, preBondingAtoms, preEdgeAtomDict, preDeleteAtoms):
    # Order mappedIDList by preAtomID
    mappedIDList = natsorted(mappedIDList, key=lambda x: x[0])

    # Save data
    outputData = output_map(mappedIDList, preBondingAtoms, preEdgeAtomDict, preDeleteAtoms)
    save_text_file('automap.data', outputData)


@contextlib.contextmanager
def restore_dir():
    startDir = os.getcwd()
    try:
        yield
    finally:
        os.chdir(startDir)

# InitialMappedIDList needs reverting to original number system and then converting to the new cutdown number system for the extended edge atom molecules
# New error message for map and tidy/build this into AutoMapper main wrap function
# Add edge atom symmetry test now that "Edge_Atom_Symmetry" test case is working
# Create more complicated test cases for a molecule with a type change 1 away from the edge atom and one that requires 3rd neighbours to differentiate

# Validation Checks:
# No repeated IDs / No IDs unassigned
# Ambiguous groups maintained pre and post aside from moved atoms
# Check that atom assigned from missingatoms is bound to the atoms that it expects to be bound to