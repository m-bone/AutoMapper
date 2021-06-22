from PathSearch import map_from_path
from LammpsToMoleculePartial import lammps_to_molecule_partial

DIR = '/home/matt/Documents/Oct20-Dec20/Bonding_Test/Generic_PU/Reaction'
EBT = ['H', 'H', 'C', 'C', 'C', 'C', 'N', 'N', 'O', 'O', 'O', 'O']

# Initial molecule creation
preRenumberedAtomIDDict = lammps_to_molecule_partial(DIR, 'cleanedpre_reaction.data', 'testpre-', EBT, ['1', '36'])
postRenumberedAtomIDDict = lammps_to_molecule_partial(DIR, 'cleanedpost_reaction.data', 'testpost-', EBT, ['1', '36'])

# Initial map creation
initialMappedIDList, extendPreEdgeAtomDict = map_from_path(DIR, 'testpre-molecule.data', 'testpost-molecule.data', EBT, debug=False)

# Convert edge atom numbers to original number system
extendPostEdgeAtomDict = {}
for key, value in extendPreEdgeAtomDict.items():
    for atomPair in initialMappedIDList:
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

lammps_to_molecule_partial(DIR, 'cleanedpre_reaction.data', 'testpre-', EBT, ['1', '36'], extendEdgeAtomsDict=revertedExtendPreEdgeAtomDict)
lammps_to_molecule_partial(DIR, 'cleanedpost_reaction.data', 'testpost-', EBT, ['1', '36'], extendEdgeAtomsDict=revertedExtendPostEdgeAtomDict)

map_from_path(DIR, 'testpre-molecule.data', 'testpost-molecule.data', EBT, debug=False, mappedIDList=[])

# InitialMappedIDList needs reverting to original number system and then converting to the new cutdown number system for the extended edge atom molecules
# New error message for map and tidy/build this into AutoMapper main wrap function
# Add edge atom symmetry test now that "Edge_Atom_Symmetry" test case is working
# Create more complicated test cases for a molecule with a type change 1 away from the edge atom and one that requires 3rd neighbours to differentiate