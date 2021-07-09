from collections import deque
from copy import deepcopy
import os
import contextlib
from natsort import natsorted

from PathSearch import map_from_path
from LammpsToMolecule import lammps_to_molecule
from LammpsToMoleculePartial import lammps_to_molecule_partial
from LammpsTreatmentFuncs import save_text_file
from LammpsSearchFuncs import element_atomID_dict
from AtomObjectBuilder import build_atom_objects

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
        lammps_to_molecule(directory, preDataFileName, preMoleculeFileName, preBondingAtoms, deleteAtoms=preDeleteAtoms)
    
    with restore_dir():
        lammps_to_molecule(directory, postDataFileName, postMoleculeFileName, postBondingAtoms, deleteAtoms=postDeleteAtoms)

    # Initial map creation
    with restore_dir():
        mappedIDList, extendPreEdgeAtomDict, commentPreBondingAtoms, commentPreEdgeAtomDict, commentPreDeleteAtoms  = map_from_path(directory, preMoleculeFileName, postMoleculeFileName, elementsByType, debug=debug)

    # Cut map down to smallest possible partial structure
    # Key considerations - watch atom type changes close to edges; if it's in pre it's in post
    os.chdir(directory)
    preElementDict = element_atomID_dict(preMoleculeFileName, elementsByType)
    postElementDict = element_atomID_dict(postMoleculeFileName, elementsByType)

    preAtomObjectDict, preBondingAtoms, preEdgeAtomDict, preDeleteAtoms = build_atom_objects(preMoleculeFileName, preElementDict)
    postAtomObjectDict, postBondingAtoms, postEdgeAtomDict, postDeleteAtoms = build_atom_objects(postMoleculeFileName, postElementDict) 

    # Determine if bonding atom is part of a cycle, and if so what atoms make up the cycle and their neighbours 
    prePreservedAtomIDs = is_cyclic(preAtomObjectDict, preBondingAtoms)
    postPreservedAtomIDs = is_cyclic(postAtomObjectDict, postBondingAtoms)
    
    def is_ring_opening(prePreservedAtomIDs, postPreservedAtomIDs, mappedIDList):
        '''Determine if a reaction is a ring opening polymerisation'''
        def find_mapped_pair(preAtom):
            for pair in mappedIDList:
                if pair[0] == preAtom:
                    return pair[1]

        for preBondingAtom, prePreservedIDSet in prePreservedAtomIDs.items():
            
            # If preBondingAtom is cyclic, get the post bonding atom
            if prePreservedIDSet is not None:
                postBondingAtom = find_mapped_pair(preBondingAtom)

                # If post bonding atom is not cyclic, ring opening polymerisation is presumed
                if postPreservedAtomIDs[postBondingAtom] is None:
                    print('ITS A RING OPENER!')
                    # Return needs to be cut of dict with bondingAtom keys and set values

    is_ring_opening(prePreservedAtomIDs, postPreservedAtomIDs, mappedIDList)

    # Returns mappedIDList for other functions to use e.g. testing
    # return mappedIDList

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

def bfs(graph, startAtom, endAtom):
        # Adapted from https://stackoverflow.com/questions/8922060/how-to-trace-the-path-in-a-breadth-first-search
        # List to track if an atomID has already been seen
        discovered = {key: False for key in graph.keys()}

        discovered[startAtom] = True

        # Break link between start atom and target atom - stops search going backwards
        newGraph = deepcopy(graph)
        newGraph[startAtom].remove(endAtom)

        # Iterate through paths whilst keeping a record of all paths
        queue = []

        queue.append([startAtom])

        while queue:
            path = queue.pop(0)
            
            # Get latest path element
            node = path[-1]

            if node == endAtom:
                return path

            for neighbour in newGraph.get(node, []):
                # Prevents path getting stuck in a loop
                if discovered[neighbour]:
                    continue
                
                # Increase the path by next neighbour and add to queue
                discovered[neighbour] = True
                newPath = list(path)
                newPath.append(neighbour)
                queue.append(newPath)
        
        # If here then no path was found
        return None

def is_cyclic(atomObjectDict, bondingAtoms):
        # Create dictionary of adjacent bonds - IMPROVEMENT: Remove H from adjacent bonds since they can't go anywhere
        moleculeGraph = {atom.atomID: atom.firstNeighbourIDs for atom in atomObjectDict.values()}
        preservedAtomIDs = {} # With respect to each bonding atom

        for bondingAtom in bondingAtoms:
            # Get starting neighbours
            startNeighbours = atomObjectDict[bondingAtom].firstNeighbourIDs

            # Setup preservedAtomIDs
            preservedAtomIDs[bondingAtom] = None
            
            # Iterate through neighbours until a cycle is found
            for startAtom in startNeighbours:
                cyclicPath = bfs(moleculeGraph, startAtom, bondingAtom)

                if cyclicPath is not None:
                    print(f'Cycle found: {[cyclicPath]}. Started with {startAtom}')
                    break
            
                # If no path is found then the bonding atom is not cyclic
                print(f'No paths found for start atom {startAtom}')

            # Found path will be converted to a set of atomIDs that need to be preserved
            if cyclicPath is not None:
                preservedIDsSet = set()
                for atomID in cyclicPath:
                    preservedIDsSet.add(atomID)
                    preservedIDsSet.update(atomObjectDict[atomID].firstNeighbourIDs)
                
                preservedAtomIDs[bondingAtom] = preservedIDsSet

        return preservedAtomIDs  

# InitialMappedIDList needs reverting to original number system and then converting to the new cutdown number system for the extended edge atom molecules
# New error message for map and tidy/build this into AutoMapper main wrap function
# Add edge atom symmetry test now that "Edge_Atom_Symmetry" test case is working
# Create more complicated test cases for a molecule with a type change 1 away from the edge atom and one that requires 3rd neighbours to differentiate

# Validation Checks:
# No repeated IDs / No IDs unassigned
# Ambiguous groups maintained pre and post aside from moved atoms
# Check that atom assigned from missingatoms is bound to the atoms that it expects to be bound to