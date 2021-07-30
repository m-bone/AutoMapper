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
# The main molecule and map generating code for AutoMapper. map_processor will
# create pre- and post-bond molecule files from regular LAMMPS input files. It
# will then create a full map file before deciding if a partial structure can
# be created. If a partial structure is required, the molecule files and map
# will be cut down and renumbered accordingly. Finally, two molecule files and
# a map file named "automap.data" will be generated.
##############################################################################

import os
import logging
import contextlib
from natsort import natsorted
from copy import deepcopy

from PathSearch import map_from_path
from LammpsToMolecule import lammps_to_molecule
from LammpsTreatmentFuncs import save_text_file
from LammpsSearchFuncs import element_atomID_dict
from AtomObjectBuilder import build_atom_objects

def map_processor(directory, preDataFileName, postDataFileName, preMoleculeFileName, postMoleculeFileName, preBondingAtoms, postBondingAtoms, deleteAtoms, elementsByType, debug=False):
    # Set log level
    if debug:
        logging.basicConfig(level='DEBUG')
    else:
        logging.basicConfig(level='INFO')
    
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
        mappedIDList = map_from_path(directory, preMoleculeFileName, postMoleculeFileName, elementsByType, debug, preBondingAtoms, preDeleteAtoms, postBondingAtoms, postDeleteAtoms)

    # Cut map down to smallest possible partial structure
    with restore_dir():
        # Load data from just created files
        os.chdir(directory)
        preElementDict = element_atomID_dict(preMoleculeFileName, elementsByType)
        postElementDict = element_atomID_dict(postMoleculeFileName, elementsByType)

        preAtomObjectDict = build_atom_objects(preMoleculeFileName, preElementDict, preBondingAtoms)
        postAtomObjectDict = build_atom_objects(postMoleculeFileName, postElementDict, postBondingAtoms) 

    # Determine if bonding atom is part of a cycle, and if so what atoms make up the cycle and their neighbours 
    prePreservedAtomIDs = is_cyclic(preAtomObjectDict, preBondingAtoms, 'Pre-bond')
    postPreservedAtomIDs = is_cyclic(postAtomObjectDict, postBondingAtoms, 'Post-bond')
    
    # Determine atoms to be kept if reaction is a ring opening
    prePartialAtomsSet, postPartialAtomsSet = is_ring_opening(prePreservedAtomIDs, postPreservedAtomIDs, mappedIDList)

    # Keep atoms up to 4 bonds away from the bonding atoms
    prePartialAtomsSet = keep_all_neighbours(preAtomObjectDict, preBondingAtoms, prePartialAtomsSet)
    postPartialAtomsSet = keep_all_neighbours(postAtomObjectDict, postBondingAtoms, postPartialAtomsSet)

    # Keep delete atoms - IMPROVEMENT: Need to detect byproducts that are not deleted
    if preDeleteAtoms is not None:
        prePartialAtomsSet.update(preDeleteAtoms)
        postPartialAtomsSet.update(postDeleteAtoms)

    # Find initial pre-bond edge atoms
    preEdgeAtoms = find_edge_atoms(preAtomObjectDict, prePartialAtomsSet)

    # Check the edges aren't too close to atoms that change type
    preExtendEdgeDict = verify_edge_atoms(preEdgeAtoms, mappedIDList, preAtomObjectDict, postAtomObjectDict)
    mappedIDList, prePartialAtomsSet, postPartialAtomsSet = extend_edge_atoms(preExtendEdgeDict, mappedIDList, preAtomObjectDict, postAtomObjectDict, prePartialAtomsSet, postPartialAtomsSet)
    
    # Refind edge atoms after potential extension
    preEdgeAtoms = find_edge_atoms(preAtomObjectDict, prePartialAtomsSet)

    # Check for and get byproduct atoms that aren't deleteIDs
    postAtomByproducts = get_byproducts(postAtomObjectDict, postBondingAtoms)
    if postAtomByproducts is not None:
        logging.debug(f'Byproducts found. Byproducts are {postAtomByproducts}')
        postPartialAtomsSet.update(postAtomByproducts)

    # Order mappedIDList by preAtomID
    mappedIDList = natsorted(mappedIDList, key=lambda x: x[0])

    # Create empty partialMappedIDList to fill the return
    partialMappedIDList = []

    # Renumber map if the partial structure has a different length to the full structure
    # If they're equal then just output the map, no changes needed
    if len(prePartialAtomsSet) != len(preAtomObjectDict):
        logging.debug(f'Creating a partial map.')
        # Build a partial map and get the renumbering dictionaries
        mappedIDList, preRenumberdAtomDict, postRenumberedAtomDict, partialMappedIDList = create_partial_map(mappedIDList, prePartialAtomsSet, postPartialAtomsSet)

        # Renumber key features for molecule creation and output
        preBondingAtoms = renumber(preBondingAtoms, preRenumberdAtomDict)
        preDeleteAtoms = renumber(preDeleteAtoms, preRenumberdAtomDict)

        postBondingAtoms = renumber(postBondingAtoms, postRenumberedAtomDict)
        postDeleteAtoms = renumber(postDeleteAtoms, postRenumberedAtomDict)

        # Renumber edge atoms if any given
        if preEdgeAtoms is not None:
            preEdgeAtoms = renumber(preEdgeAtoms, preRenumberdAtomDict)

        # Rebuild molecule files with partial structure
        with restore_dir():
            lammps_to_molecule(directory, preDataFileName, preMoleculeFileName, preBondingAtoms, deleteAtoms=preDeleteAtoms, validIDSet=prePartialAtomsSet, renumberedAtomDict=preRenumberdAtomDict)

        with restore_dir():
            lammps_to_molecule(directory, postDataFileName, postMoleculeFileName, postBondingAtoms, deleteAtoms=postDeleteAtoms, validIDSet=postPartialAtomsSet, renumberedAtomDict=postRenumberedAtomDict)

    # Output the map file
    with restore_dir():
        os.chdir(directory)
        outputData = output_map(mappedIDList, preBondingAtoms, preEdgeAtoms, preDeleteAtoms)
        save_text_file('automap.data', outputData)

    # Returns mappedIDList for other functions to use e.g. testing
    return [mappedIDList, partialMappedIDList]

def output_map(mappedIDList, preBondingAtoms, preEdgeAtoms, preDeleteAtoms):
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
    if preEdgeAtoms is not None:
        edgeIDCount.extend([[str(len(preEdgeAtoms)) + ' edgeIDs']])
        edgeAtoms.extend([['EdgeIDs', '\n']])
        for atom in preEdgeAtoms:
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

# Utility for moving to a different os path and then returning to the original directory
@contextlib.contextmanager
def restore_dir():
    startDir = os.getcwd()
    try:
        yield
    finally:
        os.chdir(startDir)

def bfs(graph, startAtom, endAtom, breakLink=False):
    # Adapted from https://stackoverflow.com/questions/8922060/how-to-trace-the-path-in-a-breadth-first-search
    # List to track if an atomID has already been seen
    discovered = {key: False for key in graph.keys()}

    discovered[startAtom] = True

    newGraph = deepcopy(graph)
     # Break link between start atom and target atom if present - stops search going backwards when searching for cycles
    if breakLink:
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

def is_cyclic(atomObjectDict, bondingAtoms, reactionType):
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
            cyclicPath = bfs(moleculeGraph, startAtom, bondingAtom, breakLink=True)

            if cyclicPath is not None:
                logging.debug(f'Cycle found: {cyclicPath}. Started with {startAtom} for the {reactionType} reaction.')
                break

        # Found path will be converted to a set of atomIDs that need to be preserved
        if cyclicPath is not None:
            preservedIDsSet = set()
            for atomID in cyclicPath:
                preservedIDsSet.add(atomID)
                preservedIDsSet.update(atomObjectDict[atomID].firstNeighbourIDs)
            
            preservedAtomIDs[bondingAtom] = preservedIDsSet

    return preservedAtomIDs  

def find_mapped_pair(preAtom, mappedIDList):
    # Get the post atom from the map for a given preAtom
    for pair in mappedIDList:
        if pair[0] == preAtom:
            return pair[1]

def is_ring_opening(prePreservedAtomIDs, postPreservedAtomIDs, mappedIDList):
    '''
    Determine if a reaction is a ring opening polymerisation.
    Return two dicts of bonding atoms keys and preserved atom sets.
    '''

    # Init dicts for storing rings
    preCyclicAtomsSet = set()
    postCyclicAtomsSet = set()

    for preBondingAtom, prePreservedIDSet in prePreservedAtomIDs.items():
        
        # If preBondingAtom is cyclic (not None), get the post bonding atom
        if prePreservedIDSet is not None:
            postBondingAtom = find_mapped_pair(preBondingAtom, mappedIDList)

            # If post bonding atom is not cyclic, ring opening polymerisation is presumed
            if postPreservedAtomIDs[postBondingAtom] is None:
                logging.debug(f'Reaction has been determined as ring opening.')

                # Store pre bond data
                preCyclicAtomsSet.add(preBondingAtom)
                preCyclicAtomsSet.update(prePreservedIDSet)

                # Get post bond data from map and store
                # This is done as the postPreservedAtomsIDs for a ring opening reaction will be None
                postCyclicAtomsSet.add(postBondingAtom)
                for preCyclicAtom in prePreservedIDSet:
                    postCyclicAtom = find_mapped_pair(preCyclicAtom, mappedIDList)
                    postCyclicAtomsSet.add(postCyclicAtom)

    return preCyclicAtomsSet, postCyclicAtomsSet

def keep_all_neighbours(atomObjectDict, bondingAtoms, partialAtomSet):
    for bondingAtom in bondingAtoms:
        # Add bonding atom to the set in case it isn't there already
        partialAtomSet.add(bondingAtom)

        # Get atom object
        atomObject = atomObjectDict[bondingAtom]

        # Add neighbourIDs to the partial set
        partialAtomSet.update(atomObject.firstNeighbourIDs)
        partialAtomSet.update(atomObject.secondNeighbourIDs)
        partialAtomSet.update(atomObject.thirdNeighbourIDs)

    return partialAtomSet

def create_partial_map(mappedIDList, prePartialAtomsSet, postPartialAtomsSet):
    # Remove all the IDs that aren't in the pre and post partial atom sets
    partialMappedIDList = []
    for pair in mappedIDList:
        if pair[0] in prePartialAtomsSet and pair[1] in postPartialAtomsSet:
            partialMappedIDList.append(pair)
        
        # If something is present one set it must be present in the other set
        if pair[0] in prePartialAtomsSet and pair[1] not in postPartialAtomsSet:
            print(f'Warning: Pre atom {pair[0]} is present but post atom {pair[1]} missing')
        if pair[0] not in prePartialAtomsSet and pair[1] in postPartialAtomsSet:
            print(f'Warning: Pre atom {pair[0]} is missing but post atom {pair[1]} present')

        # Debug tools
        # if pair[0] not in prePartialAtomsSet:
        #     print(f'Pre atom {pair[0]} not in partial atoms')
        # if pair[1] not in postPartialAtomsSet:
        #     print(f'Post atom {pair[1]} not in partial atoms')

    preRenumberedAtomDict = {}
    postRenumberedAtomDict = {}
    renumberedMappedIDList = []
    # Simply renumber the pre and post atoms based on their position in the ID list
    for index, pair in enumerate(partialMappedIDList, start=1):
        preRenumberedAtomDict[pair[0]] = str(index)
        postRenumberedAtomDict[pair[1]] = str(index)
        renumberedMappedIDList.append([str(index), str(index)])

    # Assert that the same number of atoms is in pre and post dicts. Lammps will fail otherwise
    assert len(preRenumberedAtomDict) == len(postRenumberedAtomDict), 'Different numbers of atoms have been found in the pre- and post-bond partial structures.\n Please check your input files, raise an issue on Github if the problem persists.'

    return renumberedMappedIDList, preRenumberedAtomDict, postRenumberedAtomDict, partialMappedIDList

def renumber(inputList, renumberedAtomDict):
    '''
    Take list of numbers and conversion dict as input.
    Output list of converted numbers in same order. Return None if input is None
    '''
    if inputList is None:
        return None

    outputList = []
    for value in inputList:
        outputList.append(renumberedAtomDict[value])

    return outputList

def find_edge_atoms(atomObjectDict, partialAtomSet):
    edgeAtoms = []
    for atom in partialAtomSet:
        # Skip H atoms as they can't be edge atoms
        if atomObjectDict[atom].element == 'H':
            continue
        
        # Iterate through first neighbours and check they're all in the partial atom set
        for neighbour in atomObjectDict[atom].firstNeighbourIDs:
            # If one neighbour is not in the partial atom set, atom must be an edge
            if neighbour not in partialAtomSet:
                edgeAtoms.append(atom)
                break

    if len(edgeAtoms) > 0:
        return edgeAtoms
    else: # If no edge atoms in molecule then other functions expect None
        return None

def verify_edge_atoms(preEdgeAtoms, mappedIDList, preAtomObjectDict, postAtomObjectDict):
    # If no edge atoms are given, return an empty extend list
    if preEdgeAtoms is None:
        return {}

    # Convert mappedIDList to mappedIDDict
    mappedIDDict = {pair[0]: pair[1] for pair in mappedIDList}

    # Compare if pre and post atom types are the same, return True if they are not
    def compare_atom_type(preAtom):
        preAtomType = preAtomObjectDict[preAtom].atomType
        pairAtom = mappedIDDict[preAtom]
        postAtomType = postAtomObjectDict[pairAtom].atomType

        if preAtomType != postAtomType:
            return True
        else:
            return False

    # Check for atom type changes too close to edge atoms
    extendDistanceDict = {}
    for edgeAtom in preEdgeAtoms:
        # Edge atom
        stopSearch = compare_atom_type(edgeAtom)

        if stopSearch:
            extendDistanceDict[edgeAtom] = 3
            continue

        # First neighbours
        preEdgeAtomObject = preAtomObjectDict[edgeAtom]
        for firstNeighbour in preEdgeAtomObject.firstNeighbourIDs:
            
            stopSearch = compare_atom_type(firstNeighbour)
            
            if stopSearch:
                extendDistanceDict[edgeAtom] = 2
                break

        # Second neighbours
        if stopSearch: continue # Prevents second neighbours running and overwriting the result from first neighbours
        for secondNeighbour in preEdgeAtomObject.secondNeighbourIDs:
            stopSearch = compare_atom_type(secondNeighbour)
            
            if stopSearch:
                extendDistanceDict[edgeAtom] = 1
                break

    return extendDistanceDict

def extend_edge_atoms(extendEdgeDict, mappedIDList, preAtomObjectDict, postAtomObjectDict, prePartialAtomsSet, postPartialAtomsSet):
    # Output extended mappedIDList, pre and post partial atom sets. Rerun find_edge_atoms to get new edges.
    # If an edge is in this list, it at least needs to be extended by one
    additionalPreAtoms = []
    additionalPostAtoms = []

    for preEdge, extendDist in extendEdgeDict.items():
        # For pre-bond
        additionalPreAtoms.extend(preAtomObjectDict[preEdge].firstNeighbourIDs)

        # For post-bond
        postEdge = find_mapped_pair(preEdge, mappedIDList)
        additionalPostAtoms.extend(postAtomObjectDict[postEdge].firstNeighbourIDs)

        # When further away neighbours are required
        if extendDist == 2:
            additionalPreAtoms.extend(preAtomObjectDict[preEdge].secondNeighbourIDs)
            additionalPostAtoms.extend(postAtomObjectDict[postEdge].secondNeighbourIDs)
        
        if extendDist == 3:
            additionalPreAtoms.extend(preAtomObjectDict[preEdge].thirdNeighbourIDs)
            additionalPostAtoms.extend(postAtomObjectDict[postEdge].thirdNeighbourIDs)            
        
    # Expand the mappedIDList
    for preAtom in additionalPreAtoms:
        # Could be done using the additionalPostAtoms list but this is easier
        postAtom = find_mapped_pair(preAtom, mappedIDList)
        if [preAtom, postAtom] not in mappedIDList: # Prevents repeat mapped pairs being added
            mappedIDList.append([preAtom, postAtom])

    # Update the partial atom sets
    prePartialAtomsSet.update(additionalPreAtoms)
    postPartialAtomsSet.update(additionalPostAtoms)

    return mappedIDList, prePartialAtomsSet, postPartialAtomsSet

def get_byproducts(postAtomObjectDict, postBondingAtoms):
    # Determine if there is a path from each post atom to a bonding atom in the post structure
    # If no path is present then atom must be from a byproduct that's not being deleted
    byproducts = []
    targetBondingAtom = postBondingAtoms[0] # Only need one atom, as it will be bound to the other one
    moleculeGraph = {atom.atomID: atom.firstNeighbourIDs for atom in postAtomObjectDict.values()}
    for startAtom in postAtomObjectDict.keys():
        pathToBondingAtom = bfs(moleculeGraph, startAtom, targetBondingAtom, breakLink=False)

        if pathToBondingAtom is None:
            byproducts.append(startAtom)

    if len(byproducts) > 0:
        return byproducts
    else:
        return None
