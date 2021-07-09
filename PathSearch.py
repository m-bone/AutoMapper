##############################################################################
# Developed by: Matthew Bone
# Last Updated: 16/06/2021
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
# This is the core mapping code for solving molecule maps using a path search
# approach. Requires two molecule files and a few user provided parameters.
# Designed to work with molecule files created with LammpsToMolecule or 
# LammpsToMoleculePartial as bonding, edge and delete atom values are scraped 
# from the header of molecule files.
##############################################################################

import os
import logging

from LammpsSearchFuncs import element_atomID_dict
from AtomObjectBuilder import build_atom_objects, compare_symmetric_atoms
from QueueFuncs import Queue, queue_edge_atoms, queue_bond_atoms, run_queue

def map_delete_atoms(preDeleteAtoms, postDeleteAtoms, mappedIDList):
    # If delete atoms provided, add them to the mappedIDList. No purpose to including them in the queue
    if preDeleteAtoms is not None:
        assert postDeleteAtoms is not None, 'Delete atoms found in pre-bond file but not in post-bond file.'
        assert len(preDeleteAtoms) == len(postDeleteAtoms), 'Pre-bond and post-bond files have different numbers of delete atoms.'
        for index, preAtom in enumerate(preDeleteAtoms):
            mappedIDList.append([preAtom, postDeleteAtoms[index]])
            logging.debug(f'Pre: {preAtom}, Post: {postDeleteAtoms[index]} found with specified delete atom')

def get_missing_atom_objects(missingAtomList, atomObjectDict):
    missingAtomObjects = []
    for atom in missingAtomList:
        atomObject = atomObjectDict[atom]
        missingAtomObjects.append(atomObject)

    return missingAtomObjects

def map_missing_atoms(missingPreAtomObjects, missingPostAtomObjects, mappedIDList, queue):
    missingCheckCounter = 1
    while missingCheckCounter < 4 and len(missingPostAtomObjects) > 0:
        mappedPreAtomIndex = []
        for preIndex, preAtom in enumerate(missingPreAtomObjects):
            missingPostAtomElements = [atom.element for atom in missingPostAtomObjects]
            elementOccurence = missingPostAtomElements.count(preAtom.element)

            if elementOccurence == 0:
                print(f"Couldn't find a match in post missing atom for {preAtom.atomID}")

            elif elementOccurence == 1:
                postIndex = missingPostAtomElements.index(preAtom.element)
                logging.debug(f'Pre: {preAtom.atomID}, Post: {missingPostAtomObjects[postIndex].atomID} found with missing atoms single element occurence')
                match_missing(preAtom, postIndex, missingPostAtomObjects, mappedIDList, queue, preIndex, mappedPreAtomIndex)
                
            elif elementOccurence > 1:
                if preAtom.element == 'H':
                    postHydrogenIndexList = [index for index, element in enumerate(missingPostAtomElements) if element == 'H']
                    postIndex = postHydrogenIndexList.pop()
                    logging.debug(f'Pre: {preAtom.atomID}, Post: {missingPostAtomObjects[postIndex].atomID} found with missing atoms hydrogen symmetry inference')
                    match_missing(preAtom, postIndex, missingPostAtomObjects, mappedIDList, queue, preIndex, mappedPreAtomIndex)
                else:
                    potentialPostAtomObjects = [atomObject for atomObject in missingPostAtomObjects if atomObject.element == preAtom.element]
                    postIndex = compare_symmetric_atoms(potentialPostAtomObjects, preAtom, 'index')
                    if postIndex is not None:
                        match_missing(preAtom, postIndex, missingPostAtomObjects, mappedIDList, queue, preIndex, mappedPreAtomIndex)
                        logging.debug(f'The above atomID pair was found with missing atoms symmetry comparison')                                    

        # Refresh missingPreAtomObjects so that it doesn't print needless error messages on subsequent loops
        for index in sorted(mappedPreAtomIndex, reverse=True):
            del missingPreAtomObjects[index]        
        
        missingCheckCounter += 1

def match_missing(preAtom, postAtomMissingIndex, missingPostAtomObjects, mapList, queue, preIndex, mappedPreAtomIndex):
    # Get post atom
    postAtom = missingPostAtomObjects[postAtomMissingIndex]

    mapList.append([preAtom.atomID, postAtom.atomID])
    
    if preAtom.element != 'H':
        queue.add([[preAtom, postAtom]]) # This circumvents add_to_queue()

    missingPostAtomObjects.pop(postAtomMissingIndex)
    mappedPreAtomIndex.append(preIndex)

def update_missing_list(missingAtomList, mappedIDList, mapIndex):
    mappedAtoms = [pair[mapIndex] for pair in mappedIDList]
    # Update missingAtomList to remove atoms that have been matched
    newMissingAtomList = [atom for atom in missingAtomList if atom not in mappedAtoms]

    return newMissingAtomList

def verify_edge_atoms(preAtomObjectDict, preEdgeAtomDict, postAtomObjectDict, mappedIDList):
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
    for preEdgeAtom in list(preEdgeAtomDict.keys()):
        # Edge atom
        stopSearch = compare_atom_type(preEdgeAtom)

        if stopSearch:
            extendDistanceDict[preEdgeAtom] = 3
            continue

        # First neighbours
        preEdgeAtomObject = preAtomObjectDict[preEdgeAtom]
        for firstNeighbour in preEdgeAtomObject.firstNeighbourIDs:
            
            stopSearch = compare_atom_type(firstNeighbour)
            
            if stopSearch:
                extendDistanceDict[preEdgeAtom] = 2
                break

        # Second neighbours
        if stopSearch: continue # Prevents second neighbours running and overwriting the result from first neighbours
        for secondNeighbour in preEdgeAtomObject.secondNeighbourIDs:
            stopSearch = compare_atom_type(secondNeighbour)
            
            if stopSearch:
                extendDistanceDict[preEdgeAtom] = 1
                break

        # If no extension required
        if stopSearch == False:
            extendDistanceDict[preEdgeAtom] = 0

    return extendDistanceDict

def map_from_path(directory, preFileName, postFileName, elementsByType, debug):
    # Set log level
    if debug:
        logging.basicConfig(level='DEBUG')
    else:
        logging.basicConfig(level='INFO')

    # Build atomID to element dict
    os.chdir(directory)
    preElementDict = element_atomID_dict(preFileName, elementsByType)
    postElementDict = element_atomID_dict(postFileName, elementsByType)
    elementDictList = [preElementDict, postElementDict]

    # Generate atom class objects list
    preAtomObjectDict, preBondingAtoms, preEdgeAtomDict, preDeleteAtoms = build_atom_objects(preFileName, preElementDict)
    postAtomObjectDict, postBondingAtoms, postEdgeAtomDict, postDeleteAtoms = build_atom_objects(postFileName, postElementDict)

    # Assert the same number of atoms are in pre and post - maps have the same number of atoms in
    assert len(preAtomObjectDict) == len(postAtomObjectDict), f'Different numbers of atoms in pre- and post-bond files. Pre: {len(preAtomObjectDict)}, Post: {len(postAtomObjectDict)}'

    # Initialise lists
    missingPreAtomList = []
    missingPostAtomList = []
    mappedIDList = []

    # Initialise queue
    queue = Queue()

    # Populate queue with bonding atoms and update mappedIDList
    queue_bond_atoms(preAtomObjectDict, preBondingAtoms, postAtomObjectDict, postBondingAtoms, mappedIDList, queue)

    # Populate queue with unique edge atoms if present
    queue_edge_atoms(preAtomObjectDict, preEdgeAtomDict, postAtomObjectDict, postEdgeAtomDict, mappedIDList, queue)

    # Add delete atoms to the mappedIDList
    map_delete_atoms(preDeleteAtoms, postDeleteAtoms, mappedIDList)

    # Search through queue creating new maps based on all elements in a given path
    run_queue(queue, mappedIDList, preAtomObjectDict, postAtomObjectDict, missingPreAtomList, missingPostAtomList, elementDictList)

    # Update missingPreAtoms to check if the missing atom search loop is needed
    missingPreAtomList = update_missing_list(missingPreAtomList, mappedIDList, 0)

    # If missing pre atoms are present, map missing atoms and rerun the queue until success or timeout
    timeoutCounter = 1
    while len(missingPreAtomList) > 0 and timeoutCounter < 6:
        # Update missing atom lists
        missingPreAtomList = update_missing_list(missingPreAtomList, mappedIDList, 0)
        missingPostAtomList = update_missing_list(missingPostAtomList, mappedIDList, 1)

        # Get pre atom objects
        missingPreAtomObjects = get_missing_atom_objects(missingPreAtomList, preAtomObjectDict)

        # Add any post atoms that aren't in the map or already in the missing atoms
        mappedPostAtoms = [pair[1] for pair in mappedIDList]
        totalPostAtomList = list(postAtomObjectDict.keys())
        unfoundMissingPostAtoms = [atomID for atomID in totalPostAtomList if atomID not in mappedPostAtoms and atomID not in missingPostAtomList]
        missingPostAtomList.extend(unfoundMissingPostAtoms)

        # Get post atom objects
        missingPostAtomObjects = get_missing_atom_objects(missingPostAtomList, postAtomObjectDict)

        map_missing_atoms(missingPreAtomObjects, missingPostAtomObjects, mappedIDList, queue)

        # Refresh missingAtomLists
        missingPreAtomList = update_missing_list(missingPreAtomList, mappedIDList, 0)
        missingPostAtomList = update_missing_list(missingPostAtomList, mappedIDList, 1)

        # Rerun the queue based on atom pairs added to queue from missingAtoms
        run_queue(queue, mappedIDList, preAtomObjectDict, postAtomObjectDict, missingPreAtomList, missingPostAtomList, elementDictList)
        logging.debug(f'missingPreAtoms after loop {timeoutCounter}: {missingPreAtomList}') 

        timeoutCounter += 1

    if preEdgeAtomDict is not None:
        extendDistanceDict = verify_edge_atoms(preAtomObjectDict, preEdgeAtomDict, postAtomObjectDict, mappedIDList)
    else:
        extendDistanceDict = {'Dummy': 0}
    logging.debug(f'Extended Distance: {extendDistanceDict}')

    return mappedIDList, extendDistanceDict, preBondingAtoms, preEdgeAtomDict, preDeleteAtoms
