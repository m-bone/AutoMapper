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
# This is the primary mapping code that utilises a custom path search to
# determine similar atoms in pre- and post-bond molecule files. This is only
# called from MapProcessor; it can no longer be used to create a map independently.
##############################################################################

import os
import logging
import sys

from LammpsSearchFuncs import element_atomID_dict
from AtomObjectBuilder import build_atom_objects, compare_symmetric_atoms
from QueueFuncs import Queue, queue_bond_atoms, run_queue

def map_delete_atoms(preDeleteAtoms, postDeleteAtoms, mappedIDList):
    # If delete atoms provided, add them to the mappedIDList. No purpose to including them in the queue
    if preDeleteAtoms is not None:
        assert postDeleteAtoms is not None, 'Delete atoms found in pre-bond file but not in post-bond file.'
        assert len(preDeleteAtoms) == len(postDeleteAtoms), 'Pre-bond and post-bond files have different numbers of delete atoms.'
        for index, preAtom in enumerate(preDeleteAtoms):
            mappedIDList.append([preAtom, postDeleteAtoms[index]])
            logging.debug(f'Pre: {preAtom}, Post: {postDeleteAtoms[index]} found with user specified delete atom')

def get_missing_atom_objects(missingAtomList, atomObjectDict):
    missingAtomObjects = []
    for atom in missingAtomList:
        atomObject = atomObjectDict[atom]
        missingAtomObjects.append(atomObject)

    return missingAtomObjects

def map_missing_atoms(missingPreAtomObjects, missingPostAtomObjects, mappedIDList, queue, allowInference):
    missingCheckCounter = 1
    while missingCheckCounter < 4 and len(missingPostAtomObjects) > 0:
        mappedPreAtomIndex = []
        for preIndex, preAtom in enumerate(missingPreAtomObjects):
            missingPostAtomElements = [atom.element for atom in missingPostAtomObjects]
            elementOccurence = missingPostAtomElements.count(preAtom.element)

            if elementOccurence == 0:
                print(f"Error: Couldn't find a match in post missing atom for {preAtom.atomID}. Please try again or map the atom manually.")

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
                    postIndex = compare_symmetric_atoms(potentialPostAtomObjects, preAtom, 'index', allowInference=allowInference)
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

def map_from_path(directory, preFileName, postFileName, elementsByType, debug, preBondingAtoms, preDeleteAtoms, postBondingAtoms, postDeleteAtoms, createAtoms):
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
    preAtomObjectDict = build_atom_objects(preFileName, preElementDict, preBondingAtoms)
    postAtomObjectDict = build_atom_objects(postFileName, postElementDict, postBondingAtoms, createAtoms=createAtoms)

    # Assert the same number of atoms are in pre and post - maps have the same number of atoms in unless create atoms are included
    if createAtoms is None:
        assert len(preAtomObjectDict) == len(postAtomObjectDict), f'Different numbers of atoms in pre- and post-bond files. Pre: {len(preAtomObjectDict)}, Post: {len(postAtomObjectDict)}'

    # Initialise lists
    missingPreAtomList = []
    missingPostAtomList = []
    mappedIDList = []

    # Initialise queue
    queue = Queue()

    # Populate queue with bonding atoms and update mappedIDList
    queue_bond_atoms(preAtomObjectDict, preBondingAtoms, postAtomObjectDict, postBondingAtoms, mappedIDList, queue)

    # Add delete atoms to the mappedIDList
    map_delete_atoms(preDeleteAtoms, postDeleteAtoms, mappedIDList)

    # Search through queue creating new maps based on all elements in a given path
    run_queue(queue, mappedIDList, preAtomObjectDict, postAtomObjectDict, missingPreAtomList, missingPostAtomList, elementDictList)

    # Update missingPreAtoms to check if the missing atom search loop is needed
    missingPreAtomList = update_missing_list(missingPreAtomList, mappedIDList, 0)

    # If missing pre atoms are present, map missing atoms and rerun the queue until success or timeout
    timeoutCounter = 1
    # Disable inference for the first search
    inference = False
    while len(missingPreAtomList) > 0 and timeoutCounter < 11:
        # Update missing atom lists
        missingPreAtomList = update_missing_list(missingPreAtomList, mappedIDList, 0)
        missingPostAtomList = update_missing_list(missingPostAtomList, mappedIDList, 1)

        # Get pre atom objects and record length of missing
        missingPreAtomObjects = get_missing_atom_objects(missingPreAtomList, preAtomObjectDict)
        missingPreAtomCount = len(missingPostAtomList)

        # Add any post atoms that aren't in the map or already in the missing atoms
        mappedPostAtoms = [pair[1] for pair in mappedIDList]
        totalPostAtomList = list(postAtomObjectDict.keys())
        unfoundMissingPostAtoms = [atomID for atomID in totalPostAtomList if atomID not in mappedPostAtoms and atomID not in missingPostAtomList]
        missingPostAtomList.extend(unfoundMissingPostAtoms)

        # Get post atom objects
        missingPostAtomObjects = get_missing_atom_objects(missingPostAtomList, postAtomObjectDict)

        map_missing_atoms(missingPreAtomObjects, missingPostAtomObjects, mappedIDList, queue, inference)

        # Refresh missingAtomLists
        missingPreAtomList = update_missing_list(missingPreAtomList, mappedIDList, 0)
        missingPostAtomList = update_missing_list(missingPostAtomList, mappedIDList, 1)

        # Rerun the queue based on atom pairs added to queue from missingAtoms
        run_queue(queue, mappedIDList, preAtomObjectDict, postAtomObjectDict, missingPreAtomList, missingPostAtomList, elementDictList)
        logging.debug(f'missingPreAtoms after loop {timeoutCounter}: {missingPreAtomList}') 

        # Enable inference if no new missing atoms were solved this loop
        if missingPreAtomCount == len(missingPreAtomList):
            inference = True
        else: # Disable inference if missing atoms were solved as other atoms may now be found without inference
            inference = False

        timeoutCounter += 1

    if len(missingPreAtomList) > 0:
        print('Error: Missing Atom Search timed out. Atoms will be missing from the map. Please raise an issue on GitHub if the problem persists.')
        sys.exit()

    return mappedIDList
