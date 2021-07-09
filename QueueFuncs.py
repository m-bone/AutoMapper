import logging
from natsort import natsorted
from collections import Counter, deque

# Classes and functions for search
class Queue:
    def __init__(self):
        self.elements = deque()
    
    def empty(self) -> bool:
        return not self.elements
    
    def add(self, x):
        for pair in x:
            self.elements.append(pair)

    def get(self):
        return self.elements.popleft()


def add_to_queue(queue, queueAtoms, preAtomObjectDict, postAtomObjectDict):
    queueAtomObjects = []
    for pair in queueAtoms:
        preAtom = preAtomObjectDict[pair[0]]
        postAtom = postAtomObjectDict[pair[1]]
        queueAtomObjects.append([preAtom, postAtom])
    queue.add(queueAtomObjects)

def queue_bond_atoms(preAtomObjectDict, preBondingAtoms, postAtomObjectDict, postBondingAtoms, mappedIDList, queue):
    # Loop through bonding atoms, getting atom objects and adding them to queue and mapped list
    for index, preBondAtom in enumerate(preBondingAtoms):
        preAtomObject = preAtomObjectDict[preBondAtom]
        postAtomObject = postAtomObjectDict[postBondingAtoms[index]]
        queue.add([[preAtomObject, postAtomObject]])
        mappedIDList.append([preBondAtom, postBondingAtoms[index]])
        logging.debug(f'Pre: {preBondAtom}, Post: {postBondingAtoms[index]} found with user specified bond atom')

def queue_edge_atoms(preAtomObjectDict, preEdgeAtomDict, postAtomObjectDict, postEdgeAtomDict, mappedIDList, queue):
    # Skip function call if no edge atoms present
    if preEdgeAtomDict is None:
        return

    # Check edge atom fingerprints are the same
    assert natsorted(list(preEdgeAtomDict.values())) == natsorted(list(postEdgeAtomDict.values())), 'Pre and post edge atom fingerprints do not match.'
    
    # Find unique fingerprints 
    countEdgeAtomFingerprints = Counter(list(preEdgeAtomDict.values()))
    uniqueEdgeAtomFingerprints = [fingerprint for fingerprint in preEdgeAtomDict.values() if countEdgeAtomFingerprints[fingerprint] == 1]
    for fingerprint in uniqueEdgeAtomFingerprints:
        # Find index of fingerprint values, then get key at that index
        preEdgeAtomID = list(preEdgeAtomDict.keys())[list(preEdgeAtomDict.values()).index(fingerprint)]
        postEdgeAtomID = list(postEdgeAtomDict.keys())[list(postEdgeAtomDict.values()).index(fingerprint)]
        
        preAtomObject = preAtomObjectDict[preEdgeAtomID]
        postAtomObject = postAtomObjectDict[postEdgeAtomID]
        
        queue.add([[preAtomObject, postAtomObject]])
        mappedIDList.append([preEdgeAtomID, postEdgeAtomID])
        logging.debug(f'Pre: {preEdgeAtomID}, Post: {postEdgeAtomID} found with specified edge atom')

def run_queue(queue, mappedIDList, preAtomObjectDict, postAtomObjectDict, missingPreAtomList, missingPostAtomList, elementDictList):
    while not queue.empty():
        currentAtoms = queue.get()
        for mainIndex, atom in enumerate(currentAtoms):
            atom.check_mapped(mappedIDList, mainIndex, elementDictList[mainIndex])
        
        newMap, missingPreAtoms, missingPostAtoms, queueAtoms = currentAtoms[0].map_elements(currentAtoms[1], preAtomObjectDict, postAtomObjectDict)

        # Convert queue atoms to atom class objects and add to queue
        add_to_queue(queue, queueAtoms, preAtomObjectDict, postAtomObjectDict)

        # Extend missing lists
        missingPreAtomList.extend(missingPreAtoms)
        missingPostAtomList.extend(missingPostAtoms)

        # Add new pairs to mapped ID list
        mappedIDList.extend(newMap)