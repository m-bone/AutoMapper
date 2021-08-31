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
# Contains key atom object creation and manipulation tools. The Atom object
# class and builder are the building blocks for map creation.
##############################################################################

import logging
from collections import Counter

from LammpsSearchFuncs import get_data, find_sections, get_neighbours, get_additional_neighbours
from LammpsTreatmentFuncs import clean_data

def build_atom_objects(fileName, elementDict, bondingAtoms):
    # Load molecule file
    with open(fileName, 'r') as f:
        lines = f.readlines()

    # Clean data and get coords and bonds
    data = clean_data(lines)
    sections = find_sections(data)
    types = get_data('Types', data, sections)

    atomIDs = [row[0] for row in types]
    bonds = get_data('Bonds', data, sections)
    
    # Build neighbours dict
    neighboursDict = get_neighbours(atomIDs, bonds, bondingAtoms)

    def get_elements(neighbourIDs, elementDict):
        return [elementDict[atomID]for atomID in neighbourIDs]

    atomObjectDict = {}
    for index, atomID in enumerate(atomIDs):
        atomType = types[index][1]

        # Establish all neighbours
        neighbours = neighboursDict[atomID]
        secondNeighbours = get_additional_neighbours(neighboursDict, atomID, neighbours, bondingAtoms)
        thirdNeighbours = get_additional_neighbours(neighboursDict, atomID, secondNeighbours, bondingAtoms)

        neighbourElements = get_elements(neighbours, elementDict)
        secondNeighbourElements = get_elements(secondNeighbours, elementDict)
        thirdNeighbourElements = get_elements(thirdNeighbours, elementDict)

        # Check if atom is a bonding atom, return boolean
        if atomID in bondingAtoms:
            bondingAtom = True
        else:
            bondingAtom = False

        atom = Atom(atomID, atomType, elementDict[atomID], bondingAtom, neighbours, secondNeighbours, thirdNeighbours, neighbourElements, secondNeighbourElements, thirdNeighbourElements)
        atomObjectDict[atomID] = atom
    
    return atomObjectDict

def compare_symmetric_atoms(postNeighbourAtomObjectList, preNeighbourAtom, outputType, allowInference=True):
    # Neighbour comparison - no inference
    def compare_neighbours(neighbourLevel):
        neighbourComparison = [getattr(atomObject, neighbourLevel) for atomObject in postNeighbourAtomObjectList]
        neighbourFingerprint = [''.join(sorted(elements)) for elements in neighbourComparison] # sorted to get alphabetical fingerprints

        # Remove duplicate fingerprints
        countFingerprints = Counter(neighbourFingerprint)
        tuppledFingerprints = [(index, fingerprint) for index, fingerprint in enumerate(neighbourFingerprint) if countFingerprints[fingerprint] == 1]

        # If any of the fingerprints are empty (i.e. the atom has no Xneighbours) return None
        for _, fingerprint in tuppledFingerprints:
            if fingerprint == '':
                return None

        # Any of the potential post neighbours matches the pre atom fingerprint, return the post neighbour
        for index, fingerprint in tuppledFingerprints:
            if ''.join(sorted(getattr(preNeighbourAtom, neighbourLevel))) == fingerprint:
                logging.debug(f'Pre: {preNeighbourAtom.atomID}, Post: {postNeighbourAtomObjectList[index].atomID} found with {neighbourLevel}')
                if outputType == 'index':
                    return index
                elif outputType == 'atomID':
                    return postNeighbourAtomObjectList[index].atomID
                else:
                    print('Invalid output type specified for compare_symmetric_atoms')

    # First neighbour comparison
    symmetryResult = compare_neighbours('firstNeighbourElements')

    # Second neighbour comparison
    if symmetryResult is None:
        symmetryResult = compare_neighbours('secondNeighbourElements')

    # Third neighbour comparison
    if symmetryResult is None:
        symmetryResult = compare_neighbours('thirdNeighbourElements')

    # If it makes it through all these, guess assignment and warn user about this
    if symmetryResult is not None:
        return symmetryResult
    else:
        if allowInference: # Only if inference is turned on
            # Find all potential choices by breaking the postNeighbourAtomList down into atoms that match the preAtom element
            possibleChoices = []
            for index, postNeighbourAtom in enumerate(postNeighbourAtomObjectList):
                if postNeighbourAtom.element == preNeighbourAtom.element:
                    possibleChoices.append((index, postNeighbourAtom.atomID))

            # Let the user know that an inference has been made     
            logging.debug(f'Pre: {preNeighbourAtom.atomID}, Post: {possibleChoices[0][1]} found with symmetry inference')
            print(
                f'Note: Pre-bond atomID {preNeighbourAtom.atomID} has been assigned by inference to post-bond atomID {possibleChoices[0][1]}. The potential choices were {[atom[1] for atom in possibleChoices]}. Please check this is correct.'
            )
            if outputType == 'index':
                    return possibleChoices[0][0]
            elif outputType == 'atomID':
                return possibleChoices[0][1]
            else:
                print('Invalid output type specified for compare_symmetric_atoms')


class Atom():
    def __init__(self, atomID, atomType, element, bondingAtom, neighbourIDs, secondNeighbourIDs, thirdNeighbourIDs, neighbourElements, secondNeighbourElements, thirdNeighbourElements):
        self.atomID = atomID
        self.atomType = atomType
        self.element = element
        self.bondingAtom = bondingAtom

        # Neighbours
        self.mappedNeighbourIDs = neighbourIDs # This is changed according to mapping
        self.firstNeighbourIDs = neighbourIDs.copy() # This is fixed throughout mapping process
        self.secondNeighbourIDs = secondNeighbourIDs
        self.thirdNeighbourIDs = thirdNeighbourIDs

        self.mappedNeighbourElements = neighbourElements # This is changed according to mapping
        self.firstNeighbourElements = neighbourElements.copy() # This is fixed throughout mapping process
        self.secondNeighbourElements = secondNeighbourElements
        self.thirdNeighbourElements = thirdNeighbourElements

    def check_mapped(self, mappedIDs, searchIndex, elementDict):
        """Update neighbourIDs.

        Updates neighbourIDs by removing IDs that have already been mapped.
        This will be called before all neighbour mapping attempts to stop atoms
        being mapped multiple times.

        Args:
            mappedIDs: The total list of mappedIDs at this point in the mapping. This
                will contain pre- and post-atomIDs
            searchIndex: Determines whether to use pre- or post-atomIDs

        Returns:
            Updates existing class variable self.NeighbourIDs
        """
        searchIndexMappedIDs = [row[searchIndex] for row in mappedIDs]
        
        self.mappedNeighbourIDs = [ID for ID in self.mappedNeighbourIDs if ID not in searchIndexMappedIDs]
        self.mappedNeighbourElements = [elementDict[atomID]for atomID in self.mappedNeighbourIDs]


    def map_elements(self, atomObject, preAtomObjectDict, postAtomObjectDict):
        """Map preAtom IDs to postAtomIDs by comparing neighbouring element symbols.

        Compares the occurence of string chemical element symbols of a preAtom's neighbours 
        to the known postAtom's neighbours. Creates a list of new mappedIDs and any missing IDs 
        based on the neighbours of the known pre and postAtom pair.
        Relies on compare_symmetric atoms to handle non-H atoms with >1 occurence.

        Args:
            self: This is a preAtom that has succesfully been mapped
            atomObject: The known postAtom that has already been mapped to the preAtom
            preAtomObjectDict: A dictionary of all preAtoms in the molecule
            postAtomObjectDict: A dictionary of all the postAtoms in the molecule

        Returns:
            A partial mappedIDlist, partial missing pre and postAtom lists and additional atoms for the queue
        """

        # Output variables
        mapList = []
        missingPreAtoms = []
        queueAtoms = []

        def allowed_maps(preAtom, postAtom):
            '''Populate missing atoms prematurely instead of making misleading
            maps based on element occurence'''
            # Checks if elements appear the same number of times in pre and post atoms
            # If they don't, mapping is not allowed to take place and atoms are moved to missing lists
            preElementOccurences = Counter(preAtom.mappedNeighbourElements)
            postElementOccurences = Counter(postAtom.mappedNeighbourElements)

            allowedMapDict = {}
            for element, count in preElementOccurences.items():
                if count == postElementOccurences[element]:
                    allowedMapDict[element] = True
                else:
                    allowedMapDict[element] = False
                
            # Force all H to be be True as hydrogen can be mapped by inference
            if 'H' in allowedMapDict:
                allowedMapDict['H'] = True

            return allowedMapDict

        allowedMapDict = allowed_maps(self, atomObject)

        # Match Function
        def matchNeighbour(preAtom, postAtom, preAtomIndex, postAtomIndex, mapList, queueList):
            # Append pre and post atomIDs to map
            mapList.append([preAtom.mappedNeighbourIDs[preAtomIndex], postAtom.mappedNeighbourIDs[postAtomIndex]])
            
            # Add all non-hydrogen atom atomIDs to queue
            if preAtom.mappedNeighbourElements[preAtomIndex] != 'H':
                queueList.append([preAtom.mappedNeighbourIDs[preAtomIndex], postAtom.mappedNeighbourIDs[postAtomIndex]])

            # Remove post atomID from mappedID and mappedElement atom object values
            postAtom.mappedNeighbourIDs.pop(postAtomIndex)
            postAtom.mappedNeighbourElements.pop(postAtomIndex)

        # Loop through neighbours for preAtom and compare to neighbours of postAtom
        for preIndex, neighbour in enumerate(self.mappedNeighbourElements):
            elementOccurence = atomObject.mappedNeighbourElements.count(neighbour)

            # Check if maps with the neighbour element are allowed, if not add current element to missing list
            if allowedMapDict[neighbour] == False:
                missingPreAtoms.append(self.mappedNeighbourIDs[preIndex])
                continue

            # If no match in post atom list it is a missingPreAtom
            if elementOccurence == 0:
                missingPreAtoms.append(self.mappedNeighbourIDs[preIndex])
            
            # Assign atomIDs if there is only one matching element - could this go wrong if an element moves and an identical element takes its place?
            elif elementOccurence == 1:
                postIndex = atomObject.mappedNeighbourElements.index(neighbour)
                logging.debug(f'Pre: {self.mappedNeighbourIDs[preIndex]}, Post: {atomObject.mappedNeighbourIDs[postIndex]} found with single element occurence')
                matchNeighbour(self, atomObject, preIndex, postIndex, mapList, queueAtoms)

            # More than one matching element requires additional considerations
            elif elementOccurence > 1:
                if neighbour == 'H': # H can be handled simply as all H are equivalent to each other in this case - ignores chirality
                    postHydrogenIndexList = [index for index, element in enumerate(atomObject.mappedNeighbourElements) if element == 'H']
                    postIndex = postHydrogenIndexList.pop()
                    logging.debug(f'Pre: {self.mappedNeighbourIDs[preIndex]}, Post: {atomObject.mappedNeighbourIDs[postIndex]} found with hydrogen symmetry inference')
                    matchNeighbour(self, atomObject, preIndex, postIndex, mapList, queueAtoms)
                    
                else:
                    # Get neighbour post atoms objects
                    postNeighbourIndices = [index for index, val in enumerate(atomObject.mappedNeighbourElements) if val == neighbour]
                    postNeighbourAtomIDs = [atomObject.mappedNeighbourIDs[i] for i in postNeighbourIndices]
                    postNeighbourAtomObjects = [postAtomObjectDict[atomID] for atomID in postNeighbourAtomIDs] 

                    # Get possible pre atom object
                    preNeighbourAtomObject = preAtomObjectDict[self.mappedNeighbourIDs[preIndex]]

                    # Find the post atom ID for the current pre atom
                    postNeighbourAtomID = compare_symmetric_atoms(postNeighbourAtomObjects, preNeighbourAtomObject, 'atomID')
                    if postNeighbourAtomID is not None:
                        postIndex = atomObject.mappedNeighbourIDs.index(postNeighbourAtomID)
                        matchNeighbour(self, atomObject, preIndex, postIndex, mapList, queueAtoms)
                    else:
                        # If no post atom found, add pre atom missing atom list
                        print(f'Could not find the symmetric pair for preAtom {self.mappedNeighbourIDs[preIndex]}')
                        missingPreAtoms.append(self.mappedNeighbourIDs[preIndex])                                             


        # Search mapList for missingPostAtoms
        mappedPostAtomList = [row[1] for row in mapList]
        missingPostAtoms = [neighbour for neighbour in atomObject.mappedNeighbourIDs if neighbour not in mappedPostAtomList]

        return mapList, missingPreAtoms, missingPostAtoms, queueAtoms
