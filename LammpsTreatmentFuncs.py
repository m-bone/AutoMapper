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
# A range of functions designed to pre-process, search, and save LAMMPS
# 'read_data' input files. Built for AutoMapper processing but can be applied
# to a wide range of problems
##############################################################################

import re # For clean_data, clean_settings
from collections import Counter # For refine_data
from operator import itemgetter # For refine_data
from natsort import natsorted # For refine_data

# Function maybe moved to general function file later
def clean_data(lines):
    # Remove blank lines 
    lines = [line for line in lines if line != '\n']

    # Remove comments - negative lookbehind means label comments in masses are kept e.g # C_3
    lines = [re.sub(r'(?<!\d\s\s)#(.*)', '' ,line) for line in lines]

    # Remove 5 spaces moltemplate/LAMMPS puts in front of header file
    lines = [re.sub(r'\s{5}', '' ,line) for line in lines]

    # Remove newline terminators
    lines = [re.sub(r'\n', '', line) for line in lines]

    # Remove empty strings in list caused by comments being removed
    lines = [line for line in lines if line != '']

    # Remove trailing whitespace
    lines = [re.sub(r'\s+$', '', line) for line in lines]

    return lines

def clean_settings(lines):
    # Remove newline terminators
    lines = [re.sub(r'\n', '', line) for line in lines]

    # Remove tabs
    lines = [re.sub(r'\t', '', line) for line in lines]
    
    # Replace multiple whitespaces with one
    lines = [re.sub(r'\s{2,}', ' ', line) for line in lines]

    return lines

def refine_data(data, searchIndex: list, IDset=None, newAtomIDs=None):
    '''
    Search multiple indices for matching atomID value.
    If match is found keep that list row in the data.
    Only output that row if the row appears len(searchIndex) times
    This means the data row contains valid atomIDs in all possible ID positions
    '''
    # If IDSet is not given, skip refining
    if IDset is None:
        return data

    # Convenience to convert int values to list
    if type(searchIndex) is not list:
        searchIndex = [searchIndex]

    # Adds rows that include a valid atom ID to list
    possibleValidData = []
    for atomID in IDset:
        for row in data:
            for index in searchIndex:
                if row[index] == atomID:
                    possibleValidData.append(row)

    # Lammps IDs found in above search
    # E.g. for angles, if the angle '25' appears 3 times, it means that the components for angle '25' are all in the IDSet
    possibleValidIDs = [val[0] for val in possibleValidData]
    IDCount = dict(Counter(possibleValidIDs))
    # If ID counter is the same size as the search index, ID is valid and gets added to data
    validIDs = [key for key in IDCount.keys() if IDCount[key] == len(searchIndex)]        
    validData = []
    for row in possibleValidData:
        if row[0] in validIDs:
            validData.append(row)
            # Stops duplicate IDs being refound in the future
            validIDIndex = validIDs.index(row[0])
            del validIDs[validIDIndex]

    # Re-sort validData by ID, use natsort as values are str not int
    validData = natsorted(validData, key=itemgetter(0))
    
    # If dict provided then use it to update the atom IDs from old to new
    for rowInd, row in enumerate(validData, start=1):
        if searchIndex != [0]: # Don't run this for the 'atoms' section
            row[0] = str(rowInd) # Reset the LAMMPS ID e.g. bond number
        for index in searchIndex:
            row[index] = newAtomIDs[row[index]]

    return validData

def add_section_keyword(sectionName, data):
    # Don't add keyword if list is empty - empty list means no section in file
    if len(data) == 0:
        return data

    # Add keyword name to start of list
    data.insert(0, '\n')
    data.insert(1, [sectionName])
    data.insert(2, '\n')

    return data

def save_text_file(fileName, dataSource):
    # Save to text file
    with open(fileName, 'w') as f:
        for item in dataSource:
            line = " ".join(item)
            if line != '\n':
                f.write("%s\n" % line)
            else:
                f.write(line)

# Create comment string with bond atoms and edge atoms
def format_comment(IDlist, comment):
    atomList = list(IDlist)
    atomList.insert(0, comment)
    atomString = [' '.join(atomList)] # Has to be list of lists to pass through later code

    return atomString

# Convert edge atom fingerprints from atom IDs to element strings
def edge_atom_fingerprint_strings(edgeAtomFingerprintDict, elementsByTypeDict):
    edgeElementFingerprintDict = {}
    for key, atomList in edgeAtomFingerprintDict.items():
        cutList = [elementsByTypeDict[atom] for atom in atomList]

        # Sort list alphabetically
        cutList = natsorted(cutList)

        # Join list to single string
        cutList = ''.join(cutList)

        edgeElementFingerprintDict[key] = cutList

    return edgeElementFingerprintDict
