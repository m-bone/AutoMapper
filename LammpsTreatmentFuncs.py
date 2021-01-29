import re

# Function maybe moved to general function file later
def clean_data(lines):
    '''Remove stuff - done for neatness, to pass futher checks, and makes indexing sequential'''
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

def find_sections(lines):
    # Find index of section keywords - isalpha works as no spaces, newlines or punc in section keywords
    sectionIndexList = [lines.index(line) for line in lines if line.isalpha()]

    # Add end of file as last index
    sectionIndexList.append(len(lines))

    return sectionIndexList

# Get data
def get_data(sectionName, lines, sectionIndexList):
    # Checks that section name is existing in LAMMPS data
    try:
        startIndex = lines.index(sectionName)
    except ValueError:
        # If doesn't exist, return empty list that can be added as normal to main list later
        data = []
        return data

    endIndex = sectionIndexList[sectionIndexList.index(startIndex) + 1]
    
    data = lines[startIndex+1:endIndex] # +1 means sectionName doesn't get included
    data = [val.split() for val in data]

    return data

def get_coeff(coeffName, settingsData):
    # Inputs pre-split data
    # Return all lines that include coeffName in the [0] index
    coeffs = [line for line in settingsData if line[0] == coeffName]
    
    return coeffs

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