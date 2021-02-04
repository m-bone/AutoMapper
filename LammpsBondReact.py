import sys
from LammpsUnifiedCleaner import file_unifier

# Load input from system
directory = sys.argv[1]
coeffsFileName = sys.argv[2]
dataFiles = sys.argv[3:]
print(f'DataFiles List: {dataFiles}')

file_unifier(directory, coeffsFileName, dataFiles)
