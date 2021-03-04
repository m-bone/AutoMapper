import os
import numpy as np

from LammpsSearchFuncs import get_data, find_sections
from LammpsTreatmentFuncs import clean_data


def charge_matrix(directory, fileName):
    os.chdir(directory)

    # Load molecule file
    with open(fileName, 'r') as f:
        lines = f.readlines()

    # Clean data and get charge
    data = clean_data(lines)
    sections = find_sections(data)
    charge = get_data('Charges', data, sections)

    # Create float charge list
    chargeList = [float(val[1]) for val in charge]

    # Convert np matrix
    chargeMatrix = np.array(chargeList)

    return chargeMatrix


pre_charge_matrix = charge_matrix('/home/matt/Documents/Oct20-Dec20/Bonding_Test/DGEBA_DETDA/Reaction/', 'new_start_molecule.data')
print()