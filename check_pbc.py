import numpy as np


def check_pbc_applied(coordinates, bonds, box):

    pbc = []
    min_val = min(box[0][0:1])

    for i in range(len(bonds)):
        bond_diff = coordinates[bonds[i][0]-1, :] - coordinates[bonds[i][1]-1, :]
        bond_length = np.sqrt(np.dot(bond_diff, bond_diff))
        if bond_length > min_val/2.0:
            if np.abs(bond_diff[0]) > box[0][0]/2.0:
                if bond_diff[0] > 0.0:
                    pbc.append([True, "x"]) # true if bond extends over 2 periodic images. Check also which direction.
                else:
                    pbc.append([True, "-x"])

            elif np.abs(bond_diff[1]) > box[0][1]/2.0:
                if bond_diff[1] > 0.0:
                    pbc.append([True, "y"])
                else:
                    pbc.append([True, "-y"])

            elif np.abs(bond_diff[0]) > box[0][0]/2.0 and np.abs(bond_diff[1]) > box[0][1]/2.0:
                if bond_diff[0] > 0.0 and bond_diff[1] > 0.0:
                    pbc.append([True, "xy"])
                elif bond_diff[0] > 0.0 > bond_diff[1]:
                    pbc.append([True, "x-y"])
                elif bond_diff[0] < 0.0 < bond_diff[1]:
                    pbc.append([True, "-xy"])
                else:
                    pbc.append([True, "-x-y"])

        else:
            pbc.append([False, " "])

    return pbc
