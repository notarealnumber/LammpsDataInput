import numpy as np


def calc_bond_length(coord1, coord2):
    """
    Calculates the bond length between atom 1 at coord1 and atom 2
    at coord2.

    Bonds:
    at1---at2
    at1 => coord1
    at2 => coord2

    Angles:
    at1   at3
      \   /
       \ /
       at2
    at1, at3 => coord2
    at2 => coord1
    Dihedrals, impropers:
    at1   at3
      \   / \
       \ /   \
       at2   at4
    First 3 atoms:
    at1, at3 => coord2
    at2 => coord1
    Second part, atoms 3 and 4:
    at4 => coord2
    at3 => coord1

    Args:
        coord1 (list): The atom to be subtracted.
        coord2 (list): The other atom that forms the bond with coord1-atom.
    """

    bond_diff = [coord2[0] - coord1[0],
                 coord2[1] - coord1[1],
                 coord2[2] - coord1[2]]

    bond_length = np.sqrt(bond_diff[0]**2 + bond_diff[1]**2 + bond_diff[2]**2)

    return bond_length


def define_bonds_phi(at1, at2, at3, at4, bonds):

    bonds_involved = []
    for i in range(len(bonds)):
        if bonds[i] == [at1, at2]:
            bonds_involved.append([at1, at2])
        elif bonds[i] == [at1, at3]:
            bonds_involved.append([at1, at3])
        elif bonds[i] == [at1, at4]:
            bonds_involved.append([at1, at4])
        elif bonds[i] == [at2, at1]:
            bonds_involved.append([at2, at1])
        elif bonds[i] == [at2, at3]:
            bonds_involved.append([at2, at3])
        elif bonds[i] == [at2, at4]:
            bonds_involved.append([at2, at4])
        elif bonds[i] == [at3, at1]:
            bonds_involved.append([at3, at1])
        elif bonds[i] == [at3, at2]:
            bonds_involved.append([at3, at2])
        elif bonds[i] == [at3, at4]:
            bonds_involved.append([at3, at4])
        elif bonds[i] == [at4, at1]:
            bonds_involved.append([at4, at1])
        elif bonds[i] == [at4, at2]:
            bonds_involved.append([at4, at2])
        elif bonds[i] == [at4, at3]:
            bonds_involved.append([at4, at3])

    return bonds_involved


def get_new_index(coords, ind_unk, ind_knwn, current_cell, nx, ny, nat, box):

    ind_extended = 0
    for iadd in range(nx * ny):
        # if iadd == current_cell:
        #     continue
        success = False
        for box_y in range(-1, 2):
            for box_x in range(-1, 2):
                tempcoord = [
                    coords[ind_unk + iadd * nat - 1][3] + box[0, 0] * nx * box_x,
                    coords[ind_unk + iadd * nat - 1][4] + box[0, 1] * ny * box_y,
                    coords[ind_unk + iadd * nat - 1][5]
                ]
                check_length = calc_bond_length(coords[ind_knwn][3:6], tempcoord)
                if check_length < 5.0:
                    ind_extended = ind_unk + iadd * nat
                    success = True
                if success:
                    if ind_extended == 0:
                        print(box_x, box_y, tempcoord)
                    return ind_extended


def extend_coords(coords, elements, ff_type, charges, reprod, box):
    nat = len(coords)
    extended_coords = []
    nx, ny = reprod[0][:]

    print(nx*box[0, 0], ny*box[0, 1], box[0, 2])
    print(box[0, 0], box[0, 1], box[0, 2])

    for iy in range(reprod[0][1]):
        for ix in range(reprod[0][0]):
            for k in range(nat):
                extended_coords.append([elements[k], ff_type[k], charges[k],
                                        coords[k, 0] + box[0, 0] * ix,
                                        coords[k, 1] + box[0, 1] * iy,
                                        coords[k, 2]])

    return extended_coords


def extend_bonds(bonds, coords, box, reprod, nat, nbonds):
    """
    New version checks if a given bond distance is smaller than half the cell
    in a given direction. If yes, bond will be used. If not, the corresponding atom has
    to be searched.

    The first atom index in bonds will be set into the current ix, iy cell.
    That means that only the second atom of the psf-bond will be in a different
    cell.

    -----------------------------------
    |                |                |
    |                |                |
    |--at2      at1--|--at(2+nat*((ix+iy*nx)+1))
    |                |                |
    |                |                |
    | cell0          | cell(ix+iy*nat)|
    -----------------------------------
    """

    extended_bonds = []
    nx, ny = reprod[0][:]

    imgd = 0
    celled = 0

    for iy in range(ny):  # y-direction
        for ix in range(nx):  # x direction
            addind = nat * (ix + iy * nx)
            for nb in range(nbonds):
                ind_orig1 = bonds[nb][0]
                ind_orig2 = bonds[nb][1]

                bond1 = ind_orig1 + addind - 1
                bond2 = ind_orig2 + addind - 1

                bond_length = calc_bond_length(coords[bond1][3:6], coords[bond2][3:6])

                if bond_length < 5.0:
                    celled += 1
                    extended_bonds.append([[bond1 + 1, bond2 + 1],
                                           [coords[bond1][1].lower(),
                                            coords[bond2][1].lower()]])
                    continue
                else:
                    imgd += 1
                    new_ind = get_new_index(coords, ind_orig1, bond2, ix + iy * nx, nx, ny, nat, box)
                    extended_bonds.append([[new_ind, bond2 + 1],
                                           [coords[new_ind - 1][1].lower(),
                                            coords[bond2][1].lower()]])
                    continue

    print(imgd, celled, imgd + celled, len(extended_bonds), nx*ny*len(bonds))

    return extended_bonds


def extend_thetas(angles, coords, box, reprod, nat, nangles):
    """
    iangle1       iangle3
         \        /
          \      /
           \    /
          iangle2
    """

    extended_angles = []
    nx, ny = reprod[0][:]

    cntr1 = cntr2 = cntr3 = cntr4 = 0

    for iy in range(ny):  # y-direction
        for ix in range(nx):  # x direction
            addind = nat * (ix + iy * nx)
            for nang in range(nangles):
                # ind_orig saves the index as read from the psf. The atom index starts with 1.
                ind_orig1 = angles[nang][0]
                ind_orig2 = angles[nang][1]
                ind_orig3 = angles[nang][2]

                # -1 in order to get the correct element in coords.
                iangle1 = ind_orig1 + addind - 1
                iangle2 = ind_orig2 + addind - 1
                iangle3 = ind_orig3 + addind - 1

                bond_length1 = calc_bond_length(coords[iangle2][3:6], coords[iangle1][3:6])
                bond_length2 = calc_bond_length(coords[iangle2][3:6], coords[iangle3][3:6])

                if bond_length1 < 5.0 and bond_length2 < 5.0:
                    cntr1 += 1
                    extended_angles.append([[iangle1 + 1, iangle2 + 1, iangle3 + 1],
                                           [coords[iangle1][1].lower(),
                                            coords[iangle2][1].lower(),
                                            coords[iangle3][1].lower()]])
                    continue

                elif bond_length1 > 5.0 and bond_length2 > 5.0:
                    cntr4 += 1
                    new_ind1 = get_new_index(coords, ind_orig1, iangle2, ix + iy * nx, nx, ny, nat, box)
                    new_ind2 = get_new_index(coords, ind_orig3, iangle2, ix + iy * nx, nx, ny, nat, box)
                    extended_angles.append([[new_ind1, iangle2 + 1, new_ind2],
                                           [coords[new_ind1 - 1][1].lower(),
                                            coords[iangle2][1].lower(),
                                            coords[iangle3][1].lower()]])
                    continue

                elif bond_length1 > 5.0:
                    cntr2 += 1
                    new_ind1 = get_new_index(coords, ind_orig1, iangle2, ix + iy * nx, nx, ny, nat, box)
                    extended_angles.append([[new_ind1, iangle2 + 1, iangle3 + 1],
                                           [coords[new_ind1 - 1][1].lower(),
                                            coords[iangle2][1].lower(),
                                            coords[iangle3][1].lower()]])
                    continue

                elif bond_length2 > 5.0:
                    cntr3 += 1
                    new_ind2 = get_new_index(coords, ind_orig3, iangle2, ix + iy * nx, nx, ny, nat, box)
                    extended_angles.append([[iangle1 + 1, iangle2 + 1, new_ind2],
                                           [coords[iangle1][1].lower(),
                                            coords[iangle2][1].lower(),
                                            coords[new_ind2 - 1][1].lower()]])
                    continue

    print(cntr1, cntr2, cntr3, cntr1 + cntr2 + cntr3, len(extended_angles), nx*ny*len(angles))

    return extended_angles


def extend_dihedrals(dihedrals, coords, box, reprod, nat, ndihed):
    """
    The function works for dihedral torsion angles.

    at1   at3
      \   / \
       \ /   \
       at2   at4

    The index of at2 will not change, i.e. it is the reference atom for a given dihedral
    torsion angle. at1 and at2 form one plane, at3 and at4 the second plane.

    Args:
        dihedrals (list): contains the dihedral torsion angles as read from the psf file.
        coords (list): Contains element, force field type, and coordinates for the extended
            system.
        box (array.py): Contains the cell parameters.
        reprod (list): How many times the system is extended in each direction.
        nat (int): Number of atoms in the original structure.
        ndihed (int): The number of dihedral torsion angles in the original system.
    """

    extended_dihedrals = []
    nx, ny = reprod[0][:]
    cntr1 = cntr2 = cntr3 = cntr4 = cntr5 = cntr6 = 0

    for iy in range(ny):  # y-direction
        for ix in range(nx):  # x direction
            addind = nat * (ix + iy * nx)
            for nd in range(ndihed):
                ind_orig1 = dihedrals[nd][0]
                ind_orig2 = dihedrals[nd][1]
                ind_orig3 = dihedrals[nd][2]
                ind_orig4 = dihedrals[nd][3]

                idhd1 = ind_orig1 + addind - 1
                idhd2 = ind_orig2 + addind - 1
                idhd3 = ind_orig3 + addind - 1
                idhd4 = ind_orig4 + addind - 1
                # idhdX has to be reduced by 1 in order to get the correct coordinates from coords.

                bond_length1 = calc_bond_length(coords[idhd1][3:6],
                                                coords[idhd2][3:6])
                bond_length2 = calc_bond_length(coords[idhd3][3:6],
                                                coords[idhd2][3:6])
                bond_length3 = calc_bond_length(coords[idhd4][3:6],
                                                coords[idhd3][3:6])

                if bond_length1 < 5.0 and bond_length2 < 5.0 and bond_length3 < 5.0:
                    # b1, b2, b3 ok
                    cntr1 += 1
                    extended_dihedrals.append([[idhd1 + 1, idhd2 + 1, idhd3 + 1, idhd4 + 1],
                                               [coords[idhd1][1].lower(),
                                                coords[idhd2][1].lower(),
                                                coords[idhd3][1].lower(),
                                                coords[idhd4][1].lower()]])
                    continue

                elif bond_length1 > 5.0 and bond_length2 > 5.0 > bond_length3:
                    cntr5 += 1
                    # b3 ok; b1, b2 not. Atom 4 has to be changed as well.
                    new_ind1 = get_new_index(coords, ind_orig1, idhd2, ix + iy * nx, nx, ny, nat, box)
                    new_ind2 = get_new_index(coords, ind_orig3, idhd2, ix + iy * nx, nx, ny, nat, box)
                    new_ind3 = get_new_index(coords, ind_orig4, new_ind2 - 1, ix + iy * nx, nx, ny, nat, box)
                    extended_dihedrals.append([[new_ind1, idhd2 + 1, new_ind2, new_ind3],
                                               [coords[new_ind1 - 1][1].lower(),
                                                coords[idhd2][1].lower(),
                                                coords[new_ind2 - 1][1].lower(),
                                                coords[new_ind3 - 1][1].lower()]])
                    continue

                elif bond_length1 > 5.0 and bond_length3 > 5.0 > bond_length2:
                    cntr4 += 1
                    # b2 ok; b1 and b3 not. indices of atoms 1 and 4 have to be changed.
                    new_ind1 = get_new_index(coords, ind_orig1, idhd2, ix + iy * nx, nx, ny, nat, box)
                    new_ind2 = get_new_index(coords, ind_orig4, idhd3, ix + iy * nx, nx, ny, nat, box)
                    extended_dihedrals.append([[new_ind1, idhd2 + 1, idhd3 + 1, new_ind2],
                                               [coords[new_ind1 - 1][1].lower(),
                                                coords[idhd2][1].lower(),
                                                coords[idhd3][1].lower(),
                                                coords[new_ind2 - 1][1].lower()]])
                    continue

                elif bond_length2 > 5.0 and bond_length3 > 5.0 > bond_length1:
                    cntr6 += 1
                    # b1 ok; b2 and b3 not. indices of atoms 3 and 4 have to be changed.
                    # Care has to be taken, since index of atom 4 will depend on the new
                    # index of atom 3.
                    new_ind1 = get_new_index(coords, ind_orig3, idhd2, ix + iy * nx, nx, ny, nat, box)
                    new_ind2 = get_new_index(coords, ind_orig4, new_ind1 - 1, ix + iy * nx, nx, ny, nat, box)
                    if new_ind2 is None or new_ind1 is None:
                        print(new_ind1, new_ind2, ind_orig1, ind_orig2, ind_orig3, ind_orig4)
                    extended_dihedrals.append([[new_ind1, idhd2 + 1, idhd3 + 1, new_ind2],
                                               [coords[idhd1][1].lower(),
                                                coords[idhd2][1].lower(),
                                                coords[new_ind1 - 1][1].lower(),
                                                coords[new_ind2 - 1][1].lower()]])
                    continue

                elif bond_length1 > 5.0:
                    # b2, b3 ok; b1 not. Atom 1 has to be changed.
                    cntr2 += 1
                    new_ind = get_new_index(coords, ind_orig1, idhd2, ix + iy * nx, nx, ny, nat, box)
                    extended_dihedrals.append([[new_ind, idhd2 + 1, idhd3 + 1, idhd4 + 1],
                                               [coords[new_ind - 1][1].lower(),
                                                coords[idhd2][1].lower(),
                                                coords[idhd3][1].lower(),
                                                coords[idhd4][1].lower()]])
                    continue

                elif bond_length2 > 5.0:
                    cntr3 += 1
                    # b1, b3 ok; b2 not. Indices of atoms 3  and 4 have to be changed since atom 4 index
                    # depends on the new index on atom 3
                    new_ind1 = get_new_index(coords, ind_orig3, idhd2, ix + iy * nx, nx, ny, nat, box)
                    new_ind2 = get_new_index(coords, ind_orig4, new_ind1 - 1, ix + iy * nx, nx, ny, nat, box)
                    extended_dihedrals.append([[idhd1 + 1, idhd2 + 1, new_ind1, new_ind2],
                                               [coords[idhd1][1].lower(),
                                                coords[idhd2][1].lower(),
                                                coords[new_ind1 - 1][1].lower(),
                                                coords[new_ind2 - 1][1].lower()]])
                    continue

                elif bond_length3 > 5.0:
                    cntr3 += 1
                    # b1, b2 ok; b3 not. Atom 4 has to be changed.
                    new_ind = get_new_index(coords, ind_orig4, idhd3, ix + iy * nx, nx, ny, nat, box)
                    extended_dihedrals.append([[idhd1 + 1, idhd2 + 1, idhd3 + 1, new_ind],
                                               [coords[idhd1][1].lower(),
                                                coords[idhd2][1].lower(),
                                                coords[idhd3][1].lower(),
                                                coords[new_ind - 1][1].lower()]])
                    continue

    print(cntr1, cntr2, cntr3, cntr4, cntr1 + cntr2 + cntr3 + cntr4, len(extended_dihedrals), nx*ny*len(dihedrals))

    return extended_dihedrals


def extend_impropers(impropers, coords, box, reprod, nat, nimprop):
    """
    The function works for improper torsion angles.

    at1   at3
      \   /
       \ /
       at2---at4

    The index of at2 will not change, i.e. it is the reference atom for
    a given improper torsion angle.

    Args:
        impropers (list): contains the improper torsion angles as read from the psf file.
        coords (list): Contains element, force field type, and coordinates for the extended
            system.
        box (array.py): Contains the cell parameters.
        reprod (list): How many times the system is extended in each direction.
        nat (int): Number of atoms in the original structure.
        nimprop (int): The number of dihedral/improper angles in the original system.
    """

    extended_improps = []
    nx, ny = reprod[0][:]
    cntr1 = cntr2 = cntr3 = cntr4 = 0

    for iy in range(ny):  # y-direction
        for ix in range(nx):  # x direction
            addind = nat * (ix + iy * nx)
            for nd in range(nimprop):
                ind_orig1 = impropers[nd][0]
                ind_orig2 = impropers[nd][1]
                ind_orig3 = impropers[nd][2]
                ind_orig4 = impropers[nd][3]

                impr1 = ind_orig1 + addind - 1
                impr2 = ind_orig2 + addind - 1
                impr3 = ind_orig3 + addind - 1
                impr4 = ind_orig4 + addind - 1
                # idhedX has to be reduced by 1 in order to get the correct coordinates from coords.

                bond_length1 = calc_bond_length(coords[impr1][3:6],
                                                coords[impr2][3:6])
                bond_length2 = calc_bond_length(coords[impr3][3:6],
                                                coords[impr2][3:6])
                bond_length3 = calc_bond_length(coords[impr4][3:6],
                                                coords[impr2][3:6])
                # print(bond_length1, bond_length2, bond_length3)

                if bond_length1 < 5.0 and bond_length2 < 5.0 and bond_length3 < 5.0:
                    # b1, b2, b3 ok
                    cntr1 += 1
                    extended_improps.append([[impr1 + 1, impr2 + 1, impr3 + 1, impr4 + 1],
                                             [coords[impr1][1].lower(),
                                              coords[impr2][1].lower(),
                                              coords[impr3][1].lower(),
                                              coords[impr4][1].lower()]])
                    continue

                if bond_length1 > 5.0:
                    # b2, b3 ok; b1 not. Atom 1 has to be changed.
                    cntr2 += 1
                    new_ind = get_new_index(coords, ind_orig1, impr2, ix + iy * nx, nx, ny, nat, box)
                    extended_improps.append([[new_ind, impr2 + 1, impr3 + 1, impr4 + 1],
                                             [coords[new_ind - 1][1].lower(),
                                              coords[impr2][1].lower(),
                                              coords[impr3][1].lower(),
                                              coords[impr4][1].lower()]])

                if bond_length2 > 5.0:
                    cntr3 += 1
                    # b1, b3 ok; b2 not. Atom 3 has to be changed.
                    new_ind = get_new_index(coords, ind_orig3, impr2, ix + iy * nx, nx, ny, nat, box)
                    extended_improps.append([[impr1 + 1, impr2 + 1, new_ind, impr4 + 1],
                                             [coords[impr1][1].lower(),
                                              coords[impr2][1].lower(),
                                              coords[new_ind - 1][1].lower(),
                                              coords[impr4][1].lower()]])

                if bond_length3 > 5.0:
                    cntr4 += 1
                    # b1, b2 ok; b3 not. index of atom 4 has to be changed.
                    new_ind = get_new_index(coords, ind_orig4, impr2, ix + iy * nx, nx, ny, nat, box)
                    extended_improps.append([[impr1 + 1, impr2 + 1, impr3 + 1, new_ind],
                                             [coords[impr1][1].lower(),
                                              coords[impr2][1].lower(),
                                              coords[impr3][1].lower(),
                                              coords[new_ind - 1][1].lower()]])

    return extended_improps
