import numpy as np
from multiprocessing.dummy import Pool as ThreadPool


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
        coord1 (float): The atom to be subtracted.
        coord2 (float): The other atom that forms the bond with coord1-atom.
    """

    bond_diff = [coord2[0] - coord1[0],
                 coord2[1] - coord1[1],
                 coord2[2] - coord1[2]]

    bond_length = np.sqrt(bond_diff[0]**2 + bond_diff[1]**2 + bond_diff[2]**2)

    return bond_length


def extend_coords(coords, elements, ff_type, reprod, box):
    nat = len(coords)
    extended_coords = []
    # extended_coords = np.full((reprod[0][0] * reprod[0][1] * nat, 3), 0.0)

    print(3*box[0, 0], 6*box[0, 1], box[0, 2])
    print(box[0, 0], box[0, 1], box[0, 2])

    for iy in range(reprod[0][1]):
        for ix in range(reprod[0][0]):
            for k in range(nat):
                extended_coords.append([elements[k], ff_type[k],
                                        coords[k, 0] + box[0, 0] * ix,
                                        coords[k, 1] + box[0, 1] * iy,
                                        coords[k, 2]])

    return extended_coords


def extend_bonds(bonds, coords, reprod, nat, nbonds):
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

    imaged = 0
    imgd = 0
    celled = 0

    for iy in range(ny):  # y-direction
        for ix in range(nx):  # x direction
            for nb in range(nbonds):
                # -1 in order to get the correct element in coords.
                bond1 = bonds[nb][0] + nat * (ix + iy * nx) - 1
                bond2 = bonds[nb][1] + nat * (ix + iy * nx) - 1
                bond_length = calc_bond_length(coords[bond1][2:5], coords[bond2][2:5])
                if bond_length < 5.0:
                    celled += 1
                    extended_bonds.append([[bond1 + 1, bond2 + 1],
                                           [coords[bond1][1].lower(), coords[bond2][1].lower()]])
                    continue
                else:
                    imgd += 1
                    for iat in range(len(coords)):
                        check_length = calc_bond_length(coords[bond1][2:5], coords[iat][2:5])
                        if check_length < 5.0 and iat != bond1:
                            imaged += 1
                            extended_bonds.append([[bond1 + 1, iat + 1],
                                                   [coords[bond1][1].lower(), coords[iat][1].lower()]])
                            break

    print(len(extended_bonds))
    print(imaged, imgd, celled, imaged+celled, nx*ny*len(bonds))

    return extended_bonds


def extend_thetas(angles, coords, reprod, nat, nangles):
    """
    iangle1       iangle3
         \        /
          \      /
           \    /
          iangle2
    """

    extended_angles = []
    nx, ny = reprod[0][:]

    ok12 = 0
    ok1no2 = 0
    no1ok2 = 0
    no12 = 0

    for iy in range(ny):  # y-direction
        for ix in range(nx):  # x direction
            for nang in range(nangles):
                # -1 in order to get the correct element in coords.
                iangle1 = angles[nang][0] + nat * (ix + iy * nx) - 1
                iangle2 = angles[nang][1] + nat * (ix + iy * nx) - 1 # The middle atom
                iangle3 = angles[nang][2] + nat * (ix + iy * nx) - 1

                bond_length1 = calc_bond_length(coords[iangle2][2:5], coords[iangle1][2:5])
                bond_length2 = calc_bond_length(coords[iangle2][2:5], coords[iangle3][2:5])

                if bond_length1 < 5.0 and bond_length2 < 5.0:
                    ok12 += 1
                    extended_angles.append([[iangle1 + 1, iangle2 + 1, iangle3 + 1],
                                           [coords[iangle1][1].lower(),
                                            coords[iangle2][1].lower(),
                                            coords[iangle3][1].lower()]])
                    continue
                elif bond_length1 > 5.0 and bond_length2 > 5.0:
                    chk = 0
                    bond = []
                    for iat in range(len(coords)):
                        check_length = calc_bond_length(coords[iangle2][2:5], coords[iat][2:5])
                        if check_length < 5.0 and iat != iangle2:
                            chk += 100
                            bond.append(iat)
                            no12 += 1
                        elif chk == 200:
                            extended_angles.append([[bond[0] + 1, iangle1 + 1, bond[1] + 1],
                                                   [coords[bond[0]][1].lower(),
                                                    coords[iangle1][1].lower(),
                                                    coords[bond[1]][1].lower()]])
                            break
                elif bond_length1 > 5.0 > bond_length2:
                    for iat in range(len(coords)):
                        check_length = calc_bond_length(coords[iangle2][2:5], coords[iat][2:5])
                        if check_length < 5.0 and iat != iangle2 and iat != iangle3:
                            no1ok2 += 1
                            extended_angles.append([[iat + 1, iangle2 + 1, iangle3 + 1],
                                                   [coords[iat][1].lower(),
                                                    coords[iangle2][1].lower(),
                                                    coords[iangle3][1].lower()]])
                            break
                elif bond_length2 > 5.0 > bond_length1:
                    for iat in range(len(coords)):
                        check_length = calc_bond_length(coords[iangle2][2:5], coords[iat][2:5])
                        if check_length < 5.0 and iat != iangle2 and iat != iangle1:
                            ok1no2 += 1
                            extended_angles.append([[iangle1 + 1, iangle2 + 1, iat + 1],
                                                   [coords[iangle1][1].lower(),
                                                    coords[iangle2][1].lower(),
                                                    coords[iat][1].lower()]])
                            break

    print(len(extended_angles))
    print(ok12, ok1no2, no1ok2, no12, ok12 + ok1no2 + no1ok2 + no12, 3*6*7220)

    return extended_angles


def extend_dihedrals(dihedral, coords, box, reprod, nat, ndihed):
    """
    The function works for dihedral and improper angles.

    at1   at3
      \   / \
       \ /   \
       at2   at4

    The index of at2 will not change, i.e. it is the reference atom for
    a given dihedral, improper angle.

    Args:
        addind (int): addind will be added to the original numbers that make the
            dihedral/improper angle as read from the psf file. It has to be reduced by
            one in order to get the correct coordinates from coords.
        coords (list): Contains element, force field type, and coordinates for the extended
            system.
        reprod (int): How many times the system is extended in each direction.
        ndihed (int): The number of dihedral/improper angles in the original system.
    """

    extended_dihed = []
    nx, ny = reprod[0][:]
    cntr1 = cntr2 = cntr3 = cntr4 = cntr5 = 0

    for iy in range(ny):  # y-direction
        for ix in range(nx):  # x direction
            addind = nat * (ix + iy * nx)
            for nd in range(ndihed):
                ind_orig1 = dihedral[nd][0]
                ind_orig2 = dihedral[nd][1]
                ind_orig3 = dihedral[nd][2]
                ind_orig4 = dihedral[nd][3]
                # Create these variables as they are of use later.
                idhed1 = ind_orig1 + addind - 1
                idhed2 = ind_orig2 + addind - 1
                idhed3 = ind_orig3 + addind - 1
                idhed4 = ind_orig4 + addind - 1
                # idhedX has to be reduced by 1 in order to get the correct coordinates from coords.

                bond_length1 = calc_bond_length(coords[idhed2][2:5], coords[idhed1][2:5])
                bond_length2 = calc_bond_length(coords[idhed2][2:5], coords[idhed3][2:5])
                bond_length3 = calc_bond_length(coords[idhed3][2:5], coords[idhed4][2:5])

                if bond_length1 < 5.0 and bond_length2 < 5.0 and bond_length3 < 5.0:
                    # b1, b2, b3 ok
                    cntr1 += 1
                    extended_dihed.append([[idhed1 + 1, idhed2 + 1, idhed3 + 1, idhed4 + 1],
                                           [coords[idhed1][1].lower(),
                                            coords[idhed2][1].lower(),
                                            coords[idhed3][1].lower(),
                                            coords[idhed4][1].lower()]])
                    continue

                if bond_length1 > 5.0 > bond_length2 and bond_length3 < 5.0:
                    # b2, b3 ok; b1 not. Atom 1 has to be changed.
                    for iadd in range(nx * ny):
                        if iadd == ix + iy * nx:
                            continue
                        success = False
                        for box_y in range(-1, 2):
                            for box_x in range(-1, 2):
                                tempcoord = [
                                    coords[ind_orig1 + iadd * nat - 1][2] + box[0, 0] * nx * box_x,
                                    coords[ind_orig1 + iadd * nat - 1][3] + box[0, 1] * ny * box_y,
                                    coords[ind_orig1 + iadd * nat - 1][4]
                                ]
                                check_length = calc_bond_length(coords[idhed2][2:5], tempcoord)
                                if check_length < 5.0:
                                    cntr2 += 1
                                    extended_dihed.append([
                                        [ind_orig1 + iadd * nat, idhed2 + 1, idhed3 + 1, idhed4 + 1],
                                        [coords[ind_orig1 + iadd * nat - 1][1].lower(),
                                         coords[idhed2][1].lower(),
                                         coords[idhed3][1].lower(),
                                         coords[idhed4][1].lower()]])
                                    success = True
                                if success:
                                    break
                            if success:
                                break
                        if success:
                            break

                    continue

                if bond_length1 < 5.0 and bond_length2 < 5.0 < bond_length3:
                    # b1, b2 ok; b3 not. index of atom 4 has to be changed.
                    for iat in range(len(coords)):
                        check_length = calc_bond_length(coords[idhed3][2:5], coords[iat][2:5])
                        if check_length < 5.0 and iat != idhed3 and iat != idhed4 and iat != idhed2:
                            cntr3 += 1
                            extended_dihed.append([[idhed1 + 1, idhed2 + 1, idhed3 + 1, iat + 1],
                                                   [coords[idhed1][1].lower(),
                                                    coords[idhed2][1].lower(),
                                                    coords[idhed3][1].lower(),
                                                    coords[iat][1].lower()]])
                            break
                    continue

                if bond_length1 < 5.0 < bond_length2:
                    # if b1 ok but b2 not then b3 has to be checked automatically. Search the correct
                    # index of atom 3 first, then atom 4's index, which depends on atom 3's index.
                    for iat in range(len(coords)):
                        check_length = calc_bond_length(coords[idhed2][2:5], coords[iat][2:5])
                        if check_length < 5.0 and iat != idhed2 and iat != idhed1:
                            cntr4 += 1
                            ntemp = iat
                            break
                    for iat in range(len(coords)):
                        check_length = calc_bond_length(coords[ntemp][2:5], coords[iat][2:5])
                        if check_length < 5.0 and iat != idhed2 and iat != idhed1:
                            cntr4 += 1
                            extended_dihed.append([[idhed1 + 1, idhed2 + 1, ntemp + 1, iat + 1],
                                                   [coords[idhed1][1].lower(),
                                                    coords[idhed2][1].lower(),
                                                    coords[ntemp][1].lower(),
                                                    coords[iat][1].lower()]])
                            break
                    continue

                if bond_length1 > 5.0 > bond_length2 and bond_length3 > 5.0:
                    # b2 ok; b1, b3 not. The indices of atoms 1 and 4 have to be changed.
                    # Fix b1 first.
                    for iat in range(len(coords)):
                        check_length = calc_bond_length(coords[idhed2][2:5], coords[iat][2:5])
                        if check_length < 5.0 and iat != idhed2 and iat != idhed3:
                            cntr5 += 1
                            ntemp = iat
                            break
                    # Then fix b3.
                    for iat in range(len(coords)):
                        check_length = calc_bond_length(coords[idhed3][2:5], coords[iat][2:5])
                        if check_length < 5.0 and iat != idhed3 and iat != idhed2:
                            cntr5 += 1
                            extended_dihed.append([[idhed1 + 1, idhed2 + 1, ntemp + 1, iat + 1],
                                                   [coords[ntemp][1].lower(),
                                                    coords[idhed2][1].lower(),
                                                    coords[idhed3][1].lower(),
                                                    coords[iat][1].lower()]])
                            break
                    continue

            # if len(extended_dihed)%18 != 0:
            #     print(len(extended_dihed)%18, len(extended_dihed), ix, iy)
            # print(len(extended_dihed), ix, iy)

    print(cntr1, cntr2, cntr3, cntr4, cntr5,
          cntr1 + cntr2 + cntr3 + cntr4 + cntr5, len(extended_dihed), nx*ny*len(dihedral))

    return extended_dihed
