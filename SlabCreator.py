import sys
from read_all import *
from extend_system import *
from write_lammps_data_file import *


def main():

    single_mol_info = []

    periodic_or_not, \
        single_mol_file, \
        complete_sys_file, \
        reprod, \
        box, \
        ff_file, \
        ff_name = getInfos()

    if periodic_or_not == "mol":
        print("Molecuuule")

        bonds, ff_type, charges, masses, angles, dihedrals, \
            impropers, donors, acceptors, nonbonds = readpsf(single_mol_file)
        elements, coords = readpdb(complete_sys_file)

        # These variables give the number of atoms, angles, bonds etc of the original structure.
        # nat = len(elements)
        # nbonds = len(bonds)
        # nangles = len(angles)
        # ndihed = len(dihedrals)
        # nimprp = len(impropers)
        # ndon = len(donors)
        # nacc = len(acceptors)
        # nnonb = len(nonbonds)

        for i in range(len(charges)):
            single_mol_info.append([elements[i], ff_type[i], charges[i], masses[i]])
            print(single_mol_info[-1])

        temp_type = sorted(list(set(ff_type)))
        for i in range(len(temp_type)):
            temp_type[i] = temp_type[i].lower()

        types = []
        for i in range(len(temp_type)):
            check = 0
            for j in range(len(single_mol_info)):
                if temp_type[i] in single_mol_info[j]:
                    types.append(single_mol_info[j])
                    check += 1
                    break
            if check == 1:
                continue

        print(types)

        write_lammps_file_nonperiodic(single_mol_info, types,
                                      coords, bonds, angles, dihedrals, impropers,
                                      donors, acceptors, nonbonds,
                                      complete_sys_file, box)

    else:
        elements, coords = readpdb(complete_sys_file)
        bonds, ff_type, charges, masses, angles, dihedrals, \
            impropers, donors, acceptors, nonbonds = readpsf(complete_sys_file)

        # These variables give the number of atoms, angles, bonds etc of the original structure.
        nat = len(elements)
        nbonds = len(bonds)
        nangles = len(angles)
        ndihed = len(dihedrals)
        nimprp = len(impropers)
        ndon = len(donors)
        nacc = len(acceptors)
        nnonb = len(nonbonds)

        temp_type = sorted(list(set(ff_type)))
        for i in range(len(temp_type)):
            temp_type[i] = temp_type[i].lower()
        print(temp_type)

        coords_ext = extend_coords(coords, elements, ff_type, charges, reprod, box)

        bonds_ext = extend_bonds(bonds, coords_ext, box, reprod, nat, nbonds)
        angles_ext = extend_thetas(angles, coords_ext, box, reprod, nat, nangles)
        dihedrals_ext = extend_dihedrals(dihedrals, coords_ext, box, reprod, nat, ndihed)
        # impropers_ext = extend_impropers(impropers, coords_ext, box, reprod, nat, nimprp)

        write_lammps_file_periodic(coords_ext, bonds_ext, angles_ext, dihedrals_ext, ff_type,
                                   charges, masses, complete_sys_file, reprod, box, temp_type)


def getInfos():
    """ Small function that gets information about the system to be replicated (for periodic cases)
    or for molecular systems (non-periodic cases). Read the comments in SlabCreator.input for details.
    """

    reprod = []
    pbc_box = []

    for i in range(9):
        sys.stdin.readline()
    slab_or_mol = sys.stdin.readline().strip()
    # print(slab_or_mol)
    for i in range(4):
        sys.stdin.readline()
    filename_single_mol = sys.stdin.readline().strip()
    for i in range(3):
        sys.stdin.readline()
    filename_total_sys = sys.stdin.readline().strip()
    sys.stdin.readline()
    sys.stdin.readline()
    templist = sys.stdin.readline().strip(" ").split()
    reprod.append([int(templist[0]), int(templist[1])])
    # print(reprod)
    sys.stdin.readline()
    sys.stdin.readline()
    sys.stdin.readline()
    tempbox = sys.stdin.readline().strip(" ").split()
    pbc_box.append([float(tempbox[0]), float(tempbox[1]), float(tempbox[2])])
    # print(box)
    sys.stdin.readline()
    ff_file = sys.stdin.readline().strip()
    # print(forcefield)
    sys.stdin.readline()
    ff_name = sys.stdin.readline().strip()

    box = np.array(pbc_box)

    return slab_or_mol, filename_single_mol, filename_total_sys, reprod, box, ff_file, ff_name


main()
