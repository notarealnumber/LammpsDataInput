from read_all import *
from extend_system import *
from write_lammps_data_file import *


def main():

    filebase, forcefield, reprod, pbc_box = getInfos()

    elements, coords = readpdb(filebase)
    bonds, ff_type, charges, masses, angles, dihedrals, \
        impropers, donors, acceptors, nonbonds = readpsf(filebase)

    types = sorted(list(set(ff_type)))
    for i in range(len(types)):
        types[i] = types[i].lower()
    print(types)

    # get_parms(forcefield, types)

    # These variables give the number of atoms, angles, bonds etc of the original structure.
    nat = len(elements)
    nbonds = len(bonds)
    nangles = len(angles)
    ndihed = len(dihedrals)
    nimprp = len(impropers)
    ndon = len(donors)
    nacc = len(acceptors)
    nnonb = len(nonbonds)

    coords_ext = extend_coords(coords, elements, ff_type, charges, reprod, pbc_box)

    bonds_ext = extend_bonds(bonds, coords_ext, pbc_box, reprod, nat, nbonds)
    angles_ext = extend_thetas(angles, coords_ext, pbc_box, reprod, nat, nangles)
    dihedrals_ext = extend_dihedrals(dihedrals, coords_ext, pbc_box, reprod, nat, ndihed)
    # impropers_ext = extend_impropers(impropers, coords_ext, pbc_box, reprod, nat, nimprp)

    write_lammps_file(coords_ext, bonds_ext, angles_ext, dihedrals_ext, ff_type,
                      charges, masses, filebase, reprod, pbc_box, types)


def getInfos():
    """ Small function that asks the user for filename, how many frames
    the trajectory contains, and the number of atoms a molecule has.

    :return:
    filename
    nsteps
    at_per_molec
    """
    reprod = []
    box = []
    reprod.append([15, 15])
    box.append([40.3148, 41.4320, 53.6079])

    # print(" ")
    # print("####################################")
    # basename = input("Enter base name (without extension): ")
    basename = "silica_Q3_amorph_4_7OH_0pct_ion"
    # print("####################################")
    # ffname = input("Enter base name (without extension): ")
    ffname = "cvff_interface_v1_5.frc"
    #
    # print("####################################")
    # print(" ")
    # reprod.append(int(input("How many times should the slab be replicated in x: ")))
    # print("####################################")
    # reprod.append(int(input("How many times should the slab be replicated in y: ")))
    # print(" ")
    # print("####################################")
    # print(" ")

    box = np.array(box)

    return basename, ffname, reprod, box


main()
