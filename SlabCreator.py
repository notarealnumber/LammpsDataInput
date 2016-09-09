from read_all import *
from check_pbc import *
from extend_system import *


def main():
    print("heyho")
    filebase, reprod, pbc_box = getInfos()

    elements, coords = readpdb(filebase)
    bonds, ff_type, charges, angles, dihedrals, \
    impropers, donors, acceptors, nonbonds = readpsf(filebase)

    ## These variables give the number of atoms, angles, bonds etc of the original structure.
    nat = len(elements)
    nbonds = len(bonds)
    nangles = len(angles)
    ndihed = len(dihedrals)
    nimprp = len(impropers)
    ndon = len(donors)
    nacc = len(acceptors)
    nnonb = len(nonbonds)

    coords_ext = extend_coords(coords, elements, ff_type, reprod, pbc_box)

    # extend_bonds(bonds, coords_ext, reprod, nat, nbonds)
    # angles_ext = extend_thetas(angles, coords_ext, reprod, nat, nangles)
    # dihedrals_ext = extend_dihedrals(dihedrals, coords_ext, reprod, nat, ndihed)
    impropers_ext = extend_dihedrals(impropers, coords_ext, pbc_box, reprod, nat, nimprp)

    # print(angles_ext[4558])


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
    reprod.append([3, 6])
    box.append([40.3148, 41.4320, 53.6079])

    # print(" ")
    # print("####################################")
    # basename = input("Enter base name (without extension): ")
    basename = "silica_Q3_amorph_4_7OH_0pct_ion"
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

    return basename, reprod, box


main()