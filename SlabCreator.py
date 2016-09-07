from read_all import *
from check_pbc import *
from extend_system import *

def main():
    print("heyho")
    filebase, reprod, pbc_box = getInfos()

    elements, coords = readpdb(filebase)
    bonds, ff_type, charges, angles, dihedrals, \
    impropers, donors, acceptors, nonbonds = readpsf(filebase)

    pbc_check = check_pbc_applied(coords, bonds, pbc_box)
    coords_ext = multiply_coords(coords, reprod, pbc_box)

    elements_ext = []

    for frame in range(reprod[0][0]*reprod[0][1]):
        for elmt in range(len(elements)):
            elements_ext.append(elements[elmt])

    # print("ccc", len(elements_ext), len(coords_ext))
    # thefile = "extended.xyz"
    # xyzfile = open(thefile, "w")
    # print(len(coords_ext), file=xyzfile)
    # print("Extended coordinates", file=xyzfile)
    # for i in range(len(coords_ext)):
    #     print(elements_ext[i], coords_ext[i, 0], coords_ext[i, 1], coords_ext[i, 2], file=xyzfile)

    bonds_extended = multiply_bonds(bonds, pbc_check, reprod, len(coords))
    print(angles[354])
    print(bonds_extended.index([angles[354][1], angles[354][0]]))
    print(bonds.index([angles[354][1], angles[354][0]]))
    print(bonds_extended.index([angles[354][1], angles[354][2]]))
    print(bonds.index([angles[354][1], angles[354][2]]))
    print(reprod[0][0] * reprod[0][1] * len(bonds), len(bonds_extended))


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