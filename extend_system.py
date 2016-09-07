import numpy as np


def multiply_coords(coords, reprod, box):
    nat = len(coords)
    extended_coords = np.full((reprod[0][0] * reprod[0][1] * nat, 3), 0.0)

    print(3*box[0, 0], 6*box[0, 1], box[0, 2])

    mlt = 0
    for iy in range(reprod[0][1]):
        for ix in range(reprod[0][0]):
            for k in range(nat):
                extended_coords[k + mlt * nat, 0] = coords[k, 0] + box[0, 0] * ix
                extended_coords[k + mlt * nat, 1] = coords[k, 1] + box[0, 1] * iy
                extended_coords[k + mlt * nat, 2] = coords[k, 2]
            mlt += 1

    return extended_coords


def multiply_bonds(bonds, pbc_check, reprod, nat):
    """
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

    nbonds = len(bonds)
    extended_bonds = []
    nx, ny = reprod[0][:]

    for iy in range(ny):  # y-direction
        for ix in range(nx):  # x direction
            for nb in range(nbonds):
                if pbc_check[nb][0]:

                    if pbc_check[nb][1] == "x":
                        if ix == nx-1:
                            extended_bonds.append([[ix + iy * nx],
                                                   [bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * nx * iy],
                                                   [bonds[nb][0], bonds[nb][1]]])
                        else:
                            extended_bonds.append([[bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * ((ix + iy * nx) + 1)],
                                                   [bonds[nb][0], bonds[nb][1]]])

                    if pbc_check[nb][1] == "-x":
                        if ix == 0:
                            extended_bonds.append([[bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * nx * (iy + 1)],
                                                   [bonds[nb][0], bonds[nb][1]]])
                        else:
                            extended_bonds.append([[bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * ((ix + iy * nx) - 1)],
                                                   [bonds[nb][0], bonds[nb][1]]])

                    elif pbc_check[nb][1] == "y":
                        if iy == ny-1:
                            extended_bonds.append([[bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * ix],
                                                   [bonds[nb][0], bonds[nb][1]]])
                        else:
                            extended_bonds.append([[bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * nx * (iy + 1)],
                                                   [bonds[nb][0], bonds[nb][1]]])

                    elif pbc_check[nb][1] == "-y":
                        if iy == 0:
                            extended_bonds.append([[bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * (ix + (ny - 1) * nx)],
                                                   [bonds[nb][0], bonds[nb][1]]])
                        else:
                            extended_bonds.append([[bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * (ix + (iy - 1) * nx)],
                                                   [bonds[nb][0], bonds[nb][1]]])

                    elif pbc_check[nb][1] == "xy":
                        if ix == nx-1 and iy != ny-1:
                            extended_bonds.append([[bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * ((ix + iy * nx) + 1)],
                                                   [bonds[nb][0], bonds[nb][1]]])
                        elif ix != nx-1 and iy == ny-1:
                            extended_bonds.append([[bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * (ix + 1)],
                                                   [bonds[nb][0], bonds[nb][1]]])
                        elif ix == nx-1 and iy == ny-1:
                            extended_bonds.append([[bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1]],
                                                   [bonds[nb][0], bonds[nb][1]]])
                        else:
                            extended_bonds.append([[bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * (ix - iy + (iy + 1) * (nx + 1))],
                                                   [bonds[nb][0], bonds[nb][1]]])
                        print(extended_bonds[nb])

                    elif pbc_check[nb][1] == "-x-y":
                        if ix == 0 and iy == 0:
                            extended_bonds.append([[bonds[nb][0],
                                                   bonds[nb][1] + nat * (nx * ny - 1)],
                                                   [bonds[nb][0], bonds[nb][1]]])
                        elif ix != 0 and iy == 0:
                            extended_bonds.append([[bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * ((ny - 1) * nx + (ix - 1))],
                                                   [bonds[nb][0], bonds[nb][1]]])
                        elif ix == 0 and iy != 0:
                            extended_bonds.append([[bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * ((ix + iy * nx) - 1)],
                                                   [bonds[nb][0], bonds[nb][1]]])
                        else:
                            extended_bonds.append([[bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * (ix - 1 + (iy - 1) * nx)],
                                                   [bonds[nb][0], bonds[nb][1]]])
                        print(extended_bonds[nb])

                    elif pbc_check[nb][1] == "x-y":
                        if ix < nx-1 and iy == 0:
                            extended_bonds.append([[bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * (nx * (ny - 1) + 1 + ix)],
                                                   [bonds[nb][0], bonds[nb][1]]])
                        if ix == nx-1 and iy == 0:
                            extended_bonds.append([[bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * nx * (ny - 1)],
                                                   [bonds[nb][0], bonds[nb][1]]])
                        if ix == nx-1 and iy > 0:
                            extended_bonds.append([[bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * nx * (iy - 1)],
                                                   [bonds[nb][0], bonds[nb][1]]])
                        else:
                            extended_bonds.append([[bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * (ix + nx * (iy - 1) + 1)],
                                                   [bonds[nb][0], bonds[nb][1]]])
                        print(extended_bonds[nb])

                    elif pbc_check[nb][1] == "-xy":
                        print("say hehooo")

                else:
                    extended_bonds.append([[bonds[nb][0] + nat * (ix + iy * nx),
                                           bonds[nb][1] + nat * (ix + iy * nx)],
                                           [bonds[nb][0], bonds[nb][1]]])

    return extended_bonds
