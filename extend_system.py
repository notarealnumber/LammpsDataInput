import numpy as np


def multiply_coords(coords, reprod, box):
    nat = len(coords)
    extended_coords = np.full((reprod[0][0] * reprod[0][1] * nat, 3), 0.0)

    mlt = 0
    for i in range(reprod[0][0]):
        for j in range(reprod[0][1]):
            for k in range(nat):
                extended_coords[k*mlt, 0] = coords[k, 0] + box[0][0] * i
                extended_coords[k*mlt, 1] = coords[k, 1] + box[0][1] * j
                extended_coords[k*mlt, 2] = coords[k, 2]
            mlt += 1

    return extended_coords


def multiply_bonds(bonds, pbc_check, reprod, nat):
    """
    The first atom index in bonds will be set into the current ix, iy cell.
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
                            extended_bonds.append([bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * nx * iy])
                        else:
                            extended_bonds.append([bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * ((ix + iy * nx) + 1)])

                    if pbc_check[nb][1] == "-x":
                        if ix == 0:
                            extended_bonds.append([bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * nx * (iy + 1)])
                        else:
                            extended_bonds.append([bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * ((ix + iy * nx) - 1)])

                    elif pbc_check[nb][1] == "y":
                        if iy == ny-1:
                            extended_bonds.append([bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * ix])
                        else:
                            extended_bonds.append([bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * nx * (iy + 1)])

                    elif pbc_check[nb][1] == "xy":
                        if ix == nx-1 and iy != ny-1:
                            extended_bonds.append([bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * ((ix + iy * nx) + 1)])
                        elif ix != nx-1 and iy == ny-1:
                            extended_bonds.append([bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * (ix + 1)])
                        elif ix == nx-1 and iy == ny-1:
                            extended_bonds.append([bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1]])
                        else:
                            extended_bonds.append([bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * (ix - iy + (iy + 1) * (nx + 1))])

                    elif pbc_check[nb][1] == "-x-y":
                        if ix == 0 and iy == 0:
                            extended_bonds.append([bonds[nb][0],
                                                   bonds[nb][1] + (nx * ny - 1) * nat])
                        elif ix != 0 and iy == 0:
                            extended_bonds.append([bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * ((ny - 1) * nx + (ix - 1))])
                        elif ix == 0 and iy != 0:
                            extended_bonds.append([bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * ((ix + iy * nx) - 1)])
                        else:
                            extended_bonds.append([bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + nat * (ix - 1 + (iy - 1) * nx)])

                    elif pbc_check[nb][1] == "x-y":
                        if ix != nx-1 and iy == 0:
                            extended_bonds.append([bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + (nx * (ny - 1) + 1 + ix) * nat])
                        if ix == nx-1 and iy == 0:
                            extended_bonds.append([bonds[nb][0] + nat * (ix + iy * nx),
                                                   bonds[nb][1] + (nx * (ny - 1) + 1 + ix) * nat])




                else:
                    extended_bonds.append([bonds[nb][0] + nat * (ix + iy * nx),
                                           bonds[nb][1] + nat * (ix + iy * nx)])
