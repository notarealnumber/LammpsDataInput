import numpy as np


def get_parms(ff_file, ff_types):
    """
    Reads the file containing the CVFF parameters.
    """

    with open(ff_file, "r") as parms:
        for line in parms:
            theline = line.strip(" ").split()
            if "#quadratic_bond" in theline:
                for i in range(5):
                    parms.readline()
                finish = True
                while finish:
                    thebond = parms.readline().strip(" ").split()
                    if thebond == " ":
                        print(thebond)
                        finish = False
                    elif thebond[2] in ff_types and thebond[3] in ff_types:
                        print(thebond)

                print(theline)
                exit()


def readpdb(basefile):
    """
    Reads a pdb file generated by Fateme S. Emami as included in the Interface FF 1.5 database.
    """
    pdbfile = str(basefile) + ".pdb"
    print(pdbfile)
    elements = []
    coords = []

    with open(pdbfile, "r") as ins:
        for line in ins:
            templist = line.strip(" ").split()
            if templist[0] == "ATOM":
                temp_coords = []
                for i in range(len(templist)):
                    try:
                        temp_coords.append(float(templist[i]))
                    except:
                        pass
                elements.append(templist[2])
                x = temp_coords[2]
                y = temp_coords[3]
                z = temp_coords[4]
                coords.append([x, y, z])

    coords = np.array(coords)
    return elements, coords


def readxyz(basefile):
    """
    Reads a generic xyz file.
    """
    xyzfile = str(basefile) + ".xyz"
    print(xyzfile)
    elements = []
    coords = []

    with open(xyzfile, "r") as ins:
        ins.readline()
        ins.readline()
        for line in ins:
            templist = line.strip(" ").split()
            elements.append(templist[0])
            x = float(templist[1])
            y = float(templist[2])
            z = float(templist[3])
            coords.append([x, y, z])

    coords = np.array(coords)
    print(coords[0])
    return elements, coords


def readpsf(basefile):
    """
    Reads a psf file generated by Fateme S. Emami as included in the Interface FF 1.5 database.
    """
    psffile = str(basefile) + ".psf"
    print(psffile)

    bonds = []
    ff_type = []
    charges = []
    masses = []
    angles = []
    dihedrals = []
    impropers = []
    donors = []
    acceptors = []
    nonbonds = []

    with open(psffile, "r") as psf:
        for line in psf:
            templist = line.strip(" ").split()
            # print(templist)
            if '!NATOM' in templist:
                nat = int(templist[0])
                print("Number of atoms", nat)
                for i in range(nat):
                    atomlist = psf.readline().strip(" ").split()
                    ff_type.append(atomlist[5])
                    charges.append(float(atomlist[6]))
                    masses.append(float(atomlist[7]))

            elif '!NBOND:' in templist or '!NBOND' in templist:
                nbonds = int(templist[0])
                if nbonds == 0:
                    print("The structure does not contain any bonds.")
                else:
                    print("Number of bonds", nbonds)
                    bondlist = psf.readline().strip(" ").split()
                    bonds_per_line = int(len(bondlist)/2)
                    lines2read = int(np.ceil(nbonds/bonds_per_line))
                    for nb in range(bonds_per_line):
                        bond1 = int(bondlist[nb*2])
                        bond2 = int(bondlist[(nb*2)+1])
                        bonds.append([bond1, bond2])

                    for i in range(lines2read-2):
                        bondlist = psf.readline().strip(" ").split()
                        for nb in range(bonds_per_line):
                            bond1 = int(bondlist[nb*2])
                            bond2 = int(bondlist[(nb*2)+1])
                            bonds.append([bond1, bond2])

                    bondlist = psf.readline().strip(" ").split()
                    bonds_per_line = int(len(bondlist)/2)
                    print("The last line contains", bonds_per_line, "bonds.")
                    for nb in range(bonds_per_line):
                        bond1 = int(bondlist[nb*2])
                        bond2 = int(bondlist[(nb*2)+1])
                        bonds.append([bond1, bond2])

            elif '!NTHETA:' in templist or '!NTHETA' in templist:
                nangle = int(templist[0])
                if nangle == 0:
                    print("The structure does not contain any bonding angles.")
                else:
                    print("Number of angles", nangle)
                    thetalist = psf.readline().strip(" ").split()
                    thetas_per_line = int(len(thetalist)/3)
                    lines2read = int(np.ceil(nangle/thetas_per_line))
                    for nth in range(thetas_per_line):
                        angle1 = int(thetalist[nth*3])
                        angle2 = int(thetalist[(nth*3)+1])
                        angle3 = int(thetalist[(nth*3)+2])
                        angles.append([angle1, angle2, angle3])

                    for i in range(lines2read-2):
                        # thetalist.pop()
                        thetalist = psf.readline().strip(" ").split()
                        # thetalist.append(psf.readline().strip(" ").split())
                        for nth in range(thetas_per_line):
                            angle1 = int(thetalist[nth*3])
                            angle2 = int(thetalist[(nth*3)+1])
                            angle3 = int(thetalist[(nth*3)+2])
                            angles.append([angle1, angle2, angle3])
                            # angles.append(thetalist[0][nth*3:(nth*3)+3])

                    # thetalist.pop()
                    thetalist = psf.readline().strip(" ").split()
                    # thetalist.insert(0, psf.readline().strip(" ").split())
                    thetas_per_line = int(len(thetalist)/3)
                    print("The last line contains", thetas_per_line, "angles.")
                    for nth in range(thetas_per_line):
                        angle1 = int(thetalist[nth*3])
                        angle2 = int(thetalist[(nth*3)+1])
                        angle3 = int(thetalist[(nth*3)+2])
                        angles.append([angle1, angle2, angle3])
                        # angles.append(thetalist[0][nth*3:(nth*3)+3])

            elif '!NPHI:' in templist or '!NPHI' in templist:
                nphi = int(templist[0])
                if nphi == 0:
                    print("The structure does not contain any dihedrals.")
                else:
                    print("Number of dihedrals", nphi)
                    philist = psf.readline().strip(" ").split()
                    phis_per_line = int(len(philist)/4)
                    lines2read = int(np.ceil(nphi / phis_per_line))
                    for nph in range(phis_per_line):
                        phi1 = int(philist[nph*4])
                        phi2 = int(philist[nph*4+1])
                        phi3 = int(philist[nph*4+2])
                        phi4 = int(philist[nph*4+3])
                        dihedrals.append([phi1, phi2, phi3, phi4])

                    for i in range(lines2read-2):
                        philist = psf.readline().strip(" ").split()
                        for nph in range(phis_per_line):
                            phi1 = int(philist[nph*4])
                            phi2 = int(philist[nph*4+1])
                            phi3 = int(philist[nph*4+2])
                            phi4 = int(philist[nph*4+3])
                            dihedrals.append([phi1, phi2, phi3, phi4])

                    philist = psf.readline().strip(" ").split()
                    phis_per_line = int(len(philist)/4)
                    print("The last line contains", phis_per_line, "dihedrals.")
                    for nph in range(phis_per_line):
                        phi1 = int(philist[nph*4])
                        phi2 = int(philist[nph*4+1])
                        phi3 = int(philist[nph*4+2])
                        phi4 = int(philist[nph*4+3])
                        dihedrals.append([phi1, phi2, phi3, phi4])

            elif '!NIMPHI:' in templist or '!NIMPHI' in templist:
                nimpr = int(templist[0])
                if nimpr == 0:
                    print("The structure does not contain any impropers.")
                else:
                    print("Number of impropers:", nimpr)
                    imprlist = psf.readline().strip(" ").split()
                    imprs_per_line = int(len(imprlist)/4)
                    lines2read = int(np.ceil(nimpr / imprs_per_line))
                    for nimp in range(imprs_per_line):
                        impr1 = int(imprlist[nimp*4])
                        impr2 = int(imprlist[nimp*4+1])
                        impr3 = int(imprlist[nimp*4+2])
                        impr4 = int(imprlist[nimp*4+3])
                        impropers.append([impr1, impr2, impr3, impr4])

                    for i in range(lines2read-2):
                        imprlist = psf.readline().strip(" ").split()
                        for nimp in range(imprs_per_line):
                            impr1 = int(imprlist[nimp*4])
                            impr2 = int(imprlist[nimp*4+1])
                            impr3 = int(imprlist[nimp*4+2])
                            impr4 = int(imprlist[nimp*4+3])
                            impropers.append([impr1, impr2, impr3, impr4])

                    imprlist = psf.readline().strip(" ").split()
                    imprs_per_line = int(len(imprlist[0][:])/4)
                    print("The last line contains", imprs_per_line, "impropers.")
                    for nimp in range(imprs_per_line):
                        impr1 = int(imprlist[nimp*4])
                        impr2 = int(imprlist[nimp*4+1])
                        impr3 = int(imprlist[nimp*4+2])
                        impr4 = int(imprlist[nimp*4+3])
                        impropers.append([impr1, impr2, impr3, impr4])

            elif '!NDON:' in templist or '!NDON' in templist:
                ndon = int(templist[0])
                if ndon == 0:
                    print("The structure does not contain any H-bond donors.")

            elif '!NACC:' in templist or '!NACC' in templist:
                nacc = int(templist[0])
                if nacc == 0:
                    print("The structure does not contain any H-bond acceptors.")

            elif '!NNB:' in templist or '!NNB' in templist:
                nnonb = int(templist[0])
                if nnonb == 0:
                    print("The structure does not contain any non-bonding interactions.")

    print("We have:", len(bonds), "bonds;", len(angles), "angles;", len(dihedrals), "dihedrals.")
    return bonds, ff_type, charges, masses, angles, dihedrals, impropers, donors, acceptors, nonbonds