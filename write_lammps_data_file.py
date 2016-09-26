def write_lammps_file_periodic(coords, bonds, angles, dihedrals, ff_type,
                               charges, masses, filebase, reprod, box, types):

    nat = len(coords)
    nbonds = len(bonds)
    nangles = len(angles)
    ndihedrals = len(dihedrals)
    # nimpropers = len(impropers)
    nx, ny = reprod[0][:]

    output_file = filebase + ".data"
    print(output_file)
    outname = open(output_file, "w")

    # Comment line.
    print("LAMMPS data file. SlabCreator v0.1.", file=outname, end="\n")
    print("", file=outname)
    # How many atoms, bonds, etc
    print("   ", nat, " atoms", file=outname)
    print("   ", nbonds, " bonds", file=outname)
    print("   ", nangles, " angles", file=outname)
    print("   ", ndihedrals, " dihedrals", file=outname)
    print("   0 impropers", file=outname)
    print("", file=outname)
    # How many bond types, angle types, etc
    print("   4 atom types", file=outname)
    print("   3 bond types", file=outname)
    print("   5 angle types", file=outname)
    print("   4 dihedral types", file=outname)
    print("", file=outname)
    print("   26 extra dihedral per atom", file=outname)
    print("", file=outname)

    # The box dimension
    print("    0.0", "{:10.6f}".format(nx * box[0, 0]), "xlo xhi", sep="  ", file=outname)
    print("    0.0", "{:10.6f}".format(ny * box[0, 1]), "ylo yhi", sep="  ", file=outname)
    print("    0.0", "{:10.6f}".format(box[0, 2]), "zlo zhi", sep="  ", file=outname)
    print("", file=outname)
    print("Masses", file=outname)
    print("", file=outname)
    print("   1  28.086000 # sc4", file=outname)
    print("   2  15.999400 # oc23", file=outname)
    print("   3  15.999400 # oc24", file=outname)
    print("   4   1.008000 # hoy", file=outname)
    print("", file=outname)
    print("Pair Coeffs # lj/cut/coul/long", file=outname)
    print("", file=outname)
    print("   1   0.0930001390   4.15 # sc4", file=outname)
    print("   2   0.0540006496   3.47 # oc23", file=outname)
    print("   3   0.1220010554   3.47 # oc24", file=outname)
    print("   4   0.0149825188   1.085 # hoy", file=outname)
    print("", file=outname)
    print("Bond Coeffs # harmonic", file=outname)
    print("", file=outname)
    print("  1   285.0000     1.6800  #  sc4-oc23", file=outname)
    print("  2   285.0000     1.6800  #  sc4-oc24", file=outname)
    print("  3   495.0000     0.9450  #  oc24-hoy", file=outname)
    print("", file=outname)
    print("Angle Coeffs # harmonic", file=outname)
    print("", file=outname)
    print("  1   100.0000   109.5000  #  oc23-sc4-oc23", file=outname)
    print("  2   100.0000   109.5000  #  oc23-sc4-oc24", file=outname)
    print("  3   100.0000   149.0000  #  sc4-oc23-sc4", file=outname)
    print("  4    50.0000   115.0000  #  sc4-oc24-hoy", file=outname)
    print("  5   100.0000   109.5000  #  oc24-sc4-oc24", file=outname)
    print("", file=outname)
    print("Dihedral Coeffs # harmonic", file=outname)
    print("", file=outname)
    print("  1     0.0000   1   2 # oc23-sc4-oc23-sc4", file=outname)
    print("  2     0.0000   1   2 # oc24-sc4-oc23-sc4", file=outname)
    print("  3     0.0000   1   2 # oc23-sc4-oc24-hoy", file=outname)
    print("  4     0.0000   1   2 # oc24-sc4-oc24-hoy", file=outname)
    print("", file=outname)
    # print("Improper Coeffs # harmonic", file=outname)
    # print("", file=outname)
    print("Atoms # full", file=outname)
    print("", file=outname)
    for i in range(len(coords)):
        current_type = 0
        if coords[i][1] == "SC4":
            current_type = 1
        elif coords[i][1] == "OC23":
            current_type = 2
        elif coords[i][1] == "OC24":
            current_type = 3
        elif coords[i][1] == "HOY":
            current_type = 4

        print("{0:>9}".format(i + 1), "  1", "{0:>4}".format(current_type),
              "{:10.6f}".format(coords[i][2]),
              "{:17.10f}".format(coords[i][3]),
              "{:17.10f}".format(coords[i][4]),
              "{:17.10f}".format(coords[i][5]),
              "  0   0   0 # ", coords[i][1].lower(), sep=" ", file=outname)
    print("", file=outname)
    print("Bonds", file=outname)
    print("", file=outname)
    n = 1
    for nb in bonds:
        current_bond = 0
        bond_type = nb[1][0] + nb[1][1]
        if bond_type == "sc4oc23" or bond_type == "oc23sc4":
            # "sc4" in nb[1] and "oc23" in nb[1]:
            current_bond = 1
        elif bond_type == "sc4oc24" or bond_type == "oc24sc4":
            # "sc4" in nb[1] and "oc24" in nb[1]:
            current_bond = 2
        # elif "hoy" in nb[1] and "oc23" in nb[1]:
        #     current_bond = 3
        elif bond_type == "hoyoc24" or bond_type == "oc24hoy":
        # elif "hoy" in nb[1] and "oc24" in nb[1]:
            current_bond = 3
        print("{0:>9}".format(n), "{0:>4}".format(current_bond),
              "{0:>9}".format(nb[0][0]),
              "{0:>9}".format(nb[0][1]),
              file=outname)
        if current_bond == 0:
            print(nb)
        n += 1
    print("", file=outname)
    print("Angles", file=outname)
    print("", file=outname)
    n = 1
    for nang in angles:
        current_angle = 0
        ang_type = nang[1][0] + "-" + nang[1][1] + "-" + nang[1][2]
        if ang_type == "oc23-sc4-oc23":
            current_angle = 1
        elif ang_type == "oc23-sc4-oc24":
            current_angle = 2
        elif ang_type == "sc4-oc23-sc4":
            current_angle = 3
        elif ang_type == "sc4-oc24-hoy":
            current_angle = 4
        elif ang_type == "oc24-sc4-oc24":
            current_angle = 5
        if current_angle == 0:
            print(nang)
        print("{0:>9}".format(n), "{0:>4}".format(current_angle),
              "{0:>9}".format(nang[0][0]),
              "{0:>9}".format(nang[0][1]),
              "{0:>9}".format(nang[0][2]),
              file=outname)
        n += 1
    print("", file=outname)
    print("Dihedrals", file=outname)
    print("", file=outname)
    n = 1
    for ndhd in dihedrals:
        current_dhd = 0
        dhd_type = ndhd[1][0] + "-" + ndhd[1][1] + "-" + ndhd[1][2] + "-" + ndhd[1][3]
        if dhd_type == "oc23-sc4-oc23-sc4" or dhd_type == "sc4-oc23-sc4-oc23":
            current_dhd = 1
        elif dhd_type == "oc24-sc4-oc23-sc4" or dhd_type == "sc4-oc23-sc4-oc24":
            current_dhd = 2
        elif dhd_type == "oc23-sc4-oc24-hoy" or dhd_type == "hoy-oc24-sc4-oc23":
            current_dhd = 3
        elif dhd_type == "oc24-sc4-oc24-hoy" or dhd_type == "hoy-oc24-sc4-oc24":
            current_dhd = 4
        print("{0:>9}".format(n), "{0:>4}".format(current_dhd),
              "{0:>9}".format(ndhd[0][0]),
              "{0:>9}".format(ndhd[0][1]),
              "{0:>9}".format(ndhd[0][2]),
              "{0:>9}".format(ndhd[0][3]),
              file=outname)
        n += 1
    print("", file=outname)


def write_lammps_file_nonperiodic(single_mol_info, types,
                                  coords, bonds, angles, dihedrals, impropers,
                                  donors, acceptors, nonbonds,
                                  filename, box):

    nat_per_mol = len(single_mol_info)
    nat = len(coords)
    nmol = int(nat / nat_per_mol)

    nbonds = len(bonds)
    nangles = len(angles)
    ndihedrals = len(dihedrals)
    nimpropers = len(impropers)

    output_file = filename + ".data"
    print(output_file)
    outname = open(output_file, "w")

    bond_type = []
    for ibond in range(nbonds):
        b1 = bonds[ibond][0]
        b2 = bonds[ibond][1]
        if single_mol_info[b1-1][1] == "c2" and single_mol_info[b2-1][1] == "c2":
            bond_type.append(1)
        elif single_mol_info[b1-1][1] == "c2" and single_mol_info[b2-1][1] == "h" \
                or single_mol_info[b1-1][1] == "h" and single_mol_info[b2-1][1] == "c2":
            bond_type.append(2)
        elif single_mol_info[b1-1][1] == "c2" and single_mol_info[b2-1][1] == "c3" \
                or single_mol_info[b1-1][1] == "c3" and single_mol_info[b2-1][1] == "c2":
            bond_type.append(3)
        elif single_mol_info[b1-1][1] == "c3" and single_mol_info[b2-1][1] == "h" \
                or single_mol_info[b1-1][1] == "h" and single_mol_info[b2-1][1] == "c3":
            bond_type.append(4)

    angle_type = []
    for iangle in range(nangles):
        a1 = angles[iangle][0]
        a2 = angles[iangle][1]
        a3 = angles[iangle][2]

        ang1 = single_mol_info[a1 - 1][1]
        ang2 = single_mol_info[a2 - 1][1]
        ang3 = single_mol_info[a3 - 1][1]
        ang = ang1 + ang2 + ang3

        if ang == "c2c2c2":
            angle_type.append(1)
            continue
        elif ang == "c2c2h" or ang == "hc2c2":
            angle_type.append(2)
            continue
        elif ang == "hc2h":
            angle_type.append(3)
            continue
        elif ang == "c2c2c3" or ang == "c3c2c2":
            angle_type.append(4)
            continue
        elif ang == "c3c2h" or ang == "hc2c3":
            angle_type.append(5)
            continue
        elif ang == "c2c3h" or ang == "hc3c2":
            angle_type.append(6)
            continue
        elif ang == "hc3h":
            angle_type.append(7)
            continue

    dihedral_type = []
    for idihed in range(ndihedrals):
        d1 = dihedrals[idihed][0]
        d2 = dihedrals[idihed][1]
        d3 = dihedrals[idihed][2]
        d4 = dihedrals[idihed][3]

        dihed1 = single_mol_info[d1 - 1][1]
        dihed2 = single_mol_info[d2 - 1][1]
        dihed3 = single_mol_info[d3 - 1][1]
        dihed4 = single_mol_info[d4 - 1][1]
        dihed = dihed1 + dihed2 + dihed3 + dihed4
        # print(dihed)
        if dihed == "c2c2c2c2":
            dihedral_type.append(1)
        elif dihed == "c2c2c2h" or dihed == "hc2c2c2":
            dihedral_type.append(2)
        elif dihed == "hc2c2h":
            dihedral_type.append(3)
        elif dihed == "c2c2c2c3" or dihed == "c3c2c2c2":
            dihedral_type.append(4)
        elif dihed == "c3c2c2h" or dihed == "hc2c2c3":
            dihedral_type.append(5)
        elif dihed == "c2c2c3h" or dihed == "hc3c2c2":
            dihedral_type.append(6)
        elif dihed == "hc2c3h" or dihed == "hc3c2h":
            dihedral_type.append(7)

    # Comment line.
    print("LAMMPS data file. SlabCreator v0.1.", file=outname, end="\n")
    print("", file=outname)
    # How many atoms, bonds, etc
    print("   ", nat, " atoms", file=outname)
    print("   ", nbonds * nmol, " bonds", file=outname)
    print("   ", nangles * nmol, " angles", file=outname)
    print("   ", ndihedrals * nmol, " dihedrals", file=outname)
    print("   ", nimpropers * nmol, " impropers", file=outname)
    print("", file=outname)
    # How many bond types, angle types, etc
    print("   ", len(types), "atom types", file=outname)
    print("    4 bond types", file=outname)
    print("    7 angle types", file=outname)
    print("    7 dihedral types", file=outname)
    print("", file=outname)

    # The box dimension
    print("    0.0", "{:10.6f}".format(box[0, 0]), "xlo xhi", sep="  ", file=outname)
    print("    0.0", "{:10.6f}".format(box[0, 1]), "ylo yhi", sep="  ", file=outname)
    print("    0.0", "{:10.6f}".format(box[0, 2]), "zlo zhi", sep="  ", file=outname)
    print("", file=outname)
    print("Masses", file=outname)
    print("", file=outname)
    for i in range(len(types)):
        print("   ", i + 1, types[i][3], " # ", types[i][1], file=outname)

    print("", file=outname)
    print("Pair Coeffs # lj/cut/coul/cut", file=outname)
    print("", file=outname)
    print("   1   0.0389999952   3.8754094636 # c2", file=outname)
    print("   2   0.0389999952   3.8754094636 # c3", file=outname)
    print("   3   0.0380000011   2.4499714540 # h", file=outname)
    print("", file=outname)
    print("Bond Coeffs # harmonic", file=outname)
    print("", file=outname)
    print("  1   322.7158     1.5260 # c2-c2", file=outname)
    print("  2   340.6175     1.1050 # c2-h", file=outname)
    print("  3   322.7158     1.5260 # c2-c3", file=outname)
    print("  4   340.6175     1.1050 # c3-h", file=outname)
    print("", file=outname)
    print("Angle Coeffs # harmonic", file=outname)
    print("", file=outname)
    print("  1   46.6000   110.5000 # c2-c2-c2", file=outname)
    print("  2   44.4000   110.0000 # c2-c2-h", file=outname)
    print("  3   39.5000   106.4000 # h-c2-h", file=outname)
    print("  4   46.6000   110.5000 # c2-c2-c3", file=outname)
    print("  5   44.4000   110.0000 # c3-c2-h", file=outname)
    print("  6   44.4000   110.0000 # c2-c3-h", file=outname)
    print("  7   39.5000   106.4000 # h-c3-h", file=outname)
    print("", file=outname)
    print("Dihedral Coeffs # harmonic", file=outname)
    print("", file=outname)
    print("  1     0.1581   1   3 # c2-c2-c2-c2", file=outname)
    print("  2     0.1581   1   3 # c2-c2-c2-h", file=outname)
    print("  3     0.1581   1   3 # h-c2-c2-h", file=outname)
    print("  4     0.1581   1   3 # c2-c2-c2-c3", file=outname)
    print("  5     0.1581   1   3 # c3-c2-c2-h", file=outname)
    print("  6     0.1581   1   3 # c2-c2-c3-h", file=outname)
    print("  7     0.1581   1   3 # h-c2-c3-h", file=outname)
    print("", file=outname)
    # print("Improper Coeffs # harmonic", file=outname)
    # print("", file=outname)
    print("Atoms # full", file=outname)
    print("", file=outname)
    for imol in range(nmol):
        for iat in range(nat_per_mol):
            ntype = 0
            for itype in (itype for itype, x in enumerate(types) if x == single_mol_info[iat]):
                ntype = itype + 1
            print("{:10d}".format(imol * nat_per_mol + iat + 1),
                  "{:6d}".format(imol + 1),
                  "{:4d}".format(ntype),
                  "{:9.5f}".format(single_mol_info[iat][2]),
                  "{:17.10f}".format(coords[imol * nat_per_mol + iat][0]),
                  "{:17.10f}".format(coords[imol * nat_per_mol + iat][1]),
                  "{:17.10f}".format(coords[imol * nat_per_mol + iat][2]),
                  "  0   0   0 # ", single_mol_info[iat][1],
                  file=outname)

    print("", file=outname)
    print("Bonds", file=outname)
    print("", file=outname)
    for imol in range(nmol):
        for ibond in range(nbonds):
            print("{:10d}".format(imol * nbonds + ibond),
                  "{:5d}".format(bond_type[ibond]),
                  "{:10d}".format(bonds[ibond][0] + imol * nat_per_mol),
                  "{:10d}".format(bonds[ibond][1] + imol * nat_per_mol),
                  file=outname)
    print("", file=outname)
    print("Angles", file=outname)
    print("", file=outname)
    for imol in range(nmol):
        for iangle in range(nangles):
            print("{:10d}".format(imol * nangles + iangle),
                  "{:5d}".format(angle_type[iangle]),
                  "{:10d}".format(angles[iangle][0] + imol * nat_per_mol),
                  "{:10d}".format(angles[iangle][1] + imol * nat_per_mol),
                  "{:10d}".format(angles[iangle][2] + imol * nat_per_mol),
                  file=outname)

    print("", file=outname)
    print("Dihedrals", file=outname)
    print("", file=outname)
    for imol in range(nmol):
        for idihed in range(ndihedrals):
            print("{:10d}".format(imol * ndihedrals + idihed),
                  "{:5d}".format(dihedral_type[idihed]),
                  "{:10d}".format(dihedrals[idihed][0] + imol * nat_per_mol),
                  "{:10d}".format(dihedrals[idihed][1] + imol * nat_per_mol),
                  "{:10d}".format(dihedrals[idihed][2] + imol * nat_per_mol),
                  "{:10d}".format(dihedrals[idihed][3] + imol * nat_per_mol),
                  file=outname)
    print("", file=outname)