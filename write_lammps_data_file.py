def write_lammps_file(coords, bonds, angles, dihedrals, ff_type,
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
    print("    0  impropers", file=outname)
    # print("   ", nimpropers, " impropers", file=outname)
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
    print("   1   0.0930001390   3.6972284612 # sc4", file=outname)
    print("   2   0.0540006496   3.0914157928 # oc23", file=outname)
    print("   3   0.1220010554   3.0914166375 # oc24", file=outname)
    print("   4   0.0149825188   0.9666679938 # hoy", file=outname)
    print("", file=outname)
    print("Bond Coeffs # harmonic", file=outname)
    print("", file=outname)
    print("  1   285.0000     1.6800  #  sc4-oc23", file=outname)
    print("  2   285.0000     1.6800  #  sc4-oc24", file=outname)
    print("  3   495.0000     0.9450  #  oc23-hoy", file=outname)
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
        if "sc4" in nb[1] and "oc23" in nb[1]:
            current_bond = 1
        elif "sc4" in nb[1] and "oc24" in nb[1]:
            current_bond = 2
        elif "hoy" in nb[1] and "oc23" in nb[1]:
            current_bond = 3
        elif "hoy" in nb[1] and "oc24" in nb[1]:
            current_bond = 3
        print("{0:>9}".format(n), "{0:>4}".format(current_bond),
              "{0:>9}".format(nb[0][0]),
              "{0:>9}".format(nb[0][1]),
              file=outname)
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
        if current_angle == 0: print(nang)
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
        if dhd_type == "oc23-sc4-oc23-sc4":
            current_dhd = 1
        elif dhd_type == "oc24-sc4-oc23-sc4":
            current_dhd = 2
        elif dhd_type == "oc23-sc4-oc24-hoy":
            current_dhd = 3
        elif dhd_type == "oc24-sc4-oc24-hoy":
            current_dhd = 4
        print("{0:>9}".format(n), "{0:>4}".format(current_dhd),
              "{0:>9}".format(ndhd[0][0]),
              "{0:>9}".format(ndhd[0][1]),
              "{0:>9}".format(ndhd[0][2]),
              "{0:>9}".format(ndhd[0][3]),
              file=outname)
        n += 1
    print("", file=outname)


