def write_lammps_file(coords, bonds, angles, dihedrals, impropers,
                      ff_type, charges, masses, filebase, reprod, box, types):

    nat = len(coords)
    nbonds = len(bonds)
    nangles = len(angles)
    ndihedrals = len(dihedrals)
    nimpropers = len(impropers)
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
    print("   ", nimpropers, " impropers", file=outname)
    print("", file=outname)
    # How many bond types, angle types, etc

    # The box dimension
    print("    0.0", nx * box[0, 0], "xlo  xhi", sep="  ", file=outname)
    print("    0.0", ny * box[0, 1], "ylo  yhi", sep="  ", file=outname)
    print("    0.0", box[0, 2], "zlo  zhi", sep="  ", file=outname)
    print("", file=outname)
    print("Masses", file=outname)
    print("", file=outname)
    print("Pair Coeffs # lj/cut/coul/long", file=outname)
    print("", file=outname)
    n = 1
    for typ in types:
        print(n, " # ", typ, file=outname)
        n += 1
    print("", file=outname)
    print("Bond Coeffs # harmonic", file=outname)
    print("", file=outname)
    print("  1   #  sc4-oc23", file=outname)
    print("  2   #  sc4-oc24", file=outname)
    print("  3   #  oc23-hoy", file=outname)
    print("  4   #  oc24-hoy", file=outname)
    print("", file=outname)
    print("Angle Coeffs # harmonic", file=outname)
    print("", file=outname)
    print("  1   #  oc23-sc4-oc23", file=outname)
    print("  2   #  oc24-sc4-oc24", file=outname)
    print("  3   #  oc23-sc4-oc24", file=outname)
    print("  4   #  sc4-oc23-sc4", file=outname)
    print("  5   #  sc4-oc24-hoy", file=outname)
    print("  6   #  sc4-oc23-hoy", file=outname)
    print("", file=outname)
    print("Dihedral Coeffs # harmonic", file=outname)
    print("", file=outname)
    print("Improper Coeffs # harmonic", file=outname)
    print("", file=outname)
    print("Atoms # full", file=outname)
    print("", file=outname)
    for i in range(len(coords)):
        current_type = 0
        if coords[i][1] == "SC4":
            current_type = types.index("sc4") + 1
        elif coords[i][1] == "OC23":
            current_type = types.index("oc23") + 1
        elif coords[i][1] == "OC24":
            current_type = types.index("oc24") + 1
        elif coords[i][1] == "HOY":
            current_type = types.index("hoy") + 1

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
            current_bond = 4
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
        ang_type = nang[0] + "-" + nang[1] + "-" + nang[2]
        if ang_type == "oc23-sc4-oc23":
            current_angle = 1
        elif ang_type == "oc24-sc4-oc24":
            current_angle = 2
        elif ang_type == "oc23-sc4-oc24" or ang_type == "oc24-sc4-oc23":
            current_angle = 3
        elif ang_type == "sc4-oc23-sc4":
            current_angle = 4
        elif ang_type == "sc4-oc24-hoy":
            current_angle = 5
        elif ang_type == "sc4-oc23-hoy":
            current_angle = 6


