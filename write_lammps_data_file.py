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
    print("Bond Coeffs # harmonic", file=outname)
    print("", file=outname)
    print("Angle Coeffs # harmonic", file=outname)
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
            current_type = 1
        elif coords[i][1] == "OC23":
            current_type = 2
        elif coords[i][1] == "OC24":
            current_type = 3
        elif coords[i][1] == "HOY":
            current_type = 4

        print("{0:>9}".format(i + 1), "  1", "{0:>4}".format(current_type),
              "{:17.10f}".format(coords[i][2]),
              "{:17.10f}".format(coords[i][3]),
              "{:17.10f}".format(coords[i][4]),
              "  0   0   0 # ", coords[i][1].lower(), sep=" ", file=outname)
    print("", file=outname)

