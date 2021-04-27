"""
Written by Alice Walker and Ruibin Liang

This script converts from Amber restart format to AIMS Geometry.dat files. It assumes your restart file has velocities.

To run, just give it the paths to your prmtop, restart data and desired output data as below.

main('test_data/parm.prmtop','test_data/','test_data/outs/')

"""

import argparse
from pathlib import Path

import numpy as np


def read_prmtop(prmpath: str):

    """This function parses a prmtop to extract atomic types/names and masses.

    Params:
      prmpath: the path to the prmtop file

    Returns:
       namelist: the atom names as a vector of strings.
       atommass: atomic masses as a vector of floats.

    """

    namelist = []
    atommass = []
    with open(prmpath) as f:
        for line in f:
            if "ATOM_NAME" in line:
                line = next(f)
                for line in f:
                    if "CHARGE" in line:
                        break
                    else:
                        line = line.strip("\n")
                        for c in range(0, len(line), 4):
                            namelist = namelist + [line[c : c + 3]]
            if "MASS" in line:
                line = next(f)
                for line in f:
                    if "ATOM_TYPE_INDEX" in line:
                        break
                    else:
                        a = line.split()
                        for x in a:
                            atommass.append(float(x))
    return namelist, atommass


def read_rst(rstpath: str):

    """This function parses a restart file to extract coordinates and velocities.

    Params:
      rstpath: the path to the restart file(s) to be processed. Can handle instances in subdirectories.

    Returns:
      atomcoord: the xyz coordinates as an N by 3 array.
      atomvel: the velocities as an N by 3 array.

    """

    atomcoord = np.array([])
    atomvel = np.array([])
    atomnum = 0
    with rstpath.open() as f:
        f.readline()
        line = f.readline()
        atomnum = int(line.split()[0])
        atomcoord = np.zeros((atomnum, 3))
        atomvel = np.zeros((atomnum, 3))
        for j in range(0, int((atomnum + 1) / 2)):
            line = f.readline()
            a = line.split()
            atomcoord[j * 2] = np.array([float(a[0]), float(a[1]), float(a[2])])
            if len(a) == 6:
                atomcoord[j * 2 + 1] = np.array([float(a[3]), float(a[4]), float(a[5])])
        for j in range(0, int((atomnum + 1) / 2)):
            line = f.readline()
            a = line.split()
            atomvel[j * 2] = np.array([float(a[0]), float(a[1]), float(a[2])])
            if len(a) == 6:
                atomvel[j * 2 + 1] = np.array([float(a[3]), float(a[4]), float(a[5])])
    return atomcoord, atomvel, atomnum


def output_geom(
    namelist: np.array,
    atommass: np.array,
    atomcoord: np.array,
    atomvel: np.array,
    atomnum: int,
    outpath: str,
    filename: str,
):

    """This function takes output from read_rst and read_prmtop as input and creates a Geometry.dat file, with appropriate numerical conversions to Bohr coordinates and a.u. time, to be used for AIMS simulations.

    Params:
      namelist: the vector of atom names as strings
      atommass: the vector of atomic masses as floats
      atomcoord: the 2d N by 3 array of atomic xyz coordinates
      atomvel: the 2d N by 3 array of atomic velocities
      atomnum: the total number of atoms in the system
      outpath: the directory to put all output files/directories
      filename: the name of the specific restart file to be converted

    Returns:
      A Geometry.dat file in the path outpath/filename-Geometry.dat

    """

    # Initialize new output directory if it doesn't exist

    if not outpath.exists():
        Path.mkdir(outpath)

    # Create directory structure from original files

    newdir = filename.parent
    finalpath = outpath.joinpath(newdir.relative_to(newdir.parent))
    Path.mkdir(finalpath, parents=True, exist_ok=True)

    with open(finalpath.joinpath("Geometry.dat"), "w+") as f:
        f.write("UNITS=BOHR\n")
        f.write(str(atomnum) + "\n")
        for j in range(0, atomnum):
            if namelist[j] == "Cl-":
                f.write(
                    "{0:<10s}".format(namelist[j][0:2])
                    + "".join(
                        "{0:>22.16f}".format(x / 0.52917724924) for x in atomcoord[j]
                    )
                    + "\n"
                )
            else:
                f.write(
                    "{0:<10s}".format(namelist[j][0])
                    + "".join(
                        "{0:>22.16f}".format(x / 0.52917724924) for x in atomcoord[j]
                    )
                    + "\n"
                )
        f.write("# momenta\n")
        for j in range(0, atomnum):
            for x in atomvel[j]:
                p = x * atommass[j] * 9.3499611711905144e-04 * 1.8228884855409500e03
                f.write(
                    "{0:>22.16f}".format(p)
                )  # time unit in rst7: 1/20.455 ps   length unit in rst7: angstrom    mass unit in amber: atomic mass unit
            f.write("\n")


def convert_geometries(
    prmpath: str,
    rstpath: str,
    outpath: str,
):

    """This function creates a directory structure and output to convert a set of geometries from Amber prmtop/rst format to FMS90 Geometry.dat format.

    Params:
      prmpath: the path to the prmtop file
      rstpath: the path to the restart file
      outpath: the path to the output files

    Returns:
      All converted geometries in the path outpath/filename.
    """

    findme = "*/*.rst*"

    Prmpath = Path(prmpath)
    Rstpath = Path(rstpath)
    Outpath = Path(outpath)

    # Get initial information from parmtop
    names, masses = read_prmtop(Prmpath)

    # Find all restart files in test data
    filenames = list(Rstpath.glob(findme))

    # Create new dats
    for filename in filenames:
        coords, vel, numats = read_rst(filename)
        output_geom(names, masses, coords, vel, numats, Outpath, filename)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-p", "--prm", dest="prmpath", help="/path/to/input/prmtop", required=True
    )
    parser.add_argument(
        "-r", "--rst", dest="rstpath", help="/path/to/input/restarts", required=True
    )
    parser.add_argument(
        "-o", "--out", dest="outpath", help="/path/to/outputs/", required=True
    )
    args = parser.parse_args()
    convert_geometries(args.prmpath, args.rstpath, args.outpath)
