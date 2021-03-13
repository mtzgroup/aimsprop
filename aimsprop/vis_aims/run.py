import glob
import os
import re
import sys

import numpy as np
import pkg_resources

import aimsprop as ai


def write_vmd(
    colors,  # array of rgb values
    trajectories,  # array of filenames
    states,
    opt_in={},  # dictionary of options
):

    """Writes vmd script to visualize trajectories

    Params:
        colors, np.ndarray (nstates,3), array of rgb values for each state
        trajectories, list (ntrajectories), list of xyzs filenames
        states, dictionary (ntrajectories), trajectory index to state label (0 based)
        opt_in, dictionary (options to override defaults)

    """

    header_path = pkg_resources.resource_filename(__name__, "ref/header.vmd")
    footer_path = pkg_resources.resource_filename(__name__, "ref/footer.vmd")
    molecule_path = pkg_resources.resource_filename(__name__, "ref/molecule.vmd")

    # Set User Options
    options = {
        "fn_header": header_path,
        "fn_footer": footer_path,
        "fn_mol": molecule_path,
        "fn_vmd": "vis_fms.vmd",
        "render": True,
        "exit_vmd": True,
        "opacities": True,
    }
    for key, val in list(opt_in.items()):
        options[key] = val

    fn_header = options["fn_header"]
    fn_footer = options["fn_footer"]
    fn_mol = options["fn_mol"]
    fn_vmd = options["fn_vmd"]
    render = options["render"]
    opacities = options["opacities"]
    exit_vmd = options["exit_vmd"]  # whether to render and exit vmd

    # Lines from vmd Templates
    header_lines = open(fn_header).readlines()
    footer_lines = open(fn_footer).readlines()
    mol_lines = open(fn_mol).readlines()

    # Number of trajectories
    ntrajs = len(trajectories)

    # Write the header
    # Write new material for each opacity
    # default to 1.0/ntrajs?
    fh = open(fn_vmd, "w")
    for lh in header_lines:
        mobj = re.match(r"^  # Add materials!\s*", lh)
        mobj1 = re.match(
            r"^  set mlist { Opaque Transparent BrushedMetal Diffuse Ghost Glass1 Glass2 Glass3 Glossy HardPlastic MetallicPastel Steel Translucent Edgy EdgyShiny EdgyGlass Goodsell AOShiny AOChalky AOEdgy BlownGlass GlassBubble RTChrome }",
            lh,
        )
        if mobj:
            fh.write(lh)
            for n in range(ntrajs):
                fh.write("  material change ambient Transp%d 0.000000\n" % (n + 1))
                fh.write("  material change diffuse Transp%d 0.920000\n" % (n + 1))
                fh.write("  material change specular Transp%d 0.290000\n" % (n + 1))
                fh.write("  material change shininess Transp%d 0.470000\n" % (n + 1))
                fh.write("  material change mirror Transp%d 0.000000\n" % (n + 1))
                fh.write("  material change opacity Transp%d 0.500000\n" % (n + 1))
                fh.write("  material change outline Transp%d 1.510000\n" % (n + 1))
        if mobj1:
            fh.write(
                "  set mlist { Opaque Transparent BrushedMetal Diffuse Ghost Glass1 Glass2 Glass3 Glossy HardPlastic MetallicPastel Steel Translucent Edgy EdgyShiny EdgyGlass Goodsell AOShiny AOChalky AOEdgy BlownGlass GlassBubble RTChrome "
            )
            for n in range(ntrajs):
                fh.write("Transp%d " % (n + 1))
            fh.write("}\n")
        else:
            fh.write(lh)

    # Looping through trajectories
    j = 1.0 / ntrajs
    for n in range(ntrajs):
        # Writing vmd segment for each trajecotry
        for ml in mol_lines:
            # Read in filename
            mobj1 = re.sub(r"filename", trajectories[n], ml)
            # Assign color ID number
            mobj2 = re.sub(
                r"mol color ColorID \d+", r"mol color ColorID %d" % (32 - n), mobj1
            )
            # Assign material to representation
            mobj3 = re.sub(
                r"mol material Transp\d+", r"mol material Transp%d" % (n + 1), mobj2
            )
            # Reassign RGB value to color ID
            mobj4 = re.sub(
                r"color change rgb \d+ \d\.\d+ \d\.\d+ \d\.\d+",
                r"color change rgb %d %6.5f %6.5f %6.5f"
                % (
                    32 - n,
                    colors[states[n]][0],
                    colors[states[n]][1],
                    colors[states[n]][2],
                ),
                mobj3,
            )
            fh.write("%s" % (mobj4))
        fh.write("# Finished traj %d\n" % (n + 1))

    # Writing footer for vmd file
    for lf in footer_lines:
        fh.write(lf)
    if exit_vmd:
        fh.write("source render.tcl\n")
        fh.write("exit")


def write_render(
    opt_in={},
):

    """Renders each frame of trajectory in vmd

    Params:
        opt_in, dictionary (options to override defaults)

    """

    options = {
        # visualization options
        "opacities": [],
        "end_frame": 1000,
        "increment_frame": 10,
        "sleep": 1,
        "ntrajs": 0,
        "display_size": None,
        "translation": None,
        "axis_on": False,
        "scale": 0.25,
    }
    for key, val in list(opt_in.items()):
        options[key] = val
    opacities = options["opacities"]
    end_frame = options["end_frame"]
    increment_frame = options["increment_frame"]
    ntrajs = options["ntrajs"]
    sleep = options["sleep"]
    display_size = options["display_size"]  # (1000, 600)
    translation = options["translation"]  # [-0.5, 0.0, 0.0]
    axis_on = options["axis_on"]
    scale = options["scale"]

    # TODO: Create options for:
    # Length of trajectory
    # Interval of frames (steps)
    # Translations and rotations
    # image scale and size
    # Render sleep settings

    # Number of trajectories
    fh = open("render.tcl", "w")

    # Writing Opacities based on amplitude for all trajectories
    # Setting the amplitude to zero for the time before Spawning
    # For groundstate trajectories that dropped, setting the amplitude
    # to it's last amplitude until the end of the trajectory
    for n in range(ntrajs):
        fh.write("set opa%d [list   " % (n + 1))
        for val in opacities[n]:
            fh.write("%4.3f   " % (val))
        fh.write("]\n")

    # Writing remainder of display settings
    if not axis_on:
        fh.write("axes location off\n")
    # TODO: include custom image size settings!
    fh.write("scale to %1.3f\n" % scale)
    if display_size is not None:
        fh.write("display resize %d %d\n" % (display_size[0], display_size[1]))
    fh.write("display culling off\n")
    fh.write("animate goto start\n")
    if translation is not None:
        fh.write(
            "translate 1 %1.3f %1.3f %1.3f\n"
            % (translation[0], translation[1], transpation[2])
        )

    # Writing commands to loop through the trajectories
    # TODO: include frames and increment settings
    fh.write(
        "for {set frame 0} {$frame < %d} {incr frame %d} {\n"
        % (end_frame, increment_frame)
    )
    fh.write("    animate goto $frame\n")

    # Writing commands to reset opacity
    for n in range(ntrajs):
        fh.write("    set var%d [lindex $opa%d $frame]\n" % (n + 1, n + 1))
        fh.write("    material change opacity Transp%d $var%d\n" % (n + 1, n + 1))

    # Writing commands to render tga files
    fh.write('    set fname_tga snap.[format "%04d" $frame].tga\n')
    fh.write("    render TachyonInternal snapshots/$fname_tga\n")
    fh.write("    draw delete all\n")
    # TODO: custom sleep settings
    fh.write("    sleep %1.2f\n" % sleep)
    fh.write("}\n")


# TODO: add more options
def run_vmd_render():

    """ Loads VMD and renders images, also composites into gif format """

    os.system("mkdir snapshots")
    os.system("module load VMD/1.9.2")
    print("Running VMD")
    os.system("vmd -dispdev text -e vis_fms.vmd")


def get_mp4():
    print("Writing mp4")
    print("ffmpeg must be installed - not necessarily capable of running on fire!")
    os.system(
        "ffmpeg -r 4 -i snapshots/snap.%04d.tga -c:v libx264 -pix_fmt yuv420p -r 10 snap.mp4"
    )


def get_gif():

    print("Writing gif")
    os.system("convert -delay 20 -loop 1 snapshots/snap*tga snap.gif")
