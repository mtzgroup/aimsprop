from pathlib import Path

from scripts.convert_to_geomdat import convert_geometries


def test_convert_geometries(tmp_path):
    # Get input and answer files, sorted by name
    test_data_dir = Path(__file__).parent / "test_data" / "scripts"
    prmtop_file = test_data_dir / "parm.prmtop"
    rst7_files = sorted(list(test_data_dir.glob("**/coors_*.rst7")))
    correct_geom_files = sorted(list(test_data_dir.glob("**/Geometry_*.dat")))

    # Exercise the convert_geometries function
    for rst7 in rst7_files:
        convert_geometries(prmtop_file, rst7, tmp_path)

    # Collect the created files
    outputs = sorted(list(tmp_path.glob("**/*.dat")))

    # Check that the created files match the correct Geometry.dat files
    for i, output in enumerate(outputs):
        with open(output) as f:
            created = f.read()
        with open(correct_geom_files[i]) as f:
            correct = f.read()

        # Assert that each created files matches the corresponding answer key file
        assert created == correct
