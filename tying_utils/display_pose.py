"""Use py3Dmol to display a pose in a Jupyter notebook."""
from typing import Optional, Dict
import py3Dmol
import tempfile


def pdb_string_from_pose(input_pose):
    with tempfile.NamedTemporaryFile(suffix='.pdb', delete=True) as temp_file:
        temp_file_path = temp_file.name

        input_pose.dump_pdb(temp_file_path)

        with open(temp_file_path, 'r') as f:
            return f.read()

def display_pose(input_pose, show=True, show_axes=True, residue_colors : Optional[Dict] = None):
    """Display a pose in a Jupyter notebook using py3Dmol."""
    pdb_str = pdb_string_from_pose(input_pose)

    view = py3Dmol.view(width=500, height=500)
    view.addModel(pdb_str, "pdb")
    view.setStyle({"cartoon": {"color": "white"}})
    view.addStyle({"hetflag": True}, {"stick": {}})

    # Define axis length
    axis_length = 50.0

    if show_axes:

        # X-axis (red)
        view.addLine({"start": {"x": -axis_length, "y": 0, "z": 0},
                    "end": {"x": axis_length, "y": 0, "z": 0},
                    "color": "red"})

        # Y-axis (green)
        view.addLine({
            "start": {"x": 0, "y": -axis_length, "z": 0},
            "end": {"x": 0, "y": axis_length, "z": 0},
            "color": "green"})

        # Z-axis (blue)
        view.addLine(
            {"start": {"x": 0, "y": 0, "z": -axis_length},
            "end": {"x": 0, "y": 0, "z": axis_length},
            "color": "blue"})


    if residue_colors:
        # Apply different colors to specific residues
        for res_id, color in residue_colors.items():
            view.setStyle({"resi": str(res_id)}, {"cartoon": {"color": color}})


    #view.setStyle({"cartoon": {"color": "spectrum"}})
    view.zoomTo()

    if show:
        view.show()

    return view
