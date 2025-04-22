from typing import Dict, List
import json
import os
import numpy as np

from pyrosetta import Pose, pose_from_pdb
from pyrosetta.rosetta.core.pose import append_pose_to_pose

from rprotein_utils.rfdiffusion_utils import split_chains_by_residue_distance

DEFAULT_DESIRED_ATOMS = ['N', 'CA', 'C', 'O']

CHAIN_NAMES = 'ABCDEFGHIJ'

def get_full_sequence_from_pose( pose: Pose):
    return ''.join([
        pose.residue(i).name1()
        for i in range(1, pose.total_residue()+1)
    ])


def get_pose_details_by_chain(input_pose, desired_atoms=DEFAULT_DESIRED_ATOMS):
    chains = {}

    for residue_index in range(1, input_pose.total_residue()+1):
        residue = input_pose.residue(residue_index)

        chain_num = residue.chain()
        chain_name = CHAIN_NAMES[chain_num-1]
        if chain_name not in chains:
            chains[chain_name] = {
                'xyz': {},
                'sequence': ''
            }

        chains[chain_name]['sequence'] += residue.name1()


        for desired_atom in desired_atoms:
            atom_chain_name = f'{desired_atom}_chain_{chain_name}'

            atom = residue.atom(desired_atom)
            x,y,z = atom.xyz()

            if atom_chain_name not in chains[chain_name]['xyz']:
                chains[chain_name]['xyz'][atom_chain_name] = []

            chains[chain_name]['xyz'][atom_chain_name].append([x,y,z])

    return chains


def make_protein_mpnn_pdb_input(pose_name, pose : Pose) -> Dict:
    """Return a list of dictionaries containing pose details that will be used as
    input to ProteinMPNN.

    In the ProteinMPNN examples, this is the "parsed_pdbs.jsonl" files.

    ProteinMPNN expects a "jsonl" which is a list of json dictionaries in a flat file.
    The format of each record is:
    {
        'name': str,                                        # The name of the pose
        'num_of_chains': int,                               # The number of chains in the pose
        'seq': str,                                         # The full sequence of all chains in the pose
        'seq_chain_<chain_name>': str,                      # The sequence of each chain
        'coords_chain_<chain_name>': list[np.array[3,1]]    # The coordinates of each chain
    }
    """
    full_sequence = get_full_sequence_from_pose( pose )
    pose_details = get_pose_details_by_chain(pose)

    pose_record = {
        'name': pose_name,
        'num_of_chains': len(pose_details),
        'seq': full_sequence                 # The full sequence of all chains in the pose
    }

    for chain_name, chain_details in pose_details.items():
        pose_record[f'seq_chain_{chain_name}'] = chain_details['sequence']
        pose_record[f'coords_chain_{chain_name}'] = chain_details['xyz']

    return pose_record
