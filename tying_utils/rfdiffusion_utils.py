from typing import Dict
from pyrosetta import Pose, pose_from_pdb


def split_chains_by_residue_distance(input_pose) -> Dict[str, Pose]:
    """Given a pose that has been output by rfdiffusion into multiple poses representing each chain.

    RFDiffusion outputs PDB files with only a single chain. This function iterates through the residues
    and splits into separate chains with one chain per Pose object. It splits based on the distance between
    consecutive Ca atom positions. If the distance is greater than 10 Angstroms, it is considered a new chain.

    Returns a dictionary of chain name to Pose objects.
        Includes an additional key 'all' for the combined pose with all residues.
    """

    all_residue_pose = Pose()

    chain_poses = {}
    chain_names = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    cur_chain_index = 0
    cur_chain_pose = Pose()

    for residue_index in range(1, input_pose.total_residue()):
        r1 = input_pose.residue(residue_index)
        r2 = input_pose.residue(residue_index+1)

        cur_chain_pose.append_residue_by_jump(r1, cur_chain_pose.size())

        all_residue_pose.append_residue_by_jump(
            r1, all_residue_pose.size())


        # If we have reached the last residue, we don't need to calculate the distance
        # to the next residue
        if residue_index == input_pose.total_residue():
            continue

        ca1 = r1.xyz('CA')
        ca2 = r2.xyz('CA')
        distance = ca1.distance(ca2)

        if distance > 4:
            print('>4: ',distance, residue_index+1, residue_index+2)

        # If the 2nd residue is more than 10 Angstroms away from the 1st,
        # we start a new chain
        if distance > 10:
            #print('BIGGY', distance)
            #print(distance, residue_index+1, residue_index+2)

            chain_poses[chain_names[cur_chain_index]] = cur_chain_pose
            cur_chain_index += 1
            cur_chain_pose = Pose()

    # Append the last residue!
    cur_chain_pose.append_residue_by_jump(r2, cur_chain_pose.size())
    all_residue_pose.append_residue_by_jump(
        r2, all_residue_pose.size())

    chain_poses[chain_names[cur_chain_index]] = cur_chain_pose
    chain_poses['all'] = all_residue_pose
    return chain_poses


#def get_design_name_and_sequence_from_rf_diffusion_output(pdb_dir_path: str) -> List
def get_name_and_sequence_from_pdb(pdb_path: str) -> Dict[str, str]:
    """Given a PDB path, extract the design name and sequence from the PDB file.

    The design name is extracted from the filename, and the sequence is extracted from the PDB file.
    """
    design_num = pdb_path.split('/')[-1][:-4].replace('_', '')
    design_name = f'D{design_num}'

    rfdiffusion_output_pose = pose_from_pdb(pdb_path)

    chains_poses = split_chains_by_residue_distance(rfdiffusion_output_pose)

    sequences = {}
    for chain_name, chain_pose in chains_poses.items():
        if chain_name == 'all':
            continue

        chain_sequence = ''.join([
            chain_pose.residue(i).name1()
            for i in range(1, chain_pose.total_residue() + 1)
        ])
        sequences[chain_name] = chain_sequence

    return {
        'design_name': design_name,
        'chain_sequences': sequences
    }
