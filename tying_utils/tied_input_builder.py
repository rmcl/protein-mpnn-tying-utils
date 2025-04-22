import os
import json
from typing import Dict, List, Tuple
from pyrosetta import Pose

class ProteinMPNNInputFileBuilder:
    """A class to build up tied and fixed input files for ProteinMPNN.

    If you pass a directory that already exists and has these files, this class
    will load the existing records and append to them.

    ProteinMPNN expects three files:
    - parsed_pdbs.jsonl - basically design pdb files converted to "json list" files.
    - tied_pdbs.jsonl - a dictionary of tied residues by chain.
    - fixed_pdbs.jsonl - a dictionary of fixed residues by chain.

    """

    def __init__(self, output_dir_path : str):
        self.output_dir_path = output_dir_path

        self._existing_pdb_names = set()
        self._added_pdbs = []

        if not os.path.exists(self.output_dir_path):
            os.makedirs(self.output_dir_path)
        else:
            self.load_existing_records()


    @property
    def tied_output_file_path(self):
        return os.path.join(self.output_dir_path, 'tied_pdbs.jsonl')

    @property
    def fixed_output_file_path(self):
        return os.path.join(self.output_dir_path, 'fixed_pdbs.jsonl')

    @property
    def output_parsed_pdb_path(self):
        return os.path.join(self.output_dir_path, 'parsed_pdbs.jsonl')

    def load_existing_records(self):
        """Load existing records from the tied_pdbs, and ."""
        self._all_tied_records = load_and_merge_existing_json_records(
            self.tied_output_file_path)

        self._all_fixed_records = load_and_merge_existing_json_records(
            self.fixed_output_file_path)

        try:
            with open(self.output_parsed_pdb_path, 'r') as f:
                for line in f.readlines():
                    if line.strip() == '':
                        continue

                    parsed_pdb = json.loads(line)
                    self._existing_pdb_names.add(parsed_pdb['name'])
        except FileNotFoundError:
            pass


    def store(self):
        """Store the tied and fixed records to the appropriate files."""
        with open(self.tied_output_file_path, 'w+') as f:
            f.write(json.dumps(self._all_tied_records) + '\n')

        with open(self.fixed_output_file_path, 'w+') as f:
            f.write(json.dumps(self._all_fixed_records) + '\n')

        with open(self.output_parsed_pdb_path, 'a+') as f:
            for parsed_pdb in self._added_pdbs:
                f.write(json.dumps(parsed_pdb) + '\n')

            # Clear the added pdbs after writing to file
            self._added_pdbs = []


    def add_tied_and_fixed_pdb(
        self,
        pdb_name : str,
        parsed_pdb : Dict,
        tied_residues_by_chain,
        fixed_residues_by_chain,
    ):
        """Add a new PDB to the records."""

        if pdb_name in self._existing_pdb_names:
            raise ValueError(f"PDB {pdb_name} already exists in records!")

        self._all_tied_records[pdb_name] = tied_residues_by_chain
        self._all_fixed_records[pdb_name] = fixed_residues_by_chain

        self._added_pdbs.append(parsed_pdb)


def load_and_merge_existing_json_records( jsonl_path: str) -> Dict:
    """Load existing json records from a jsonl file and merge them into a single dictionary."""
    all_existing_records = {}
    try:
        with open(jsonl_path, 'r') as f:
            for line in f.readlines():
                if line.strip() == '':
                    continue

                existing_records = json.loads(line)
                all_existing_records.update(existing_records)
    except FileNotFoundError:
        pass

    return all_existing_records


def create_fixed_and_tied_residue_sets(combined_pose : Pose, designed_residues : List[int], tied_chains : List[Tuple[str]]) -> Dict:
    """Given a Pose object create a dictionary of fixed and tied residues.

    This method expects a pose with two or more chains. The last chain is assumed
    to be the combined chain with all residues from other chains.

    Returns a dictionary containing two
        fixed_residues_by_chain:
            A dictionary containing the fixed residues by chain. These are the residues
            which ProteinMPNN should not modify
            {"A": [1, 2, 3], "B": [1, 2, 3]}

        tied_residues_by_chain: A list of dictionaries containing the tied residues
            [{"A": [1], "C": [1]}, {"A": [2], "C": [2]}]

    """

    CHAIN_NAMES = 'ABCDEFGH'

    #assert tied_chains, [
    #    ('A', 'C'),
    #    ('B', 'D'),
    #]
    all_tied_chains = []
    known_tied_chains = set()
    for row in tied_chains:
        # If the weight is not specified, set it to 1.0
        if len(row) == 2:
            chain_1 = row[0]
            chain_2 = row[1]
            weight = 1.0
        elif len(row) == 3:
            chain_1 = row[0]
            chain_2 = row[1]
            weight = row[2]

        known_tied_chains.add(chain_1)
        known_tied_chains.add(chain_2)

        all_tied_chains.append((row[0], row[1], weight))

    #print('TOTAL COMBINED LENGTH', combined_pose.total_residue())

    # Calulate the length of each chain
    chain_lengths = {}
    chain_starts = {}
    for i in range(1, combined_pose.total_residue()+1):

        chain = combined_pose.residue(i).chain()
        chain_name = CHAIN_NAMES[chain-1]

        if chain_name not in chain_lengths:
            chain_lengths[chain_name] = 0
            chain_starts[chain_name] = i

        chain_lengths[chain_name] += 1



    print( 'CHAIN LENGTHS', chain_lengths)
    print( 'CHAIN STARTS', chain_starts)

    tied_residues_by_chain = []
    for chain_1, chain_2, weight in all_tied_chains:

        #print( 'TIED CHAINS', chain_1, chain_2)
        if chain_1 not in chain_lengths or chain_2 not in chain_lengths:
            raise ValueError(f"Invalid chain names: {chain_1}, {chain_2}")

        total_chain_residue_len = chain_lengths[chain_1]
        chain_start_pos = chain_starts[chain_1]

        assert total_chain_residue_len == chain_lengths[chain_2], f"Chains {chain_1} ({total_chain_residue_len}) and {chain_2} ({chain_lengths[chain_2]}) must have the same length"

        for i in range(1, total_chain_residue_len + 1):

            # Input is a residue position of original design where residue indexes are all 1-relative to start of pose
            chain_pos = i + chain_start_pos

            #if chain_pos in designed_residue_set:
            #    continue

            tied_residues_by_chain.append({
                chain_1: [[i], [weight]],
                chain_2: [[i], [weight]]
            })

    designed_residues_set = set(designed_residues)

    fixed_residues_by_chain = {}
    """
    NOT USING THIS RIGHT NOW AND IT BROKE SO REVISIT IF NEEDED
    for i in range(1, combined_pose.total_residue()+1):

        if i not in designed_residues_set:
            chain_num = combined_pose.residue(i).chain()
            chain_name = CHAIN_NAMES[chain_num-1]

            # LETS SKIP THE CHAINS THAT ARE TIED TO OTHER CHAINS.
            # ONLY ITERATE OVER TIED CHAINS ONCE
            # also tied represents chains that are only the target of tying
            print('KNOWN TIED', known_tied_chains, chain_name)
            if chain_name not in known_tied_chains:
                continue

            if chain_name not in fixed_residues_by_chain:
                fixed_residues_by_chain[chain_name] = []

            # Adjust I to be relative to the start of the chain
            actual_chain_pos = i - chain_starts[chain_name]

            fixed_residues_by_chain[chain_name].append(actual_chain_pos)

            print(chain_name, actual_chain_pos)

            # TODO: SOMETHING MAY BE WRONG HERE! WE HAVEN'T BEEN USING FIXED RESIDUES SO IGNORE FOR NOW!
            #raise Exception('AAHH MAY BE PROBLEM HERE')
            #fixed_residues_by_chain[tied_chain_name].append(actual_chain_pos)

    """

    return {
        'fixed_residues_by_chain': fixed_residues_by_chain,
        'tied_residues_by_chain': tied_residues_by_chain
    }
