def create_protein_mpnn_run_command(
    protein_mpnn_install_dir_path : str,
    input_dir_path : str,
    output_dir_path : str,
    num_sequence_per_target : int = 2,
    sampling_temp : float = 0.2,
    batch_size : int = 1,
    seed : int = None,
    python_cmd_or_path : str = 'python'
):
    """Create the shell command to run ProteinMPNN."""

    protein_mpnn_args = {
        'jsonl_path': f'{input_dir_path}/parsed_pdbs.jsonl',
        'tied_positions_jsonl': f'{input_dir_path}/tied_pdbs.jsonl',
        'fixed_positions_jsonl': f'{input_dir_path}/fixed_pdbs.jsonl',

        'out_folder': output_dir_path,
        'num_seq_per_target': f'{num_sequence_per_target}',
        'sampling_temp': str(sampling_temp),
        'batch_size': str(batch_size),
    }

    if seed:
        protein_mpnn_args['seed'] = int(seed)

    py_file_path = f'{protein_mpnn_install_dir_path}/protein_mpnn_run.py'

    cmd = [
        python_cmd_or_path,
        py_file_path
    ]
    for arg_name, arg_value in protein_mpnn_args.items():
        cmd.append(f'--{arg_name}')
        cmd.append(arg_value)

    return ' '.join(cmd)
