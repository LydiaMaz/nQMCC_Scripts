import os

from wavefunction_input import WavefunctionInput, VariationalWavefunctionInput, ProductWavefunctionInput, DeuteronWavefunctionInput

# Class to group control parameters together
class ControlGroup:
    def __init__(self, owner, param_names):
        self.owner = owner
        self.param_names = param_names

    def write_output(self, stream=None):
        if stream is None:
            stream = print
        values = []
        labels = []
        for param in self.param_names:
            val = self.owner.get(param, '')
            if isinstance(val, WavefunctionInput): # if it's a wavefunction input, skip it
                continue
            values.append(str(val))
            labels.append(param)
        if not values:  # skip empty groups
            return
        values_str = "  ".join(values)
        label_str = "  ".join(labels)
        spacing = max(80, len(values_str) + 2)
        line = f"{values_str:<{spacing}}{label_str}"
        stream(line.rstrip() + "\n")

# Class to read and write the control file
class Control:

    int_values = {
        'read_write_walk', 
        'ke', 'lemp',
        'rng_seed',
        'num_blocks', 'block_size',
        'burn_in', 'num_moves_between',
        'npts', 'one_over_dx', 'fd_factor',
        'lastp_sample_type',
        'stop_after_group',
        'nlopt_method',  'num_opt_walks',  'num_opt_evaluations', 'walkers_per_node'
    }

    strings = {
        'basis', 'walk_file', 'constants_file', '2b_file', '3b_file', 'group_file', 'nortab_file', 'optimization_input', 'optimized_deck', 'scratch_dir'
    }

    boolean_values = {
        'bra_eq_ket', 'sample_L2'
    }


    def __init__(self):
        self.Name = ""
        self.parameters = {}
        self.sections = [] # stores ControlGroup or WavefunctionInput objects

    # helper functions
    def parse_group(self, group, tokens):
        for key, value in zip(group, tokens):
            if key in self.int_values:
                parsed = int(value)
            elif key in self.strings:
                parsed = str(value)
            elif key in self.boolean_values:
                val = value.strip().lower()
                parsed = val
            else:
                parsed = float(value)
            self.parameters[key] = parsed

    def set_param(self, tokens):
        if len(tokens) % 2 != 0:            
            tokens.append('blank')                # pad with a blank string for consistency ('ance reg_flag gam e1 v1 e2 v2' line has 8 values but 7 param names)
        values = tokens[:len(tokens) // 2]        # first half are values
        keys = tokens[len(tokens) // 2:]          # second half are parameter names
        self.parse_group(keys, values)
        self.sections.append(ControlGroup(self.parameters, keys))


    def read_control(self, control_file):
        with open(control_file, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]
            self.Name = os.path.splitext(os.path.basename(control_file))[0]

        # Phase 1: Read initial control lines
        index = 0
        for _ in range(2):                                         
            tokens = lines[index].split()
            self.set_param(tokens)
            index += 1

        # Phase 2: Handle wavefunction(s)
        def read_wavefunction_block(start_line, label='wf'):
            wf_type = lines[start_line].split()[0]
            wf_type_key = f"{label}_type"
            self.parameters[wf_type_key] = wf_type
            self.sections.append(ControlGroup(self.parameters, [wf_type_key]))
            wf = create_wavefunction_input(wf_type)
            end_line = wf.read_wavefunction_input(lines, start_line + 1)

            self.parameters[label] = wf  # store wavefunction object directly
            self.sections.append(ControlGroup(self.parameters, [label]))  # keep section info for writing
            self.sections.append(wf)
            return end_line

        if self.parameters.get('bra_eq_ket') == '.false.':
            index = read_wavefunction_block(index, label='ket')
            index = read_wavefunction_block(index, label='bra')
        else:
            index = read_wavefunction_block(index, label='wf')

        # Phase 3: Read the rest of control parameters
        while index < len(lines):
            tokens = lines[index].split()
            self.set_param(tokens)
            index += 1

    # Write the control file in-place
    def write_control(self, stream=None):
        if stream is None:
            filename = f"{self.Name}.ctrl" 
            with open(filename, 'w') as f:
                for section in self.sections:
                    section.write_output(stream=f.write)
        else:
            for section in self.sections:
                section.write_output(stream=stream)

    # Update paths in the control file
    def update_paths(self, input_dir, working_dir):
        # Update global string paths
        for key in self.strings:
            val = self.parameters.get(key)

            if not isinstance(val, str):
                continue

            stripped = val.strip("'\"")

            # Special case: optimized_deck should go to working_dir
            if key == 'optimized_deck':
                self.parameters[key] = f"'{os.path.join(working_dir, stripped)}'"
                continue  # Skip rest of loop to avoid overwriting

            # Skip basis (as per your logic)
            if key == 'basis':
                continue

            self.parameters[key] = f"'{os.path.join(input_dir, stripped)}'"

        # Update wavefunction (wf, bra, ket) paths
        for label in ['wf', 'bra', 'ket']:
            wf = self.parameters.get(label)
            if isinstance(wf, WavefunctionInput):
                wf.resolve_paths(input_dir)



# helper function to create wavefunction input based on type
def create_wavefunction_input(wf_type: str) -> WavefunctionInput:
    if wf_type == "'variational'":
        return VariationalWavefunctionInput()
    elif wf_type == "'product'":
        return ProductWavefunctionInput()
    elif wf_type == "'deuteron'":
        return DeuteronWavefunctionInput()
    else:
        raise ValueError(f"Unknown wavefunction type: {wf_type}")