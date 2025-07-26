from spatial_symmetry import read_spatial_symmetries
from parameters import Parameters
import os


class Deck:

    int_values = {
        'PI', 'Pi', 'lwf', 'lsc', 'l3bc', 'lopc', 'lcut',
        'lfsp0', 'e_or_v', 'lscat', 'lmu1', 'Lmu2', 'LL', 'lnodes',
        'lfsp', 'e_v', 'lscat', 'm1', 'm2', 'L', 'nodes'
    }

    arrays = {
        'esep': 4,
        'eta': 2, 'zeta': 2,
        'fscal': 8, 'ac': 8, 'aa': 8, 'ar': 8, 'alpha': 8, 'beta': 8, 'gamma': 8
    }
    
    def __init__(self):
        self.parameters = {}
        self.parameter_groups = [] # list of lists for related parameters
        self.spatial_symmetries = False # flag to indicate if spatial symmetries are present
        
        # Initialize arrays
        for attr, size in self.arrays.items():
            self.parameters[attr] = [0.0 for _ in range(size)]
    
    # helper functions
    def parse_group(self, group, tokens):
        for param, value in zip(group, tokens):
            parsed = int(value) if param in self.int_values else float(value)
            self.parameters[param] = parsed
        self.parameter_groups.append(group)
  
    def set_array_from_line(self, tokens, array_name):
        self.parameters[array_name] = [float(x) for x in tokens[:self.arrays[array_name]]]
        self.parameter_groups.append([array_name])

    def set_param(self, tokens):
        if '-' in tokens[-1]: # some dk files have range formatting (qsss1-2)
            values = list(map(float, tokens[:-1]))
            label = tokens[-1]
            prefix = ''.join(filter(str.isalpha, label))
            param_names = [f"{prefix}{i}" for i in range(1,3)]  # e.g. qsss1, qsss2
        else:
            values = tokens[:len(tokens) // 2]        # first half are values
            param_names = tokens[len(tokens) // 2:]   # second half are parameter names

        self.parse_group(param_names, values)
      


    def read_deck(self, filepath, parameters):

        with open(filepath, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]
            self.parameters['Name'] = lines[0].rsplit(' ', 1)[0].strip()            # Extract name from first line
            self.parameter_groups.append(['Name'])

        for i in range(1,19):                                         # Logic always the same for the first 18 lines
            tokens = lines[i].strip().split()
            if tokens[-1] in self.arrays:
                self.set_array_from_line(tokens, tokens[-1])
            else:
                self.set_param(tokens)
        
        # Logic changes at this line based on .params file
        current_line = 19

        NPPART = int(parameters.NPPART)
        NDPART = int(parameters.NSDPART)

        nppart_lines = min(NPPART, 3)                   # max 3 correlation lines (spu, ppu, etc.)
        ndpart_lines = sum([3, 2, 1][:min(NDPART, 3)])  # up to 6 lines total for sd/pp/dd correlation
        extra_lines = nppart_lines + ndpart_lines

        for i in range(current_line, current_line + extra_lines):
            tokens = lines[i].strip().split()
            self.set_param(tokens)
        
        # Read the spatial symmetries
        if (NPPART != 0):
            self.spatial_symmetries = True
            symmetries = read_spatial_symmetries(
                lines, current_line + extra_lines,
                nbeta=int(parameters.NBETA),
                nppart=NPPART,
                ndpart=NDPART,
                phi_type=int(parameters.PHI_TYPE)
            )
            self.parameters.update(symmetries)
    

    
    def write_deck(self, filename, extension=None, stream=None):
        close_stream = False

        # Determine output stream
        if stream is None:
            if extension is None:
                filename = f"{filename}.dk"
            else:
                filename = f"{filename}.{extension}"
            stream = open(filename, 'w')
            close_stream = True

        def write_line(line):
            if callable(stream):  # e.g., print
                stream(line)
            else:
                stream.write(line)

        # Write all parameter groups
        for group in self.parameter_groups:
            values = []
            labels = []

            for param in group:
                val = self.parameters.get(param)
                if isinstance(val, list):
                    values.extend(f"{x:.5f}" for x in val)
                    labels.append(param)
                else:
                    values.append(f"{val:<8}")
                    labels.append(param)

            values_str = " ".join(values)
            label_str = " ".join(labels)
            spacing = max(80, len(values_str) + 2)
            line = f"{values_str:<{spacing}}{label_str}"
            write_line(line.rstrip() + "\n")

        # Write spatial symmetries
        if self.spatial_symmetries is True:
            for sym_obj in self.parameters.values():
                if hasattr(sym_obj, 'correlation_groups'):
                    for group in sym_obj.correlation_groups:
                        values = [getattr(sym_obj, name) for name in group]
                        line = format_group_line(values)
                        comment = "  " + " ".join(group)
                        write_line(f"{line:<80}{comment}\n")

        if close_stream:
            stream.close()


            
# Helper function to format group lines
def format_group_line(values):
    formatted = []
    for val in values:
        if isinstance(val, float):
            formatted.append(f"{val:>8.5f}")
        elif isinstance(val, int):
            formatted.append(f"{val:>4d}")
        elif isinstance(val, str):
            formatted.append(f"{val:<10}")
        else:
            formatted.append(str(val))
    return " ".join(formatted)

# Helper function to read the deck and parameters
def read_params_and_deck(param_file, deck_file):
    params = Parameters()
    params.read_parameters(param_file)

    deck = Deck()
    deck.read_deck(deck_file, params)

    deck_name = os.path.splitext(os.path.basename(deck_file))[0]

    return deck, deck_name