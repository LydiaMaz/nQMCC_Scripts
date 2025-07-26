import os

from deck import Deck, read_params_and_deck

# Function to zero out all float parameters in a deck and its spatial symmetry objects
# Set as % of current value or to a fixed value based on correlation groups

def zero_var_params(param_file, deck_file, spatial_sym_name, correlation_groups):
    deck, deck_name = read_params_and_deck(param_file, deck_file)

    # Create a flat dictionary: {param_name: {'mode': ..., 'value': ...}}
    correlation_instructions = {}
    for group in correlation_groups:
        for param in group['params']:
            correlation_instructions[param] = {
                'mode': group['mode'],
                'value': group['value']
            }

    correlation_set = set(correlation_instructions.keys())

    # Zero all regular float parameters (non-spatial symmetry)
    for key, value in deck.parameters.items():
        if not hasattr(value, '__dict__'):  # regular param
            if key in correlation_set:
                print(key,value)
                instr = correlation_instructions[key]
                if isinstance(value, list):
                    for i in range(len(value)):
                        if isinstance(value[i], float):
                            if instr['mode'] == 'scale':
                                deck.parameters[key][i] = instr['value'] * value
                            elif instr['mode'] == 'set':
                                deck.parameters[key][i] = instr['value']
                elif isinstance(value, float):
                    if instr['mode'] == 'scale':
                        deck.parameters[key] = instr['value'] * value
                    elif instr['mode'] == 'set':
                        deck.parameters[key] = instr['value']
            else:
                if isinstance(value, list):
                    for i in range(len(value)):
                        if isinstance(value[i], float):
                            value[i] = 0.0
                elif isinstance(value, float):
                    deck.parameters[key] = 0.0
        elif hasattr(value, '__dict__'):
            if key == spatial_sym_name:
                for attr_name, attr_value in vars(value).items():
                    if isinstance(attr_value, float):
                        if attr_name in correlation_set:
                            instr = correlation_instructions[attr_name]
                            if attr_value != 0.0:
                                if instr['mode'] == 'scale':
                                    setattr(value, attr_name, instr['value'] * attr_value)
                                elif instr['mode'] == 'set':
                                    setattr(value, attr_name, instr['value'])
                            else:
                                setattr(value, attr_name, 0.0)
                        else:
                            setattr(value, attr_name, 0.0)
            else:
                # Zero all float attributes of other spatial symmetry objects
                for attr_name, attr_value in vars(value).items():
                    if isinstance(attr_value, float):
                        setattr(value, attr_name, 0.0)
    return deck, deck_name

"""
"""
def save_opt_file(bscat, deck, base_dir="./opt", suffix=""):
    os.makedirs(base_dir, exist_ok=True)
    deck_base = f"{bscat:.4f}"+suffix
    opt_path = os.path.join(base_dir, deck_base)
    deck.write_deck(opt_path, extension="opt")  # specify .opt file
    return opt_path + ".opt"


