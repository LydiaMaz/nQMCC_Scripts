
import os
import sys

import deck as Deck
from deck import read_params_and_deck

# Function to zero out all float parameters in a deck and its spatial symmetry objects
# If correlation_group is provided, only those params will be set to default value if they are non-zero
def zero_var_params(param_file, deck_file, spatial_sym_name, correlation_group, default=None):
    deck, deck_name = read_params_and_deck(param_file, deck_file)
    correlation_group = set(correlation_group)  # for faster lookup

    # Zero all regular float parameters (non-spatial symmetry)
    for key, value in deck.parameters.items():
        if not hasattr(value, '__dict__'):  # regular param
            if isinstance(value, list):
                for i in range(len(value)):
                    if isinstance(value[i], float):
                        value[i] = 0.0
            elif isinstance(value, float):
                deck.parameters[key] = 0.0

    if deck.spatial_symmetries is True:
        for key, sym_obj in deck.parameters.items():
            if hasattr(sym_obj, '__dict__'):
                # spatial symmetry object
                if key == spatial_sym_name:
                    # For input spatial symmetry object, set according to correlation group and non-zero logic
                    for attr_name, attr_value in vars(sym_obj).items():
                        if isinstance(attr_value, float):
                            if attr_name in correlation_group:
                                if attr_value != 0.0:
                                    setattr(sym_obj, attr_name, default * attr_value)
                                else:
                                    setattr(sym_obj, attr_name, 0.0)
                            else:
                                setattr(sym_obj, attr_name, 0.0)
                else:
                    # For other spatial symmetry objects, set all float attributes to zero
                    for attr_name, attr_value in vars(sym_obj).items():
                        if isinstance(attr_value, float):
                            setattr(sym_obj, attr_name, 0.0)

    # save_opt_file(deck, deck_file, base_dir="./opt")
    return deck, deck_name






def save_opt_file(deck, working_deck_path, base_dir="./opt"):
    os.makedirs(base_dir, exist_ok=True)

    deck_base = os.path.splitext(os.path.basename(working_deck_path))[0]
    opt_path = os.path.join(base_dir, deck_base)

    deck.write_deck(opt_path, extension="opt")  # specify .opt file

    return opt_path + ".opt"

# zero_var_params('nuclei/params/li7.params', 'nuclei/decks/li7.dk', '4P[21]', ['spu', 'spv', 'spr', 'spa', 'spb', 'spc', 'spk', 'spl'], default=0.25)
# zero_var_params('nuclei/params/he4n.params', 'nuclei/decks/he5_1hp_av18.dk', '1S[0]', ['spu', 'spv', 'spr', 'spa', 'spb', 'spc', 'spk', 'spl'], default=0.2)
