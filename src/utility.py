import os
import shlex
import shutil
from control import Control

class Utility:
    def __init__(self, filename, run_cmd=None):
        self.util_file = filename
        self.read_util(filename)
        if run_cmd is not None:
            self.run_cmd = shlex.split(run_cmd)
            self.write_util()

    def read_util(self, filename):
        if not filename.endswith(".util"):
            raise ValueError("Utility file must end with '.util'")

        with open(filename, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]

        for line in lines:
            tokens = line.split()
            key = tokens[-1].lower()
            value = tokens[:-1]
            if key == "run_cmd":
                value_str = " ".join(tokens[:-1])
                setattr(self, key, shlex.split(value_str))  # properly parsed into list
            else:
                setattr(self, key, " ".join(value))  # store as string

    def write_util(self, filename=None):
        if filename is None:
            if hasattr(self, 'util_file'):
                filename = self.util_file
            else:
                raise ValueError("No filename provided and no previously read file to overwrite.")

        with open(filename, 'w') as f:
            for key, value in self.__dict__.items():
                if key in ['target_control', 'scattering_control', 'util_file']:
                    continue
                if isinstance(value, list):
                    value_str = " ".join(value)
                else:
                    value_str = str(value)
                f.write(f"{value_str:<40} {key}\n")

    def copy_control_files(self):
        print(f"Target copy destination directory: {self.working_dir}")

        target_copy_path = os.path.join(self.working_dir, os.path.basename(self.target_file))
        scattering_copy_path = os.path.join(self.working_dir, os.path.basename(self.scattering_file))

        shutil.copy(self.target_file, target_copy_path)
        shutil.copy(self.scattering_file, scattering_copy_path)

        # Load control files from the copied locations
        target_control = Control()
        target_control.read_control(target_copy_path)

        scattering_control = Control()
        scattering_control.read_control(scattering_copy_path)

        # Update and rewrite paths
        target_control.update_paths(self.input_dir, self.working_dir)
        target_control.write_control()

        scattering_control.update_paths(self.input_dir, self.working_dir)
        scattering_control.write_control()

        setattr(self, 'target_control', target_control)
        setattr(self, 'scattering_control', scattering_control)

