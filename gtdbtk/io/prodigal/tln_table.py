from typing import Optional

from gtdbtk.exceptions import GTDBTkExit


class TlnTableFile(object):

    def __init__(self, path: str,
                 best_tln_table: Optional[int] = None,
                 coding_density_4: Optional[float] = None,
                 coding_density_11: Optional[float] = None):
        self.path = path
        self.best_tln_table = best_tln_table
        self.coding_density_4 = coding_density_4
        self.coding_density_11 = coding_density_11

    def read(self):
        with open(self.path, 'r') as fh:
            for line in fh.readlines():
                idx, val = line.strip().split('\t')
                if idx == 'best_translation_table':
                    try:
                        self.best_tln_table = int(val)
                    except ValueError:
                        raise GTDBTkExit(f'Invalid translation table: {val} for {self.path}')
                elif idx == 'coding_density_4':
                    try:
                        self.coding_density_4 = float(val)
                    except ValueError:
                        raise GTDBTkExit(f'Invalid coding density: {val} for {self.path}')
                elif idx == 'coding_density_11':
                    try:
                        self.coding_density_11 = float(val)
                    except ValueError:
                        raise GTDBTkExit(f'Invalid coding density: {val} for {self.path}')

    def write(self):
        with open(self.path, 'w') as fh:
            fh.write(f'best_translation_table\t{self.best_tln_table}\n')
            fh.write(f'coding_density_4\t{self.coding_density_4}\n')
            fh.write(f'coding_density_11\t{self.coding_density_11}\n')
