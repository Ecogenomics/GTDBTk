import logging

from gtdbtk.exceptions import GTDBTkExit


class Batchfile(object):

    def __init__(self, path: str):
        self.path = path
        self.genome_path, self.genome_tln = self.read(path)

    @staticmethod
    def read(path):
        logger = logging.getLogger('timestamp')
        genomes, tln_tables = dict(), dict()
        seen_paths = set()
        warnings = list()
        with open(path) as fh:
            for line_no, line in enumerate(fh):
                line_split = line.strip().split("\t")
                if line_split[0] == '':
                    continue  # blank line

                if len(line_split) not in {2, 3}:
                    raise GTDBTkExit('Batch file must contain either 2 '
                                     'columns (detect translation table), '
                                     'or 3 (specify translation table).')

                if len(line_split) == 2:
                    genome_file, genome_id = line_split
                elif len(line_split) == 3:
                    genome_file, genome_id, tln_table = line_split
                    if tln_table not in {'4', '11'}:
                        raise GTDBTkExit('Specified translation table must '
                                         'be either 4, or 11.')
                    tln_tables[genome_id] = int(tln_table)

                if genome_file is None or genome_file == '':
                    warnings.append(f'Missing genome path on line {line_no + 1}.')
                elif genome_id is None or genome_id == '':
                    warnings.append(f'Missing genome ID on line {line_no + 1}.')
                elif genome_id in genomes:
                    warnings.append(f'Genome ID {genome_id} appears multiple times.')
                if genome_file in seen_paths:
                    logger.warning(f'Genome file appears multiple times: {genome_file}')

                # All good, record the value.
                genomes[genome_id] = genome_file
                seen_paths.add(genome_file)

        # Check if any warnings were raised.
        if len(warnings) > 0:
            warning_str = '\n'.join(warnings)
            raise GTDBTkExit(f'Please check the format of your batchfile, '
                             f'the following errors were found: {warning_str}')

        return genomes, tln_tables
