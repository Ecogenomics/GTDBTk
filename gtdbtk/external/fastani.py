###############################################################################
#                                                                             #
#    This program is free software: you can redistribute it and/or modify     #
#    it under the terms of the GNU General Public License as published by     #
#    the Free Software Foundation, either version 3 of the License, or        #
#    (at your option) any later version.                                      #
#                                                                             #
#    This program is distributed in the hope that it will be useful,          #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of           #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the            #
#    GNU General Public License for more details.                             #
#                                                                             #
#    You should have received a copy of the GNU General Public License        #
#    along with this program. If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################

import multiprocessing
import os
import shutil
import tempfile

import gtdbtk.config.config as Config
from gtdbtk.biolib_lite.common import make_sure_path_exists, remove_extension
from gtdbtk.exceptions import FastANIException
from gtdbtk.tools import splitchunks


class FastANI(object):
    """Python wrapper for FastANI (https://github.com/ParBLiSS/FastANI)"""

    def __init__(self, cpus):
        """Instantiate the class.

        Parameters
        ----------
        cpus : int
            The number of CPUs to use.

        """
        self.cpus = cpus

    def run(self, fastani_verification, genomes):
        """Using the instance defined number of CPUs, run FastANI."""

        manager = multiprocessing.Manager()
        out_q = manager.dict()
        procs = []
        nprocs = self.cpus

        for item in splitchunks(fastani_verification, nprocs):
            p = multiprocessing.Process(
                target=self._fastaniWorker,
                args=(item, genomes, out_q))
            procs.append(p)
            p.start()

        # Wait for all worker processes to finish
        for p in procs:
            p.join()
            if p.exitcode == 1:
                raise FastANIException("A process returned a non-zero exit code.")

        all_fastani_dict = dict(out_q)

        return all_fastani_dict

    def _fastaniWorker(self, sublist_genomes, genomes, out_q):
        """Multi thread worker to calculate FastANI"""
        try:
            for userleaf, potential_nodes in sublist_genomes.items():
                dict_parser_distance = self._calculate_fastani_distance(
                    userleaf, potential_nodes, genomes)
                for k, v in dict_parser_distance.items():
                    if k in out_q:
                        raise Exception("{} not in output.".format(k))
                    out_q[k] = v
            return True
        except Exception as error:
            print(error)
            FastANIException(error)

    def _calculate_fastani_distance(self, user_leaf, list_leaf, genomes):
        """ Calculate the FastANI distance between all user genomes and the reference to classfy them at the species level

        Parameters
        ----------
        user_leaf : User genome node
        list_leaf : Dictionary of nodes including one or many user genomes and one reference genome.
        genomes : Dictionary of user genomes d[genome_id] -> FASTA file

        Returns
        -------
        dictionary
            dict_results[user_g]={ref_genome1:{"af":af,"ani":ani},ref_genome2:{"af":af,"ani":ani}}
        """
        try:
            self.tmp_output_dir = tempfile.mkdtemp()
            make_sure_path_exists(self.tmp_output_dir)
            dict_parser_distance = {}

            # we first calculate the user genome vs the reference
            # we write the two input files for fastani, the query file and
            # reference file
            path_user_list = os.path.join(self.tmp_output_dir, 'query_list.txt')
            with open(path_user_list, 'w') as f:
                f.write('{0}\n'.format(genomes.get(user_leaf.taxon.label)))

            leafnodes = list_leaf.get("potential_g")
            for node in leafnodes:
                leafnode = node[0]
                shortleaf = leafnode.taxon.label
                path_ref_list = os.path.join(self.tmp_output_dir, 'ref_{}.txt'.format(shortleaf))
                if leafnode.taxon.label.startswith('GB_') or leafnode.taxon.label.startswith('RS_'):
                    shortleaf = leafnode.taxon.label[3:]
                with open(path_ref_list, 'w') as f:
                    f.write('{}\n'.format(os.path.join(
                        Config.FASTANI_GENOMES, shortleaf + Config.FASTANI_GENOMES_EXT)))
                # run fastANI
                if not os.path.isfile(path_user_list) or not os.path.isfile(path_ref_list):
                    raise FastANIException

                path_results = os.path.join(self.tmp_output_dir, 'results_{}_UvsRef.tab'.format(shortleaf))
                path_error = os.path.join(self.tmp_output_dir, 'error_{}.log'.format(shortleaf))

                cmd = 'fastANI --ql {0} --rl {1} -o {2} > /dev/null 2>{3}'.format(path_user_list,
                                                                                  path_ref_list,
                                                                                  path_results,
                                                                                  path_error)

                os.system(cmd)

                if not os.path.isfile(path_results):
                    errstr = 'FastANI has stopped:\n'
                    if os.path.isfile(path_error):
                        with open(path_error) as debug:
                            for line in debug:
                                finalline = line
                            errstr += finalline
                    raise ValueError(errstr)

                dict_parser_distance = self._parse_fastani_results(path_results, dict_parser_distance)

                # We then calculate the reference vs user genome
                path_results_reverse = os.path.join(self.tmp_output_dir, 'results_{}_RefvsU.tab'.format(shortleaf))
                cmd_reverse = 'fastANI --ql {0} --rl {1} -o {2} > /dev/null 2>{3}'.format(path_ref_list,
                                                                                          path_user_list,
                                                                                          path_results_reverse,
                                                                                          path_error)
                os.system(cmd_reverse)
                if not os.path.isfile(path_results_reverse):
                    errstr = 'FastANI has stopped:\n'
                    if os.path.isfile(path_error):
                        with open(path_error) as debug:
                            for line in debug:
                                finalline = line
                            errstr += finalline
                    raise ValueError(errstr)
                dict_parser_distance = self._parse_fastani_results_reverse(path_results_reverse, dict_parser_distance)

            shutil.rmtree(self.tmp_output_dir)
            return dict_parser_distance

        except ValueError as error:
            if os.path.exists(self.tmp_output_dir):
                shutil.rmtree(self.tmp_output_dir)
            raise error
        except Exception as error:
            if os.path.exists(self.tmp_output_dir):
                shutil.rmtree(self.tmp_output_dir)
            raise error

    def _parse_fastani_results(self, fastout_file, dict_results):
        """ Parse the fastani output file

        Parameters
        ----------
        fastout_file : str
            fastani output file.


        Returns
        -------
        dictionary
            dict_results[user_g]={ref_genome1:{"af":af,"ani":ani},ref_genome2:{"af":af,"ani":ani}}
        """
        with open(fastout_file, 'r') as fastfile:
            for line in fastfile:
                info = line.strip().split()
                ref_genome = os.path.basename(info[1]).replace(
                    Config.FASTANI_GENOMES_EXT, "")
                user_g = remove_extension(os.path.basename(info[0]))
                ani = float(info[2])
                af = round(float(info[3]) / float(info[4]), 2)
                if user_g in dict_results:
                    dict_results[user_g][ref_genome] = {"ani": ani, 'af': af}
                else:
                    dict_results[user_g] = {ref_genome: {"ani": ani, "af": af}}
        return dict_results

    def _parse_fastani_results_reverse(self, fastout_file, dict_parser_distance):
        # TODO: Merge _parse_fastani_results and _parse_fastani_results_reverse
        """ Parse the fastani output file for the reverse comparison and pick the best ANI and AF


        Parameters
        ----------
        fastout_file : fastani output file.
        dict_parser_distance: dictionaryof user genomes vs list of refrence genomes with ANI and AF


        Returns
        -------
        dictionary
            dict_parser_distance[user_g]={ref_genome1:{"af":af,"ani":ani},ref_genome2:{"af":af,"ani":ani}}
        """
        with open(fastout_file) as fastfile:
            for line in fastfile:
                info = line.strip().split()
                ref_genome = os.path.basename(info[0]).replace(
                    Config.FASTANI_GENOMES_EXT, "")
                user_g = remove_extension(os.path.basename(info[1]))
                ani = float(info[2])
                af = round(float(info[3]) / float(info[4]), 2)
                if user_g in dict_parser_distance:
                    if ref_genome in dict_parser_distance.get(user_g):
                        if dict_parser_distance.get(user_g).get(ref_genome).get('ani') < ani:
                            dict_parser_distance[user_g][ref_genome]["ani"] = ani
                        if dict_parser_distance.get(user_g).get(ref_genome).get('af') < af:
                            dict_parser_distance[user_g][ref_genome]["af"] = af
                    else:
                        dict_parser_distance[user_g][ref_genome] = {"ani": ani, 'af': af}
                else:
                    dict_parser_distance[user_g] = {ref_genome: {"ani": ani, "af": af}}
        return dict_parser_distance
