import logging
import os

from gtdbtk.ani_rep import ANIRep, ANISummaryFile
from gtdbtk.biolib_lite.common import make_sure_path_exists, canonical_gid
from gtdbtk.biolib_lite.seq_io import read_fasta

from gtdbtk.biolib_lite.taxonomy import Taxonomy
from gtdbtk.classify import Classify
from gtdbtk.config.output import DIR_ANISCREEN

from gtdbtk.files.gtdb_radii import GTDBRadiiFile
from gtdbtk.config.common import CONFIG

class ANIScreener(object):
    """Computes a list of genomes to a list of representatives."""

    def __init__(self, cpus,af_threshold=None):
        """Instantiate the ANI rep class.

        Parameters
        ----------
        cpus : int
            The maximum number of CPUs available to this workflow.
        """
        self.logger = logging.getLogger('timestamp')
        self.cpus = cpus
        self.af_threshold = af_threshold if af_threshold else CONFIG.AF_THRESHOLD
        self.gtdb_radii = GTDBRadiiFile()

    def run_aniscreen(self,genomes,out_dir,prefix):

        # If prescreen is set to True, then we will first run all genomes against a skani database
        # of all genomes in the reference package.
        # All genomes classified with skani will be removed from the input genomes list for the
        # rest of the pipeline.



        ani_rep = ANIRep(self.cpus)

        #we create a copy of the genomes dictionary to avoid modifying it
        genomes_copy = genomes.copy()
        # we remove all empty files from genomes.
        for k in list(genomes_copy.keys()):
            if os.path.getsize(genomes_copy[k]) == 0:
                self.logger.warning(f'Genome {k} file is invalid for skani. It will be removed from the sketch step.')
                del genomes_copy[k]




        skani_results = ani_rep.run_skani(genomes, prefix)

        taxonomy = Taxonomy().read(CONFIG.TAXONOMY_FILE, canonical_ids=True)

        skani_classified_user_genomes = self.sort_skani_ani_screen(
             skani_results,taxonomy)

        #We write the results in 2 different files for each domain
        reports = {}
        if skani_classified_user_genomes:
            for domain,results in skani_classified_user_genomes.items():
                # we create the directory if it does not exist;
                # It is not created when using skani as the comparison is done in the temp dir and the sketch
                # is in memory
                make_sure_path_exists(os.path.join(out_dir,DIR_ANISCREEN))
                ani_summary_file = ANISummaryFile(os.path.join(out_dir,DIR_ANISCREEN),prefix,results,taxonomy,domain)
                ani_summary_file.write(ani_screen_step=True)
                reports[domain] = os.path.join(out_dir,DIR_ANISCREEN,prefix + '.' + domain + '.ani_summary.tsv')
        len_skani_classified_bac120 = len(skani_classified_user_genomes['bac120']) \
            if 'bac120' in skani_classified_user_genomes else 0

        len_skani_classified_ar53 = len(skani_classified_user_genomes['ar53']) \
            if 'ar53' in skani_classified_user_genomes else 0

        self.logger.info(f'{len_skani_classified_ar53 + len_skani_classified_bac120} genome(s) have '
                         f'been classified using the ANI pre-screening step.')

        return skani_classified_user_genomes,reports

    def sort_skani_ani_screen(self,skani_results,taxonomy,bac_ar_diff=None):
        """ When run skani on all genomes before using pplacer, we need to sort those results and store them for
        a later use

        Parameters
        ----------
        skani_results : dict
            The results of the skani run
        taxonomy : dict
            The taxonomy of the reference genomes
        """
        classified_user_genomes = {}

        # sort the dictionary by ani then af
        for gid in skani_results.keys():
            thresh_results = [(ref_gid, hit) for (ref_gid, hit) in skani_results[gid].items() if
                              hit['af'] >= self.af_threshold and hit['ani'] >= self.gtdb_radii.get_rep_ani(
                                  canonical_gid(ref_gid))]
            all_results = [(ref_gid, hit) for (ref_gid, hit) in skani_results[gid].items()]
            closest = sorted(thresh_results, key=lambda x: (-x[1]['ani'], -x[1]['af']))
            all_closest = sorted(all_results, key=lambda x: (-x[1]['ani'], -x[1]['af']))
            if len(closest) > 0:
                ref_gid, hit = closest[0]
                hit_taxonomy = taxonomy[canonical_gid(ref_gid)]
                if len(all_results) > 1:
                    other_ref = '; '.join(Classify.formatnote(
                        all_results,taxonomy,Classify.parse_radius_file(), [ref_gid]))
                    if len(other_ref) > 0:
                        hit['other_related_refs'] = other_ref

                if hit_taxonomy[0] == 'd__Bacteria':
                    classified_user_genomes.setdefault('bac120', {})[gid]={ref_gid:hit}

                elif hit_taxonomy[0] == 'd__Archaea':
                    classified_user_genomes.setdefault('ar53', {})[gid]={ref_gid:hit}


        return classified_user_genomes