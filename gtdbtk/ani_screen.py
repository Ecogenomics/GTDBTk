import logging
import os

from gtdbtk.ani_rep import ANIRep, ANISummaryFile
from gtdbtk.biolib_lite.common import make_sure_path_exists, canonical_gid

from gtdbtk.biolib_lite.taxonomy import Taxonomy
from gtdbtk.classify import Classify
from gtdbtk.config.output import DIR_ANISCREEN

from gtdbtk.config.config import (TAXONOMY_FILE, MASH_SKETCH_FILE, AF_THRESHOLD)
from gtdbtk.files.gtdb_radii import GTDBRadiiFile


class ANIScreener(object):
    """Computes a list of genomes to a list of representatives."""

    def __init__(self, cpus):
        """Instantiate the ANI rep class.

        Parameters
        ----------
        cpus : int
            The maximum number of CPUs available to this workflow.
        """
        self.logger = logging.getLogger('timestamp')
        self.cpus = cpus
        self.gtdb_radii = GTDBRadiiFile()

    def run_aniscreen(self,genomes, no_mash,out_dir,prefix, mash_k, mash_v, mash_s, mash_max_dist, mash_db):

        # If prescreen is set to True, then we will first run all genomes against a mash database
        # of all genomes in the reference package. The next step will be to classify those genomes with
        # FastANI.
        # All genomes classified with FastANI will be removed from the input genomes list for the
        # rest of the pipeline.
        mash_classified_user_genomes = {}
        #if mash_db finishes with a backslash, it should be considered a directory
        if mash_db.endswith('/'):
            make_sure_path_exists(mash_db)
        if os.path.isdir(mash_db):
            mash_db = os.path.join(mash_db, MASH_SKETCH_FILE)

        #we set mash_d == mash_max_dist to avoid user to run mash with impossible values
        mash_d = mash_max_dist

        ani_rep = ANIRep(self.cpus)
        # we store all the mash information in the classify directory
        fastani_results = ani_rep.run_mash_fastani(genomes, no_mash, mash_d, os.path.join(out_dir, DIR_ANISCREEN),
                                                    prefix, mash_k, mash_v, mash_s, mash_max_dist, mash_db)

        taxonomy = Taxonomy().read(TAXONOMY_FILE, canonical_ids=True)

        mash_classified_user_genomes = self.sort_fastani_ani_screen(
             fastani_results,taxonomy)

        #We write the results in 2 different files for each domain
        reports = {}
        if mash_classified_user_genomes:
            for domain,results in mash_classified_user_genomes.items():
                ani_summary_file = ANISummaryFile(os.path.join(out_dir,DIR_ANISCREEN),prefix,results,taxonomy,domain)
                ani_summary_file.write(ani_screen_step=True)
                reports[domain] = os.path.join(out_dir,DIR_ANISCREEN,prefix + '.' + domain + '.ani_summary.tsv')
        len_mash_classified_bac120 = len(mash_classified_user_genomes['bac120']) \
            if 'bac120' in mash_classified_user_genomes else 0

        len_mash_classified_ar53 = len(mash_classified_user_genomes['ar53']) \
            if 'ar53' in mash_classified_user_genomes else 0

        self.logger.info(f'{len_mash_classified_ar53 + len_mash_classified_bac120} genome(s) have '
                         f'been classified using the ANI pre-screening step.')

        return mash_classified_user_genomes,reports

    def sort_fastani_ani_screen(self,fastani_results,taxonomy,bac_ar_diff=None):
        """ When run mash/FastANI on all genomes before using pplacer, we need to sort those results and store them for
        a later use

        Parameters
        ----------
        fastani_results : dict
            The results of the FastANI run
        taxonomy : dict
            The taxonomy of the reference genomes
        """
        classified_user_genomes = {}

        # sort the dictionary by ani then af
        for gid in fastani_results.keys():
            thresh_results = [(ref_gid, hit) for (ref_gid, hit) in fastani_results[gid].items() if
                              hit['af'] >= AF_THRESHOLD and hit['ani'] >= self.gtdb_radii.get_rep_ani(
                                  canonical_gid(ref_gid))]
            all_results = [(ref_gid, hit) for (ref_gid, hit) in fastani_results[gid].items()]
            closest = sorted(thresh_results, key=lambda x: (-x[1]['ani'], -x[1]['af']))
            all_closest = sorted(all_results, key=lambda x: (-x[1]['ani'], -x[1]['af']))
            if len(closest) > 0:
                ref_gid, hit = closest[0]
                hit_taxonomy = taxonomy[canonical_gid(ref_gid)]
                if len(all_closest) > 1:
                    other_ref = '; '.join(Classify.formatnote(
                        closest,taxonomy,Classify.parse_radius_file(), [ref_gid]))
                    if len(other_ref) > 0:
                        hit['other_related_refs'] = other_ref

                if hit_taxonomy[0] == 'd__Bacteria':
                    classified_user_genomes.setdefault('bac120', {})[gid]={ref_gid:hit}

                elif hit_taxonomy[0] == 'd__Archaea':
                    classified_user_genomes.setdefault('ar53', {})[gid]={ref_gid:hit}


        return classified_user_genomes