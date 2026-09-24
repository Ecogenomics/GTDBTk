import logging
import os

from gtdbtk.ani_rep import ANIRep, ANISummaryFile, SkaniRawHitsFile
from gtdbtk.biolib_lite.common import make_sure_path_exists, canonical_gid
from gtdbtk.biolib_lite.seq_io import read_fasta

from gtdbtk.biolib_lite.taxonomy import Taxonomy
from gtdbtk.classify import Classify
from gtdbtk.config.output import DIR_ANISCREEN, PATH_ANISCREEN_UNASSIGNED_HITS

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

        skani_results = ani_rep.run_skani(genomes, prefix,output_dir=out_dir)

        taxonomy = Taxonomy().read(CONFIG.TAXONOMY_FILE, canonical_ids=True)

        skani_classified_user_genomes = self.sort_skani_ani_screen(
             skani_results,taxonomy)

        # Keep the raw hits of genomes NOT classified here: classify needs them for
        # below-radius / low-AF reporting (Issue #717) since it does not rerun skani.
        # Written before ANISummaryFile.write(), which rounds AF in place.
        # Always written (even if empty) so a stale file from a previous run is never reused.
        classified_gids = {gid for d in skani_classified_user_genomes.values() for gid in d}
        unassigned_hits_file = os.path.join(out_dir, PATH_ANISCREEN_UNASSIGNED_HITS.format(prefix=prefix))
        make_sure_path_exists(os.path.dirname(unassigned_hits_file))
        unassigned_gids = set(skani_results) - classified_gids
        SkaniRawHitsFile.write(unassigned_hits_file, skani_results, unassigned_gids)
        self.logger.info(f'skani hits of {len(unassigned_gids)} genome(s) not assigned by the ANI screen '
                         f'saved to: {unassigned_hits_file}')

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

        return skani_classified_user_genomes,reports,unassigned_hits_file

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
        species_radius = Classify.parse_radius_file()   # read once, not once per genome

        for gid in skani_results.keys():
            # Same rule as Classify._sort_skani_results_pre_pplacer (standalone classify):
            # take the closest AF-passing hit, then check only that hit against its own radius.
            # A lower hit within its own radius does NOT rescue the genome.
            af_pass = sorted([(ref_gid, hit) for (ref_gid, hit) in skani_results[gid].items()
                              if hit['af'] >= self.af_threshold],
                             key=lambda x: (-x[1]['ani'], -x[1]['af']))
            if not af_pass:
                continue
            ref_gid, hit = af_pass[0]
            if hit['ani'] < self.gtdb_radii.get_rep_ani(canonical_gid(ref_gid)):
                continue

            hit_taxonomy = taxonomy[canonical_gid(ref_gid)]
            # Issue #717: other references are capped at the AF threshold, closest-first
            # (formatnote's top-50 cap therefore keeps the 50 closest).
            if len(af_pass) > 1:
                other_ref = '; '.join(Classify.formatnote(
                    af_pass, taxonomy, species_radius, [ref_gid]))
                if len(other_ref) > 0:
                    hit['other_related_refs'] = other_ref

            if hit_taxonomy[0] == 'd__Bacteria':
                classified_user_genomes.setdefault('bac120', {})[gid] = {ref_gid: hit}
            elif hit_taxonomy[0] == 'd__Archaea':
                classified_user_genomes.setdefault('ar53', {})[gid] = {ref_gid: hit}

        return classified_user_genomes