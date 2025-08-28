import json
import os
from typing import List, Optional, Dict

from pydantic import BaseModel


class Steps(BaseModel):
    """ Information about a specific step of GTDB-Tk"""
    name: str
    output_dir: Optional[str]
    starts_at: Optional[str]
    ends_at: Optional[str]
    duration: Optional[str]
    status: Optional[str]

    def is_complete(self) -> bool:
        if self.status == "completed":
            return True
        else:
            return False

class ANIScreenStep(Steps):
    """ Information about the ANI screen step of GTDB-Tk"""
    name: str = "ANI screen"
    genome_dir: Optional[str]
    batchfile: Optional[str]
    min_af: Optional[str]
    output_files: Optional[Dict]

class IdentifyStep(Steps):
    """ Information about the Identify step of GTDB-Tk"""
    name: str = "identify"
    genes: Optional[bool]
    extension: Optional[str]
    write_single_copy_genes: Optional[bool]
    genome_dir: Optional[str]
    batchfile: Optional[str]
    output_files: Optional[Dict]

class AlignStep(Steps):
    """ Information about the Align step of GTDB-Tk"""
    name: str = "align"
    identify_dir: Optional[str]
    skip_gtdb_refs: Optional[bool]
    taxa_filter: Optional[str]
    min_perc_aa: Optional[str]
    custom_msa_filters: Optional[bool]
    skip_trimming: Optional[bool]
    rnd_seed: Optional[str]
    cols_per_gene: Optional[str]
    min_consensus: Optional[str]
    max_consensus: Optional[str]
    min_perc_taxa: Optional[str]
    outgroup_taxon: Optional[str]
    output_files: Optional[Dict]

class ClassifyStep(Steps):
    """ Information about the Classify step of GTDB-Tk"""
    name: str = "classify"
    align_dir: Optional[str]
    genome_dir: Optional[str]
    batchfile: Optional[str]
    scratch_dir: Optional[str]
    debug_option: Optional[bool]
    full_tree: Optional[bool]
    skip_ani_screen: Optional[bool]
    output_files: Optional[Dict]

class InferStep(Steps):
    """ Information about the Infer step of GTDB-Tk"""
    name: str = "infer"
    msa_file: Optional[str]
    prot_model: Optional[str]
    no_support: Optional[bool]
    gamma: Optional[bool]

class RootStep(Steps):
    """ Information about the Root step of GTDB-Tk"""
    name: str = "root"
    input_tree: Optional[str]
    gtdbtk_classification_file: Optional[str]
    custom_taxonomy_file: Optional[str]
    outgroup_taxon: Optional[str]

class DecorateStep(Steps):
    """ Information about the Decorate step of GTDB-Tk"""
    name: str = "decorate"
    input_tree: Optional[str]
    gtdbtk_classification_file: Optional[str]
    custom_taxonomy_file: Optional[str]
    suffix: Optional[str]
    output_files: Optional[Dict]



class StageLogger(object):
    class __StageLogger(BaseModel):

        version: Optional[str]
        command_line: Optional[str]
        database_version: Optional[str]
        database_path: Optional[str]
        steps: List[Steps]
        output_dir: Optional[str]
        path: Optional[str]

        def __str__(self):
            return repr(self) + self.val

        def write(self):
            with open(self.path, "w") as f:
                f.write(self.json(indent=4))

        def has_stage(self, stage_object: object) -> bool:
            processed_steps = [x for x in self.steps if isinstance(x, stage_object)]
            if len(processed_steps) > 0:
                return True
            else:
                return False

        def get_stage(self, stage_object: object) -> object:
            processed_steps = [x for x in self.steps if isinstance(x, stage_object)]
            if len(processed_steps) > 0:
                return processed_steps[0]
            else:
                return None

        # create reset steps function, which takes an optional list of steps to keep
        # if no list is provided, all steps are removed

        def reset_steps(self, keep_steps: List[str] = None):
            if keep_steps is None:
                self.steps = []
            else:
                self.steps = [x for x in self.steps if x.name in keep_steps]

        def read_existing_steps(self):
            with open(self.path, "r") as f:
                data = json.load(f)
                #stage_logger = self.parse_obj(data)
                steps = []
                for step in data.get('steps'):
                    if step.get('name') == "ANI screen":
                        step_object = ANIScreenStep(**step)
                    elif step.get('name') == "identify":
                        step_object = IdentifyStep(**step)
                    elif step.get('name') == "align":
                        step_object = AlignStep(**step)
                    elif step.get('name') == "classify":
                        step_object = ClassifyStep(**step)
                    elif step.get('name') == "infer":
                        step_object = InferStep(**step)
                    else:
                        raise Exception(f"Unknown step name {step.name}")
                    steps.append(step_object)
                self.steps = steps

    instance = None
    def __new__(cls): # __new__ always a classmethod
        if not StageLogger.instance:
            StageLogger.instance = StageLogger.__StageLogger(steps=[])
        return StageLogger.instance
    def __getattr__(self, name):
        return getattr(self.instance, name)
    def __setattr__(self, name):
        return setattr(self.instance, name)





# class StageLoggerFile:
#     def __init__(self,output_dir: str):
#         self.output_dir = output_dir
#
#     # We generate the StageLogger object
#     def setStageLogger(self, version:str, command_line: str,
#                  database_version: str, database_path: str):
#         self.path = os.path.join(self.output_dir, "gtdbtk.json")
#         self.stage_logger = StageLogger(version = version,
#                                         command_line=command_line,
#                                         database_version=database_version,
#                                         database_path=database_path,
#                                         steps=[])
#     def write(self):
#         with open(self.path, "w") as f:
#             f.write(self.stage_logger.json(indent=4))
#
#     def read(self):
#         with open(self.path, "r") as f:
#             self.stage_logger = StageLogger.parse_raw(f.read())
#             steps = []
#             for step in self.stage_logger.steps:
#                 if step.name == "ANI screen":
#                     step_object = ANIScreenStep(**step.dict())
#                 elif step.name == "identify":
#                     step_object = IdentifyStep(**step.dict())
#                 elif step.name == "align":
#                     step_object = AlignStep(**step.dict())
#                 elif step.name == "classify":
#                     step_object = ClassifyStep(**step.dict())
#                 elif step.name == "infer":
#                     step_object = InferStep(**step.dict())
#                 else:
#                     raise Exception(f"Unknown step name {step.name}")
#                 steps.append(step_object)
#             self.stage_logger.steps = steps