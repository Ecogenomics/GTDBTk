import json
from datetime import datetime

from dataclasses import dataclass, field, asdict
from typing import List, Optional, Dict, ClassVar



def _convert(obj):
    if isinstance(obj, datetime):
        return obj.isoformat()
    if isinstance(obj, list):
        return [_convert(o) for o in obj]
    if isinstance(obj, dict):
        return {k: _convert(v) for k, v in obj.items()}
    return obj

@dataclass
class Steps:
    """Information about a specific step of GTDB-Tk"""
    name: str
    output_dir: Optional[str] = None
    starts_at: Optional[str] = None
    ends_at: Optional[str] = None
    duration: Optional[str] = None
    status: Optional[str] = None

    def is_complete(self) -> bool:
        return self.status == "completed"


@dataclass
class ANIScreenStep(Steps):
    """Information about the ANI screen step of GTDB-Tk"""
    name: str = field(default="ANI screen", init=False)
    genome_dir: Optional[str] = None
    batchfile: Optional[str] = None
    min_af: Optional[str] = None
    output_files: Optional[Dict] = None


@dataclass
class IdentifyStep(Steps):
    """Information about the Identify step of GTDB-Tk"""
    name: str = field(default="identify", init=False)
    genes: Optional[bool] = None
    extension: Optional[str] = None
    write_single_copy_genes: Optional[bool] = None
    genome_dir: Optional[str] = None
    batchfile: Optional[str] = None
    output_files: Optional[Dict] = None


@dataclass
class AlignStep(Steps):
    """Information about the Align step of GTDB-Tk"""
    name: str = field(default="align", init=False)
    identify_dir: Optional[str] = None
    skip_gtdb_refs: Optional[bool] = None
    taxa_filter: Optional[str] = None
    min_perc_aa: Optional[str] = None
    custom_msa_filters: Optional[bool] = None
    skip_trimming: Optional[bool] = None
    rnd_seed: Optional[str] = None
    cols_per_gene: Optional[str] = None
    min_consensus: Optional[str] = None
    max_consensus: Optional[str] = None
    min_perc_taxa: Optional[str] = None
    outgroup_taxon: Optional[str] = None
    output_files: Optional[Dict] = None


@dataclass
class ClassifyStep(Steps):
    """Information about the Classify step of GTDB-Tk"""
    name: str = field(default="classify", init=False)
    align_dir: Optional[str] = None
    genome_dir: Optional[str] = None
    batchfile: Optional[str] = None
    scratch_dir: Optional[str] = None
    debug_option: Optional[bool] = None
    full_tree: Optional[bool] = None
    skip_ani_screen: Optional[bool] = None
    output_files: Optional[Dict] = None


@dataclass
class InferStep(Steps):
    """Information about the Infer step of GTDB-Tk"""
    name: str = field(default="infer", init=False)
    msa_file: Optional[str] = None
    prot_model: Optional[str] = None
    no_support: Optional[bool] = None
    gamma: Optional[bool] = None


@dataclass
class RootStep(Steps):
    """Information about the Root step of GTDB-Tk"""
    name: str = field(default="root", init=False)
    input_tree: Optional[str] = None
    gtdbtk_classification_file: Optional[str] = None
    custom_taxonomy_file: Optional[str] = None
    outgroup_taxon: Optional[str] = None


@dataclass
class DecorateStep(Steps):
    """Information about the Decorate step of GTDB-Tk"""
    name: str = field(default="decorate", init=False)
    input_tree: Optional[str] = None
    gtdbtk_classification_file: Optional[str] = None
    custom_taxonomy_file: Optional[str] = None
    suffix: Optional[str] = None
    output_files: Optional[Dict] = None





@dataclass
class _StageLoggerImpl:
    """Internal implementation of StageLogger"""
    version: Optional[str] = None
    command_line: Optional[str] = None
    database_version: Optional[str] = None
    database_path: Optional[str] = None
    steps: List[Steps] = field(default_factory=list)
    output_dir: Optional[str] = None
    path: Optional[str] = None

    def __str__(self):
        return repr(self)

    def write(self):
        if not self.path:
            # Silently skip writing if no path is set
            # This is expected for commands like check_install
            return
        with open(self.path, "w") as f:
            json.dump(_convert(asdict(self)), f, indent=4)

    def has_stage(self, stage_class: type) -> bool:
        return any(isinstance(x, stage_class) for x in self.steps)

    def get_stage(self, stage_class: type) -> Optional[Steps]:
        for x in self.steps:
            if isinstance(x, stage_class):
                return x
        return None

    def reset_steps(self, keep_steps: Optional[List[str]] = None):
        if keep_steps is None:
            self.steps = []
        else:
            self.steps = [x for x in self.steps if x.name in keep_steps]

    def read_existing_steps(self):
        if not self.path:
            raise ValueError("Path not set for StageLogger")

        with open(self.path, "r") as f:
            data = json.load(f)
            steps = []
            for step in data.get('steps', []):
                name = step.get('name')
                if name == "ANI screen":
                    step_object = ANIScreenStep(**step)
                elif name == "identify":
                    step_object = IdentifyStep(**step)
                elif name == "align":
                    step_object = AlignStep(**step)
                elif name == "classify":
                    step_object = ClassifyStep(**step)
                elif name == "infer":
                    step_object = InferStep(**step)
                elif name == "root":
                    step_object = RootStep(**step)
                elif name == "decorate":
                    step_object = DecorateStep(**step)
                else:
                    raise Exception(f"Unknown step name {name}")
                steps.append(step_object)
            self.steps = steps


class StageLogger:
    """Singleton wrapper for StageLogger"""
    _instance = None

    def __new__(cls):
        if not cls._instance:
            cls._instance = super().__new__(cls)
            cls._instance._impl = _StageLoggerImpl(steps=[])
        return cls._instance

    @classmethod
    def instance(cls):
        return cls._instance

    def __getattr__(self, name):
        return getattr(self._impl, name)

    def __setattr__(self, name, value):
        if name == '_impl':
            super().__setattr__(name, value)
        else:
            setattr(self._impl, name, value)