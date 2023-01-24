import os
from typing import List, Optional, Dict

from pydantic import BaseModel


class Steps(BaseModel):
    """ Information about a specific step of GTDB-Tk"""
    name: str
    input: Optional[str]
    output_dir: Optional[str]
    output_files: Optional[List[str]]
    starts_at: Optional[str]
    ends_at: Optional[str]
    duration: Optional[str]
    status: Optional[str]

class ANI_Screen_Step(Steps):
    """ Information about the ANI screen step of GTDB-Tk"""
    name: str = "ANI screen"
    output_files: Optional[Dict[str, str]]

class Identify_Step(Steps):
    """ Information about the Identify step of GTDB-Tk"""
    name: str = "identify"

class Align_Step(Steps):
    """ Information about the Align step of GTDB-Tk"""
    name: str = "align"

class Classify_Step(Steps):
    """ Information about the Classify step of GTDB-Tk"""
    name: str = "classify"

class StageLogger(BaseModel):
    version : str
    command_line: str
    database_version: str
    database_path: str
    steps: List[Steps]

class StageLoggerFile:
    def __init__(self,version:str, command_line: str,
                 database_version: str, database_path: str,
                 output_dir: str):
        self.output_dir = output_dir
        self.stage_logger = StageLogger(version = version,
                                        command_line=command_line,
                                        database_version=database_version,
                                        database_path=database_path,
                                        steps=[])
    def write(self):
        with open(os.path.join(self.output_dir, "gtdbtk.json"), "w") as f:
            f.write(self.stage_logger.json(indent=4))

    def read(self):
        with open(os.path.join(self.output_dir, "gtdbtk.json"), "r") as f:
            self.stage_logger = StageLogger.parse_raw(f.read())