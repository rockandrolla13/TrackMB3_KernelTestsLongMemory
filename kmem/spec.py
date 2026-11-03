"""Typed experiment specs (plan step 3.1, AR-ABS-001, AR-EXT-003, CR-008, CR-016, CR-017).

An ExperimentSpec describes one run as data: its cells with explicit seed indices, the masters, the options
the driver was run with, the writer format, the final number of replications per cell, and the prerequisite
files. The specs of EXP2-EXP7 (experiments/EXPn/spec.json) were built from the saved main-run manifests by
tools/build_specs.py, so every seed index and every cell name is the one the saved rows were made with.

A cell's seed index is data, never a list position: inserting a cell cannot reseed the others.
Each cell also keeps `fields`, the exact dict its driver's run_chunk consumes (as the manifest recorded it).
"""
from __future__ import annotations

import hashlib
import json
import os
from typing import Any, Literal

from pydantic import BaseModel, ConfigDict, Field, model_validator

from kmem import seeds

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EXPERIMENTS = {
    "EXP2": "experiments/EXP2-20261023",
    "EXP3": "experiments/EXP3-power-20261025",
    "EXP4": "experiments/EXP4-deseason-20261025",
    "EXP5": "experiments/EXP5-scaled-n-20261025",
    "EXP6": "experiments/EXP6-bandwidth-20261025",
    "EXP7": "experiments/EXP7-memory-kernel-20261101",
}


class Cell(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")

    name: str
    seed_index: int = Field(ge=0)
    family: Literal["ind", "mmd", "growth"]
    n: int = Field(gt=0)
    target: int = Field(ge=0, description="replications in the saved main run")
    fields: dict[str, Any] = Field(description="the cell dict the driver's run_chunk consumes")

    @model_validator(mode="after")
    def _agrees_with_fields(self):
        f = self.fields
        if f.get("name") != self.name or f.get("idx") != self.seed_index or f.get("n") != self.n:
            raise ValueError(f"cell {self.name}: name, seed_index (idx) or n disagree with fields")
        return self

    def driver_cell(self) -> dict:
        return dict(self.fields)


class Prerequisite(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")

    path: str  # relative to the repo root
    sha256: str
    made_by: str  # the command that builds it

    def verify(self) -> bool:
        with open(os.path.join(REPO, self.path), "rb") as f:
            return hashlib.sha256(f.read()).hexdigest() == self.sha256


class ExperimentSpec(BaseModel):
    model_config = ConfigDict(frozen=True, extra="forbid")

    id: str
    folder: str
    driver: str
    mode: Literal["main"] = "main"
    masters: dict[str, int]  # key -> master; a cell uses masters[master_key.format(**cell.fields)]
    master_key: str
    seed_rule: str
    options: dict[str, Any] = Field(default_factory=dict, description="run options recorded in the manifest")
    writer: Literal["str", "fmt12g"] = Field(description="CSV float format: str() or EXP7's 12 digits")
    chunk: dict[str, dict[str, int]] = Field(description="family -> n -> replications per task")
    max_workers: int
    cells: list[Cell]
    prerequisites: list[Prerequisite] = Field(default_factory=list)
    source_manifests: list[str]

    @model_validator(mode="after")
    def _consistent(self):
        names = [c.name for c in self.cells]
        if len(set(names)) != len(names):
            raise ValueError(f"{self.id}: cell names are not unique")
        idx = {}
        for c in self.cells:
            key = (self.master_for(c), c.seed_index)
            if key in idx:
                raise ValueError(f"{self.id}: cells {idx[key]} and {c.name} share master and seed index")
            idx[key] = c.name
        for key, m in self.masters.items():
            if seeds.MASTERS.get((self.id, key)) != m:
                raise ValueError(f"{self.id}: master {key}={m} is not in kmem.seeds.MASTERS")
        return self

    def master_for(self, cell: Cell) -> int:
        return self.masters[self.master_key.format(**cell.fields)]

    def cell(self, name: str) -> Cell:
        for c in self.cells:
            if c.name == name:
                return c
        raise KeyError(f"{self.id}: no cell {name!r}")

    def chunk_for(self, cell: Cell) -> int:
        return self.chunk[cell.family][str(cell.n)]

    @property
    def total_rows(self) -> int:
        per = {c.name: c.target * len(c.fields.get("adjustments", [None])) for c in self.cells}
        return sum(per.values())


def spec_path(exp: str) -> str:
    return os.path.join(REPO, EXPERIMENTS[exp], "spec.json")


def load(exp: str) -> ExperimentSpec:
    with open(spec_path(exp)) as f:
        return ExperimentSpec.model_validate(json.load(f))
