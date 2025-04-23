import collections
from dataclasses import (
    dataclass,
    field,
)

import functools

import itertools

import math

import os

import pandas as pd

import python_codon_tables as pct

from typing import (
    Final,
    Iterator,
    Self,
)


CodonTable = dict[str, dict[str, float]]
PathType = str | os.PathLike[str]


ALL_CODONS: Final[list[str]] = [
    c 
    for codon_dict in pct.get_codons_table("h_sapiens_9606").values()
    for c in codon_dict.keys()
] 
STOP_CODONS: Final[set[str]] = {"TAA", "TAG", "TGA"}


def read_fasta(filepath: PathType) -> Iterator[tuple[str, str]]:
    seq = ""
    header = ""
    with open(filepath, "r") as f:
        for line in f:
            if line.startswith(">"):
                seq = ""
                header = line[1 :].strip()
            else:
                seq += line.strip()
                yield header, seq


codon_to_aa = {
    codon: aa
    for aa, codon_dict in pct.get_codons_table("h_sapiens_9606").items()
    for codon in codon_dict.keys()
}
                

@dataclass
class CodonCollection:
    seq: str
    table: CodonTable
    _rf: int = field(
        default = 0,
    )
    codon_usage: pd.DataFrame = field(
        init = False,
    )

    def __post_init__(self) -> None:
        if not (len(self.seq) - self.rf) % 3 == 0:
            raise ValueError(
                "Sequence length must be a multiple of 3 or the reading frame must be adjusted."
            )

        if self.rf < 0 or self.rf >= 3:
            raise ValueError("Reading frame must be 0, 1, or 2.")

        self.codon_usage = self._create_codon_usage_table()

    def aa_iter(self) -> Iterator[str]:
        yield from (codon_to_aa[codon] for codon in self.codon_iter())

    def codon_usage_bias(self) -> pd.DataFrame:
        cu = self.codon_usage.copy()
        cu.replace(0.0, math.nan, inplace = True)
        for aa in cu.index:
            cu.loc[aa, :] = (
                cu.loc[aa, :]
                - self.table[aa][
                    max(self.table[aa].items(), key = lambda x: x[1])[0]
                ]
            ) ** 2
        return cu.sum(1) / cu.notna().sum(1)

    def codon_adaptability_index(self) -> float:
        cub = self.codon_usage_bias()
        weights = (
            self.table[aa][codon] / max(self.table[aa].items(), key = lambda x: x[1])[1]
            for aa, codon in zip(self.aa_iter(), self.codon_iter(), strict = True)
        )
        n_reciprocal = 1 / (len(self.seq) // 3)
        return functools.reduce(
            lambda x, y: x * y, 
            weights,
        ) ** n_reciprocal

    def set_codon_table(self, table: CodonTable) -> None:
        self.table = table

    @property
    def rf(self) -> int:
        return self._rf

    @rf.setter
    def rf(self, rf: int) -> None:
        if not (len(self.seq) - rf) % 3 == 0:
            raise ValueError(
                "Sequence length must be a multiple of 3 or the reading frame must be adjusted."
            )
        self._rf = rf

    def _create_codon_usage_table(self) -> pd.DataFrame:
        ca_matrix = pd.DataFrame(
            0,
            index = self.table.keys(),
            columns = ALL_CODONS,
        )
        for codon in self.codon_iter():
            aa = codon_to_aa[codon]
            ca_matrix.loc[aa, codon] += 1

        return (
            ca_matrix
                .apply(
                    lambda x: x / x.sum(),
                    axis = 1,
                )
                .fillna(0.0)
        )

    def codon_iter(self) -> Iterator[str]:
        yield from (
            self.seq[i : i + 3]
            for i in range(
                self.rf, len(self.seq) - self.rf, 3
            )
        )
    
