"""Builders for synthetic sequences, hits, and barcode matrices."""

import numpy as np
import pytest

from fshdetector import barcode


def random_dna(length: int, seed: int) -> str:
    generator = np.random.default_rng(seed)
    return "".join(generator.choice(list("ACGT"), size=length))


def substitute(sequence: str, rate: float, generator) -> str:
    """Replaces a fraction of the bases with other bases."""
    bases = np.array(list(sequence))
    for i in np.flatnonzero(generator.random(len(bases)) < rate):
        bases[i] = generator.choice([b for b in "ACGT" if b != bases[i]])
    return "".join(bases)


def array_matrix(units: int, sites: int = 80, seed: int = 1) -> np.ndarray:
    """An allele array whose units carry independent random alleles."""
    generator = np.random.default_rng(seed)
    return generator.integers(0, 2, size=(units, sites)).astype(np.int8)


def sample_model(*arrays: np.ndarray) -> barcode.Model:
    pool = np.concatenate(arrays)
    return barcode.Model(pool, np.full(pool.shape[1], 0.01))


@pytest.fixture(name="generator")
def fixture_generator():
    return np.random.default_rng(7)
