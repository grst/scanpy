from pathlib import Path

import anndata
import pytest

import scanpy.api as sc


HERE = Path(__file__).parent


@pytest.fixture
def adata_kbet_sim():
    return anndata.read_h5ad(HERE / '_data' / 'kbet-sim.h5ad')


def test_kbet_needs_neighbors(adata_kbet_sim):
    with pytest.raises(ValueError):
        sc.tl.kbet(adata_kbet_sim)


def test_kbet_basic(adata_kbet_sim):
    sc.pp.neighbors(adata_kbet_sim)
    acceptance = sc.tl.kbet(adata_kbet_sim, k=75)  # type: float
    # R determined: rejection=0.0436 k=75
    assert .040 < acceptance < .050
