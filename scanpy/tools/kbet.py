from typing import Union

from scipy.sparse import spmatrix
from scipy.stats import chisquare
from anndata import AnnData


def kbet(
    adata: AnnData,
    batch_key: str = 'batch',
    *,
    adjacency: spmatrix = None,
    copy: bool = False,
) -> Union[AnnData, float]:
    """kBET: k-nearest neighbour batch effect test.

    Use the heuristic :func:`sc.pp.kbet_n_neighbors` to find the ideal
    neighborhood size to pass to :func:`sc.pp.neighbors`.

    Parameters
    ----------
    adata
        Annotated data matrix.
    batch_key
        The column in :attr:`anndata.AnnData.uns` to use as batch ID.
    adjacency
        Sparse adjacency matrix of the graph, defaults to
        ``adata.uns['neighbors']['connectivities']``.
    copy
        Copy instance before computation and return a copy.
        Otherwise, perform computation in-place and return the kBET score.

    Returns
    -------
    adata
        If ``copy == True``, a copy of the input ``adata`` will be returned.
        Else, ``adata.uns['kbet']`` will contain the score.
    score
        If ``copy == False``, the kBET score is returned.
    """
    if adjacency is None:
        if 'neighbors' not in adata.uns:
            raise ValueError('No neighbors found. Provide the `adjacency` parameter or run `sc.pp.neighbors(adata)`')
        adjacency = adata.uns['neighbors']['connectivities']

    f_obs = None
    f_exp = None
    ddof = None
    scores, p_vals = chisquare(f_obs, f_exp, ddof)
    score_mean = scores.mean()

    if copy:
        ad_ret = adata.copy()
        ad_ret.uns['kbet'] = score_mean  # TODO: column for e.g. kbet-per-louvain
        return ad_ret
    else:
        return score_mean
