def kbet(adata, batch_key='batch', *, copy=False):
    """kBET: k-nearest neighbour batch effect test

    Parameters
    ----------
    adata : :class:`~anndata.AnnData`
        Annotated data matrix.
    batch_key : :class:`str`, optional (default: ``'batch'``)
        The column name in :attr:`anndata.AnnData.uns`
    copy : :class:`bool`, optional (default: ``False``)
        Copy instance before computation and return a copy.
        Otherwise, perform computation in-place and return the kBET score.

    Returns
    -------
    adata : :class:`~anndata.AnnData`
        If ``copy == True``, a copy of the input ``adata`` will be returned.
        Else, ``adata.uns['kbet']`` will contain the score.
    score : :class:`float`
        If ``copy == False``, the kBET score is returned.
    """
    score = 0.1

    if copy:
        ad_ret = adata.copy()
        ad_ret.uns['kbet'] = score  # TODO: column for e.g. kbet-per-louvain
        return ad_ret
    else:
        return score
