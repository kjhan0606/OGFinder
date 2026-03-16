"""Find candidate merge/separate pairs using spatial proximity."""

import numpy as np
from scipy.spatial import cKDTree


def find_candidate_pairs(objects, kron_radius, flux_auto, alpha=1.1,
                         radius_scale=2.0, min_flux=0.0):
    """
    Find candidate pairs of potentially overlapping sources.

    Uses cKDTree for O(N log N) neighbor search.

    Parameters
    ----------
    objects : structured array
        SEP objects with 'x', 'y', 'a', 'b' fields.
    kron_radius : np.ndarray
        Kron radii for each source (unused, kept for API compat).
    flux_auto : np.ndarray
        AUTO fluxes for filtering.
    alpha : float
        Distance factor: pair if d < alpha * (R_i + R_j).
    radius_scale : float
        Unused (kept for API compat).
    min_flux : float
        Minimum flux to consider a source.

    Returns
    -------
    pairs : list of (int, int)
        List of candidate pair indices.
    """
    n = len(objects)
    if n < 2:
        return []

    x = objects['x']
    y = objects['y']
    a = np.maximum(objects['a'], 0.5)
    b = np.maximum(objects['b'], 0.5)

    # Effective radii: mean of semi-major and semi-minor axes
    radii = 0.5 * (a + b)
    radii = np.maximum(radii, 1.0)

    # Flux filter
    valid = flux_auto > min_flux
    indices = np.where(valid)[0]
    if len(indices) < 2:
        return []

    # Build KD-tree
    coords = np.column_stack([x[indices], y[indices]])
    tree = cKDTree(coords)

    # Maximum search radius
    max_radius = alpha * 2 * radii[indices].max()

    pairs = []
    seen = set()

    for k, idx_i in enumerate(indices):
        r_i = radii[idx_i]
        search_r = alpha * (r_i + radii[indices].max())
        neighbors = tree.query_ball_point(coords[k], r=search_r)

        for m in neighbors:
            if m <= k:
                continue
            idx_j = indices[m]
            r_j = radii[idx_j]
            d = np.sqrt((x[idx_i] - x[idx_j])**2 + (y[idx_i] - y[idx_j])**2)
            if d < alpha * (r_i + r_j):
                pair = (min(idx_i, idx_j), max(idx_i, idx_j))
                if pair not in seen:
                    seen.add(pair)
                    pairs.append(pair)

    return pairs


def build_merge_groups(pairs, n_sources):
    """
    Build connected components from candidate pairs.

    Parameters
    ----------
    pairs : list of (int, int)
    n_sources : int

    Returns
    -------
    groups : list of list of int
        Each group is a list of source indices that should be considered together.
    """
    if not pairs:
        return []

    # Union-Find
    parent = list(range(n_sources))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb

    for a, b in pairs:
        union(a, b)

    # Group by root
    groups_dict = {}
    involved = set()
    for a, b in pairs:
        involved.add(a)
        involved.add(b)

    for idx in involved:
        root = find(idx)
        if root not in groups_dict:
            groups_dict[root] = []
        groups_dict[root].append(idx)

    return [sorted(g) for g in groups_dict.values() if len(g) >= 2]
