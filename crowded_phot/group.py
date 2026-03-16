"""Group overlapping sources for simultaneous fitting."""

import numpy as np


def group_sources(sources, radius_scale=2.0):
    """Group sources by proximity (Kron radius overlap).

    Parameters
    ----------
    sources : list of dicts with x, y, kron_radius, a
    radius_scale : float, grouping radius multiplier

    Returns
    -------
    groups : list of lists (indices into sources)
    """
    n = len(sources)
    if n == 0:
        return []

    # Compute grouping radii
    radii = np.array([
        radius_scale * s.get('kron_radius', 3.5) * max(s.get('a', 2.0), 1.0)
        for s in sources
    ])

    # Build adjacency
    visited = [False] * n
    groups = []

    for i in range(n):
        if visited[i]:
            continue

        group = []
        stack = [i]
        while stack:
            idx = stack.pop()
            if visited[idx]:
                continue
            visited[idx] = True
            group.append(idx)

            xi = sources[idx]['x']
            yi = sources[idx]['y']
            ri = radii[idx]

            for j in range(n):
                if visited[j]:
                    continue
                xj = sources[j]['x']
                yj = sources[j]['y']
                rj = radii[j]
                dist = np.sqrt((xi - xj)**2 + (yi - yj)**2)
                if dist < ri + rj:
                    stack.append(j)

        groups.append(group)

    return groups
