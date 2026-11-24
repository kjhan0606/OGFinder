#!/usr/bin/env python3
"""Compare concentration-associated injection support across amplitudes."""
import argparse
import json
from pathlib import Path

import numpy as np

from audit_s23b_icl_satellite_field_omission import DATA, read_tsv, sha256, write_tsv

DIRECTORIES={
    0.5: DATA/'s23b_injection_support_amp_05_2026-11-21_01',
    1.0: DATA/'s23b_injection_support_selection_2026-11-21_01',
    2.0: DATA/'s23b_injection_support_amp_20_2026-11-21_01',
}


def run(output):
    if output.exists():
        raise FileExistsError(output)
    tables={}
    for amp,directory in DIRECTORIES.items():
        summary=json.loads((directory/'summary.json').read_text())
        path=directory/'matched_support_ledger.tsv'
        if sha256(path)!=summary['output_sha256'][path.name]:
            raise ValueError('Changed ledger')
        rows=read_tsv(path)
        mapping={(r['target_id'],r['filter']):r for r in rows}
        if len(mapping)!=len(rows):
            raise ValueError('Duplicate matched target')
        tables[amp]=mapping
    if not tables[0.5].keys()==tables[1.0].keys()==tables[2.0].keys():
        raise ValueError('Matched samples differ')
    results=[]
    for band_index,band in enumerate(('i','i2')):
        cells={}
        for key,base in tables[1.0].items():
            if key[1]!=band or base.get('quartile') not in ('low','high'):
                continue
            records=[tables[a][key] for a in (0.5,1.0,2.0)]
            if any((r['quartile'],r['concentration'],r['ra'],r['dec']) !=
                   (base['quartile'],base['concentration'],base['ra'],base['dec']) for r in records):
                raise ValueError('Quartile or coordinate mismatch')
            cell=(int(np.floor(((float(base['ra'])-84.056247280625)%360)/5)),
                  int(np.floor((float(base['dec'])+90)/5)))
            cells.setdefault(cell,[]).append((base['quartile'],[r['supported']=='True' for r in records]))
        keys=sorted(cells)
        def statistic(sample):
            low=np.asarray([v for group,v in sample if group=='low'],float)
            high=np.asarray([v for group,v in sample if group=='high'],float)
            if not len(low) or not len(high):
                return None
            return high.mean(axis=0)-low.mean(axis=0)
        source=[r for c in keys for r in cells[c]]
        point=statistic(source)
        rng=np.random.default_rng(2026112125+band_index)
        draws=[]
        for _ in range(5000):
            sample=[r for i in rng.integers(0,len(keys),len(keys)) for r in cells[keys[i]]]
            stat=statistic(sample)
            if stat is not None:
                draws.append(stat)
        if len(draws)<4500:
            raise ValueError('Too few valid draws')
        draws=np.asarray(draws)
        for index,amplitude in ((0,0.5),(2,2.0)):
            diff=draws[:,index]-draws[:,1]
            low,high=np.percentile(diff,[2.5,97.5])
            results.append(dict(filter=band,amplitude=amplitude,reference_amplitude=1.0,n_cells=len(keys),
                                n_extreme_clusters=len(source),difference_of_gaps=float(point[index]-point[1]),
                                interval_low=float(low),interval_high=float(high),valid_draws=len(draws)))
    output.mkdir(parents=True)
    table=output/'paired_amplitude_support.tsv'
    write_tsv(table,results)
    summary=dict(science_measurement=False,replicas=5000,
                 inputs={str(p/'summary.json'):sha256(p/'summary.json') for p in DIRECTORIES.values()},
                 runner_sha256=sha256(Path(__file__)),output_sha256={table.name:sha256(table)},
                 limitation='Post-hoc unadjusted extreme-quartile contrasts in the selected one-member '
                            'injection sample, with fixed quartile boundaries and marginal intervals. '
                            'This does not measure population completeness or PSF wings.')
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output.resolve())
