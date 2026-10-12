#!/usr/bin/env python3
"""Match CAMIRA member rows to the selected clusters.

The file is camira_s20a_wide_sm_member.dat.gz from the Google Drive folder
linked by the Oguri cluster_catalogs README
(https://drive.google.com/drive/folders/19dKiNs7Wdq44X3AtYy2opOf8ulHTLWeM).
It has no header. Eight columns were identified by matching the cluster
coordinate, richness, and redshift to camira_s20a_wide_sm.tbl:

    ra_cl_deg dec_cl_deg N_mem z_cl ra_deg dec_deg field7 weight

field7 is stored and is not used as a magnitude, a radius, or a redshift.
Every member row of a matched cluster is kept. There is no weight cut.
The BCG position used later is the cluster-catalog coordinate. A member row
on that coordinate is marked role=bcg.
"""
import gzip
import hashlib
import json
import os
import subprocess
import sys

import numpy as np

ROOT = os.environ.get('ASTRAFEK_DATA', '/gpfs/kjhan/Astrafek/data')
FILE_ID = '14mRh_zPdrIlTSa5_pw3iFk5lN8M3cXzY'
URL = (
    'https://drive.usercontent.google.com/download?id=%s'
    '&export=download&confirm=t' % FILE_ID
)
SHA = 'bf0c1fec5e11c6378ea69e3a65e49e9b32beecc83e5e68000bfa997e4ef2059e'
DEST_NAME = 'camira_s20a_wide_sm_member.dat.gz'
MATCH_DEG = 1e-4
EXACT_ARCSEC = 0.05


def sha256(path):
    digest = hashlib.sha256()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def ensure_catalog():
    dest = os.path.join(ROOT, 'catalog', DEST_NAME)
    os.makedirs(os.path.dirname(dest), exist_ok=True)
    if os.path.isfile(dest) and sha256(dest) == SHA:
        print('member catalog present %s' % dest)
        return dest
    tmp = dest + '.partial'
    subprocess.check_call(['curl', '-fsSL', '--retry', '2', '--max-time', '180', '-L', '-o', tmp, URL])
    found = sha256(tmp)
    if found != SHA:
        os.remove(tmp)
        raise SystemExit('member catalog sha256 is %s' % found)
    os.replace(tmp, dest)
    print('member catalog written %s' % dest)
    return dest


def read_selected():
    path = os.path.join(ROOT, 'logs', 'selected_clusters.tsv')
    with open(path, encoding='utf-8') as handle:
        header = handle.readline().rstrip('\n').split('\t')
        rows = [line.rstrip('\n').split('\t') for line in handle if line.strip()]
    need = ['id', 'ra_deg', 'dec_deg', 'z_cl', 'N_mem']
    if any(name not in header for name in need):
        raise SystemExit('selected_clusters.tsv is missing a column')
    index = {name: header.index(name) for name in header}
    return index, rows


def angular_sep_arcsec(ra1, dec1, ra2, dec2):
    dra = (ra1 - ra2) * np.cos(np.deg2rad(dec1))
    return np.hypot(dra, dec1 - dec2) * 3600.0


def main():
    path = ensure_catalog()
    with gzip.open(path, 'rt') as handle:
        members = np.loadtxt(handle)
    if members.ndim != 2 or members.shape[1] != 8:
        raise SystemExit('member catalog does not have 8 columns')
    groups = {}
    for row in members:
        key = (round(float(row[0]), 5), round(float(row[1]), 5))
        groups.setdefault(key, []).append(row)
    index, selected = read_selected()
    out_rows = []
    cluster_rows = []
    n_bcg = 0
    for item in selected:
        ra = float(item[index['ra_deg']])
        dec = float(item[index['dec_deg']])
        key = (round(ra, 5), round(dec, 5))
        block = groups.get(key)
        if block is None:
            raise SystemExit('cluster %s is absent from the member catalog' % item[index['id']])
        block = np.vstack(block)
        # The member file prints z_cl to three decimals. 0.001 still rejects a different cluster.
        if abs(float(block[0, 2]) - float(item[index['N_mem']])) > 0.001:
            raise SystemExit('cluster %s richness does not match the member file' % item[index['id']])
        if abs(float(block[0, 3]) - float(item[index['z_cl']])) > 0.001:
            raise SystemExit('cluster %s redshift does not match the member file' % item[index['id']])
        sep = angular_sep_arcsec(block[:, 4], block[:, 5], ra, dec)
        order = np.argsort(sep)
        bcg_at = None
        for pos in order:
            if sep[pos] <= EXACT_ARCSEC:
                bcg_at = int(pos)
                break
        if bcg_at is not None:
            n_bcg += 1
        for pos in order:
            role = 'bcg' if bcg_at is not None and int(pos) == bcg_at else 'member'
            out_rows.append((
                item[index['id']], role,
                float(block[pos, 4]), float(block[pos, 5]),
                float(block[pos, 7]), float(block[pos, 6]), float(sep[pos]),
            ))
        cluster_rows.append((
            item[index['id']], len(block),
            int(np.sum(block[:, 7] >= 1.0)),
            'yes' if bcg_at is not None else 'no',
            float(sep[order[0]]), float(block[:, 7].sum()),
        ))
    logs = os.path.join(ROOT, 'logs')
    os.makedirs(logs, exist_ok=True)
    member_path = os.path.join(logs, 'members.tsv')
    with open(member_path, 'w', encoding='utf-8') as handle:
        handle.write('id\trole\tra_deg\tdec_deg\tweight\tfield7\tsep_arcsec\n')
        for row in out_rows:
            handle.write('%s\t%s\t%.8f\t%.8f\t%.6g\t%.4f\t%.4f\n' % row)
    cluster_path = os.path.join(logs, 'member_clusters.tsv')
    with open(cluster_path, 'w', encoding='utf-8') as handle:
        handle.write('id\tn_members\tn_weight_eq_1\tbcg_row\tnearest_sep_arcsec\tweight_sum\n')
        for row in cluster_rows:
            handle.write('%s\t%d\t%d\t%s\t%.4f\t%.4f\n' % row)
    summary = {
        'source_file': DEST_NAME,
        'sha256': SHA,
        'drive_folder': '19dKiNs7Wdq44X3AtYy2opOf8ulHTLWeM',
        'drive_file_id': FILE_ID,
        'n_member_rows': int(members.shape[0]),
        'n_clusters_in_member_file': len(groups),
        'columns': [
            'ra_cl_deg', 'dec_cl_deg', 'N_mem', 'z_cl',
            'ra_deg', 'dec_deg', 'field7', 'weight',
        ],
        'field7_min': float(members[:, 6].min()),
        'field7_max': float(members[:, 6].max()),
        'field7_note': (
            'The member file has no header. field7 is stored and is not used '
            'as a magnitude, a radius, or a redshift.'),
        'weight_note': (
            'Every row of a matched cluster is kept. There is no weight cut. '
            'The weight sum is a column of the member file and is not N_mem.'),
        'n_selected': len(selected),
        'n_matched_clusters': len(cluster_rows),
        'n_bcg_rows': n_bcg,
        'n_member_rows_selected': len(out_rows),
        'bcg_position': 'cluster catalog coordinate',
        'science_measurement': False,
    }
    with open(os.path.join(logs, 'member_summary.json'), 'w', encoding='utf-8') as handle:
        json.dump(summary, handle, indent=1)
        handle.write('\n')
    print(
        'members clusters %d rows %d bcg_rows %d'
        % (summary['n_matched_clusters'], summary['n_member_rows_selected'], n_bcg)
    )


if __name__ == '__main__':
    sys.exit(main() or 0)
