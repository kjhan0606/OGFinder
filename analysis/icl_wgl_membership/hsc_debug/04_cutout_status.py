#!/usr/bin/env python3
"""Record whether HSC cutouts are present. This script does not download them.

Images come from the HSC archive after a login. The release page is
https://hsc-release.mtk.nao.ac.jp/doc/index.php/data/
PDR3 calls the global-sky coadd the patch product and the local-sky coadd
calexp (https://hsc-release.mtk.nao.ac.jp/doc/index.php/available-data__pdr3/).
A later stack uses the global-sky frames. One calexp file is reserved for the
refusal check and is not a stack input.
"""
import json
import os
import sys

ROOT = os.environ.get('ASTRAFEK_DATA', '/gpfs/kjhan/Astrafek/data')


def main():
    folder = os.path.join(ROOT, 'cutouts')
    os.makedirs(folder, exist_ok=True)
    names = sorted(
        name for name in os.listdir(folder)
        if name.endswith(('.fits', '.fits.gz', '.fz'))
    )
    status = {
        'n_files': len(names),
        'files': names,
        'downloaded_by_this_script': False,
        'science_measurement': False,
        'status': 'empty' if not names else 'present',
    }
    path = os.path.join(ROOT, 'logs', 'cutout_status.json')
    with open(path, 'w', encoding='utf-8') as handle:
        json.dump(status, handle, indent=1)
        handle.write('\n')
    print('cutouts %s n %d' % (status['status'], status['n_files']))


if __name__ == '__main__':
    sys.exit(main())
