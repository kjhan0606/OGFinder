#!/usr/bin/env python3
"""Record why this debug does not write a stack mass.

Oguri et al. 2018 (arXiv:1701.00818), section 3.2, compare the catalog with
M_200m. The sentence that N_mem > 15 roughly corresponds to
M_200m ≳ 10^14 h^{-1} Msun is a sample threshold, not a per-cluster M_200c.
Murata et al. 2019 (arXiv:1904.07524) define halo mass inside R_200m.
Chiu et al. 2019 (arXiv:1909.02042) calibrate richness to M_500.
The scaled stack aligns on M_200c. This script does not convert those masses
and does not write clusters.json.
"""
import json
import os
import sys

ROOT = os.environ.get('ASTRAFEK_DATA', '/gpfs/kjhan/Astrafek/data')


def main():
    summary_path = os.path.join(ROOT, 'logs', 'selection_summary.json')
    with open(summary_path, encoding='utf-8') as handle:
        summary = json.load(handle)
    gate = {
        'mass_column_in_catalog': False,
        'alignment_mass': 'M_200c',
        'oguri2018_mass_definition': 'M_200m',
        'oguri2018_arxiv': '1701.00818',
        'murata2019_mass_definition': 'M_200m',
        'murata2019_arxiv': '1904.07524',
        'chiu2019_mass_definition': 'M_500',
        'chiu2019_arxiv': '1909.02042',
        'clusters_json_written': False,
        'stack_executed': False,
        'science_measurement': False,
        'status': 'stopped',
    }
    logs = os.path.join(ROOT, 'logs')
    with open(os.path.join(logs, 'mass_gate.json'), 'w', encoding='utf-8') as handle:
        json.dump(gate, handle, indent=1)
        handle.write('\n')
    text = """# ICL 스택 디버그 상태

작업 디렉터리: `{root}`
프로그램: `OGFinder/plugins/icl/icl_stack.py`
이 실행은 파이프라인 디버그이다. `science_measurement` 는 거짓이다.
절차는 `scripts/run_debug.sh` 가 재실행한다.

## 목록과 선택

- 목록은 `catalog/camira_s20a_wide_sm.tbl` 이다. 데이터 행은 {n_rows}개이고, BCG 분광 적색이동이 있는 행은 {n_spec}개이다.
- 분광 행에서 `|z_spec-z_cl|/(1+z_cl)` 의 중앙값은 {p50:.5f}, 68퍼센타일은 {p68:.5f}이다.
- 선택 창은 그 분수가 0.01보다 작은 행이다. Oguri et al. 2018 (arXiv:1701.00818) 은 `Δz/(1+z)` 의 산란이 대부분의 구간에서 0.01보다 좋다고 적는다.
- 적색이동 구간은 {z_lo:.2f} 이상 {z_hi:.2f} 미만이다. 0.10부터 1.20까지 0.05 폭의 정확한 경계 가운데, 이 창을 통과한 행이 가장 많은 구간이다. 그 행은 {n_window}개이다.
- 그 구간의 `N_mem` 90퍼센타일은 {nmem:.4f}이다. 이 값 이상인 {n_rich}개 가운데 richness 가 큰 순서로, 서로 0.3도보다 멀리 떨어진 4개를 골랐다. 최소 각거리는 {sep:.4f}도이다.
- 고른 식별자는 {ids} 이다. 표는 `logs/selected_clusters.tsv` 이다.

## 질량 관문

Oguri et al. 2018 (arXiv:1701.00818) 3.2절의 비교 질량은 `M_200m` 이다. `N_mem > 15` 가 `M_200m` 약 `10^14 h^{{-1}} Msun` 이상에 대응한다는 문장은 표본 문턱이다. 은하단마다의 `M_200c` 식이 아니다.

Murata et al. 2019 (arXiv:1904.07524) 의 헤일로 질량은 `R_200m` 안의 질량이다. Chiu et al. 2019 (arXiv:1909.02042) 의 관계는 `M_500` 이다.

`u` 정렬에 쓰는 질량은 `M_200c` 이다. 이 스크립트는 다른 정의를 `M_200c` 로 환산하지 않는다. `clusters.json` 은 만들지 않았고 `icl_stack.py` 는 실행하지 않았다.

## 영상

cutout 상태는 `logs/cutout_status.json` 에 있다. 영상이 없으면 스택 입력 파일을 만들지 않는다.
""".format(
        root=ROOT,
        n_rows=summary['n_rows'],
        n_spec=summary['n_spec'],
        p50=summary['frac_p50'],
        p68=summary['frac_p68'],
        z_lo=summary['z_lo'],
        z_hi=summary['z_hi'],
        n_window=summary['n_window'],
        nmem=summary['nmem_p90'],
        n_rich=summary['n_rich'],
        sep=summary['min_sep_deg'],
        ids=' '.join(str(i) for i in summary['ids']),
    )
    with open(os.path.join(logs, 'debug_status.md'), 'w', encoding='utf-8') as handle:
        handle.write(text)
    print('mass gate %s' % gate['status'])


if __name__ == '__main__':
    sys.exit(main())
