# S23B ICL–위성 표본의 공간 후보 장부 — 2026-10-22

## 결과

현재 S20A 기반 800개 target 표본은 `prestack_core.spatial_labels`의 5도 RA/Dec 격자에서 78개 셀을 차지한다. 이미 ICL–위성 외곽구조 분석에 들어간 94개 cluster/filter 행은 i에서 26행·13셀, i2에서 68행·33셀이다. 로컬 S23B CAMIRA 후보 1,333개는 총 100셀에 놓이고, 현재 800개 target의 좌표와 1초각 안에서 일치하는 후보는 없다. 숫자형 catalog ID끼리는 서로 다른 릴리스 사이에 결합하지 않고, 좌표를 비교했다.

후보 셀 중 78개는 기존 800개 footprint 셀과 겹치며 22개는 새 셀이다. 기존 연관 분석 target 셀 기준으로 보면 후보는 i의 13개 기존 셀 전부와 겹치고 새 셀 87개, i2의 33개 기존 셀과 겹치며 새 셀 67개를 제공한다. 따라서 현재 i2 target 셀 33개와 후보 측 새 셀을 합친 **100개는 target 셀 수의 낙관적 상한**이다. 앞선 설계 민감도에서 제시한 nominal 92개 선에는 카탈로그 기하만으로는 도달 가능하지만, 30비교 보수 설계선 180개보다 최소 80개 부족하다. 이 비교는 실제 독립 공분산 블록 수가 아니며, 후보 영상 지원이 없거나 matched-control 셀이 연결요소를 합치면 usable block은 더 적다.

이미 받은 800개 target 영상의 WCS 안에 후보 중심은 없었다. 기존 798개 matched-control 영상에는 후보 중심 7개가 들어가므로, 이 일곱 위치는 새로운 과학 target로 자동 채택하지 않고 control 오염 가능성 표식으로 유지한다. 나머지 후보들의 새 target coadd, mask, depth, PSF, 위성 member 완전도, shape catalog 및 matched control은 확인되지 않았다. 이 감사는 추가 HSC 요청·자료 반출·이미지 픽셀 판독을 하지 않았다.

## 재현

- 실행기: [`build_s23b_spatial_candidate_ledger.py`](build_s23b_spatial_candidate_ledger.py)
- 테스트: [`test_build_s23b_spatial_candidate_ledger.py`](test_build_s23b_spatial_candidate_ledger.py)
- 결과: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_spatial_candidate_ledger_2026-10-22_01/`
- TSV: `candidate_ledger.tsv`, 1,333행, SHA-256 `ab8df9fce1435d1d9eca15f742f43356cbd430b2159852a14d9b67963d83b775`
- JSON: `summary.json`, SHA-256 `f6ac867881d8184d02a580e9a18540b7513a90f530e51e6cdf04205db535e8d1`
- 입력 후보표 SHA-256 `1b59ca9941849fa7970c43c6e0374db6587716fc4ae5d71b11f84a245ebe6b2c`; parent manifest `95388814062e95f6e8429e2543a5432433515c9d9607680d2e873488d2758401`; 결합표 `9a6befa87a81257f8faf1d34fa7c019a91d64484fbfa64136bd93496dfc7b58c`; 기존 WCS 감사 JSON `071377813fcfbf15f298a5d48f3e7a3f0b74f60899d88c46eafb7d5bb58ffebb`.

실행 예:

```bash
python3 analysis/icl_wgl_membership/hsc_debug/build_s23b_spatial_candidate_ledger.py \
  --output-dir /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_spatial_candidate_ledger_2026-10-22_01
```

결과는 계속 `science_measurement=false`다. 신규 영상 및 지원 여부를 먼저 확인하지 않고 이 후보를 스택에 합치지 않는다.
