# S23B 후보 확장: 지원률과 독립 공간표본 상한 (2026-11-05)

## 요약

현재 800개 표본에 catalog-only S23B 후보 상한 1,333개를 더하면 총 2,133개다. 지원률이 그대로 유지된다고 가정할 때 50–75 kpc `i2` pair는 60% 기준 약 181개, 40% 기준 약 373개로 예상된다. 이는 후보 위치의 HSC 영상 가용성과 matched-control 지원을 모르는 낙관적 계획값이다.

| Satellite coverage floor | Current i2 pairs | Current center-eligible | Observed partial r (95% block interval) | Projected pairs after all candidates | Projected center-eligible |
|---:|---:|---:|---:|---:|---:|
| 60% | 68 | 15 | 0.293 `[0.028, 0.500]` | 181.3 | 40.0 |
| 40% | 140 | 30 | 0.249 `[-0.002, 0.454]` | 373.3 | 80.0 |

40% 기준은 현재 표본의 pair 수를 늘리지만 i2 부분상관 구간이 0을 포함하고 control 선택도 달라진다. Baseline으로 채택하거나 60% 결과와 합치지 않는다.

## 공간표본 상한

30개 비교를 보정하는 Fisher-z 계획에서 현재 40% 점추정 `r=0.249`를 가정해도 80% 검정력에 필요한 독립 블록은 251개다. 사후 탐색값 `r=0.293`에 대한 단일 사전 지정 endpoint 요구량은 92개, 30개 비교 보정 요구량은 180개다.

S23B candidate ledger에는 전체 100개 5° 셀이 있지만 기존 discovery와 겹치지 않는 셀은 22개뿐이다. 따라서 S23B만으로 독립 재현에 추가할 수 있는 target 셀 상한은 22개이며, 현재 800개 표본과 합친 셀 상한 100개도 가족보정 계획의 180개보다 작다. 셀은 실제 독립 공분산 블록이 아니며, 대조장 연결은 상한을 더 줄일 수 있다.

## 결정과 한계

S23B 이미지를 모두 받아도 이 catalog pool은 confirmatory한 공간 독립성 목표를 충족할 수 없다. 다만 탐색적 stack 확장에는 도움이 될 수 있다. 우선 지원 표본의 PSF/selection 검증을 마치고, 독립 블록이 더 큰 HSC 영역 또는 다른 적격 공개 footprint를 찾은 뒤 새 영상 수신을 계획한다.

본 분석은 기존 catalog 후보·지원률·power 요약만 읽었다. 후보별 coadd, mask, depth, PSF, member completeness, shape support, matched controls는 검증하지 않았다. HSC query와 데이터 수신은 없으며 `science_measurement=false`다.

## 2026-11-06 연결 tract 성분 감사

기존에 저장된 1,333개 S23B 후보 좌표를 800개 discovery lens와 함께 HSC
patch-QA tract-center 그래프에 배정했다. 840개 tract는 2도 인접 기준으로
4개 연결 성분을 만들었다. Discovery lens는 그중 세 성분에 499, 246, 55개씩
배정됐다. S23B 후보 1,333개는 기존 세 성분에 694, 564, 72개씩 배정됐고,
discovery 표본이 차지하지 않은 네 번째 성분에는 3개뿐이다. 따라서 22개 새
5도 면적 셀은 독립된 HSC 연결 영역 22개를 뜻하지 않는다. 후보 pool은 추가
분석 지원을 늘릴 수 있지만, 기존 footprint 안의 영역별 독립 재현 표본을
충분히 보충하지 못한다. 후보는 아직 이미지나 측정이 아니며 추가 영상 수신은
승인되지 않았다.

재현기 `audit_s23b_candidate_tract_region_support.py`, 결과
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_tract_region_support_2026-11-06_01/summary.json`
(SHA-256 `e9a9b9424ff1e8ec9f7014a2ab323645d0bc492c292e6714339167b79733c703`).
Tract 중심 기반 분류는 근사이며 polygon overlap이나 공식 독립 survey-field
정의를 대신하지 않는다. `science_measurement=false`를 유지한다.

## 재현성

- Astrafek script: [`assess_s23b_support_yield_and_blocks.py`](assess_s23b_support_yield_and_blocks.py)
- Machine output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_support_yield_and_blocks_2026-11-05_01/`
- The summary records SHA-256 values for all five input summaries and the analysis script.
