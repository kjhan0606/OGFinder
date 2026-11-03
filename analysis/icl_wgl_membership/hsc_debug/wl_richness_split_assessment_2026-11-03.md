# 약중력 풍부도 분할 진단 — 2026-11-03

## 목적과 방법

기존 HSC S19A source shape 행을 사용해 CAMIRA 풍부도가 낮은 렌즈와 높은 렌즈의
cluster-only 약한렌즈 신호가 다른지 진단했다. 표본은 HSC 형상 profile을 보기 전에
부모 800개 전체에서 풍부도와 lens ID 순으로 정렬해 정확히 400개씩 나눴다. 낮은 쪽의
최대 풍부도는 20.490, 높은 쪽의 최소 풍부도는 20.498이다. 각 그룹에 동일한 HSC Y3
4-bin `n(z)`, 기존 raw B-mode 제외값 0, 기존 일곱 proper-radius bin을 적용했다.

공통 equal-area 하늘 셀 52개를 하나씩 삭제하고, lens 또는 source가 삭제 셀에 속한 모든
쌍을 제외한 뒤 두 그룹 stack과 bin별 보정을 다시 계산했다. 오류는 상위−하위 profile의
짝지은 delete-one covariance에서 얻었다. 계산은 lageunha CPU에서 수행했으며, HSC 접속,
질의, 데이터 다운로드는 없었다.

## 결과

각 부모 그룹은 400개다. shape 선택·적색편이 효율·반경을 모두 통과한 표본은 낮은 풍부도
그룹에서 렌즈 299개와 1,043,918쌍, 높은 쪽에서 렌즈 310개와 1,109,086쌍이다.

`DeltaSigma`와 상위−하위 차이의 단위는 `10^12 M_sun Mpc^-2`다.

| Proper R (Mpc) | 낮은 풍부도 | 높은 풍부도 | 높은−낮은 | 짝지은 잭나이프 오차 |
|---|---:|---:|---:|---:|
| 0.1–0.2 | 153.33 | 168.45 | 15.13 | 21.71 |
| 0.2–0.4 | 93.94 | 135.52 | 41.58 | 12.26 |
| 0.4–0.7 | 53.72 | 89.10 | 35.38 | 8.71 |
| 0.7–1.0 | 39.94 | 67.86 | 27.92 | 6.76 |
| 1.0–1.5 | 30.46 | 45.78 | 15.32 | 4.88 |
| 1.5–2.0 | 18.97 | 36.21 | 17.23 | 3.56 |
| 2.0–3.0 | 10.50 | 19.75 | 9.25 | 2.12 |

상위 풍부도 그룹의 중심 `DeltaSigma`가 일곱 반경 bin 모두에서 크다. 접선 profile의
7차원 상위−하위 차이와 잭나이프 공분산을 사용한 조건부 Gaussian 진단은
`chi2=82.66` (rank 7, 명목 `p=3.95e-15`)다. 이 p-value는 추정된 잭나이프 공분산과
공간 셀의 독립성 가정에 민감하므로 유의도나 검출로 해석하지 않는다. 접선성분 7×7
차이 공분산은 양의 고유값을 갖고 수치 조건수는 164다.

부모 sample manifest에서 두 풍부도 그룹의 렌즈 적색편이 중앙값은 각각 0.442, 0.427이고
평균은 0.435, 0.426이다. 부모 800개에서 풍부도와 적색편이의 Spearman 상관은 −0.046
(`p=0.195`)이지만, 이것만으로 그룹별 redshift, selection, source geometry 효과가 제거된
것은 아니다. 교차성분 차이는 반경별로 일관된 유의한 null 검증으로 취급하지 않았고,
이 분석은 농도·질량·후광 모델을 맞추지 않았다.

## 해석과 다음 관문

기존 표적 표본에서 풍부도에 따라 lensing profile이 달라지는 현상은 보이며, 향후 사전에
고정한 풍부도 의존성 가설을 검토할 동기를 준다. 다만 이는 무작위 중심 차감 전의
cluster-only profile이다. 실제 CAMIRA 풍부도–질량 보정, lens selection과 redshift 교란,
형상 보정 및 `n(z)` 불확실도, reduced shear, 검증된 연속 HSC 관측영역 공분산을 포함하지
않았으므로 질량 관계나 과학 측정으로 승격할 수 없다. 현재 상태는 `science_measurement=false`다.

다음 필수 단계는 승인된 범위의 무작위 대조 중심 집계/형상 자료 확보 후, 동일 source와
공간 영역을 함께 삭제하는 cluster-minus-random paired stack이다. 특히 새 random-center
count-only 쿼리는 별도 사용자 승인 대기 중이므로 이 실행에서 제출하지 않았다.

## 재현 및 출처

- 계산기: `audit_wl_richness_split.py`
- 산출물: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_richness_split_2026-11-03_01/summary.json`
- 산출물 SHA-256: `20f28c5cdec8a8a6ca9d3041bfcdad352999d1a6e4b3e809eb6bc5e6e1827a33`
- 계산기 SHA-256: `f108dddf49812afc4405df705e1663734f2d32f19b8c7504ace0aa2e4eb06c13`
- 결과는 기존의 800-lens shape-row export, count preflight/plan, sample manifest,
  HSC Y3 SACC `n(z)`를 `verify_and_load_pairs`로 재검증한 뒤 만들었다.
