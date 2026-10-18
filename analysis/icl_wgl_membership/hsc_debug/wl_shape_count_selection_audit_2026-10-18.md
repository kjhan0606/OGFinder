# HSC S19A 렌즈별 shape 후보 수의 표본 의존성 감사

평가일: 2026-10-18 (Asia/Seoul)<br>
상태: count-only 탐색 감사; `science_measurement=false`

## 목적과 범위

800개 CAMIRA 렌즈에 대해 이미 승인·완료된 32개 CAS count-only batch를 다시 읽어,
렌즈 원뿔 안의 raw S19A shape 후보 수가 렌즈 적색편이와 richness에 따라 어떻게
달라지는지 확인했다. 조회 결과는 렌즈당 정수 하나뿐이다. source 행, 좌표, shape,
photo-z 값은 새로 조회하거나 저장하지 않았다.

각 렌즈는 기하 적색편이와 고유반경 3 Mpc 원뿔을 사용했다. `count(DISTINCT object_id)`는
shape 품질, HSC Y3 tomographic bin, magnitude, resolution, B-mode, 반경 컷을 적용하기
전의 후보량이다. 같은 source가 여러 렌즈 원뿔에 있으면 렌즈마다 다시 세므로 합계는
고유 source 수가 아니라 lens-source pair 후보량이다.

## 결과

800개 렌즈의 3,069,544 pair 후보에서 원뿔 count가 0인 렌즈는 168개(21.0%)였다.
렌즈별 count의 중앙값은 4,438, 10/90 백분위는 0/6,939, 최댓값은 9,767이다.
최소 1,000개 후보인 렌즈는 610개, 3,000개는 553개, 5,000개는 292개였다. 이 숫자는
최종 분석 적격 수나 측정 S/N이 아니다.

| 렌즈 기하 적색편이 | 렌즈 수 | count=0 | 0 비율 | count 중앙값 |
|---|---:|---:|---:|---:|
| 0.29–0.35 | 147 | 22 | 15.0% | 6,473 |
| 0.35–0.40 | 117 | 21 | 17.9% | 5,307 |
| 0.40–0.45 | 204 | 42 | 20.6% | 4,573 |
| 0.45–0.50 | 162 | 28 | 17.3% | 4,354 |
| 0.50–0.57 | 170 | 55 | 32.4% | 3,601 |

count와 기하 적색편이의 Spearman 상관은 `rho=-0.427`이었다. 이는 기술적 의존성
요약이며 유의도 검정은 하지 않았다. 고적색편이 렌즈에서 각거리 원뿔이 작아지고,
관측 깊이·footprint·catalog 선택도 달라질 수 있으므로 물리적 관계로 해석하지 않는다.
Richness와 count의 Spearman 상관은 `rho=+0.048`로 단조 관계가 약했다. 이는 shape
count가 richness만으로 예측되지 않는다는 점을 보여주는 표본 계획 진단이다.

## 분석에 대한 의미

원뿔 count가 0인 렌즈가 21%이고 적색편이별 차이도 있으므로, 800개 렌즈가 균일하게
약한렌즈에 기여한다고 가정하면 안 된다. 최종 렌즈 표본은 source 품질·tomographic
bin·photo-z 기하·반경·B-mode 컷과 random-lens 공분산을 적용한 뒤 다시 정해야 한다.
그때 richness와 z, sky footprint를 함께 보고하고, 탈락 렌즈와 유효 source weight를
기록해야 한다. 현재 count-only 집계로는 최종 유효 source 수, shear, `DeltaSigma`, 또는
적정 표본력이 결정되지 않는다.

## 재현성과 산출물

```bash
python3 analysis/icl_wgl_membership/hsc_debug/audit_wl_shape_count_selection.py
```

계산기는 preflight의 manifest/hash, 32 count CSV, 800개 unique lens ID, count 합계,
원본 sample manifest와 기하 적색편이 일치를 검증한다. 재현 요약과 10개 적색편이/
richness 구간 표는
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_shape_count_selection_audit_2026-10-18_01/`
에 있다. 요약은 모든 source transfer 및 `science_measurement` 상태를 false로 기록한다.
