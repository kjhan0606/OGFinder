# S23B S16A ICL 대조 위치 count-only 감사 (2026-11-11)

## 범위와 실행

사용자가 승인한 계획 SHA-256
`832f2de18f5d512396cc84673e458932635c1f766ea91cdc04c7b47125738dfa`를
LagEunha에서 제출했다. 대상은 S16A 표적 raw count가 1,000 이상이면서 기존
S19A 427중심 count 감사 밖에 있고, 저장된 ICL 표적–대조 영상 쌍이 있는
은하단 116곳의 대조 위치다. S16A PDR1에서 짝지은 은하단 redshift 기준
3 proper-Mpc 원뿔의 `COUNT(DISTINCT object_id)`만 조회했다. 다섯 batch가
완료됐고, 결과는 `control_id`, `n_shape_sources` 두 열 116행이다.

대조 ID·정수 count 외 source ID, 좌표, shape 값, photo-z, image, PSF는
반환하지 않았다. 인증 정보는 결과 파일에 저장되지 않았다.
`science_measurement=false`다.

## 정량 결과

| ICL 필터 | 쌍 수 | 표적 count 중앙값 [Q16, Q84] | 대조 count 중앙값 [Q16, Q84] | paired 대조/표적 비율 중앙값 [Q16, Q84] | 양수 대조 수 |
|---|---:|---:|---:|---:|---:|
| 전체 | 116 | 4,778 [2,684, 6,818] | 4,727 [2,801, 6,900] | 0.988 [0.729, 1.201] | 115/116 |
| i | 98 | 4,916 [3,509, 6,877] | 4,740 [3,167, 7,248] | 0.988 [0.752, 1.189] | 97/98 |
| i2 | 18 | 2,895 [1,600, 5,942] | 4,450 [1,208, 5,319] | 0.992 [0.647, 1.721] | 18/18 |

표적 count가 1,000 이상인 116곳은 모두 양수였고, 대조는 115곳이 양수,
1곳이 0이다. 대조 원뿔 count 합은 550,390이지만 원뿔 간 중복을 제거한
고유 source 수가 아니다. 표적·대조 count 간 Spearman `rho=0.737`은 기술
통계로만 기록했으며 p-value는 계산하거나 사용하지 않았다.

## 해석과 한계

이 선택된 116쌍에서 raw S16A catalog-count의 중앙값은 서로 비슷하다. 이는
두 위치의 shape catalog availability가 큰 차이를 보이지 않는다는 사전
점검일 뿐이다. Source-redshift, shape quality, shear calibration, mask,
B-mode, radial-bin, effective weight와 estimator cuts가 적용되지 않아 usable
background source 수나 weak-lensing signal을 말해주지 않는다. count 차이와
비율도 측정값이나 검출 유의도가 아니다.

이들은 렌즈당 ICL 영상 대조 위치 한 곳씩이다. established WL 설계의 렌즈당
5개 random centers를 대체하지 않으므로 최종 random-subtracted WL stack에
바로 쓸 수 없다. 126개 고-count 표적 중 ICL 쌍 자체가 없는 10곳도 이 비교에
포함되지 않았다. S16A와 S19A source 행은 release/calibration 검증 전에
합치지 않는다.

## 다음 관문

남은 후보에서 기존 WL random-center 후보 5곳/렌즈의 적격성, 표적과의 공통
support, 위치 재사용을 로컬에서 점검한다. S16A random-center count-only 조회와
shape-row 반출은 서로 다른 요청이며 각각 정확한 계획과 별도 사용자 승인이
필요하다. 그 전까지 과학 측정은 수행하지 않는다.

## 재현 정보

- 제출기: [submit_s23b_s16a_icl_control_counts.py](submit_s23b_s16a_icl_control_counts.py)
- 로컬 감사기: [audit_s23b_s16a_icl_control_counts.py](audit_s23b_s16a_icl_control_counts.py)
- 승인 계획 경로: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_s16a_icl_control_count_plan_2026-11-11_01/`
- 승인 계획 SHA-256: `832f2de18f5d512396cc84673e458932635c1f766ea91cdc04c7b47125738dfa`
- CAS job IDs (batch order): `1074892`, `1074893`, `1074917`, `1074932`, `1074933`
- 결과 경로: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_s16a_icl_control_count_plan_2026-11-11_01/results/`
- CAS count TSV SHA-256: `b72d15093bdc3f4ddc6e03a1f85f9911b8999ed5e463e50c1eb7d6d4d4edbcd7`
- CAS count summary SHA-256: `136b86970f15f0d96118420a6ffe9ca9d510ad603b16cdfa28063fb55f4386f2`
- Local audit path: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_s16a_icl_control_count_audit_2026-11-11_01/`
- Local audit summary SHA-256: `eb12a650f0718d0fa74f5857eaa7042edd6a18d6798e7688d69fd95a7d9771b2`
- Paired target/control TSV SHA-256: `4223245a45264063e04fd05f2ed1b648fef239f88872486f754e7facc40eb79d`
- Control submitter SHA-256: `b6b3336e16df9abbd8100960e223a454aece4c6d799a10f62eace540a246fe1c`
- Local auditor SHA-256: `5f3d4eea0f26ab60fd317ee396fc79b0e10c92883d08c80aea2ed8cb75e5eeea`
- New CAS query beyond this approved plan, source-row transfer, image read: none.
