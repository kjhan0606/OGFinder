# S16A 표적–대조장 원시 지지와 보정 경로 점검 (2026-11-17)

기존 114개 렌즈·570개 무작위 중심의 승인된 count-only 결과만 수치 분석했다. 원래 선정 원장은 변경하지 않았다. 새 CAS 조회, source 행·영상·PSF·photo-z 또는 보정 데이터 수신은 없다. 보정 경로는 공개 공식 문서와 해당 문서가 연결한 소프트웨어를 읽어 정적 검토했다. `science_measurement=false`를 유지한다.

## 원시 개수와 관측 지지

원본 count/선정/표적 manifest 및 SQL 해시를 확인하고, 모든 random ID를 표적 원시 개수와 일대일 결합했다. 저장된 각거리와 실제 구면 각거리, H0=70·Omega_m=0.3에서의 proper 거리를 대조했다. 570개 무작위 위치·object ID는 전역 중복이 없고 각 렌즈는 슬롯 1–5를 정확히 유지한다.

| ICL 필터 | 렌즈 | random 중심 | 0개 중심 | 5곳 모두 양수인 렌즈 | 5곳 모두 ≥1,000 | 5곳 비율 모두 0.5–2 |
|---|---:|---:|---:|---:|---:|---:|
| all | 114 | 570 | 45 | 96 | 86 | 67 |
| i | 96 | 480 | 31 | 84 | 75 | 60 |
| i2 | 18 | 90 | 14 | 12 | 11 | 7 |

마지막 두 열은 진단용 민감도 비교이며 과학 컷이나 새로 선정된 표본 수가 아니다. 원시 개수 양수도 radial-bin별 유효 source/weight를 보장하지 않는다.

570개 개별 중심의 `N_random/N_lens` 중앙값은 **0.927**, 16–84백분위는 **0.402–1.173**이다. 렌즈별 고정 5곳 평균을 먼저 계산하면 그 비율의 중앙값은 **0.897 [0.561, 1.110]**이다. 두 통계는 계산 순서가 다르다. 0개 중심도 다섯 곳 평균에 포함했지만, 빈 중심에 전단 0을 부여하거나 이를 shear estimator로 사용하지 않았다. 개수 비율은 boost 보정이나 검출 유의도가 아니다.

| PDR3 부모 영상 필드 | 렌즈 | random 중심 | 0개 중심 | 5곳 모두 양수 |
|---|---:|---:|---:|---:|
| AUTUMN | 22 | 110 | 15 | 17 |
| HECTOMAP | 8 | 40 | 3 | 7 |
| SPRING | 84 | 420 | 27 | 72 |

기존 공식 PDR3 patch polygon에서 표적 114곳과 random 570곳 모두 부모 필드에 포함됐고 필드 간 모호성이 없었다. 모든 random은 짝지은 표적과 같은 부모 필드다. **이 분류는 S16A shape mask나 3 proper-Mpc 원뿔의 유효 면적이 아니다.** 같은 PDR3 필드·깊이·seeing 선정만으로 S16A 형상원 지지가 확보되지는 않았다.

0개 중심과 양수 중심의 PDR3 seeing 비율 중앙값은 각각 1.0002와 1.0000, 깊이 차이 중앙값은 −0.0214와 +0.0050 mag다. 이러한 기술 통계는 0개 결과의 원인을 특정하지 않는다. 마스크·S16A footprint 경계·source 선택 중 무엇이 영향을 줬는지는 총 개수만으로 분리하지 않았다. 같은 렌즈의 다섯 중심 및 겹치는 원뿔은 독립이 아니므로 p-value·공분산은 계산하지 않았다.

![원시 개수와 하늘 분포](/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_s16a_wl_random_support_audit_2026-11-17_01/raw_count_support.png)

## 미완전한 18개 그룹

| 렌즈 | ICL 필터 | 부모 필드 | 양수 random 수/5 | 표적 원시 개수 |
|---|---|---|---:|---:|
| S23B_10185 | i | AUTUMN | 4 | 1463 |
| S23B_10600 | i | AUTUMN | 0 | 4640 |
| S23B_2195 | i | AUTUMN | 2 | 1655 |
| S23B_2417 | i | AUTUMN | 0 | 1878 |
| S23B_2548 | i | AUTUMN | 4 | 5777 |
| S23B_3111 | i | SPRING | 2 | 3598 |
| S23B_5959 | i | SPRING | 4 | 4348 |
| S23B_6125 | i | SPRING | 3 | 4622 |
| S23B_6290 | i | SPRING | 0 | 1536 |
| S23B_8173 | i2 | SPRING | 4 | 1597 |
| S23B_8180 | i2 | SPRING | 4 | 1122 |
| S23B_8206 | i2 | SPRING | 2 | 1773 |
| S23B_8310 | i2 | SPRING | 4 | 1130 |
| S23B_8627 | i | SPRING | 4 | 4379 |
| S23B_9013 | i | SPRING | 2 | 1795 |
| S23B_9183 | i | SPRING | 4 | 3412 |
| S23B_9307 | i2 | SPRING | 0 | 1601 |
| S23B_9809 | i2 | HECTOMAP | 2 | 2617 |

모두 0인 그룹은 S23B_10600, S23B_2417, S23B_6290, S23B_9307이다. 일부가 0인 14개 그룹과 합쳐 18개다. 이 그룹을 삭제·대체하지 않았고, 기존의 query 전 제외 2개 렌즈도 별도 원장에 그대로 보존했다.

## S16A 보정 경로

공식 S16A 안내에는 S16A distortion, shape error/RMS, weight, multiplicative/additive bias가 제공되며, Y1 ensemble estimator와 source-redshift 관련 보정을 따라야 한다. Sirius mask 및 선택·responsivity 보정 조건도 해당 릴리스에 맞춰 확인해야 한다. [HSC S16A 공식 안내](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s16a-shape-catalog-pdr2/)

공식 Y1 calibration 소프트웨어는 보정에 사용한 선택과 가중치의 적용 범위를 명시한다. 이번에는 코드를 import하거나 실행하지 않았다. 제공된 catalog factor를 받을 경우에도 값과 선택의 적용성을 확인해야 한다. [HSC 연결 Y1 calibration 저장소](https://github.com/PrincetonUniversity/hsc-y1-shear-calib)

현지 `s19a_shape_stack.py`는 S19A/Y3의 선택 계수와 tomographic n(z) 경로를 사용한다. S16A 열 이름만 바꾸어 이 경로에 넣는 것은 검증되지 않았다. 보정 데이터와 source selection을 같은 릴리스에서 검증하는 어댑터가 필요하다. [HSC S19A 공식 안내](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s19a-shape-catalog-pdr3/)

S16A에 맞춘 COSMOS 재가중 보정 경로가 공식적으로 공개돼 있다. 그 자료는 S17A 측광을 이용하며 SOM weight와 shape weight, 실제 적용한 source cut의 일치를 확인해야 한다. 현 파이프라인에서의 선택 일치 검증은 아직 완료하지 않았고 자료를 새로 받지 않았다. [S17A Wide-COSMOS 안내](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s17a-wide-cosmos/)

따라서 전체 S16A 반출이나 기존 S19A와의 결합 측정에 앞서, 소규모 source-row QA로 문서의 열·좌표 단위·보정값·radial 지지·source 재사용부터 확인하도록 계획했다. 이 QA에도 photo-z·전체 선택/PSF 필드가 없어 완전한 보정이나 과학 측정을 검증할 수는 없다.

## 준비한 QA 파일럿 — 미제출

두 ICL 필터와 원시 지지 유형(5곳 양수/일부 양수/전부 0)의 각 층에서 표적 원시 개수·ID 순으로 정렬해 floor(n/4), floor(3n/4) 순위를 골랐다. 단일 원소 층은 한 곳만 골랐다. 이는 진단 목적의 결정적 11렌즈 파일럿이며 대표 표본이나 최종 과학 표본이 아니다.

| 렌즈 | 진단 층 | 부모 필드 | 양수 random 수/5 |
|---|---|---|---:|
| S23B_3012 | i:five_positive | SPRING | 5 |
| S23B_9131 | i:five_positive | SPRING | 5 |
| S23B_9013 | i:partial_positive | SPRING | 2 |
| S23B_8627 | i:partial_positive | SPRING | 4 |
| S23B_6290 | i:all_zero | SPRING | 0 |
| S23B_10600 | i:all_zero | AUTUMN | 0 |
| S23B_9828 | i2:five_positive | HECTOMAP | 5 |
| S23B_3620 | i2:five_positive | SPRING | 5 |
| S23B_8310 | i2:partial_positive | SPRING | 4 |
| S23B_8206 | i2:partial_positive | SPRING | 2 |
| S23B_9307 | i2:all_zero | SPRING | 0 |

- 전체 원장: 표적 11곳 + random 55곳 = 66곳. 이미 0으로 확인된 random 23곳은 보존하되 SQL에서는 제외한다.
- 요청 제안: S16A `pdr1`, 3 proper-Mpc, 양수 43곳(표적 11 + random 32), 순차 5 batch.
- 반환 12열: `center_id`, `object_id`, `source_ra`, `source_dec`, `e1`, `e2`, `sigma_e`, `rms_e`, `weight`, `m`, `c1`, `c2`.
- 예상 원뿔–source 행: **169,535**. 원뿔 간 중복을 제거한 고유 source 수가 아니다.
- 예상 CSV: **27,125,600–61,032,600 bytes (25.9–58.2 MiB)**. 추정치이며 수신 상한은 **134,217,728 bytes (128 MiB)**다.
- 영상·PSF·photo-z·PDF·추가 selection/PSF-moment 필드는 요청하지 않는다. 전단·질량·boost·ICL–WL 상관을 측정하지 않는다.
- 계획 SHA-256: `7af70c2be81f83e61bc6f74caf98b66a6309ae1b0caed0b9b7d1a32f42bd3f10`.
- 계획과 SQL을 원본 count/좌표/파일 해시에 독립 대조했다. 새 SQL의 실제 서버 스키마 검증은 아직 하지 않았다. 오류가 나면 중단하며 승인 범위를 자동 수정하지 않는다.
- **아직 제출하지 않았다. 사용자께서 정한 다운로드별 사전 승인 조건에 따라, 정확한 이 계획의 별도 승인을 기다린다.**

## 재현과 검증

- 분석: `audit_s23b_s16a_wl_random_support.py` (기본값 또는 `--output NEW_DIR`).
- 계획: `prepare_s23b_s16a_shape_qa_pilot.py` (SQL만 작성하며 제출 기능 없음).
- 테스트: 진단 13개, 파일럿 7개 통과. 0개 유지, 잘못된 count/좌표/슬롯/중복 거부, 구면 RA 래핑, 결정적 층화, 12열 제한을 확인했다.
- 분석 결과: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_s16a_wl_random_support_audit_2026-11-17_01/`.
- 분석 summary SHA-256: `56c463d81b7c099dec6380aa1feb84efacd2625d6fbf16e8c189eec7719dea79`.
- 미제출 계획: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_s16a_shape_qa_pilot_plan_2026-11-17_01/`.
- 독립 감사 코드·계획 검증·보정 준비도: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_s16a_wl_random_support_review_2026-11-17_01/`.
- 보정 준비도 SHA-256: `aa64d5d85178983d6457a9f5756e297a3c7e6b5538e28f7b1bee8a3c017da101`.
- 원래 source/영상 자료는 Git에 추가하지 않았다. 코드·인수인계 변경은 현재 미커밋 상태다.
