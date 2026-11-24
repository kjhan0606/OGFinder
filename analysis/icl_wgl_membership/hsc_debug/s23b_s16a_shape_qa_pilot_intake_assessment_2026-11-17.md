# S16A 파일럿 수신·입력 QA 평가 (2026-11-17)

승인된 11렌즈 파일럿의 수신과 로컬 QA를 완료했다. 수신 자료는 입력·구현 검증에 사용할 수 있다. Y1 전체 source 선택과 배경 source-redshift 보정, 관측영역·공분산이 검증되기 전에는 전단·질량 또는 ICL–WL 과학 측정으로 사용하지 않는다. `science_measurement=false`다.

## 승인 범위와 실제 수신

원래 승인 범위는 11개 표적과 짝지은 random 55곳의 원장이며, 그중 이전 count가 양수인 43곳만 요청했다. 이미 0으로 확인된 23곳은 새 SQL에서 제외하고 원장에 명시적으로 유지했다. 3 proper-Mpc, 12열, 순차 5배치와 128 MiB 누적 CSV 수신 상한을 유지했다.

| 항목 | 결과 |
|---|---:|
| 완료 배치 / 요청 중심 | 5 / 43 |
| 원장의 전체 중심 / 이전 count-only 0 중심 | 66 / 23 |
| 원뿔–source 행 | 169,535 |
| 고유 catalog object | 115,993 |
| 실제 CSV 수신량 | 29,274,595 bytes (27.9 MiB) |
| 수신 상한 | 134,217,728 bytes (128 MiB) |
| 유한 수치 열 + 양수 weight 행 | 169,535 / 169,535 |

첫 원본 요청은 `ishape_hsm_regauss_derived_weight`라는 열을 찾지 못해 HTTP 406으로 거절됐으며, CAS job 생성과 CSV 수신은 없었다. 앞선 계획에서 개요 문서의 이름을 사용했으나 실제 공개 스키마는 `ishape_hsm_regauss_derived_shape_weight`였다. 모든 요청 열·타입·degree 단위를 공개 schema-browser JSON에서 다시 대조했다. [HSC 실제 스키마](https://hsc-release.mtk.nao.ac.jp/schema/#pdr1.s16a_wide.weaklensing_hsm_regauss)

출력 alias `weight`의 원본 식만 정정한 v2를 보존했다. 각 SQL에서 이 치환을 되돌리면 원본 SQL과 정확히 같고, 중심·좌표·원뿔·join·반환 열·행 수·상한은 불변임을 검사했다. 이미 승인된 같은 의미의 자료 범위로 계속했으며, 새로운 해시에 대한 사용자 승인을 별도로 받았다고 기록하지 않았다. 원본 계획과 거절 기록은 그대로 남겼다.

- 승인된 원래 자료 범위 해시: `7af70c2be81f83e61bc6f74caf98b66a6309ae1b0caed0b9b7d1a32f42bd3f10`.
- 실제 실행한 식별자 정정 계획 해시: `1e3ac9451e519761e0f773915de79728d5443808efedc74155963a79c3ae4b30`.
- 수신 파일에는 source ID·좌표·shape와 weight/bias만 있다. 영상·PSF·photo-z·PDF·추가 photometry/selection/PSF-moment 필드는 수신하지 않았다.
- 인증 정보는 파일·체크포인트에 저장하지 않았다. 초기 스키마 거절 후 source CSV를 재수신한 시도는 없고, 확정 수신량·보수적 예산 차감량·완료 CSV 크기가 모두 일치한다.

## 정량 입력 검증

5개 CSV의 정확한 12열, 중심 ID, 각 중심의 이전 count와 행 수, 고유 `(center_id, object_id)` 쌍, 파일 SHA-256과 크기를 검증했다. Object ID는 64-bit 부동소수점으로 바꾸지 않았다. 실제 source 좌표 전부가 degree 좌표의 해당 원뿔 안에 들어가는지 확인했다(허용 수치 여유 0.01 arcsec). 알 수 없는 중심·중복 쌍·개수 불일치가 없었다.

공개 Y1 코드가 정의한 `weight × (sigma_e² + rms_e²) = 1` 관계를 169,535행에서 점검했다. 상대 잔차 중앙값은 4.00×10⁻⁸, 최댓값은 2.53×10⁻⁷이며 진단 허용오차 5×10⁻⁶ 초과는 0행이다. 이는 저장된 보정 열 사이의 산술 일관성으로, shear calibration의 물리적 정확도나 불확실성 한계가 아니다. [공개 Y1 calibration utility](https://github.com/PrincetonUniversity/hsc-y1-shear-calib/blob/main/utilities.py)

| 열 | 최소 | 중앙값 | 최대 |
|---|---:|---:|---:|
| sigma_e | 0.00206347 | 0.107927 | 0.632325 |
| rms_e | 0.329917 | 0.409356 | 0.453399 |
| weight | 1.80083 | 5.35982 | 9.18481 |
| m | -1.01742 | -0.0911513 | 0.0321612 |
| c1 | -0.015952 | -0.000131646 | 0.00732727 |
| c2 | -0.0125631 | -7.52769e-05 | 0.0103854 |

기본 유한값 조건을 통과해도 보정 적용 범위가 검증된 것은 아니다. `m`의 극단값을 보고 사후 QC 플래그를 추가했으며, 원본 자료를 자르거나 값을 수정하지 않았다.

| 사후 진단 플래그 | source-pair 행 | 고유 object |
|---|---:|---:|
| m ≤ −1 | 2 | 2 |
| distortion 크기 ≥ 2 | 0 | 0 |
| derived sigma_e > 0.4 | 547 | 427 |

이 문턱들은 전체 Y1 source cut이 아니다. 특히 `m≤−1`은 ensemble m 보정에서 개별 object를 자동 무효로 만드는 규칙으로 취급하지 않았다. 받은 `derived sigma_e`도 source-selection recipe의 raw HSM sigma와 다른 열이므로 0.4 문턱을 공식 raw-sigma 컷으로 대신 쓰지 않았다. SNR·resolution·photometry·flag가 없어 극단값의 발생 원인이나 최종 채택 여부를 결정할 수 없다.

## source 재사용과 반경별 지지

115,993개 고유 object 중 34,482개가 여러 중심에서 나타났으며, 한 object가 최대 5개 중심에 쓰였다. 재사용은 이 파일럿에서 같은 렌즈 그룹의 random 중심 사이에만 있었다. 서로 다른 렌즈 그룹 사이의 재사용과 lens–random 사이의 공유는 0이었다. 같은 object의 받은 열 값이 서로 다른 경우도 0이었다. 이러한 결과가 더 큰 표본의 독립성이나 오염 부재를 증명하지는 않는다. 원뿔 간 source 재사용을 향후 공분산에 반영해야 한다.

아래는 표적과 random 5곳 각각에 적어도 한 행이 있는 그룹 수다. 유한 보정값·양수 weight 조건에서도 수는 같았다. 면적이나 유효 background-source 선택을 검증한 결과는 아니다.

| proper 반경 [Mpc] | 공통 원시 지지 그룹 / 11 | 수치 조건 지지 그룹 / 11 |
|---|---:|---:|
| 0–0.1 | 2 | 2 |
| 0.1–0.2 | 3 | 3 |
| 0.2–0.4 | 3 | 3 |
| 0.4–0.7 | 3 | 3 |
| 0.7–1 | 4 | 4 |
| 1–1.5 | 4 | 4 |
| 1.5–2 | 4 | 4 |
| 2–3 | 4 | 4 |

동일한 그룹을 모든 반경에 공통으로 요구하면 0.1–3 Mpc에서 S23B_3012, S23B_3620, S23B_9131의 **3개**다. 중앙 0–0.1 Mpc까지 포함하면 S23B_3012와 S23B_9131의 2개다. 이는 새 과학 표본 선정이 아니다. 이 파일럿은 5곳 양수·일부 양수·모두 0인 원시 지지 층을 의도적으로 골랐으므로 3/11을 전체 114렌즈나 1,043표적의 유지율로 해석하지 않는다. 빈 annulus의 원인도 총 개수와 이 12열만으로 Poisson·mask·footprint 중 하나로 지정하지 않았다.

![S16A 파일럿 가중치와 공통 반경 지지](/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_s16a_shape_qa_pilot_intake_2026-11-17_01/pilot_intake_qa.png)

## 현재 사용 범위와 다음 관문

입력 무결성·좌표·기본 수치 일관성은 확인했다. 이 자료로 S16A 전용 intake/선택 어댑터를 개발하고 source 재사용·지지를 진단할 수 있다. 완전한 Y1 선택·selection bias·photo-z/배경 source 보정·관측 영역과 공분산이 남아 있어 전단·DeltaSigma·질량·boost·ICL–WL 상관을 측정하지 않았다. 기존 S19A/Y3 컷·계수·tomographic n(z)도 적용하거나 결합하지 않았다.

다음 제안은 기존 자료와 공식 스키마로 **Y1 선택·보정 어댑터의 입력 규약 및 최소 추가 source-quality/photo-z 반출 계획**을 설계하는 것이다. 새 데이터 수신은 정확한 계획별 별도 승인 대상이다. 요청된 단계 경계에 따라 현재는 승인 대기하며 후속 분석·표본 확대·조회는 시작하지 않는다.

## 재현 기록

- 수신 코드: `submit_s23b_s16a_shape_qa_pilot.py`. 현재 기본 계획은 정정된 `_02`이며, 이미 완료된 요청을 다시 실행하면 인증이나 네트워크 없이 로컬 결과만 검증한다. 이를 실제 경로에서 금지된 network/credential stub으로 확인했다.
- 입력 QA: `audit_s23b_s16a_shape_qa_pilot.py` (`--output NEW_DIR`; 덮어쓰기 거부).
- 사후 tail QC: `audit_s23b_s16a_shape_qa_tails.py` (원본 불변, 사후 플래그만 기록).
- 테스트: 수신 20개, 입력 QA 7개, tail 1개 통과. 실제 수신 5배치와 전체 QA·독립 CSV 대조도 완료했다.
- 원본 CSV·수신 체크포인트: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_s16a_shape_qa_pilot_plan_2026-11-17_02/results/`.
- QA 표·PNG/PDF·요약: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_s16a_shape_qa_pilot_intake_2026-11-17_01/`.
- 실행 코드 스냅샷·독립 검증·tail QC·실행 환경: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_s16a_shape_qa_pilot_evidence_2026-11-17_01/`.
- 실제 수신 당시의 코드와 완료 후 기본 경로를 수정한 코드는 서로 구분해 보존했다. QA 및 수신 원본 CSV는 변경하지 않았다.
- 수신 요약 SHA-256: `373b1bca841d733bd4269421dd00fb2b6afb89295cc2f39dbab3d8424bd45c6c`.
- QA 요약 SHA-256: `fd04520cc46254c09c2dc3e011e1e376f8f0169d33412057edc635db119fa40d`.
- Tail 요약 SHA-256: `438d885bacffc24341d88c6cb5aebf1f672347dea8bbdf19d336eb6d5255d907`.
- 독립 대조 SHA-256: `f369a9edacad0e6b48f011c97efd3ebe62c54443d72477cdacf2591bfcf0f3e0`.
- 자료는 `/gpfs`에만 있고 Git에 추가하지 않았다. 이번 코드·기록 변경은 미커밋 상태다.
