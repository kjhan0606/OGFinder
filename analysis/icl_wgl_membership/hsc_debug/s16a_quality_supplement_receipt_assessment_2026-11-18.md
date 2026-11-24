# S16A 보충 자료 수신·로컬 QA — 2026-11-18

승인된 추가 조회·수신과 로컬 진단을 완료했다. **입력·object 선택 진단에는 사용할 수 있다. 보정된 과학용 WL 측정으로는 아직 사용하지 않는다.** Photo-z 결측이 공간적으로 불균일하고, Y1 보정·사용 가능 면적·공분산 관문이 남아 있다. `science_measurement=false`다.

## 승인 범위와 수신 무결성

사용자는 계획 `44b9d4df4904dac295f08a426abb20aa536ae85fc2c91239ac599f39deb412a1`의 조회·수신과 수신 후 로컬 진단을 명시적으로 승인했다. 기존 S16A 파일럿 object ID만 사용했고 계획·SQL 변경 없이 첫 배치부터 성공했다. LagEunha에서 CPU 수동 프로세스로 실행했으며 Slurm/GPU는 사용하지 않았다.

| 항목 | 결과 |
|---|---:|
| 완료 SQL/CSV | 24 / 24 |
| 고유 object / 반환 열 | 115,993 / 31 |
| CSV bytes | 45,335,621 (약 43.2 MiB) |
| 확정 수신량 / 예산 차감량 | 둘 다 45,335,621 bytes |
| 누적 CSV 상한 | 134,217,728 bytes (128 MiB) |
| 결합한 원래 source–center 쌍 | 169,535 |
| 유지한 중심 / 이전 count-only zero 중심 | 66 / 23 |

모든 배치에서 정확한 31열, 요청 ID당 한 행, 누락·추가·중복 ID 부재와 파일 해시를 확인했다. 큰 정수 ID를 float로 바꾸지 않았다. NULL은 보존했고 원래 source–center 원장의 순서·멤버십·값도 유지했다. 수신 재시도·실패·스키마 정정은 없었다. 인증은 숨김 입력만 사용했고 저장하지 않았다. 완료된 제출 명령을 다시 실행했을 때도 network/credential 금지 hook 아래 로컬 검증만 수행했다.

새 중심·object·좌표·shape/weight/m/c 재수신, 영상·PSF·전체 PDF·COSMOS calibration 파일·S19A/Y3 자료는 이번 범위에 없다. [계획과 입력 규약](s16a_y1_input_contract_and_supplement_plan_2026-11-18.md)을 그대로 유지했다.

## Object 조건 진단

다음은 Y1 Table-4 object 조건의 로컬 재현 결과이며, 보정된 배경 소스 선택이 아니다. Photo-z/risk를 선택 컷으로 쓰지 않았다.

| 판정 | 고유 object | source–center 쌍 |
|---|---:|---:|
| 통과 | 115,978 (99.987%) | 169,518 |
| 미통과 | 4 | 4 |
| 판정 불가 | 11 | 13 |
| 전체 | 115,993 | 169,535 |

차이는 전부 **grzy 중 최소 두 밴드의 S/N ≥ 5** 조건에서 발생했다. 다른 품질 플래그, i S/N, resolution, distortion, raw HSM sigma, extendedness, 보정 i 등급과 blendedness 조건은 모두 통과했다. 4개 미통과 object의 실제 밴드 S/N과 11개 판정 불가의 입력은 별도 표에 보존했다. 문턱을 바꾸거나 결측을 채워 통과시키지 않았다. 이 작은 차이만으로 공개 카탈로그 오류를 단정하지 않으며, upstream 선택 provenance는 추후 규약 검토 대상이다.

기존 tail QC는 기본 object 컷과 구분해야 한다.

| 이전 사후 tail | 고유 object | object 조건 통과 / 미통과 |
|---|---:|---:|
| m ≤ −1 | 2 | 1 / 1 |
| derived sigma_e > 0.4 | 427 | 427 / 0 |
| 위 두 조건의 합집합 | 429 | 428 / 1 |

Derived sigma_e는 raw HSM sigma와 다른 양이다. 실제로 위 427개는 raw-sigma 조건도 통과한다. 두 m tail은 resolution 약 0.9823/0.9814이고, 그중 하나는 모든 object 조건을 통과했다. 따라서 기본 컷만으로 calibration 극단값 문제가 해소됐다고 주장하지 않는다. m≤−1을 임의로 삭제하거나 source별 1+m으로 나누지도 않았다. Ensemble 보정과 calibration 적용 범위는 별도 검증 대상이다.

## Point photo-z 완비성

전체 115,993개 중 **109,321개(94.2479%)**에 유한한 Ephor_AB best/risk가 있고, **6,672개**는 없다. Source–center 쌍 기준 결측은 8,483행이다. Object 조건 통과 표본에서도 6,671개에 point-z/risk가 없다. 유한 값에서 best-z 중앙값은 0.82854, risk 중앙값은 0.15164다. 음수 best-z 또는 [0,1] 밖의 유한 risk는 없었다. 이 분포는 n(z) 보정이나 배경 효율을 검증한 결과가 아니다.

| ICL 영상 분류 | 역할 | 고유 object | 유한 best/risk 비율 |
|---|---|---:|---:|
| i | lens | 23,203 | 99.35% |
| i | random | 38,700 | 99.68% |
| i2 | lens | 13,826 | 87.95% |
| i2 | random | 40,264 | 88.25% |

분류는 원장의 **ICL 영상 filter label**이며, S16A shape/photo-z의 원인 변수나 독립된 통계 실험으로 해석하지 않는다. Source가 있는 43개 중심의 입력 완비율은 **47.79–100%**다. 최저는 S23B_3620에 짝지은 `RC-cbf934142d614f978176`: 5,509행 중 2,633행에만 best/risk가 있다. 전체 평균만 보고 photo-z 없는 행을 임의로 버리면 공간 선택 함수가 달라질 수 있다. 결측의 원인이나 이를 보정하는 방법을 이번 진단만으로 정하지 않았다. 23개 zero 중심의 완비율은 0%로 채우지 않고 정의 불가로 보존했다.

## 반경 지지와 source 공유

Object 컷과 유한 calibration·양수 shape weight를 요구한 뒤에도, 표적과 고정 random 5곳 모두에 최소 한 source가 있는 그룹 수는 이전 raw 지지와 같다.

| proper 반경 [Mpc] | raw 공통 그룹 | object 컷 이후 공통 그룹 |
|---|---:|---:|
| 0–0.1 | 2 | 2 |
| 0.1–0.2 | 3 | 3 |
| 0.2–0.4 | 3 | 3 |
| 0.4–0.7 | 3 | 3 |
| 0.7–1 | 4 | 4 |
| 1–1.5 | 4 | 4 |
| 1.5–2 | 4 | 4 |
| 2–3 | 4 | 4 |

모든 0.1–3 Mpc bin에 공통인 그룹은 S23B_3012, S23B_3620, S23B_9131의 3개다. 중앙 bin까지 포함하면 3012/9131 두 개이며, 0.7–3 Mpc는 S23B_9828을 더한 네 개다. 이는 양수·일부 양수·전부 zero 층을 의도적으로 뽑은 11렌즈 파일럿의 진단으로, 전체 후보 표본의 유지율이 아니다. Known-zero 중심에 zero shear를 부여하지 않았다.

Object 조건 통과 표본에서도 34,480개 object가 여러 중심에 쓰이고, 고유 object보다 쌍이 53,540개 더 많다. 한 object가 최대 5중심에 나타난다. 원래 파일럿에서 확인한 random 간 source 공유를 무시한 독립 표본 공분산은 정당화되지 않는다. 보고한 유효 수는 catalog shape weight의 기술 통계이며 lensing 효율을 반영한 수가 아니다.

![S16A object 선택과 공통 반경 지지](/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s16a_quality_supplement_qa_2026-11-18_01/selection_qa.png)

## 사용 판단과 다음 단계

입력 무결성, object 조건, tail과 반경 지지 및 point-z 결측 진단에는 사용할 수 있다. 기존 숫자만 보고 미확인 컷을 대신 적용하던 위험을 줄였지만, 과학용 source-z/희석 보정·Y1 selection bias/PSF null test·redshift responsivity·S16A usable area/Sirius-mask geometry·공분산을 해결하지는 않았다. Source 수신 성공이나 99.987% object 통과율을 전단·질량·ICL–WL 결과의 과학 검증으로 해석하지 않는다.

다음 제안은 기존 자료와 공개 문서로 **공간별 photo-z 입력 결측, S16A 영역/마스크, Y1 calibration 적용 규약**을 검토하는 것이다. 추가 데이터가 필요하면 정확한 범위를 계획해 별도 승인받는다. 현 단계는 종료하고 다음 단계 승인을 기다린다. 추가 조회·표본 확대·S16A/S19A 결합이나 과학 측정은 자동으로 시작하지 않는다.

## 재현·검증 기록

- 수신기: `submit_s23b_s16a_quality_supplement.py`; 실제 사용 코드 SHA-256 `df04b1861ee0a573f2ab447a39201e8b004b103fdda7f47235c66db82316dd11`.
- QA: `audit_s23b_s16a_quality_supplement.py`; 실행 코드 SHA-256 `fdbd5993ac122b1c7960c562bd9fa314fd5b4c3d626311eeac4b524381bd86ec`.
- 새 receiver 13개와 QA 5개를 포함해 관련 **95개 테스트 통과**. 실제 수신·전체 QA·별도 검산도 통과했다. 그림을 직접 확인했고 제목·범례는 패널 내부에 있다.
- 별도 stdlib 검산은 원본 CSV/ID/쌍/선택 이유·tail·point-z 통계·528개 반경 프로파일·유효 수를 대조했다. 반경 거리에는 Astropy 대신 독립 Simpson 적분을 사용했다.
- 수신 요약 SHA-256: `bb9ee55540a3f7acdafff1f6778f190ea3b9807a071df9c821d5ff74ac92c27d`.
- QA 요약 SHA-256: `d7e474ea88260e6ae5bc1f70c086ae0b4b7753a9d0224ce5e746cca540210474`.
- 공간별 photo-z 진단 요약 SHA-256: `ec15fe062c85a15b4166640d063dc0ae02fa8a01f620d50bb7b17e08cd7c654f`.

| 산출물 | 절대 경로 |
|---|---|
| 24 CSV·체크포인트·수신 요약 | `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_s16a_quality_supplement_plan_2026-11-18_01/results/` |
| 선택·tail·반경·photo-z QA·PNG/PDF | `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s16a_quality_supplement_qa_2026-11-18_01/` |
| 승인 범위·코드 스냅샷·독립 검산·공간별 완비성·테스트/환경 기록 | `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s16a_quality_supplement_evidence_2026-11-18_01/` |

원본 및 앞선 기록을 보존했다. 데이터는 Git에 넣지 않았고 이번 단계의 커밋·푸시는 하지 않았다. 활성 수신·QA 프로세스는 정상 종료했다.
