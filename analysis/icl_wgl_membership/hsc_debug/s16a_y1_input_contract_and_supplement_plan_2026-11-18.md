# S16A/Y1 입력 규약과 최소 보충 계획 — 2026-11-18

사용자가 승인한 로컬 규약·계획 작성 단계를 완료했다. 기존 파일럿 169,535개 source–center 행의 115,993개 고유 object에 필요한 입력만 보충하는 계획이다. **새 CAS 조회·수신은 아직 실행하지 않았다.** 공개 문서와 schema-browser JSON만 확인했다. 다음 데이터 수신은 아래 정확한 계획에 대한 별도 승인 대상이다.

## 구현한 규약

`s16a_y1_source_contract.py`는 [Y1 논문 Table 4](https://arxiv.org/html/1705.06745v3)의 object 조건을 진단한다. 품질 플래그, i-band S/N ≥ 10, resolution ≥ 0.3, distortion 크기 < 2, raw sigma 범위 [0, 0.4], extinction 보정 i ≤ 24.5, blendedness < 10^−0.375 및 grzy 중 두 밴드 이상의 S/N ≥ 5를 사용한다. 유한하지 않거나 없는 입력은 통과시키지 않는다. 네 밴드 모두를 요구하지 않는다. 원래 pipeline의 cosmic-ray/bad 플래그 전파 결함도 규약에 기록했다.

이는 object 조건 재현이며, 관측 영역·별 마스크·배경 소스 보정까지 끝난 최종 과학 선택이 아니다. [공식 공개 안내](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s16a-shape-catalog-pdr2/)는 weaklensing 표가 weak-lensing flag와 Sirius mask를 적용한 부분집합이라고 설명한다. 이 사실만으로 random 중심의 사용 가능 면적이나 마스크 지도가 제공되지는 않는다.

기존 `sigma_e`는 calibration으로 얻은 오차다. 0.4 source 컷에는 별도 **raw `ishape_hsm_regauss_sigma`**가 필요하다. 기존 derived sigma_e > 0.4 또는 m ≤ −1 tail을 사후에 자동 삭제하지 않았다. 가중치와 m/c calibration은 함께 유지하고, 향후 전단은 ensemble m·responsivity·회전된 additive 성분을 검증해 사용해야 한다. 개별 source마다 1+m으로 나누는 방식은 규약에서 금지했다.

기존 12열 중 이 object 진단에 직접 필요한 입력은 e1/e2 두 개뿐이다. 빠진 28개 입력을 0으로 채우지 않아 **169,535행 모두 현재 전체 object 판정은 unknown**이다. 이는 소스가 전부 불량이라는 뜻이나 과학 선택의 유지율이 아니다.

## 정확한 추가 수신 제안

| 항목 | 제안 범위 |
|---|---|
| Release / tables | `pdr1`, `s16a_wide.meas`, `meas2`, `forced`, `photoz_ephor_ab` |
| 대상 | 기존 파일럿의 고유 object 115,993개만 |
| 새 중심 / 새 object | 0 / 0 |
| 반환 | object_id + 품질·측광 입력 28개 + photo-z 2개 = 31열 |
| 요청 / 예상 반환 행 | 순차 SQL 24개 / 115,993행 |
| 배치 | 최대 5,000 ID, 마지막 993 ID |
| SQL 파일 총크기 | 3,632,407 bytes |
| 예상 비압축 CSV | 14,895,104–59,436,416 bytes, 약 14.2–56.7 MiB |
| 누적 CSV 상한 | 134,217,728 bytes, 128 MiB; 재시도·중단 수신분 포함 |
| 계획 SHA-256 | `44b9d4df4904dac295f08a426abb20aa536ae85fc2c91239ac599f39deb412a1` |

용량은 object당 128–512 bytes 및 배치 header 여유로 계산한 추정이다. 실제 수신량은 아직 없다. 기존 좌표·shape·weight·m/c, 영상·PSF·PSF moment, 전체 photo-z PDF, COSMOS calibration 파일과 S19A/Y3 자료는 요청하지 않는다.

추가 열의 정확한 매핑은 계획의 `column_expressions`와 `schema_column_check.json`에 있다. 중요 매핑은 다음과 같다.

| 로컬 입력 | 공개 PDR1 식별자 |
|---|---|
| is_primary | `meas.idetect_is_primary` — band별 unforced flag |
| centroid_sdss_failed | `meas2.icentroid_sdss_flags` |
| raw_sigma / resolution | `meas2.ishape_hsm_regauss_sigma`, `ishape_hsm_regauss_resolution` |
| i_flux / error / mag | `meas.icmodel_flux`, `icmodel_flux_err`, `icmodel_mag` |
| grzy flux / error | `meas.[grzy]cmodel_flux`, `[grzy]cmodel_flux_err` |
| a_i | `forced.a_i` |
| photo-z / risk | `photoz_ephor_ab.photoz_best`, `photoz_risk_best` |

SQL은 exact signed-bigint ID의 `VALUES` 목록과 object_id 기준 `LEFT JOIN`만 쓴다. 선택 컷이나 `DISTINCT`, `COALESCE`로 행·결측값을 숨기지 않는다. 미대응 photometry/photo-z는 NULL로 보존한다. 반환 시 요청 ID당 정확히 한 행이어야 하며 중복·누락·추가 ID 또는 열 불일치는 중단 사유다. ID는 부동소수점으로 바꾸지 않는다.

기존 center–object 쌍 169,535행은 별도 원장으로 보존했다. 수신 뒤 `merge_supplement`가 기존 순서·값·중복 source의 여러 중심 소속을 유지해 결합한다. 원본 열 덮어쓰기와 중복 center–object 쌍은 거부한다. 이미 알려진 23개 zero 중심도 기존 66중심 원장에 유지하며 새로 조회하지 않는다.

## Photo-z와 과학 측정 관문

Ephor_AB 선택은 [S16A에 맞춰 재가중한 S17A Wide-COSMOS 안내](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s17a-wide-cosmos/)에 Ephor_ab 보정 경로가 있다는 점에 근거한다. 알고리즘 정확성이나 calibration 버전 일치를 검증했다는 뜻은 아니다. 이번 대표 photo-z/risk는 결측·분포 점검용이며, 임의의 z_lens+0.1/risk 컷을 미리 적용하지 않았다.

실제 배경 효율·dilution 보정에는 동일 방법과 동일 선택을 쓴 COSMOS 보정 표본, SOM×shape weight 등의 대응 검증이 필요하다. 이 보충 수신만으로 해결되지 않는다. 기존 S19A/Y3 n(z)·selection-bias 계수를 S16A에 옮기지 않는다. Y1 selection bias·PSF null test, redshift-dependent responsivity, usable area, 같은 컷 이후 반경 지지, source 재사용·공간 공분산도 별도 관문이다. **전단·질량·boost·ICL–WL 측정은 하지 않았으며 science_measurement=false다.**

## 검증과 재현

- 신규 29개 테스트: 컷 경계·불리언 문자열·결측·두 밴드 조건·raw/derived sigma 구분·m tail 비삭제·큰 ID·정확한 원장 결합·SQL 범위 등을 검증했다.
- 기존 support/pilot/receiver/intake/tail 테스트를 포함해 **77개 통과** (`/home/kjhan/miniconda3/bin/python3 -m pytest`). 시스템 `/usr/bin/python3`에는 pytest가 없어 기존 연구용 Miniconda 환경을 사용했다.
- 별도 stdlib 코드가 원본 CSV를 다시 읽어 모든 ID·쌍·SQL VALUES 목록·projection·join·31열 공개 schema 대응·파일 해시를 대조했다. 생성기 SQL 함수를 호출해 비교한 결과와 구분해 보존했다.
- 실제 계획 재검증 때 network/credential 접근을 금지한 hook을 적용했고 시도는 0회였다. 기존 출력 덮어쓰기와 SQL 변조 거부도 확인했다.
- 서버 SQL 호환성은 아직 미검증이다. 실행 승인을 받은 뒤 첫 배치에서 확인하며, 오류·불명확한 제출 상태이면 중단한다. 현재는 실행기/수신기를 추가로 실행하지 않았다.

계획 생성과 검증은 로컬 파일만 읽는다. 이미 생성된 계획은 다음으로 확인한다.

```bash
cd /home/kjhan/BACKUP/Astrafek/OGFinder
/home/kjhan/miniconda3/bin/python3 analysis/icl_wgl_membership/hsc_debug/prepare_s23b_s16a_quality_supplement.py --verify
```

| 산출물 | 절대 경로 |
|---|---|
| 계획·SQL·ID/쌍 원장·입력 규약·schema 근거 | `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_s16a_quality_supplement_plan_2026-11-18_01/` |
| 공개 schema 원본 | `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s16a_y1_input_contract_sources_2026-11-18_01/` |
| 별도 감사·금지 hook·테스트/실행 기록·코드 스냅샷 | `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s16a_quality_supplement_review_2026-11-18_01/` |

원본 파일럿·앞선 assessment·과거 수치·해시와 타 작업 변경을 보존했다. Git 데이터 반입이나 이번 단계의 커밋·푸시는 하지 않았다. **다음 승인 대상은 위 해시의 추가 조회·수신과 수신 후 로컬 진단이다.** 승인되면 먼저 상한·재개·ID/열 검증을 갖춘 순차 수신기를 테스트하고, 이 계획 범위만 제출한다. 과학 측정이나 표본 확대 승인을 포함하지 않는다.
