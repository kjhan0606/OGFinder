# ICL–위성 고리 제외: 확장 PSF 혼합 민감도 — 2026-11-03

## 질문과 분석 범위

확장 HSC-I 공개 PSF 날개를 위치별 compact PSF 코어에 이어 붙였을 때, 관측 segmentation에서
지원되는 위성 50–75 kpc 고리를 ICL 200–300 kpc 측정에서 제외하는 조치가 합성광 회수에
영향을 주는지 평가했다. 이는 앞서 수행한 local-core-only 공유 고리 주입검증의 wing 민감도
확장이다. 고정된 HSC-I 공개 PSF 날개는 대상별·위치별·시기별 seeing과 일치하지 않는다.
특히 i2 결과는 공개 HSC-I를 대용한 탐색 민감도다.

분석은 고정 74개 ICL–위성 연관쌍 중 compact PSF core 위치가 있는 48개 target/filter 쌍만
사용했다(12 `i`, 36 `i2`). 매 배치에 0 오류가 있었고 프레임 JSON의 SHA-256을 독립 재검증했다.
48개는 공개 74쌍 전체를 대표하는 완전 표본이 아니다. 고정 association table은 74행이며
SHA-256 `483f81540543fc46ec1cf4333bbd3e238799b493e48fe03235da3038f5714917`이다.

혼합 PSF 매개변수는 splice 반경 2.4 arcsec, blend 폭 0.6 arcsec, 날개 진폭 0.65다. 공개 FITS는
승인된 CC-BY-4.0 HSC-I `coadd/bg` 참조(`HSC_PSF_i.fits`, SHA-256
`8897af89767e5fbb8d5f9c70965126bc89409f6c65ecd8d117c7995f34be2d85`)이며 그 자체도 승인된
대상별 PSF가 아니다. diffuse, 위성 envelope, 동시 주입의 응답을 같은 저장 마스크로 재측정했고,
지원된 50–75 kpc 위성 target-pixel 합집합을 제외한 결과와 짝지었다. 별도로 단위 점원 모델이
저장된 ICL 마스크의 200–300 kpc 환형영역에 남기는 wing flux를 기록했다. 관측 FITS는 바뀌지
않았다.

수치는 프레임 간 q16–q84와 중앙값이다. 이 범위는 표본 간 산포 요약이며 모집단 신뢰구간이 아니다.
응답차는 각 프레임에서 고리 제외−기존 고정 마스크를 계산한 짝지은 차다.

## 결과

| 필터·합성 주입 | 프레임 | 기존 고정 마스크 응답 q16 / 중앙값 / q84 | 위성 고리 제외 응답 q16 / 중앙값 / q84 | 짝지은 변화 q16 / 중앙값 / q84 |
|---|---:|---:|---:|---:|
| i · diffuse → ICL | 12 | 0.22268 / 0.22538 / 0.22745 | 0.22275 / 0.22472 / 0.22787 | −0.00167 / −0.00018 / +0.00102 |
| i · 위성 envelope → ICL | 12 | 0.15678 / 0.46817 / 0.77260 | 0.14941 / 0.46687 / 0.77849 | −0.01379 / −0.00022 / +0.00274 |
| i2 · diffuse → ICL | 36 | 0.22255 / 0.22568 / 0.22927 | 0.22200 / 0.22602 / 0.22929 | −0.00110 / 0.00000 / +0.00132 |
| i2 · 위성 envelope → ICL | 36 | 0.29497 / 0.68183 / 0.91962 | 0.28831 / 0.65932 / 0.91905 | −0.06592 / 0.00000 / +0.01161 |

점원 1 nJy 입력에서 저장된 마스크를 통과해 ICL 200–300 kpc에 남은 모델 wing flux fraction은
중앙값 기준 i `0.000994 → 0.000902`, i2 `0.000970 → 0.000896`으로 작아졌다. 동시에 해당
환형영역의 유효 마스크 coverage 중앙값은 i `0.681 → 0.618`, i2 `0.652 → 0.611`으로 감소했다
(각각 2.46, 3.70 percentage points). 합계 143개의 지원 위성 고리에서 ICL 영역과 겹친
유효 target 픽셀 320,748개를 제외했다. 따라서 작은 모델 wing 감소와 함께 측정 지지 손실도
보고해야 한다.

고정 마스크 기준 diffuse 응답 변화의 q16–q84 범위는 두 대역 모두 0을 포함한다. 위성 envelope
응답 변화도 q16–q84가 0을 포함하고 중앙 변화는 i `−0.00022`, i2 `0.00000`이다. i2의 paired
최솟값은 `−0.266`으로 프레임별 영향은 균일하지 않다. 이 결과는 검토한 표본과 합성 조건에서
고리 제외가 주입 회수율의 일정한 보정량을 제공한다는 증거가 아니다. 실제 측정에 나눗셈 보정이나
응답 보정은 적용하지 않는다.

## 해석 및 한계

- 이 결과는 synthetic light transfer 민감도이며 관측 ICL 검출이나 weak-lensing 측정이 아니다.
- 공개 wing은 위치·epoch·seeing 불일치가 남는다. i2에는 필터별 공개 wing이 없어 HSC-I proxy를
  사용했으므로 결과를 i2 PSF calibration으로 해석할 수 없다.
- 선택된 48프레임만 실행했다. association 표본의 나머지 26쌍은 이 position-core 지원 교집합에 없다.
- 위성 고리 제외 영역은 기존 segmentation과 saved profile support에서 만들고 주입 간 고정했다.
  자동 검출 재실행, 실제 영상 마스크 변경, 마스크 선택의 불확실성을 보정하지 않는다.
- source Poisson noise, coadd 상관잡음의 전체 estimator covariance, matched-control 차감 분산,
  상류 sky-subtraction 전달, 별도의 확장 PSF 공간블록 불확실성은 포함하지 않는다.
- 따라서 `science_measurement=false`, `psf_approved=false`, `observational_validation_passed=false`,
  `flux_covariance_available=false`를 유지한다.

## 재현과 산출물

실행기는 [`extended_psf_injection.py`](extended_psf_injection.py)이며, 독립 해시·표본 소속·guard
flag·고정 마스크 짝지은 응답 집계는 [`summarize_icl_hybrid_shell_exclusion.py`](summarize_icl_hybrid_shell_exclusion.py)로
수행한다. 둘 다 lageunha의 이미 내려받은 로컬 파일만 사용한다. HSC CAS 요청이나 데이터 다운로드는
없었다.

- i batch: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_hybrid_shell_exclusion_association_i12_2026-11-03_01/summary.json`
  SHA-256 `ee3c1ae31b6ae2cef099702f071899e6b54a652935d8456f22893b37f476ec9e`
- i2 batch: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_hybrid_shell_exclusion_association_i2_36_2026-11-03_01/summary.json`
  SHA-256 `56e628cfcf39fbdf80d499376d158b0056bbafeab7de8a648c9e0b327b7ce7b3`
- 독립 집계기는 두 summary에서 각각 12/12, 36/36 frame JSON checksum과 74행 association membership을 확인했으며
  48개의 유일 target/filter key가 일치했다.

```bash
python analysis/icl_wgl_membership/hsc_debug/summarize_icl_hybrid_shell_exclusion.py \
  --association-table /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_outer_sky_association_reblocked_2026-11-03_01/association_rows.tsv \
  --i-summary /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_hybrid_shell_exclusion_association_i12_2026-11-03_01/summary.json \
  --i2-summary /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_hybrid_shell_exclusion_association_i2_36_2026-11-03_01/summary.json
```
