# ICL–위성 표본 공분산 범위와 공유 고리 픽셀 민감도 — 2026-11-03

## 지원표 범위 정정

이전 63쌍 공분산 실행은 ICL 지원 입력으로 `prestack_flux_profiles`의 공통-stack 표를
사용했다. 그 표는 137쌍만 담아 annulus별 지원 분석에 필요한 전체 표본이 아니었다.
`icl_annulus_support`의 200–300 kpc pair table은 coverage 40% 이상인 224개 표적/필터 쌍을
담는다. 지원 선별 코드는 이제 이 annulus별 table을 기본값으로 사용한다. 기존 63쌍
산출물은 덮어쓰지 않았고, 전체 표본 결과는 새 디렉터리에 저장했다.

위성 50–75 kpc 표적·대조 고리 coverage 각각 60% 이상, 위성당 허용 방위각 대조점 3개 이상,
ICL 고리 coverage 40% 이상을 다시 적용하면 94개 표적/필터 쌍과 235개 위성 프로파일이
지원된다(26 `i`, 68 `i2`). 고정된 74개 ICL–위성 연관쌍(21 `i`, 53 `i2`)은 모두 이 표본에
포함되고 누락은 없다. 지원 94쌍 중 나머지 20쌍은 해당 연관표본 밖이다.

새 독립 범위 검증기는 ICL support와 위성 프로파일의 조건을 별도로 다시 적용해 두 공분산
산출물의 target/filter 쌍 집합이 정확히 94개인지 확인했다. 이어 표적+matched-control
공분산 공식 및 입력계보 검증기도 94행 모두 통과했다. HSC 재접속·쿼리·다운로드는 없다.

| Filter | Pair 수 | 허용 위성 프로파일 | 대조/표적 ICL lag 분산 중앙비 | 전체 ICL sigma/대각-only | ICL–위성 픽셀 잡음 상관 중앙값 [16–84 백분위] |
|---|---:|---:|---:|---:|---:|
| i | 26 | 47 | 0.979 | 1.248 | −0.00703 [−0.01785, +0.00537] |
| i2 | 68 | 188 | 0.995 | 1.243 | 약 0 [−0.01462, +0.00730] |

9×9 유한 lag 지지에서 표적·대조 ICL 고리는 모두 별도 영상의 비중첩 위치에 있다. 이
결과는 단거리 픽셀 잡음만 다루며 장파장 sky 공분산, 실제 위치 의존 공분산, source
Poisson, PSF·마스크·photometric 오차를 포함하지 않는다.

## 동일 픽셀 사용과 ICL–위성 관계의 민감도

94쌍의 저장된 segmentation/mask를 실제 픽셀 단위로 대조했다. 고정된 연관표본 74쌍 안에는
206개 허용 위성 프로파일(41 `i`, 165 `i2`)이 있다. 이 가운데 ICL 200–300 kpc 고리와
위성 50–75 kpc 표적 추정량이 공유하는 유효 측정 픽셀은 140개 프로파일, 즉 68%에서
1개 이상이다(`i`: 29/41, `i2`: 111/165). 공유 비율 중앙값은 위성 유효 픽셀 대비 약
0.34다. 이는 노이즈 공분산이 0에 가깝다는 결과와 모순되지 않는다. 위성 추정량이
target에서 방위각 control을 빼기 때문에 signed estimator 노이즈가 상쇄될 수 있지만,
표적 위성 고리의 실제 광량은 ICL 고리 측정에도 들어갈 수 있다.

이를 직접 점검하기 위해 각 연관 은하단의 지원된 위성 50–75 kpc **표적 고리에서 실제
사용된 픽셀의 합집합**을 만들고, 같은 표적 ICL 200–300 kpc 고리에서 제외했다. matched
control의 ICL 값과 위성 프로파일은 고정했다. 원래 표적 ICL 평균은 저장된 annulus 측정과
각 74쌍에서 `1e-5 nJy arcsec^-2` 이내로 일치했다. 위성 고리 픽셀을 빼면 ICL 지원 픽셀은
중앙값 `i` 4.3%, `i2` 5.3% 감소했다. 한 쌍씩 coverage가 40% 아래로 내려가 같은 기준을
유지한 재분석에서는 `i` 20/21, `i2` 52/53쌍이 남았다.

반경별 profile은 raw, 미소광 보정값을 사용해 픽셀 공분산 단위와 맞췄고, 적색편이와
풍부도를 선형 통제한 부분상관을 5,000회 필터별 연결 5도 블록 재표집했다.

| Filter / support | Pair 수 | 원래 부분상관 r | 위성 고리 제외 r | 짝지은 변화 95% 블록 구간 |
|---|---:|---:|---:|---:|
| i, 원래 74표본의 filter subset | 21 | 0.120 | 0.082 | [−0.084, −0.003] |
| i, coverage ≥40% 공통 표본 | 20 | 0.118 | 0.083 | [−0.080, +0.004] |
| i2, 원래 74표본의 filter subset | 53 | 0.234 | 0.193 | [−0.064, −0.021] |
| i2, coverage ≥40% 공통 표본 | 52 | 0.217 | 0.176 | [−0.065, −0.018] |

`i2` 주 표본에서는 target 위성 고리 픽셀을 제외한 뒤 부분상관이 약 0.04 감소한다.
이는 관측 관계의 일부가 동일 픽셀의 재사용과 연결될 수 있음을 보여준다. 공통 coverage
표본의 위성 고리 제외 상관의 95% 블록 구간은 `[-0.034, 0.371]`로 0을 포함한다. 따라서
이 관계는 유의한 독립 물리 연관으로 확정할 수 없고, 단순 noise covariance만으로도
설명되지 않는다. `i` 표본은 특히 구간이 넓다.

원래 raw 측정의 per-pair 9×9 finite-lag 잡음 분산·공분산을 partial residual cross-products에서
빼는 별도 method-of-moments 진단은 상관을 `i` 0.120→0.122, `i2` 0.234→0.249로 조금
올렸다. 5,000 블록 재표집의 noise-corrected 95% 구간은 각각 `[-0.232, 0.555]`와
`[0.018, 0.449]`다. `i2`에서 per-pair 픽셀 잡음이 통제 후 satellite residual sum-of-squares의
약 10.2%, ICL은 약 1.7%를 차지했다. 이 진단은 lag/variance 모델의 독립 pair pixel-noise
항만 보정하며 앞의 shared-shell light 민감도에는 적용되지 않았다.

이 민감도는 작은 지원 픽셀 집합을 제거하는 고정표본 실험이지 PSF/광량 deblending 또는
최종 ICL 정의가 아니다. 위성 표적 고리 밖의 빛, 마스크 밖 PSF 날개, 재검출·완전성,
큰 규모 sky 변화, extinction/passband·보정 불확실도를 고치지 않는다. 노이즈 보정 부분상관
분석은 raw 값에 대한 별도의 대각 per-pair 측정오차 진단이며, 여기서 제외된 aperture의
새로운 오차 및 공분산은 아직 전파되지 않았다. 모든 결과는 `science_measurement=false`다.

## 재현 및 계보

- 지원 covariance selection 기본값 수정: `audit_satellite_icl_pixel_overlap.py`
- 범위 독립 검증기: `verify_satellite_icl_covariance_support_scope.py`
- noise-error 부분상관 민감도: `assess_icl_satellite_pixel_noise_bias.py`
- 동일 고리 픽셀 제외 민감도: `measure_icl_satellite_overlap_sensitivity.py`
- 94쌍 표적/satellite lag 공분산: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_icl_lag_covariance_2026-11-03_02/`
- 94쌍 target/control 공분산: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_icl_target_control_covariance_2026-11-03_02/`
- independent scope audit: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_icl_covariance_scope_verification_2026-11-03_01/verification.json`
- shared-pixel table: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_icl_pixel_overlap_2026-11-03_01/pixel_overlap.tsv`
- overlap sensitivity: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_overlap_sensitivity_2026-11-03_01/`
- overlap sensitivity summary SHA-256: `733433f2495ef4cf1592d7885a9845de57c726b6b64a2fdf93ca91514a0260e9`
- pixel-noise bias summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_pixel_noise_bias_2026-11-03_01/summary.json` (SHA-256 `268c5ee1a420e9cd97c479f8fe5ea8000bd3afaac41c7469ff3a1028bfed8396`)
- target/control covariance summary SHA-256: `90d9bac38f3bdea287fa0f5d54b75f76d8b842e80c32fe55674d0a95f9a174c9`
- scope verification SHA-256: `ce30786c1a2e394f1210e0410d597ec75c726dd34c8b416dd72690d5a24d6361`

다음 ICL 관문은 위성 outer-shell 제거 민감도를 출발점으로 삼아 PSF와 sky를 보존한 joint
injection/recovery를 수행하고, 각 지점에서 satellite-removed ICL profile과 matched-control
profile의 공간 block 공분산을 평가하는 것이다.
