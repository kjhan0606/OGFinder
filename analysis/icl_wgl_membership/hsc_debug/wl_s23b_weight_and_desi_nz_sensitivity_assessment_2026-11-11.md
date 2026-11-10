# S23B 대조 가중치 및 DESI `n(z)` 민감도 평가 (2026-11-11)

## 범위와 방법

승인된 S23B HSC S19A shape 행 80개 렌즈/400개 대조 중심, 기존 HSC Y3 SACC 기준, 앞서 승인된 DESI/HSC `n(z)` 자료만 재사용했다. 새 CAS 조회·네트워크 요청·자료 다운로드·영상 픽셀 읽기는 없었다. 렌즈 세트와 대조를 모두 지지하는 공통 표본을 반경 구간별로 고정한 뒤 pooled, 대조 중심 동일 가중, 렌즈 세트 동일 가중 추정량을 비교했다. 별도로 두 DESI scale cut과 6개 `n(z)` 곡선을 SACC 기준 프로파일과 비교했다. 질량 적합은 하지 않았고 `science_measurement=false`다.

## 가중 방식 민감도

값의 단위는 `10^13 Msun Mpc^-2`다. `공통 pooled`는 해당 반경에서 아래 N개 렌즈 세트만 사용한 기준값이다. `대조 동일 가중`은 렌즈 항은 pooled로 두고, 세트별 random profile을 같은 가중치로 평균해 뺀다. `짝 동일 가중`은 각 세트의 렌즈−대조 차이를 같은 가중치로 평균한다.

| R (proper Mpc) | 공통 pooled | 대조 동일 가중 | 짝 동일 가중 | 짝 동일 가중 변화 | 세트 수 |
|---:|---:|---:|---:|---:|---:|
| 0.1–0.2 | 23.5766 | 24.5626 | 23.8681 | +1.24% | 60 |
| 0.2–0.4 | 6.9182 | 7.0503 | 8.4897 | +22.71% | 67 |
| 0.4–0.7 | 3.1736 | 2.8952 | 2.5583 | −19.39% | 69 |
| 0.7–1.0 | 4.1110 | 4.4934 | 3.6189 | −11.97% | 73 |
| 1.0–1.5 | 2.6300 | 2.7464 | 2.6955 | +2.49% | 74 |
| 1.5–2.0 | 1.5032 | 1.1423 | 1.4923 | −0.73% | 74 |
| 2.0–3.0 | 0.7621 | 0.7119 | 0.7758 | +1.80% | 74 |

중간 반경에서 렌즈 세트 가중이 프로파일을 약 12–23% 바꾼다. 대조 중심 동일 가중만으로도 1.5–2 Mpc에서 −24.0% 차이가 난다. 따라서 기존 pooled random-subtraction 결과는 세트 구성과 aperture별 source 지지에 민감하다. 내측 첫 구간의 공통 지지는 60세트뿐이다. 이 가중 민감도 추정량에는 독립 공간 공분산이 없고 세트 간 산포는 오차막대가 아니다.

## DESI/HSC `n(z)` 민감도

모든 비교는 기존 주 표본 74개 렌즈/367개 대조 중심, source bins 2–4를 사용한다. 각 행은 7개 반경에서 SACC 기준 대비 `DeltaSigma_t` 변화의 중앙값과 범위다.

| DESI scale cut / 곡선 | 반경별 변화 중앙값 | 최솟값–최댓값 | 불확도 처리 |
|---|---:|---:|---|
| 0.3–3 `mean` | −2.95% | −3.11%–−2.33% | 5° endpoint jackknife |
| 0.3–3 `mean_mu` | −2.41% | −2.57%–−2.04% | 5° endpoint jackknife |
| 0.3–3 `lower` | −2.46% | −2.60%–−2.08% | 점별 주변 envelope, 공분산 없음 |
| 0.3–3 `upper` | −3.34% | −3.60%–−2.55% | 점별 주변 envelope, 공분산 없음 |
| 0.3–3 `lower_mu` | −2.08% | −2.30%–−1.84% | 점별 주변 envelope, 공분산 없음 |
| 0.3–3 `upper_mu` | −2.61% | −2.76%–−2.17% | 점별 주변 envelope, 공분산 없음 |
| 1–5 `mean` | −3.18% | −3.60%–−2.67% | 5° endpoint jackknife |
| 1–5 `mean_mu` | −2.75% | −3.07%–−2.36% | 5° endpoint jackknife |
| 1–5 `lower` | −2.74% | −2.94%–−2.37% | 점별 주변 envelope, 공분산 없음 |
| 1–5 `upper` | −3.60% | −4.16%–−2.92% | 점별 주변 envelope, 공분산 없음 |
| 1–5 `lower_mu` | −2.45% | −2.58%–−2.14% | 점별 주변 envelope, 공분산 없음 |
| 1–5 `upper_mu` | −3.00% | −3.48%–−2.50% | 점별 주변 envelope, 공분산 없음 |

중앙 `mean` 곡선 민감도는 약 −3%다. 이는 profile을 modest하게 이동시키지만 photo-z 불확도를 완전히 전파하지 않는다. DESI `lower`/`upper`는 점별 주변 envelope라 joint posterior 또는 신뢰구간으로 해석할 수 없다.

## 판단 및 남은 관문

이 검사는 pooled 대조 가중이 무해하다고 볼 수 없음을 보여준다. 후속 estimator는 lens-set matching과 radial-bin support를 명시하고, 0 지지 렌즈·대조 처리 규칙을 사전에 고정해야 한다. DESI `n(z)`는 SACC 기준을 몇 퍼센트 이동시키는 보정 민감도지만, 완전한 광도거리 분포 공분산을 제공하지 않는다.

5° RA/`sin(Dec)` 블록은 분리된 HSC 영역을 하나로 묶을 수 있다. 기존 cross-shear 진단도 대조 선택과 블록 폭에 따라 달라졌다. reduced shear, magnification, cluster-member dilution/boost, 대조 위치의 depth/mask/PSF/sky 일치, selection-response, 반복 source 및 survey-field covariance, injection/recovery가 해결되지 않았다. 따라서 질량 적합이나 과학 검출 주장을 진행할 수 없고 `science_measurement=false`를 유지한다.

## 재현 정보

- 실행기: [`analyze_s23b_wl_weight_and_desi_nz_sensitivity.py`](analyze_s23b_wl_weight_and_desi_nz_sensitivity.py)
- 실행기 SHA-256: `af42929b6c392a2ccf98cf11780de46f141450ec9fc43f4e06af51492094ae18`
- 기준 프로파일 summary SHA-256: `2c899b0090dd94a34f1fd6b0cb111e3d1ada696da2178f98e9fd5b659e699edf`
- 승인된 shape-row 계획 SHA-256: `7bfb998149f30a50b5ee7a537bb703ed302b7ea269a3296c35204a836e75aabc`
- QA summary SHA-256: `eb73aa0365813d880d3edf60b12e4039f2575f8a97eab8868fc918121a34d445`
- HSC Y3 SACC `n(z)` SHA-256: `a28f9e2e088e92d96d2d87083a6eeeac4c0c84c6b48e033bd3d1dc1bf0d8958f`
- DESI 자료 commit: `abc70af29910e49fce17a63ef3083160d515fb6c`
- 최종 산출물 (Git 외): `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_weight_and_desi_nz_sensitivity_2026-11-11_03/`
- `summary.json` SHA-256: `9c66b281e0920b79a5c7854e40081b873621867cfaca2a9991d5e93e031f1aa5`
- `assessment.md` SHA-256: `4e4eca2a5c1461c42ca9137c2bf48e575d50656939ddc1cccf984f8512b8de20`
- weighting TSV SHA-256: `9618569b5e4da15dafbff8d1be84e86ab29113bde45a14d6d79d1c314595d97e`
- sensitivity TSV SHA-256: `b626a5d703063b5914eac18b76276c540ba077b90a11f83b29290cbc8676d898`
- plot SHA-256: `47920fae74098ae7fd5b74c97f1f0cb39a82f3eacad492dad85c66b47bbc0b3a`

`_01`은 radial별 지지가 다른 random profile과 전체 74세트 target profile을 빼어 가중 비교에 사용하지 않는다. `_02` 수치는 같지만 설명문의 지지 세트 범위가 잘못됐다. 최종 판정은 공통 지지와 설명 문구를 바로잡은 `_03`만 사용한다.
