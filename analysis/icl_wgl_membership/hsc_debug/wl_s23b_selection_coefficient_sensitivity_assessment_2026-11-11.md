# S23B 고정 52-lens profile의 HSC 선택계수 민감도

## 범위와 결론

승인된 HSC S19A shape-row export, HSC Y3 `n(z)`, 기존 지원표와 고정 profile을 사용해 HSC edge-selection 보정계수의 오차 반응을 계산했다. 기준 profile은 pinned profile과 최대 절대차 `0.03125 Msun Mpc^-2` 이내로 재현됐다. 실행은 `lageunha` CPU에서 끝났다. 추가 CAS 질의, 자료 다운로드, 이미지 읽기는 없었다.

표본은 고정된 렌즈 세트 52개와 matched control 260개로 구성됐다. 각 렌즈는 7개 radial annulus 모두에서 positive lens support와 최소 네 개의 positive control support를 가졌다. 각 보정계수는 공식 HSC Y3 calibration utility의 인용된 1σ 오차만큼 양·음 방향으로 하나씩 바꿨다. 각 변경에서 per-source 보정과 적층 전체를 다시 계산하고, 동일한 52개 렌즈 세트의 target-minus-control 평균을 구했다. HSC의 selection-bias 공식과 계수 정의는 [공식 HSC shape-catalog 문서](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s19a-shape-catalog-pdr3/) 및 문서에서 연결된 [HSC Y3 calibration utility](https://github.com/mr-superonion/hsc-y3-calib/blob/main/utils.py)를 따랐다.

각 계수의 ±1σ 변화로 얻은 profile 반응은 같은 표본의 2.5°·5° equal-area 공간 잭나이프 표준오차보다 작았다. Tangential 반응의 최대 비율은 1.21%였고 cross 반응의 최대 비율은 0.83%였다. 결과는 각 계수를 한 번에 하나씩 바꾼 민감도다. 결합 보정 공분산이나 전체 systematic error budget으로 해석하지 않는다. `science_measurement=false`를 유지한다.

## 최대 계수 반응

표의 profile 반응 단위는 `10^13 Msun Mpc^-2`다. 마지막 두 열은 각 반응의 절댓값을 fixed-sample spatial jackknife의 대응 성분 오차로 나눈 비율이다. 첫 비율은 2.5° 셀의 tangential/cross 오차이고, 둘째 비율은 5° 셀의 tangential/cross 오차다.

| R proper (Mpc) | 최대 tangential 반응 | 최대 cross 반응 | 오차 대비 2.5° | 오차 대비 5° |
|---:|---:|---:|---:|---:|
| 0.1–0.2 | 0.0433 | 0.0140 | 0.82% / 0.21% | 0.98% / 0.25% |
| 0.2–0.4 | 0.0247 | 0.0228 | 1.04% / 0.77% | 1.21% / 0.83% |
| 0.4–0.7 | 0.0093 | 0.0034 | 0.63% / 0.22% | 0.69% / 0.21% |
| 0.7–1.0 | 0.0047 | 0.0018 | 0.31% / 0.14% | 0.38% / 0.15% |
| 1.0–1.5 | 0.0032 | 0.0007 | 0.31% / 0.08% | 0.33% / 0.09% |
| 1.5–2.0 | 0.0013 | 0.0008 | 0.19% / 0.10% | 0.23% / 0.11% |
| 2.0–3.0 | 0.0029 | 0.0008 | 0.52% / 0.13% | 0.54% / 0.13% |

Tangential 응답은 각 annulus에서 `m_mag` perturbation이 가장 컸다. Cross 응답을 가장 크게 만든 계수는 radial bin마다 달랐다. 최대값 표는 네 계수의 반응을 annulus별로 취합했다. 네 계수의 부호가 다른 반응을 하나의 수로 합산하지 않았다.

## 해석과 남은 관문

단독 선택계수 오차는 현재 fixed-sample spatial error보다 훨씬 작다. Equal-area 셀은 검증된 연속 HSC 관측 영역이 아니므로 spatial error 자체도 확정된 survey covariance가 아니다. 작은 계수 반응은 HSC 선택효과 전체가 작음을 뜻하지 않는다.

Calibration utility는 magnitude-edge 및 resolution-edge 기여를 각 marginal `m_sel`·`a_sel` 오차 안에서 독립으로 합산하지만, multiplicative와 additive 보정 사이의 covariance를 제공하지 않는다. 따라서 본 계산은 coefficient-by-coefficient 반응만 보고하며, 전체 계통오차로 합치지 않는다.

다음 weak-lensing 주요 관문은 검증된 HSC footprint 기반 spatial covariance, target/control 영상의 depth·mask·PSF·sky 일치, photo-z dilution과 `n(z)` posterior, reduced shear, magnification, cluster-member boost다. 현재 결과만으로 질량이나 검출을 보고하지 않는다.

## 재현 정보

- 실행기 [analyze_s23b_wl_selection_coefficient_sensitivity.py](analyze_s23b_wl_selection_coefficient_sensitivity.py)
- 실행 결과 경로 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_fixed52_selection_coefficient_sensitivity_2026-11-11_01/`
- 결과 `summary.json` SHA-256 `3a92c5770fe55889f05d532a2ce81a8100e79df1ad93393fbf747d2d1f14df43`
- 결과 `coefficient_sensitivity.tsv` SHA-256 `4a727b3725f52e814e3dcec96ec7b818e5a011268a6fb13e94e114ac0aece2ed`
- 실행기 SHA-256 `07da1b74c09fe50080aa343b50bb434e7f82fc3ea988a2d837191e612441af3b`
- 승인 row-plan SHA-256 `7bfb998149f30a50b5ee7a537bb703ed302b7ea269a3296c35204a836e75aabc`
- 비교한 2.5°·5° spatial covariance summary SHA-256 `acdd4da69e4d7fe3dbb52d8584b3996bd836a41196da59aecfe3fcdf3e1fc6e9`
