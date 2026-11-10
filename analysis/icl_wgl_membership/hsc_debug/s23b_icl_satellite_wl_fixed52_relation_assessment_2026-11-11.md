# S23B 고정 52개 표본 ICL–위성–약중력 탐색 분석 (2026-11-11)

## 범위와 지위

해시가 고정된 기존 Astrafek ICL–CAMIRA member profile과 weak-lensing per-lens profile만 결합했다. HSC/CAS 요청, 새 다운로드, 원본 shape-row CSV 열람, 이미지 픽셀 분석은 없었다. 결과는 모두 `science_measurement=false`이며 보정된 약중력 측정이나 확인적 검출이 아니다.

전체 교집합은 52개이고, i 필터 25개와 i2 필터 27개는 서로 다른 은하단 집합이다. ICL 두 annulus의 최소 공통 픽셀 지원을 10% 이상으로 두면 i 18개, i2 22개가 남는다. 20% 이상은 각각 5개와 7개뿐이어서 상관 적합에는 사용하지 않았다.

## 계산

ICL–위성 지표는 ICL inner−outer contrast와 CAMIRA member-candidate inner/outer surface-density ratio의 로그다. 약중력 지표는 저장된 target-minus-random ΔΣ와 cross profile이며, 7개 radial bin에서 equal-lens mean을 계산했다. 관계 표에는 redshift와 ln(richness)을 조정한 partial correlation 및 target/control 위치의 연결된 5° 공간 블록 재표집 구간을 기록했다. 95% 구간은 점별 탐색용이며 다중비교 보정이나 물리적 유의도 판정이 아니다.

고정 52개 profile estimator에는 HSC S19A per-source multiplicative/additive shear terms (m, c1, c2), radial-bin responsivity, edge-selection response, HSC Y3 SACC tomographic n(z), 그리고 target-minus-random subtraction이 이미 적용되어 있다. 이는 알려진 응답 보정의 적용 사실이며 보정 불확실성이나 잔여 편향이 충분히 검증됐다는 뜻은 아니다.

### ICL–위성 관계

| 필터 | 공통 픽셀 선택 | N | 블록 수 | 조정 partial r | 점별 95% 블록-bootstrap 구간 |
|---|---|---:|---:|---:|---:|
| i | all | 25 | 13 | 0.046 | [-0.403, 0.752] |
| i | support_ge_0p10 | 18 | 10 | -0.227 | [-0.790, 0.352] |
| i2 | all | 27 | 19 | 0.069 | [-0.335, 0.561] |
| i2 | support_ge_0p10 | 22 | 16 | -0.116 | [-0.499, 0.239] |

네 점별 구간 모두 0을 포함한다. 현재 고정 표본에서는 적색편이와 풍부도를 조정한 뒤 ICL contrast–CAMIRA 후보 밀도비의 명확한 연관성이 분해되지 않는다.

### 약중력 stack 및 관계 민감도

| 필터 | R (proper Mpc) | N | tangential ΔΣ (10¹³ M☉ Mpc⁻²) | 점별 95% 블록-bootstrap 구간 | 중앙 lens/random pair 수 |
|---|---:|---:|---:|---:|---:|
| i | 0.1–0.2 | 18 | 32.618 | [23.880, 42.533] | 12.000 / 49.500 |
| i | 0.2–0.4 | 18 | 7.640 | [1.263, 12.525] | 42.000 / 178.000 |
| i | 0.4–0.7 | 18 | 2.678 | [-1.159, 9.937] | 111.500 / 481.500 |
| i | 0.7–1.0 | 18 | 1.932 | [-3.686, 6.170] | 171.500 / 735.500 |
| i | 1.0–1.5 | 18 | 1.094 | [-3.043, 3.652] | 389.000 / 1826.500 |
| i | 1.5–2.0 | 18 | 1.194 | [-0.542, 2.829] | 553.500 / 2621.500 |
| i | 2.0–3.0 | 18 | 0.430 | [-0.574, 1.769] | 1579.000 / 7442.000 |
| i2 | 0.1–0.2 | 22 | 16.925 | [2.745, 32.659] | 12.000 / 51.000 |
| i2 | 0.2–0.4 | 22 | 6.506 | [0.215, 12.827] | 37.000 / 188.500 |
| i2 | 0.4–0.7 | 22 | 7.317 | [4.878, 9.685] | 94.500 / 480.000 |
| i2 | 0.7–1.0 | 22 | 5.036 | [2.229, 7.918] | 152.000 / 778.500 |
| i2 | 1.0–1.5 | 22 | 2.723 | [0.237, 5.719] | 359.500 / 1908.000 |
| i2 | 1.5–2.0 | 22 | 1.329 | [-0.580, 3.378] | 502.000 / 2475.000 |
| i2 | 2.0–3.0 | 22 | 1.297 | [-0.010, 2.848] | 1462.000 / 6885.500 |

공통 픽셀 10% 이상 표본의 cross-profile 14개 점별 구간은 모두 0을 포함한다. 이는 약한 null 점검일 뿐 PSF·shear 보정을 검증하지 않는다.
전체 분석에서 약중력 partial-correlation 구간 8개가 112개 점별 구간 중 0을 포함하지 않았다. 위치는 `metrics.tsv`에 기록했지만 radial bin·필터·관계의 다중 탐색을 거쳤고 독립 검증이 없어 검출로 해석하지 않는다.

공통 픽셀 10% 미만 표본을 포함한 전체 결과와 10% 이상 결과를 나란히 계산했다. 20% 이상 표본은 수가 너무 적어 표본 수만 기록했다. 약중력 profile의 렌즈별 pair 수는 작을 수 있으므로 각 bin의 pair-count 열과 stack 표를 함께 확인해야 한다.

## 한계

표본은 높은 원시 S19A 형상원 가용성 및 ICL pixel support로 선택되어 전체 S23B 은하단을 대표하지 않는다. 필터별 표본이 disjoint하므로 i와 i2 값을 서로 짝지은 필터 비교로 해석하지 않는다. CAMIRA member positions는 위성 은하의 확정 분류가 아니며, 이미지 지지율이 낮은 프로파일에는 강한 선택 민감도가 있다. 이미 적용한 HSC 응답 보정과 HSC Y3 n(z)의 불확실성, 외부 검증표본에 기반한 photo-z dilution 보정, reduced shear, magnification, cluster-member dilution/boost, random aperture와 target의 깊이/mask/PSF/sky 매칭, HSC 조사영역 covariance는 해결되지 않았다. 별도 74-lens 민감도 분석에서 DESI n(z)와 magnification 대안을 비교했지만, 그것은 이 고정 52개 표본별 보정오차를 측정하지 않는다. 블록 bootstrap은 근사 공간 민감도만 제공하며 survey-region covariance를 검증하지 않는다.

## 산출물

- `metrics.tsv`: filter/support별 ICL–위성 및 radial-bin별 ΔΣ 관계와 점별 블록-bootstrap 구간.
- `stack_profiles.tsv`: 10% 공통 지지 표본의 tangential/cross equal-lens mean profiles와 pair counts.
- `joined_lens_profiles.tsv`: 분석에 사용한 52×7 lens/bin 조인 테이블.
- `fixed52_relation.png`: filter별 ICL–위성 산점도와 weak-lensing stack.
- `summary.json`: 입력/출력 SHA-256, 실행 범위, 상태.

완전한 관계별 수치와 점별 구간은 `metrics.tsv`를 기준으로 한다. 이 탐색 분석만으로 ICL–위성 연관성 또는 halo-mass 관계를 검출했다고 주장하지 않는다.

## 재현

실행기: [analyze_s23b_icl_satellite_wl_fixed52_relation.py](analyze_s23b_icl_satellite_wl_fixed52_relation.py)
결과 디렉터리: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_satellite_wl_fixed52_relation_2026-11-11_03`
결과 요약 SHA-256: `a7a9f30107bb800986d5d1aa4bf4e2012497092170e3bc27bf246587f7840e2f`
