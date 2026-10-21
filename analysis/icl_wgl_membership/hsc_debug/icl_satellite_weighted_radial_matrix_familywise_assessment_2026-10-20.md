# CAMIRA 가중 프로필을 포함한 ICL–위성 반경 가족 민감도 — 2026-10-20

## 목적과 입력

이전 equal-satellite 분석에서 84개 반경/필터 가설을 한 가족으로 보정했지만, 위성들을 같은 비중으로 평균하는 정의 하나에 의존했다. 기존 저장 member-profile 상세표에서 각 cluster/filter/위성 반경 안의 측광값을 CAMIRA 열 8의 명시된 `weight`로 가중해 위성 평균을 다시 계산하고, 그 프로필과 equal-satellite 프로필을 동시에 검정했다. 각 클러스터는 두 방식 모두에서 동일 가중했다. CAMIRA member 원자료, hash-matched acquisition summary, 기존 800-frame satellite profile 및 ICL annulus 표를 해시 검증했다. 가중·비가중 프로필 모두 84개 반경 조합에서 표본 행과 5도 공간블록 수가 정확히 일치했다.

## 방법

6개 ICL 반경 × 7개 위성 반경 × i/i2 필터 × 두 satellite estimator를 합친 **168개 가설**을 한 가족으로 두었다. 각 셀에서 `satellite ΔSB ~ ICL ΔSB + redshift + CAMIRA richness`를 적합하고 cluster-robust CR1 t 통계를 계산했다. 420개 고유 target/control pair로 정의한 같은 5도 연결공간블록을 모든 검정에서 공유하고, 제한 귀무모형의 HC2 잔차에 Webb 6점 wild-cluster bootstrap 5,000회를 적용했다. 각 검정은 최소 10 cluster/filter 행, 8개 공간블록이어야 했다. 168개 모두 기준을 통과했다.

## 결과

관심 조합인 ICL 200–300 kpc × 위성 50–75 kpc에서 두 estimator는 같은 쌍(26 i, 68 i2)과 같은 블록 수(9 i, 25 i2)를 사용했다.

| 위성 estimator | 필터 | 쌍 / 블록 | 부분상관 r | 원시 wild-bootstrap p | 168개 Holm p | 168개 max-T p |
|---|---|---:|---:|---:|---:|---:|
| equal-satellite | i | 26 / 9 | 0.079 | 0.621 | 1.000 | 1.000 |
| equal-satellite | i2 | 68 / 25 | 0.293 | 0.030 | 1.000 | 0.999 |
| CAMIRA `weight` | i | 26 / 9 | 0.076 | 0.602 | 1.000 | 1.000 |
| CAMIRA `weight` | i2 | 68 / 25 | 0.241 | 0.055 | 1.000 | 1.000 |

i2에서 CAMIRA 가중 프로필의 점 추정은 equal-satellite 값보다 작았다(차이 약 `−0.052`). 이 표는 estimator 차이의 별도 paired 유의성 검정을 제공하지 않으며, 그 차이를 물리 효과로 판정하지 않는다. 전체 168개 원시 p값 중 17개가 0.05 아래였지만 Holm 조정값은 모두 0.05 이상(최소 1.000), max-T p값은 모두 0.05 이상(최소 0.255)이었다. member weighting으로도 다중 반경·필터 탐색 후의 연관성 검출은 나오지 않았다.

## 해석과 제한

CAMIRA weight는 알고리즘의 질량 및 공간 필터가 이미 반영돼 독립적·불편한 위성 멤버십 확률이 아니다. 필터·반경별 선택 지원이 다르고 missing-not-at-random을 교정하지 않았다. 5도 블록은 공간 공분산의 근사치며 i 대역 최저 지원 셀은 8블록이다. max-T는 완전 귀무에서의 단일 단계 진단이고 subset pivotality를 확인하지 않아 강한 family-wise 통제를 주장하지 않는다. Holm p값도 소수 cluster 블록에서의 wild-bootstrap 주변 p값 근사를 사용한다. PSF 날개·하늘 전달·위성 완전도·소광·passband·coadd 공분산이 남아 `science_measurement=false`다. 따라서 분석은 현재 HSC 표본의 estimator 민감도와 미검출을 보고할 뿐, stripping이나 ICL 형성의 물리 증거가 아니다. 신규 HSC 자료를 조회하거나 받지 않았다.

## 재현과 검증

분석기 [`analyze_icl_satellite_weighted_radial_matrix_familywise.py`](analyze_icl_satellite_weighted_radial_matrix_familywise.py)는 이전 84셀 공간/가족 분석기를 재사용하고 CAMIRA weighted member profile을 같은 표본에 연결한다. Slurm job 430809 (`a10`, GPU 1, 8 GB, 1시간 요청)는 4분 16초에 완료했고 최대 RSS는 57,896 KB, stderr는 비어 있었다. 독립 validator [`validate_icl_satellite_weighted_radial_matrix_familywise.py`](validate_icl_satellite_weighted_radial_matrix_familywise.py)가 입력·출력 해시, 168개 정확 셀, 두 estimator의 지원 일치, 5,000개 null값에서 max-T p, 원시 p에서 Holm p를 재계산해 통과했다. 관련 단위검사 8개가 통과했다.

권위 결과는 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_weighted_radial_matrix_familywise_2026-10-20_01/`이다. `summary.json` SHA-256 `0a337df80fda6fb92f178d251e09347f27012f5cc53f088fd967311b31ef8068`; 168행 association table `weighted_radial_matrix_association.tsv` SHA-256 `6e6c37a44823ab0d6ac5f1795af87f4d47cc3dda1f65b0872ea9e622828142e7`; support table SHA-256 `1411eca0e449f1d95dac24ffb6cf0eaa498032f964ca87b380b5ba1a528b5b42`; joint max-T null table SHA-256 `04386baebe37577b40b96d2c859f9202e5cdc73b5eaf01b7ed6c1cc6b53e56a5`이다. 분석기 SHA-256은 `7ec81f30a1491a7ff2f6602a8572654c738cd0513dc31024351aea2ca77e6d3d`다.

## 동일 표본의 estimator 차이 paired 민감도 — 2026-10-20

관심 조합 두 필터에서 같은 cluster/filter 행과 연결 5도 공간블록을 고정하고, equal-satellite 대비 CAMIRA 가중 부분상관 차이의 10,000회 paired block-bootstrap 95% percentile 구간을 계산했다. 이는 앞선 탐색에서 선택된 단일 endpoint의 기술적 민감도이며 별도 귀무검정이나 보정된 발견 주장이 아니다.

| 필터 | 행 / 블록 | equal r | CAMIRA-weighted r | 차이 (weighted−equal) | paired block-bootstrap 95% 구간 |
|---|---:|---:|---:|---:|---:|
| i | 26 / 9 | 0.079 | 0.076 | −0.003 | [−0.225, 0.092] |
| i2 | 68 / 25 | 0.293 | 0.241 | −0.052 | [−0.139, 0.061] |

두 차이 구간 모두 0을 포함한다. i2 bootstrap 재표본 중 84.7%에서 차이가 음수였지만, 이는 재표집 부호 빈도이며 posterior 확률이나 유의확률이 아니다. CAMIRA 가중으로 점 추정이 감소하는 방향은 보이나 해당 감소 자체는 정밀하게 분리되지 않는다. 이전 168가설 가족조정 결과(유의 셀 0)도 그대로 유지한다. `science_measurement=false`다.

재현기 [`analyze_icl_satellite_weighting_contrast.py`](analyze_icl_satellite_weighting_contrast.py), 회귀검사 3개 통과. 입력 지원·가중 프로필과 10,000회 bootstrap을 독립 재계산해 산출물과 일치함을 확인했다. 결과 디렉터리는 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_weighting_contrast_2026-10-20_01/`; `summary.json` SHA-256 `b1a07b7cd46c8a2dc8fd562663a3297ba1d4f46c3a0b25df9ab35f9b56ffc699`, `paired_weight_contrast.tsv` SHA-256 `992bf7e7883b3b21babf88d6ba73c0d54e52d1e2ed8b1a7ff2ed07b294d0a012`, 재현기 SHA-256 `142c5ef82483ed9f9cee383b75c59aaa7b22173546a3f6be91de93c45d72e5a9`다.
