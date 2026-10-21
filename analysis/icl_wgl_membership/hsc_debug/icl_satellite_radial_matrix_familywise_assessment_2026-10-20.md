# ICL–위성 반경행렬의 공간·다중검정 평가 — 2026-10-20

## 분석

앞선 지원 감사의 84개 ICL×위성 반경 조합을 모두 회귀에 넣었다. 표적과 매칭 대조장 좌표는 ICL annulus 표와 이미 검증된 `prestack_extraction` 좌표를 `control_id`로 연결해 복원했다. 420개의 고유 target/control cluster-filter 쌍에 동일한 5도 하늘 셀 기반 연결블록을 부여했다. 각 필터의 여섯 ICL 반경(50–100부터 400–500 kpc)과 일곱 위성 반경(5–10부터 50–75 kpc)을 교차했다.

각 교집합에서 `satellite ΔSB ~ ICL ΔSB + redshift + CAMIRA richness` cluster-level 선형모형의 ICL 부분상관과 cluster-robust t 통계를 계산했다. 제한 귀무모형 잔차에 HC2 조정 후, 모든 검정에 같은 공간블록별 Webb 6점 가중치를 적용한 5,000회 wild-cluster bootstrap을 사용했다. 84개 셀 모두 사전 설정한 최소 10행·8개 블록 기준을 통과했으며, 실제 셀별 지원은 i에서 21–78행/8–16블록, i2에서 57–150행/24–37블록이었다. 전체 84개 검정을 한 가족으로 보고 Holm 조정 p값과 단일 단계 max-T 완전귀무 p값을 계산했다.

## 결과

관심 조합인 ICL 200–300 kpc × 위성 50–75 kpc의 i2 부분상관은 `r=0.293` (68개, 25블록), 원시 wild-bootstrap `p=0.0288`이었다. 그러나 84개 전체 검정의 Holm 조정값은 `p=1.000`, max-T 값은 `p=0.987`로 가족 수준의 증거가 없다. 같은 조합의 i 결과는 `r=0.079` (26개, 9블록), `p_raw=0.621`이다. 표본 지지와 효과가 두 필터 사이에서 일관되지 않는다.

84개 중 미보정 `p<0.05`는 8개였으나, Holm 조정 최소값은 `0.588`, max-T 조정 최소값은 `0.156`이다. 따라서 반경을 탐색한 뒤 고른 특정 양의 계수를 검출로 해석할 수 없다. 예를 들어 i의 ICL 400–500 kpc × 위성 40–50 kpc는 `r=0.156`, `p_raw=0.0070`이지만 `p_Holm=0.588`이다. 반경군 어디에서도 보정 후 5% 수준의 연관성을 검출하지 않았다.

| 필터 | ICL (kpc) | 위성 (kpc) | 쌍 / 블록 | 부분상관 r | 원시 p | Holm p | max-T p |
|---|---:|---:|---:|---:|---:|---:|---:|
| i | 200–300 | 50–75 | 26 / 9 | 0.079 | 0.621 | 1.000 | 1.000 |
| i2 | 200–300 | 50–75 | 68 / 25 | 0.293 | 0.029 | 1.000 | 0.987 |
| i | 400–500 | 50–75 | 21 / 8 | 0.368 | 0.023 | 1.000 | 0.915 |
| i2 | 400–500 | 50–75 | 57 / 25 | −0.041 | 0.944 | 1.000 | 1.000 |

## 해석과 한계

이 결과는 사후 선택된 탐색 비교의 전체 반경 가족을 보정했을 때 현재 데이터에서 ICL–위성 연관을 검출하지 못한다는 뜻이다. 이는 연관이 없다는 증명이 아니다. 5도 공간 연결블록은 독립 관측영역을 근사하며, 일부 i 셀은 8–10블록뿐이다. max-T p값은 완전 귀무 기준 진단이고 subset pivotality가 검증되지 않아 강한 가족 오류율 통제로 주장하지 않는다. Holm 결과도 소수 블록 wild-bootstrap의 주변 p값 근사에 의존한다.

셀마다 측정 지원 표본이 달라 비무작위 결측/선택 편향을 제거하지 않았다. PSF 날개, 상류 sky 전달, 위성 검출 완전도, 소광, passband 차이, coadd 픽셀 공분산도 해결되지 않았다. 그러므로 `science_measurement=false`이며 인과적 stripping/ICL 생성 주장을 허용하지 않는다. 신규 HSC 자료 요청이나 다운로드는 없었다.

## 재현·검증

Astrafek 분석기는 [`analyze_icl_satellite_radial_matrix_familywise.py`](analyze_icl_satellite_radial_matrix_familywise.py), 회귀검사는 [`test_analyze_icl_satellite_radial_matrix_familywise.py`](test_analyze_icl_satellite_radial_matrix_familywise.py)다. 전체 84개 관측계수·공간블록 수, 다중보정, 출력 해시와 5,000개 max-T 값을 독립 재계산해 확인했다. 관련 분석기와 지원 감사 테스트 6개가 통과했고 `git diff --check`도 통과했다.

권위 산출물은 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_radial_matrix_familywise_2026-10-20_01/`이다. `summary.json` SHA-256은 `9e4bf2e8369c2c7b3de025a51eafc74538b560b9c42a9cf2fc68de49e1a8ccfb`, 결과표 `radial_matrix_association.tsv` SHA-256은 `99ae3afd7a8845de9fb668cb415dcc2b7f8edded60bba550e2d69b5765014b0d`, 5,000회 null 분포 `maxT_null.tsv` SHA-256은 `2fcfe2a763d33a1ee12255c35dc77e5623dfdcf5121f4283597629d3368f4c19`이다. 분석기 SHA-256은 `4aa446fd4b4ca38b58fb28ec92264ed8a51718ad00a369b1fdacd20f5e705717`다.
