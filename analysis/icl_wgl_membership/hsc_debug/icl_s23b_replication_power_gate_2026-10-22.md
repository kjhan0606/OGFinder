# S23B ICL–위성 독립 재현 표본력 관문 — 2026-10-22

## 과학적 판단

S23B의 현재 selection-compatible 후보 1,333개는 기존 ICL–위성 i2 분석과 다른 67개의 5도 target 셀을 차지한다. 이는 새 표본이 제공할 수 있는 **target 셀 상한**이며 아직 실사용 영상 표본도, 독립 연결 블록도 아니다. query coverage, masks, PSF, satellite completeness, matched-control assignments를 적용하면 블록 수는 줄거나 합쳐질 수 있다.

향후 표본이 들어오면 독립 재현 endpoint로 i2 ICL 200–300 kpc target-minus-control와 위성 50–75 kpc target-minus-control의 partial Pearson correlation을 redshift 및 richness로 조정해 검사하는 안을 제안한다. 검정은 양측 α=0.05 단일 endpoint로 고정하고 후보 제외·마스크·지원 기준을 S23B 결과를 열기 전에 얼려야 한다. 이 endpoint는 기존 94행에서 관찰한 뒤 정한 것이므로 독립 표본에서만 재현 검정으로 해석할 수 있다.

Fisher-z 부분상관 근사에서 단일 고정 endpoint 80% 검정력은 r=0.293일 때 92 independent blocks, r=0.25일 때 126, r=0.20일 때 196이다. S23B 단독 67-cell 상한의 nominal power는 각각 0.662, 0.520, 0.358이다. 즉 사후 추정값 r=0.293이 참이라고 해도 S23B 후보 상한만으로 독립 재현 80% 목표에 못 미치며, 최소 25개의 추가 유효 블록이 필요하다. 더 현실적인 r=0.20에서는 최소 129개가 더 필요하다. 어느 수치도 실제 독립 블록 수를 보증하지 않는다.

탐색 분석에서 30개 비교를 유지하면 r=0.293에서도 180 block이 필요하고, 기존 discovery 33개 target 셀과 S23B 상한 67개를 단순 합한 100은 80개가 모자란다. 이 33개를 S23B 독립 재현 표본에 포함할 수는 없다. 따라서 S23B metadata query가 통과하더라도 다른 독립 footprint 후보 또는 별도 사전 설계가 필요할 수 있다. 더 넓은 redshift selection을 늘려 후보 cluster 수만 증가시키는 것은 이미 S20A 카탈로그 분석에서 coarse-cell 증가가 작았으므로 해결책으로 간주하지 않는다.

## 재현과 한계

- 계산기: [`assess_s23b_replication_power.py`](assess_s23b_replication_power.py)
- 테스트: [`test_assess_s23b_replication_power.py`](test_assess_s23b_replication_power.py), 5 passed
- 입력 장부 SHA-256: `f6ac867881d8184d02a580e9a18540b7513a90f530e51e6cdf04205db535e8d1`
- 결과: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_replication_power_2026-10-22_01/`
- 요약 JSON SHA-256: `2b7a4bf4eeef17c74efb239565ca585fbf08cf547c0c173f254348cdf5f793bc`

계산은 target sky cell을 독립 block로 간주한 계획 상한에 적용한 Fisher-z 근사다. 공간 종속성, 후보 영상 미지원, 매칭 대조장에 의한 연결, 사후 효과크기 부풀림, HSC sky/PSF/마스크 계통효과를 반영하지 않아 실제 power를 과대평가할 수 있다. HSC query/download/측정은 하지 않았고 `science_measurement=false`다.
