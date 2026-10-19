# ICL–위성 연관성 측정 잡음 Monte Carlo 민감도 — 2026-10-19

## 목적과 표본

Astrafek은 HSC 표적·매칭 대조 coadd의 실제 ICL 및 위성 측정 가중치로 계산한 pixel-noise 공분산을, 그 공분산을 보유한 동일 cluster/filter의 연관성 관측값에 전파했다. 전체 반경별 연관성 표본은 94쌍(i 26, i2 68)이지만 측정 공분산이 있는 교집합은 57쌍(i 13, i2 44)이다. 따라서 이번 결과는 전체 94쌍 분석을 대체하지 않는 별도 민감도 표본이다.

실행기는 [`propagate_satellite_icl_measurement_noise.py`](propagate_satellite_icl_measurement_noise.py), Slurm 진입점은 [`run_satellite_icl_measurement_noise_mc.sbatch`](run_satellite_icl_measurement_noise_mc.sbatch)이다. Python 함수 단위시험 3개가 통과했고, Slurm job 428798은 exit 0으로 완료됐다. 기존 연관성 행과 target/control 공분산 표만 사용했으며 HSC 질의나 자료 반출은 없었다.

## 방법과 결과

적색이동과 richness를 고정한 채, 각 cluster/filter 행의 관측 ICL·위성 값을 행별 2변량 Gaussian 잡음으로 10,000회 교란하고 편상관을 다시 계산했다. diagonal 모델은 전달된 variance plane의 분산을 쓰고, finite-lag 모델은 대칭 9×9 잔차 공분산 kernel을 실제 측정 가중치에 전파한 분산·공분산을 쓴다. 두 모델을 따로 실행했으며 공간 블록 재표본화는 섞지 않았다.

| 필터 | 쌍 | 측정값 편상관 r | Diagonal 잡음 95% Monte Carlo 범위 | Finite-lag 잡음 95% Monte Carlo 범위 | Finite-lag에서 r>0인 비율 |
|---|---:|---:|---:|---:|---:|
| i2 | 44 | 0.409 | [0.337, 0.458] | [0.320, 0.465] | 100.0% |
| i | 13 | 0.164 | [0.022, 0.300] | [−0.016, 0.333] | 96.3% |

측정값 편상관은 입력 행으로 독립 재계산해 i2 `0.4088323816774757`, i `0.1637092174310721`로 확인했으며 보고서 수치와 일치한다. 전체지원 94행 분석의 i2 `r≈0.293`보다 공분산 교집합 44행의 `r≈0.409`가 크다. 이 차이는 PSF·공분산 지원 가능성에 따른 표본 선택이 결과를 바꿀 수 있음을 보여주므로, 44행 결과를 전체 은하단 관계로 확대 해석하지 않는다.

## 해석 한계

Monte Carlo 범위는 고정된 57개 측정쌍 아래의 잡음 모델 민감도이며 신뢰구간, 귀무가설 검정 p 값, 총 오차막대가 아니다. 필터·행별 Gaussian 및 pair 간 독립 잡음을 가정하고, 공간 블록 변동·위성 선택 완전도·HSC 상류 sky 전달·PSF 날개·source Poisson 잡음·소광·물리적 공변량 오차를 포함하지 않는다. 따라서 spatial-block 범위와 합쳐 총 불확실성으로 보고하지 않는다. 이 계산은 pixel-noise만으로 생길 수 있는 변화를 분리해 보여주지만 ICL–위성 검출 주장을 허용하지 않는다. `science_measurement=false`를 유지한다.

## 재현 산출물

- 결과: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_icl_measurement_noise_mc_2026-10-19_01/`
- `summary.json` SHA-256: `87ea03ffda04d60f0bc6e70f1206c147fa2d3262f46fea6e2c80b8add93e34f7`
- `measurement_noise_partial_r.tsv` SHA-256: `b7d2f0e7f42bfa0ff69a4d9b51d988b09bc5d00d452bb4ae4a9ddeb1c819ce00`
- 분석 스크립트 SHA-256: `b3433169f445622575b7391643333e37b03a6708e3cee835e88b3410feac8da3`
- 출력은 전체 94쌍의 공간 블록 민감도 및 가족단위 다중비교 검정과 함께 해석한다. 어느 분석도 현재 ICL–위성 관계의 확증 검출로 남지 않았다.
