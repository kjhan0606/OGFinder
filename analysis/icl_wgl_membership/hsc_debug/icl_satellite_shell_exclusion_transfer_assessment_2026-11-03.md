# ICL–위성 고리 제외 주입·회수 민감도 — 2026-11-03

## 표본 범위와 방법

2026-11-03에 고정한 ICL–위성 연관표본 74개 target/filter 쌍과 94개 joint-annulus
지원쌍의 교집합을 기준으로 삼았다. 승인된 위치별 HSC PSF 코어까지 있는 쌍은 정확히
48개(12 `i`, 36 `i2`)였으며, 26개 연관쌍은 이 PSF 지원 집합에 들어오지 않았다. 샘플
선택기는 independent audit에서 검증된 74행 association table을 직접 읽고 이 교집합만
실행하도록 확장했다. 넓은 57프레임 joint-support batch 중 association 밖 9개는 이번
연관성 결과에서 제외했다.

각 실제 target coadd에 diffuse 성분과 위성 envelope 성분을 따로 또는 함께 주입하고,
기존 위치별 HSC-I 모델 PSF 코어로 합성 성분을 convolution했다. 동일한 고정 segmentation,
source mask, variance를 사용해 기존 ICL 200–300 kpc 추정량과 위성 50–75 kpc target 픽셀
합집합을 제외한 ICL 추정량을 둘 다 계산했다. shell 픽셀은 이미 지원된 표적 위성 profile의
50–75 kpc 고리에서 가져왔고, 각 프레임에서 원래 74표본 규약과 같은 coverage·대조 프로필
수 조건을 적용했다. 실제 science coadd는 수정하지 않았다.

결과는 48개 선택 프레임의 기술적 중앙값 및 q16–q84 범위다. q16–q84는 모집단 신뢰구간이
아니다. 각 프레임별 shell-excluded−raw 차이도 따로 산출했다.

## 결과

회수 응답은 단위 주입 세기를 기존 ICL 측정에 전달한 비율이다.

| 필터·주입 성분 | 프레임 수 | 기존 ICL 응답 중앙값 [q16, q84] | 위성 고리 제외 응답 중앙값 [q16, q84] | 짝지은 응답차 중앙값 [q16, q84] |
|---|---:|---:|---:|---:|
| i · diffuse | 12 | 0.22564 [0.22284, 0.22764] | 0.22491 [0.22298, 0.22805] | −0.00018 [−0.00168, +0.00102] |
| i · satellite envelopes | 12 | 0.46709 [0.15592, 0.77427] | 0.46624 [0.14866, 0.78033] | +0.00012 [−0.01361, +0.00287] |
| i2 · diffuse | 36 | 0.22594 [0.22282, 0.22951] | 0.22623 [0.22222, 0.22957] | 0.00000 [−0.00111, +0.00133] |
| i2 · satellite envelopes | 36 | 0.68125 [0.29405, 0.91706] | 0.65981 [0.28765, 0.91710] | 0.00000 [−0.06558, +0.01211] |

고리 제외로 빠진 유효 ICL 기하 픽셀 비율의 중앙값은 `i` 2.46%, `i2` 3.70%였다.
고리와 겹치는 유효 픽셀이 하나도 없는 프레임은 각각 1/12와 6/36이었다. 네 필터/주입
조합 모두 짝지은 응답차 q16–q84 범위가 0을 포함했다. 결합 주입의 ICL 응답에서 두 개별
응답의 합을 뺀 중앙 잔차는 필터별로 `4e-8 nJy arcsec^-2` 이하였다. 이는 해당 고정 마스크
모델에서 주입 합성이 사실상 선형임을 보일 뿐 관측 신호의 보정 정확도를 검증하지 않는다.

## 해석과 제한

공유 위성 픽셀을 제외해도 diffuse ICL 회수는 거의 바뀌지 않았다. satellite-envelope→ICL
응답은 `i2`에서 marginal median 기준 약간 낮아졌지만, 프레임별 짝지은 차이의 중앙값은 0이고
q16–q84가 넓다. 따라서 현재 결과는 단순한 scalar shell 보정을 뒷받침하지 않는다. 위성 고리
제외로 관측 ICL–위성 부분상관이 낮아졌다는 앞선 민감도 결과와 함께 해석해야 한다. 두 분석은
서로 다른 양(측정된 상관과 합성광 전달)을 다루며 하나가 다른 하나를 보정하지 않는다.

이 실행은 74개 연관쌍 중 PSF 위치가 확보된 48개만 사용해 선택 표본이다. PSF 입력은 compact
HSC-I 모델 코어라 확장 wing을 검증하지 않으며 i2에 대한 전용 wing도 아니다. source Poisson
noise를 추가하지 않았고 coadd의 상관 noise, 상류 sky subtraction, 재검출/마스크 재선택,
matched-control 분산과 전체 공간블록 공분산도 포함하지 않았다. 그러므로 `science_measurement=false`,
`psf_approved=false`, `observational_validation_passed=false`, `flux_covariance_available=false`를
유지한다. 주입 회수율을 관측 profile에 나누어 보정하지 않는다.

## 재현과 provenance

- 단일 프레임 주입기: [`joint_icl_satellite_injection.py`](joint_icl_satellite_injection.py)
- batch selector/runner: [`run_joint_icl_satellite_batch.py`](run_joint_icl_satellite_batch.py)
- 실행에는 verified `--association-pairs`와 `--exclude-supported-satellite-shells-from-icl`를 사용했다.
- 48프레임 결과: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_shell_exclusion_association_psf48_2026-11-03_01/`
- batch summary SHA-256: `2f6c3d2ae393dd14532d640a1f9710ddee3064f8daf4feaca8f63c57aaf1821d`
- association input SHA-256: `483f81540543fc46ec1cf4333bbd3e238799b493e48fe03235da3038f5714917`
- batch는 48개 프레임 파일 checksum, verified 74행 표본 membership, 완료/오류 수, 각 필터의 프레임별 회수 응답을 독립 재집계해 확인했다.
- 별도로 보존된 57개 joint-support 결과는 association 표본 밖 프레임 9개를 포함하므로 이 association 평가에는 쓰지 않았다.
