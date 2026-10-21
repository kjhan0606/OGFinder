# ICL–위성 표본수 계획 민감도 — 2026-10-21

## 질문과 현재 지원

현재 CAMIRA parent 표본 800개 중 ICL 200–300 kpc와 위성 50–75 kpc의 공통 영상·마스크 지원을 통과한 결합 행은 i 26개, i2 68개다. 이전 168개 반경·필터·추정량 가족 검정은 두 필터 모두에서 가족 보정 후 연관성을 검출하지 못했다. 여기서는 같은 정도의 지원률을 가정할 때, 더 큰 parent 표본이 어느 정도 필요한지 계획 민감도로 계산했다.

Astrafek 재현기 [`plan_icl_sample_size.py`](plan_icl_sample_size.py)는 두 공변량(redshift, richness)을 둔 부분상관의 양측 비중심 t 검정에 대해 80% 검정력 표본수를 계산한다. 단일 사전지정 검정 `alpha=0.05`와 168개 가설에 대한 보수적인 Bonferroni `alpha=0.05/168`을 구분하고, 공간 의존 민감도용 설계효과 1.0 및 1.5를 제시한다. 설계효과는 현 표본에서 추정한 값이 아니라 입력 가정이다.

## 168개 가족 보정 시나리오

| 가정한 부분상관 | 독립 유효 행 수 (DE=1.0) | parent 렌즈 환산 i (DE=1.0) | parent 렌즈 환산 i2 (DE=1.0) | 유효 행 수 (DE=1.5) | parent 렌즈 환산 i2 (DE=1.5) |
|---:|---:|---:|---:|---:|---:|
| 0.10 | 2,849 | 87,662 | 33,518 | 4,274 | 50,283 |
| 0.20 | 536 | 16,493 | 6,306 | 804 | 9,459 |
| 0.293 | 323 | 9,939 | 3,800 | 485 | 5,706 |

Parent 표본 환산은 현재 i 26/800, i2 68/800의 joint-support 수율이 미래 표본에서도 그대로 유지된다고 가정한 단순 외삽이다. 필터별 행은 독립 관측으로 합산하지 않았다. 마지막 행의 `r=0.293`은 이미 탐색한 i2 endpoint의 점추정치에서 온 사후 가정이므로 편향될 수 있으며, 예측 효과나 확증 근거로 사용하면 안 된다. 0.10과 0.20도 관측 결과가 아닌 설계 시나리오다.

## 현재 공개 후보 풀과의 비교 (2026-10-21)

이미 승인된 S23B CAMIRA 카탈로그의 오프라인 감사를 연결했다. 기존 800개와 1-arcsec S20A 중복 제거, 현재 `0.30 <= z < 0.55` 및 BCG 분광 적색편이 일치 컷, 6-arcmin 분리 조건 뒤 추가 후보 상한은 1,333개다. 이를 전부 얻는다고 가정해도 parent 표본은 2,133개다. 현재 joint-support 수율을 그대로 적용하면 예상 i2 지원 행은 약 181개(현재 68 + 추가 113개)다. 이는 사후 시나리오 `r=0.293`, DE=1.0의 가족 보정 목표 323행보다도 작고, `r=0.20` 목표 536행에는 더 못 미친다.

이 후보는 카탈로그상 상한일 뿐이다. S23B 표는 PDR3 영상 manifest에 포함되지 않으며 후보 지점의 PDR3 coadd, 충분한 주변 footprint, 마스크·깊이·PSF 또는 매칭 대조 영상은 확인되지 않았다. 따라서 1,333개를 실제 분석 가능한 추가 은하단으로 세지 않는다. HSC 서비스 조회나 새 데이터 수신은 하지 않았다. 재현기 [`audit_icl_s23b_expansion_feasibility.py`](audit_icl_s23b_expansion_feasibility.py)의 입력 SHA 검증과 출력은 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_expansion_feasibility_2026-10-21_02/summary.json`에 저장했다. 출력 SHA-256: `e482ebd640230126b5d7d7e4438b78609b22899728765e512cc43b6bd21ab4e5`; 재현기 SHA-256: `c284eb8edc0fcc604b6c88022fed8ffe0f1c49f5c89230dc51b7664e3afdf528`.

별도로 로컬 S20A 부모 카탈로그의 기존 오프라인 확장 감사에서는 `0.20 <= z_cl < 0.70` 및 BCG 분광 redshift 조건 완화 시 2,596개가 추가 분리 후보로 남았지만, 새 5도 RA/Dec 셀은 6개뿐이었다(기존 78개에서 총 84개). 현재 i2 지원률을 가정하면 약 221행이 추가되어 총 약 289행으로 예상되며, 사후 `r=0.293`의 DE=1.0 목표 323행에도 못 미친다. 이 후보는 넓은 적색편이 범위와 cluster photo-z 선택을 사용하므로 현재 표본과 같은 모집단으로 단순 결합할 수 없다. S20A member 표는 로컬에 있지만 새 후보와의 member 연결·PDR3 영상·대조장 지원은 아직 이 확장 선택으로 검증하지 않았다. 결과와 해시는 [S20A 확장 풀 평가](camira_expansion_pool_assessment_2026-10-19.md)에 있고, 해당 member 파일은 6,139,116 bytes, SHA-256 `bf0c1fec5e11c6378ea69e3a65e49e9b32beecc83e5e68000bfa997e4ef2059e`다.

## 연구 판단

현재 800개 parent 표본은 이 joint-support 정의에서 가족 보정된 작거나 중간 정도 상관을 검출하기 위한 규모보다 작다. 하지만 이 표본수 계산만으로 더 많은 HSC 영상을 요청하는 것은 충분한 다음 단계가 아니다. 현재의 제한된 joint support는 마스크와 영상 품질에 따른 선택을 포함하고, PSF 날개·sky 전달·위성 완전도·coadd 공분산도 해결되지 않았다. 새 표본은 이 관측 계통을 동일한 기준으로 측정할 수 있을 때만 검정력을 높인다. 따라서 이 분석은 데이터 확대 우선순위를 정하는 계획 민감도이며 `science_measurement=false`를 유지한다.

## 재현

```bash
python3 analysis/icl_wgl_membership/hsc_debug/plan_icl_sample_size.py \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_sample_size_plan_2026-10-21_01/sample_size_plan.json
```

Focused tests: `python3 -m pytest -q analysis/icl_wgl_membership/hsc_debug/test_plan_icl_sample_size.py` — 3 passed. Result SHA-256: `df3d6d06adecaa310ab94df94bbebd1441262087520c16385927f285846f0b88`. Script SHA-256: `5a0a30ead5d559cc8e3ab23cbb7e6c7224f155a9d56718c2a2e714923c728724`.
