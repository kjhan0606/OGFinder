# ICL 초과광과 위치별 PSF 코어 크기의 공변동 — 2026-10-20

위치별 모델 PSF 코어가 관측 ICL annulus 측정과 함께 변하는지, 기존 공동 주입 표본 57개 프레임에서 점검했다. 각 프레임의 중심 위치 HSC-I PSF 장축 모멘트 등가 FWHM을, 해시 검증한 ICL 200–300 kpc target-minus-control 표면밝기와 결합했다. PSF 요약·파일, 57개 주입 프레임 결과, ICL annulus 입력, matched-control 좌표를 모두 재검증했다. 표본은 i 13개와 i2 44개이고, 각각 10개와 24개의 연결된 5° target/control 공간블록으로 나뉜다.

| 필터 | N | 보정 모델 | 부분상관 `r` | 블록 재표집 2.5–97.5 백분위 | 한 블록씩 제외한 `r` 범위 |
|---|---:|---|---:|---:|---:|
| i | 13 | 무보정 | 0.511 | [0.148, 0.855] | [0.378, 0.633] |
| i | 13 | z, richness, common coverage | 0.500 | [−0.026, 1.000] | [0.398, 0.627] |
| i | 13 | z, richness, target/control coverage 각각 | 0.534 | [−0.980, 1.000] | [0.405, 0.777] |
| i2 | 44 | 무보정 | 0.242 | [0.071, 0.482] | [0.206, 0.318] |
| i2 | 44 | z, richness, common coverage | 0.247 | [−0.008, 0.458] | [0.199, 0.306] |
| i2 | 44 | z, richness, target/control coverage 각각 | 0.246 | [0.002, 0.445] | [0.193, 0.303] |

따라서 이 선택 표본에서는 PSF 코어가 넓을수록 측정된 ICL 초과광도 큰 양의 공변동이 관찰됐다. 보정 모델의 leave-one-block-out 점추정은 계속 양수지만, 재표집 범위는 i에서 매우 넓고 i2에서도 경계값 부근을 포함한다. i의 표적·대조 지원률 별도 조정은 10,000회 중 51회, 공통지원 조정은 1회가 특이한 nuisance 표본 때문에 유효 상관을 만들지 못했다. 재표집 백분위는 기술적 민감도 범위이지 교정된 신뢰구간이나 검정 p값이 아니다. 이 패턴은 PSF·시야 조건이 ICL 신호와 관련될 가능성을 보여 주지만, 원인이나 보정량을 정하지 않는다.

공통지원 공변량에는 저장된 ICL annulus 마스크의 coverage를 사용했다. 공동 주입 실행은 검출·마스크를 다시 만들기 때문에 그 실행의 coverage가 저장 annulus 값과 달랐다. 재검출 minus 저장 coverage 차이의 중앙값은 i `0.016`(q16–q84 `0.000–0.120`), i2 `0.032`(`0.000–0.143`)였고 최대 차이는 i2 `0.343`이었다. 두 coverage 정의를 혼합하지 않고 산출물에 모두 보존했다.

이 검사는 57개 PSF 검증 표본이며 전체 800개 모집단이 아니다. 중심 PSF 모델은 약 7 arcsec 크기의 코어뿐이므로 확장 PSF 날개, 대규모 sky 잔차, 마스크·noise 계통효과는 포함하지 않는다. 따라서 이 상관은 ICL 검출이나 관측값 보정이 아니다. `science_measurement=false`, `psf_wings_validated=false`를 유지한다.

재현기는 [`audit_icl_psf_core_dependence.py`](audit_icl_psf_core_dependence.py), 테스트는 [`test_audit_icl_psf_core_dependence.py`](test_audit_icl_psf_core_dependence.py)다. 산출물은 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_psf_core_dependence_2026-10-20_02/summary.json` (SHA-256 `2fbbf6f89c16eea38153f19002b6e14b8b63d28b72d73c6fabd7ab000af5a805`)에 있다. 신규 HSC 자료 요청은 없었다.

영역 평균 차이 민감도는 [`icl_psf_region_fixed_effects_assessment_2026-10-20.md`](icl_psf_region_fixed_effects_assessment_2026-10-20.md)에서 추가로 확인했다. 근사 tract 영역 고정효과 후 부분상관은 i `0.484`, i2 `0.231`이며, 이 분석도 탐색적 민감도이지 검출이나 유의성 검정은 아니다.
