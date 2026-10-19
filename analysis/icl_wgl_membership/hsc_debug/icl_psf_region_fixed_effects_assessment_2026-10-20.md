# ICL–PSF 코어 공변동의 관측영역 민감도 — 2026-10-20

기존 57개 프레임 PSF–ICL 상관이 HSC 관측영역 사이의 평균 차이에서 비롯됐는지 확인했다. patch-QA 표의 840개 tract 중심을 구면거리 2° 기준으로 연결요소 네 개로 묶고, 각 target/control 좌표를 가장 가까운 tract 영역에 배정했다. 57쌍 모두 두 좌표가 같은 근사 영역에 배정됐다. 따라서 target과 control에 대한 영역 표시변수는 중복이며 한 세트만 모델에 넣었다.

| 필터 | 조정 | 영역 고정효과 전 r | 영역 고정효과 후 r | 잔여 자유도 |
|---|---|---:|---:|---:|
| i (N=13) | 없음 | 0.511 | 0.470 | 11 |
| i | z, richness, 공통 coverage | 0.500 | 0.484 | 8 |
| i | z, richness, target/control coverage 별도 | 0.534 | 0.549 | 7 |
| i2 (N=44) | 없음 | 0.242 | 0.234 | 41 |
| i2 | z, richness, 공통 coverage | 0.247 | 0.231 | 38 |
| i2 | z, richness, target/control coverage 별도 | 0.246 | 0.231 | 37 |

공통 coverage 모델에서 영역 연결 기준을 1.6°에서 3.0°까지 바꾸어도 tract 분할과 부분상관은 변하지 않았다. 이 단순 근사 영역의 평균 효과를 넣어도 양의 점추정은 유지됐다. 다만 i 표본은 영역 보정 후 자유도가 8에 불과하고, 이 분석은 불확실도 구간이나 p값을 산출하지 않았다. 영역 보정 후 양의 관계가 유의하거나 물리적 효과라고 결론 내릴 수 없다.

이는 broad field-region 평균 차이에 대한 민감도만 다룬다. tract 중심의 연결 영역은 실제 관측 footprint를 대신하지 않으며 PSF 날개, sky, 마스크, noise 또는 선택 공분산을 보정하지 않는다. 기존 제한과 같이 `science_measurement=false`, `psf_wings_validated=false`를 유지한다.

재현 코드는 [`audit_icl_psf_region_fixed_effects.py`](audit_icl_psf_region_fixed_effects.py), 테스트는 [`test_audit_icl_psf_region_fixed_effects.py`](test_audit_icl_psf_region_fixed_effects.py)다. 결과는 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_psf_core_region_sensitivity_2026-10-20_01/summary.json` (SHA-256 `967a33b5290620e28ae504a33d2fcb7a438825b1adcc4d13ac53abb5ca149d91`)에 있다. 입력 PSF–ICL 요약의 SHA-256은 `2fbbf6f89c16eea38153f19002b6e14b8b63d28b72d73c6fabd7ab000af5a805`이다. 신규 HSC 자료 요청은 하지 않았다.
