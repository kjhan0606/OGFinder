# i2 확장 PSF 프록시와 위치별 PSF 코어 주입 비교 — 2026-10-21

## 목적과 범위

이미 완료된 44개 `i2` 프레임의 공개 HSC-I 확장 PSF 주입 민감도를 같은 프레임의 위치별 HSC PSF 코어 주입과 비교했다. 새 HSC 자료나 영상은 받지 않았고, 기존 영상의 과학 측정값도 바꾸지 않았다. 두 주입 결과의 지원·입력 진폭·프로필 척도가 같은지 확인한 뒤 프레임별 전송 응답을 쌍으로 비교했다.

확장 PSF는 `i2`에 대응하는 위치·epoch별 PSF가 아니라 공개 HSC-I 모델이다. 따라서 이 비교는 필터 특화 PSF 교정이나 관측 검증이 아닌 프록시 모델 민감도다. 프레임 간 분위수는 표본의 기술 통계이며 신뢰구간이 아니다.

## 결과

각 응답은 단위 주입 표면밝기 진폭에 대한 측정/주입 비율이다. paired 열은 같은 프레임에서 `공개 HSC-I 확장 모델 − 위치별 코어` 응답차의 중앙값과 q16–q84다.

| 주입 성분 → 측정량 | 위치별 코어 median [q16, q84] | 공개 확장 프록시 median [q16, q84] | 쌍 차이 median [q16, q84] | 쌍 차이 전체 범위 |
|---|---:|---:|---:|---:|
| diffuse → ICL 200–300 kpc | 0.2245 [0.1834, 0.2519] | 0.2235 [0.1693, 0.2546] | −0.00036 [−0.00092, 0.00093] | −0.0330–0.0576 |
| satellite envelope → ICL 200–300 kpc | 0.4904 [0.1946, 0.8373] | 0.4849 [0.1975, 0.8545] | 0.00846 [−0.00773, 0.02847] | −0.0705–0.0995 |
| diffuse → satellite 50–75 kpc | 0.0156 [−0.1840, 0.1842] | 0.0174 [−0.1837, 0.1740] | −0.000002 [−0.00129, 0.00475] | −0.0617–0.6939 |
| satellite envelope → satellite 50–75 kpc | 0.4094 [−0.3839, 1.2626] | 0.4010 [−0.5457, 1.1109] | 0.00480 [−0.1914, 0.0810] | −3.5513–0.3770 |

공개 대칭·비대칭 템플릿 자체의 차이는 44개 프레임에서 거의 수치 정밀도 수준이었다(응답별 paired 중앙값 절댓값 최대 `3.7×10⁻⁹`, 최대 편차 절댓값 `9.5×10⁻⁷`). 두 템플릿은 같은 공개 별 표본을 사용하므로 독립 검증은 아니다.

## 해석과 한계

ICL 200–300 kpc의 median diffuse 회수율은 두 입력 간 비슷했지만, paired 차이의 꼬리와 위성 측정의 큰 프레임 변동은 사라지지 않았다. 특히 위성 50–75 kpc 회수 응답은 음수에서 큰 양수까지 넓게 분포한다. 이 결과는 안정적인 단일 회수계수나 PSF 보정을 정당화하지 않는다. `i2`의 ICL–위성 연관성을 과학 측정으로 승격하는 근거도 아니다.

합성광에만 PSF를 적용했고 주입원 광자 잡음, coadd 픽셀 공분산, 상류 하늘 제거 전달, 위치·epoch별 PSF 날개 및 관측 완전도를 검증하지 않았다. 따라서 `science_measurement=false`, `psf_approved=false`, `observational_validation_passed=false`, `flux_covariance_available=false`를 유지한다.

## 재현과 산출물

재현 코드는 [`audit_extended_psf_core_vs_wide.py`](audit_extended_psf_core_vs_wide.py), 테스트는 [`test_audit_extended_psf_core_vs_wide.py`](test_audit_extended_psf_core_vs_wide.py)다. 새 필터 선택 경로까지 포함한 focused 테스트 3개가 통과했다.

권위 출력은 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_extended_psf_i2_proxy_variant_comparison_2026-10-21_01/core_vs_wide_audit.json`이며 SHA-256은 `ba1287493d5d3db4d8050f601ac982a659e1ca8588c024c90d8a7814c30e0f97`이다. 분석기 SHA-256은 `8d08096c9d4babab32aa61fc6bef34d710267e65873c7cc1296ac3c36009f625`다. 비교는 44/44 프레임을 사용했다.

입력 요약 SHA-256: 위치별 코어 `df53035c6066bb4215c0d7d2aec9fef423e62a52bb6c60b89eb5be0e9d17ccc1`; 공개 대칭 프록시 `5f2e86d6c3499148b1dbb3151e261daf827de5c434b2b32da7c9af1be846c59d`; 공개 비대칭 프록시 `03bdf8f12fde0f9e84ab28ae6c39819f4129564aab0dd4e0c3bb3c1e0340194e`.

두 44프레임 입력 배치는 독립 산술 검사에서 모두 최대 저장/재계산 차이 0으로 검증됐다. 대칭 검증 JSON SHA-256은 `31b8443bbc46ff8c57c711eeb1ac886f21ca63e0056631f535c002a05dd8e5f5`, 비대칭은 `2e7c32c6cb0e68a6aa69f25d272d42c4f1d683d1a204ee48b15ed1e5e9feea2b`다. 대칭·비대칭 paired 비교 파일은 `17ab4d0db5ac73ce4b772a1ef0dc176f2796d0bdecf6f4f4b114188cec70526d`다.
