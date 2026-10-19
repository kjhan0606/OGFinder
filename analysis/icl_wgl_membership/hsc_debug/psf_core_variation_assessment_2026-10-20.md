# HSC 위치별 모델 PSF 코어 분포 — 2026-10-20

승인된 PDR3 `pdr3_wide` `HSC-I` coadd 모델 PSF 744개와 별도 승인 시험 PSF 1개를 다시 읽고, 계획 좌표·역할·파일 SHA-256을 하나씩 대조했다. 배치 744개는 은하단 중심 56곳과 위성 주입 위치 688곳이며, 시험 파일을 더해 주입 표본 57개 프레임 중심과 688개 위성 위치, 총 745개 PSF를 확인했다. 새 HSC 요청이나 다운로드는 없었다.

음수 보간 픽셀은 커널 입력에서 유지했다. 아래 크기 지표만 일관된 설명 통계를 얻기 위해 각 stamp의 양수 성분을 다시 정규화해 계산했다. 따라서 `moment FWHM`은 양수 픽셀 2차 모멘트의 가우스 등가 코어 폭이고, 적합한 seeing이나 확장 날개 폭은 아니다.

| 표본 | N | 장축 등가 FWHM 중앙값 [q16, q84] | 단축 등가 FWHM 중앙값 | 축비 중앙값 | R50 중앙값 | R90 중앙값 |
|---|---:|---:|---:|---:|---:|---:|
| 전체 | 745 | 1.084 [1.015, 1.147] arcsec | 1.043 arcsec | 0.970 | 0.363 arcsec | 0.913 arcsec |
| 프레임 중심 | 57 | 1.087 [1.010, 1.153] arcsec | 1.049 arcsec | 0.969 | 0.361 arcsec | 0.917 arcsec |
| 위성 위치 | 688 | 1.084 [1.016, 1.147] arcsec | 1.043 arcsec | 0.970 | 0.363 arcsec | 0.912 arcsec |

전체 장축 폭의 5–95백분위 범위는 0.965–1.198 arcsec이고 전체 최솟값–최댓값은 0.859–1.366 arcsec다. 중심과 위성 위치의 중앙 폭은 비슷하지만, 위치 전체에서 코어 폭이 완전히 일정하지는 않다. 따라서 기존 주입 실험처럼 각 위치의 커널을 사용하는 것이 전역 단일 코어를 쓰는 것보다 입력 관측을 더 충실히 반영한다. 이 분포만으로 스택 편향의 크기를 정하거나 보정계수를 만들 수는 없다.

대부분의 stamp에 음수 픽셀이 있으며 중앙 음수 픽셀 수는 500개다. 전체 양수 픽셀 합/원래 부호 포함 합의 중앙값은 1.0033이다. 양수 부분 기반 지표를 커널에 적용하거나 원 커널을 클리핑·재정규화하면 안 된다. 모델은 41–43픽셀, 0.168 arcsec/pixel로 폭이 약 7 arcsec뿐이므로, 수십–수백 kpc ICL 신호에 필요한 확장 PSF 날개는 여전히 검증되지 않았다. `science_measurement=false`, `psf_approved=false`를 유지한다.

재현기는 [`summarize_hsc_psf_core_variation.py`](summarize_hsc_psf_core_variation.py), 테스트는 [`test_summarize_hsc_psf_core_variation.py`](test_summarize_hsc_psf_core_variation.py)다. 요약은 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/psf_core_variation_2026-10-20_01/summary.json`에 있으며 SHA-256은 `2a1a0f63f7cbabadaf6f1a99f27d53ea9b87d8df2232305c73a89e8ce75807ad`다.
