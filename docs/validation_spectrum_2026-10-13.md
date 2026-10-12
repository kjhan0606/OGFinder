# 2026-10-13 분광 통합 해석 추가 검증

대상은 `/home/kjhan/BACKUP/Astrafek/OGFinder`의 Desktop bridge/순서도와 `/home/kjhan/BACKUP/Astrafek/ds10-web`의 native core/API/Web 순서도입니다. 합성 스펙트럼을 사용한 소프트웨어·수치 검증이며 실제 관측의 과학적 승인 검증은 아닙니다. 모든 계산과 브라우저 테스트는 lageunha의 CPU, OPENBLAS/OMP/MKL 각 1 thread로 실행했습니다. GPU와 HSC 계정/데이터는 사용하지 않았습니다.

## 새로 확인한 범위

| 검증 | 합격 조건 |
|---|---|
| 독립 합성 주입→공동 적합→네 해석 | 초기 z=0 또는 0.35, 참 z=초기 z+0.0008, flux 수치 scale=1 또는 10^-17, 총 네 조합 |
| 적색편이 회수 | 참값과 절대차 ≤2×10^-6 |
| Balmer E(B−V) | 참값과 절대차 ≤0.002 |
| 보정 선 광도 | 참값과 상대차 ≤0.5% |
| 기기 폭 제거 후 넓은 Hα FWHM | 참값과 상대차 ≤0.3% |
| 조건부 log 질량 | 참값과 절대차 ≤0.005 dex |
| 질량의 전체 공분산 오차 | 원본 적합 파라미터·독립 보조입력을 각각 유한차분한 결과와 상대차 ≤10^-5. 외부 소광 및 Balmer 소광 모두 검사 |
| 보고서 운동학 불일치 | common_redshift 변조, NaN redshift, 기준 rest wavelength=0인 세 경우 질량 숫자 없이 blocked |

회수 검증의 flux·거리·소광계수·scatter는 합성 입력으로, 연구 대상에 적용할 기본값이 아닙니다. 선 성분을 사용자가 지정한 뒤 계산하는 경로를 검증했으며 자동 선 식별을 검증한 것은 아닙니다.

## 발견한 결함과 수정

기존 정상 입력 검사는 모두 통과했지만, 새 비정상 보고서 검사의 세 경우에서는 질량이 conditional_mass_estimate로 표시되었습니다. common_redshift가 fitted shift parameter와 일치하는지 검사하지 않았고, 기준파장 0으로 생긴 NaN 오차가 `max(0, variance)`에서 0으로 처리될 수 있었습니다.

`spectrum_inference.py`에서 reference rest wavelength>0, finite redshift/guess와 z>−1, common_redshift=redshift_guess+fitted_shift/reference의 일치, 양의 sigma를 검사하도록 수정했습니다. 오차 gradient와 propagated variance가 finite가 아니면 거부합니다. 부적절한 운동학은 질량 항목을 blocked로 만들며 독립적인 다른 항목은 계속 계산할 수 있습니다. 추가 비정상 보고서 세 경우가 모두 통과했으며 당시 inference 20개 테스트는 경고를 오류로 취급한 실행에서도 통과했습니다. 이후 Balmer 공분산 유한차분 경우를 한 개 추가했습니다.

## 기존 경로 재검증

BPT 값/경계 모호성, 소광·광도·질량 값과 covariance, 낮은 S/N/flags 및 잘못된 가정의 차단을 검사합니다. API 소유권과 source binding, Desktop subprocess·산출물 checkpoint 복원, 실제 브라우저의 inline JSON 입력→공동 적합 자동 연결→네 결과→ZIP 저장→새로고침→재실행도 포함합니다.

BPT 경계와 질량 계수는 [Kauffmann 2003](https://arxiv.org/abs/astro-ph/0304239), [Kewley 2001](https://arxiv.org/abs/astro-ph/0106324), [Greene & Ho 2005 식6](https://arxiv.org/abs/astro-ph/0508335)의 정의를 대조했습니다. 중앙 calibration 계수의 조건부 계산이며 식 계수/virial factor/소광법칙/템플릿 모양 체계오차의 posterior는 계산하지 않습니다.

최종 회귀 결과는 아래에 기록합니다.

| 최종 검사 | 결과 |
|---|---|
| 서버의 분광/워크플로 회귀 | 108 passed (56.94 s) |
| Desktop workflow/bridge/checkpoint | 18 passed (3.46 s) |
| 실제 브라우저 workflow/저장·복원 | 8 passed (28.11 s) |
| 합계 | 134 passed |
| Web typecheck / build | 통과; 615334 bytes, 74 source files |
| Desktop manifest | 39 plugins / 296 steps / 894 parameters, 0 problems |

테스트 명령은 `python3 -m pytest`의 서버 `test_spectrum_*`, `test_workflow_support.py`, `test_science_workflows.py`, Desktop `plugins/workflow_designer/test_workflow.py`, Web `e2e/test_science_workflow_ui.py`입니다. Web은 npm typecheck/build 후 실제 Chromium을 사용했습니다. 과거 126개 결과는 수정 전 회귀 결과이며, 현재 최종 결과는 추가 8개 경우를 포함한 134개입니다. 실제 관측 데이터의 외부 기준 비교·대체 모형·체계오차 검증은 수행하지 않았고 science_validated=false를 유지합니다. 커밋·푸시는 하지 않았습니다.


## 커밋 대상만 추출한 체크아웃 검증

2026-10-13 커밋 준비 시 index의 파일만 lageunha `/tmp/astrafex-commit-check-fdid3h72`에 추출하여 재검증했습니다. 기존 런처·요금제·persona·회절 PSF 등 별도 변경은 제외했고, 파이프라인의 직접 의존인 native morphology/PSF stamp/counts 코드·manifest·테스트는 포함했습니다. 서버 126개(기존 108개+의존 도구 18개), Desktop 18개, 실제 브라우저 8개, 총 152개 통과했습니다. 이 체크아웃의 Web typecheck/build도 통과했습니다(593082 bytes / 73 source files). 서버 테스트 중 zero-noise 합성 PSF 사례에서 기존 noisemodel/detect의 divide RuntimeWarning 3개가 발생했으며 테스트는 통과했습니다. 앞 절의 615334 bytes는 미커밋 런처가 포함된 작업 트리 빌드 크기입니다.

커밋 대상 Desktop manifest는 39 plugins / 295 steps / 894 parameters, 0 problems입니다. 앞 절의 296 steps에는 별도 미커밋 회절 PSF 단계가 포함되어 있습니다.
