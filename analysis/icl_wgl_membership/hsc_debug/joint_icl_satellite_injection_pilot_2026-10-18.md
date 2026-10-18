# 공동 ICL–위성 외곽광 주입·회수 파일럿

**상태:** 단일 프레임 고정 마스크 전달시험; `science_measurement=false`,
`observational_validation_passed=false`, `flux_covariance_available=false`.

## 목적과 실행

기하학적 중첩 진단에서 기존 matched ICL 표본에 속한 위성 50–75 kpc 고리 162개 중
115개가 ICL 200–300 kpc 고리와 교차할 수 있었다. 그 후속으로 실제 HSC i2 프레임
cluster 149 (`z=0.37556`)의 science 배열을 메모리에서만 복사해 다음 네 조건을
측정했다: 기준, 원형 확산 성분만, 위성 외곽 성분만, 두 성분 동시 주입. 저장된 K2.5
ICL 마스크와 위성 검출·기본 마스크를 고정하고 기존 annulus 측정 함수를 재사용했다.
13개 CAMIRA 위치 대응 검출에 위성 성분을 주입했고, 그중 9개가 50–75 kpc 프로필
및 최소 3개 방위 대조점 기준을 통과했다.

확산 성분은 `A exp[(100−R)/100]`의 형태로 R=100 kpc에서 1 nJy/arcsec²가 되게
했고, 각 위성 외곽 성분은 `A exp[(50−r)/15]`로 r=50 kpc에서 위성당
1 nJy/arcsec²가 되게 했다. 모델은 nJy/arcsec²에서 해당 FITS의 native flux 단위로
변환한 뒤 메모리 영상에만 더했다. 분산은 그대로 두었으며 photon noise는 추가하지
않았다.

## 측정 전달

| 주입 | ICL 200–300 kpc 측정 변화 | 위성 50–75 kpc 차감 측정 변화 |
|---|---:|---:|
| 원형 확산 성분 1단위 | +0.22969 nJy/arcsec² | +0.01660 nJy/arcsec² |
| 13개 위성 외곽 성분, 위성당 1단위 | +0.51930 nJy/arcsec² | +0.28085 nJy/arcsec² |

ICL 추정량은 target-side raw annulus 평균을 사용했다. 합성 신호는 대상 은하단 영상에만
넣고 matched control은 바꾸지 않아, 주입 신호의 target-minus-control 변화와 같다.
위성 값은 9개 지원 위성의 각 target 프로필에서 같은 은하단 거리의 방위 대조점 평균을
뺀 뒤 위성에 동일 가중했다. 동시 주입 결과에서 두 단독 주입 변화 합을 뺀 잔차는 ICL
`2.25×10⁻⁷`, 위성 `8.25×10⁻⁷` nJy/arcsec²로 작았다. 이는 고정 마스크·고정 지원에서
측정기가 선형 합산되는지 확인한 수치 검사다.

단위 주입에 대한 이 한 프레임 응답행렬은 위성 성분이 ICL 고리에, 매끄러운 확산 성분이
위성 target-minus-control 값에 얼마나 전달되는지 보여 준다. 모집단 평균이나 편향 보정값은
아니다. 주입 형태, 한 은하단, 선택된 검출·마스크, 고정 검출, 대각 분산에 조건부인
파일럿이므로 관측 상관을 교정하거나 유의도를 바꾸는 데 사용할 수 없다.

## 감사·재현

- 권위 있는 산출물: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/joint_icl_satellite_injection_2026-10-18_02/cluster_149_summary.json`
- Slurm job **427603**은 00:00:06에 종료 코드 0으로 완료했다. 원본 FITS는 수정되지 않았다.
- 코드: [`joint_icl_satellite_injection.py`](joint_icl_satellite_injection.py)
- 단위·주입 모양 회귀시험 9개, `py_compile`, `git diff --check` 통과.
- 산출물 입력 해시는 원본 HSC 영상, ICL 마스크/짝지음 표, 위성 검출·멤버·기본 마스크,
  준비 요약, 전체 manifest를 기록한다.

첫 실행 폴더 `_2026-10-18_01/`에는 모델의 nJy/arcsec² 값을 native flux에 직접 더한
잘못된 결과가 보존되어 있다. 해당 숫자는 단위 오류이므로 사용하면 안 된다. 단위 변환을
수정한 재실행 `_2026-10-18_02/`만 위 결과의 근거다.

## 다음 검증

다음에는 i/i2와 적색편이·coverage 층을 포함하는 여러 독립 하늘 블록에서 같은 주입을
반복하고, 위성 수·반경 분포별 전달행렬과 공간 블록 분산을 추정해야 한다. 그 뒤 검출과
마스크를 다시 만드는 완전 파이프라인 주입 및 실제 source Poisson noise/PSF 처리가 필요하다.
이 단일 프레임 결과로 `science_measurement` 또는 과학용 공분산 상태를 승격하지 않는다.
