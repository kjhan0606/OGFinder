# HSC 디버그 절차

이 디렉터리의 스크립트는 ICL 스케일 스택을 실제 목록으로 준비하는 절차이다. 논문에 함께 제공한다. 실행 결과는 과학 측정이 아니다. `science_measurement` 는 거짓으로 남는다.

작업 디렉터리는 `/gpfs/kjhan/Astrafek/data` 이다. 다른 경로를 쓰려면 `ASTRAFEK_DATA` 를 둔다. OGFinder 위치는 `ASTRAFEX_OGF_ROOT` 이다.

```bash
bash scripts/run_debug.sh
```

순서:

1. `01_fetch_catalog.sh` 는 CAMIRA S20A Wide 별마스크 목록을 받는다. SHA-256 이 기록된 값과 같으면 다시 받지 않는다.
2. `02_select_clusters.py` 는 분광 적색이동 창을 통과한 가장 촘촘한 0.05 적색이동 구간의 은하단을 고르고 `logs/selected_clusters.tsv` 를 쓴다. 6각분보다 가까운 쌍은 richness 가 큰 쪽만 남긴다.
3. `03_mass_gate.py` 는 목록의 질량이 `M_200c` 가 아니라는 점(Oguri et al. 2018, arXiv:1701.00818 은 `M_200m`)을 기록한다. `clusters.json` 을 만들지 않는다.
4. `04_cutout_status.py` 는 `cutouts/` 안의 FITS 유무만 적는다. 영상을 내려받지 않는다. HSC 영상은 로그인 뒤 아카이브에서 가져온다.
5. `05_stack.py` 는 `masses.tsv` 와 `cutouts/frames.tsv` 가 있을 때 메뉴와 같은 명령으로 `icl_stack.py` 를 실행한다. 둘 중 하나가 없으면 스택을 실행하지 않는다.
6. `06_fetch_members.py` 는 CAMIRA S20A Wide 별마스크 멤버 목록을 받아 고른 은하단에 맞춘다. 파일에는 머리글이 없다. 멤버 행은 모두 남기고, 가중치 하한을 두지 않는다. BCG 위치는 은하단 목록의 좌표이다.
7. `07_fetch_cutouts.py` 는 unagi 가 적는 PDR3 컷아웃 주소로 HSC-I `type=coadd` 한 장(20각초)을 확인한다. 계정은 환경 변수 `HSC_USER`, `HSC_PASSWORD` 가 있을 때만 붙인다. 응답이 FITS 가 아니면 다운로드를 멈춘다. FITS 이면 첫 2차원 HDU 의 헤더만 남기고, global-sky 코애드로 확인되기 전에는 나머지 은하단을 받지 않는다.
8. `08_mask_stack.py` 는 `cutouts/frames.tsv` 의 sky 가 `global` 또는 `global-sky` 이고 product 가 비어 있거나 `deepCoadd` 또는 `coadd/bg` 일 때, BCG와 멤버 발자국을 가린 뒤 `warp_metric` 과 `mean_stack` 으로 kpc 평균을 만든다. 산출물 이름은 `bcg-and-members-masked` 이다. `M_200c` 를 쓰지 않고 `icl_stack.py` 를 호출하지 않는다. SEP 은 목록 좌표의 픽셀 발자국만 잡고, 맞추지 못한 검출은 표에만 남긴다.

2026-10-12 실행에서 멤버 파일은 고른 213개 은하단에 모두 맞춰졌고, 멤버 행은 9607개이다. 그 가운데 205개는 은하단 목록 좌표와 같은 멤버 행이 있다. 컷아웃 확인의 HTTP 상태는 401이고 FITS 는 받지 못했다. `cutouts/frames.tsv` 가 없어 가린 평균은 쓰이지 않았다.

`masses.tsv` 의 머리글은 `id`, `m200c_msun`, `mass_definition`, `citation` 이다. `mass_definition` 값은 `200c` 이다. 식은 그 인용에서 그대로 옮긴다.

`cutouts/frames.tsv` 의 머리글은 `id`, `filename`, `sky`, `product` 이다. `sky` 는 `global` 또는 `global-sky` 이다. `product` 는 비우거나 `deepCoadd` 이다. 중심은 astropy `world_to_pixel` 의 0-기반 픽셀이고, 픽셀 크기는 헤더의 CD 또는 CDELT 이다.

2026-10-12 인증 프로브 후속: `probe_login.py --global-sky`가 HSC-I
`coadd/bg` 요청에 HTTP 200을 받았다. 저장한 공식 매뉴얼은 이를
`deepCoadd`의 coaddition 직후, aggressive background subtraction 이전
산출물로 연결한다. 전체 FITS 헤더에서 primary의 `FLUXMAG0=63095734448.01944`
(AB 영점 27), IMAGE/MASK/VARIANCE HDU 1/2/3을 확인했다. 이번 헤더 전용
프로브는 영상 픽셀을 남기지 않았고 다른 은하단 요청도 하지 않았다.
교정·제품 출처 검토는 `/gpfs/kjhan/Astrafek/data/logs/cutout_global_review.json`이다.
이 한 건의 제품 출처 관문은 통과했으나 sky 잔류량·PSF·마스크·random·주입/회수
검증과 12개 시범 영상은 미완료이다. 이 기록은 위의 초기 401 이력을 대체하지 않는다.

2026-10-19 사용자 승인으로 기존 주입 위치 744곳의 HSC PDR3 `HSC-I`
`coadd` 모델 PSF 코어를 추가 수신하고 좌표·FITS·정규화·중심 품질을
검증했다. 절차는 `prepare_hsc_psf_picker_plan.py`와
`audit_hsc_psf_picker_batch.py`, 결과와 한계는
[`psf_picker_batch_2026-10-19.md`](psf_picker_batch_2026-10-19.md)에 있다.
PSF 코어 입력은 준비됐지만 확장 날개와 영상 노이즈는 별도로 검증해야 하며,
과학 측정은 수행하지 않았다.

같은 승인 PSF를 실제 위치별 주입에 연결한 재검출 민감도도 57/57 프레임에서
완료했다(i 13, i2 44, 오류 0). PSF 코어와 무PSF 실행의 짝지은 응답 차이는
실제 관측 coadd 픽셀 배열의 복사본에서 얻었고, 각 필터에서 중앙값이 작으며
16–84백분위가 0을 포함한다. 해당 57개 target FITS 헤더에는 `FLUXMAG0`와
variance plane은 있지만 `GAIN`·`EXPTIME`가 없어 source photon variance를
정량 추가할 수 없다. 확장 PSF 날개, 상관 coadd noise, HSC 상류 sky 전달도
별도 검증이 필요하며 science measurement로 승격하지 않았다. 이 헤더 검사는
[`audit_icl_coadd_noise_metadata.py`](audit_icl_coadd_noise_metadata.py)로 재현하고
[`icl_coadd_noise_metadata_assessment_2026-10-19.md`](icl_coadd_noise_metadata_assessment_2026-10-19.md)에 기록했다.
PSF 반응 상세·SHA-256은
[`psf_convolved_joint_injection_assessment_2026-10-19.md`](psf_convolved_joint_injection_assessment_2026-10-19.md),
재현 명령은 [`run_joint_icl_satellite_psf_redetection.sbatch`](run_joint_icl_satellite_psf_redetection.sbatch),
감사기는 [`audit_joint_icl_satellite_psf_redetection.py`](audit_joint_icl_satellite_psf_redetection.py)에 있다.

2026-10-12 수집 확장: 이후 12개 FITS 다운로드가 모두 완료됐으며, 파일 크기와
SHA-256, FITS 구조와 WCS를 재검사했다. 이제 213개 후보 전체를 받도록
`hsc_pilot.py prepare --all-selected`로 목록을 만들고 기존 `download_pilot.py`의
입력 제한을 확장했다.
기존 12개를 재사용하고 추가 201개(약 10.28 GiB)를 요청한다.
계획은 `../implementation_plan.md`의 첫 절, 실행 데이터는
`/gpfs/kjhan/Astrafek/data/sample_2026-10-12`이다.

별도로 만들었던 제출 스크립트는 Slurm 환경 전달 오류로 철회했다. 기존에
성공했던 다운로드·비밀번호 입력 경로를 그대로 사용한다. OGFinder 루트의
SSH 터미널에서 아래 명령을 실행하고 노드 배정 뒤 비밀번호를 입력한다.
`/`를 포함한 특수문자는 그대로 입력한다.

```bash
srun -p a10 -N1 -n1 --cpus-per-task=1 --mem=4G --gres=gpu:1 --time=12:00:00 --pty \
  /home/kjhan/miniconda3/bin/python3 analysis/icl_wgl_membership/hsc_debug/download_pilot.py \
  --user kjhan0606 --resume \
  --manifest /gpfs/kjhan/Astrafek/data/sample_2026-10-12/sample_manifest.json \
  --work /gpfs/kjhan/Astrafek/data/sample_2026-10-12/download_01
```

오류가 나면 완료분을 보존하고 멈추며 같은 명령으로 이어받는다. 터미널을
열어 둔 채 실행하고 중복 실행하지 않는다. 실제 진행 상태는 실행 디렉터리의
`download_01/download_status.json`, 전체 완료 목록은 그 디렉터리의
`download_manifest.json`이다. 상위 디렉터리의 이전 제출 상태는 과거 이력이다.
비밀번호는 기존 프로그램의 메모리에만 둔다. GPU 연산은 없지만 장시간 작업은
Slurm에서 실행하고 GPU를 요청해야 한다는 서버 규정을 따른다.

## 2026-10-18 CAMIRA 멤버 밀도–ICL 프로파일 예비 비교

`analyze_satellite_association_pilot.py`는 이미 저장된 50개 검증 은하단 프로필과
대응 대조장 프로필, 800개 표본 manifest의 CAMIRA 멤버 위치를 해시 검사 후
연결한다. 지원 민감도에서 이미 선택된 K2.5/40%/50–300 kpc 24쌍만 사용하고,
클러스터와 대조장의 필터가 정확히 일치하는지 다시 확인한다. 결과를 기존 자료에
덮어쓰지 않도록 새 `--work` 경로를 지정한다.

```bash
python3 analysis/icl_wgl_membership/hsc_debug/analyze_satellite_association_pilot.py \
  --data /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13 \
  --work /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_association_pilot_2026-10-18_02
```

출력 `cluster_metrics.tsv`는 은하단별 멤버 위치 수밀도와 matched-control 차감
표면밝기를 기록하고, `summary.json`은 필터별 탐색 상관 및 모든 제한사항을
기록한다. 멤버 광도·별질량·멤버십 확률은 사용하지 않았고 위성 중심 외곽광은
측정하지 않았으므로 `science_measurement=false`가 유지된다. 이 분석은 이미
생성된 작은 프로필 배열을 읽는 진단이며 800개 이미지 전체의 광도 스택을 만들지
않는다.

## 2026-10-18 전체 800개 대상 영상 프로필

`measure_full_sample_profiles.py`는 다운로드된 800개 영상과
`prestack_cluster_masks_2026-10-16_01`의 마스크를 읽어 반경별 HSC-I 프로필을
만들고, 저장된 798개 대조장 프로필과 동일한 필터·유효면적 조건으로 차감 평균과
5° 공동 하늘 블록 공분산을 계산한다. 전체 이미지 입력이므로 Slurm 작업으로만
실행하며, 기존 대조장 프로필을 재사용한다. 새 FITS 다운로드는 없다.

```bash
sbatch analysis/icl_wgl_membership/hsc_debug/run_measure_full_sample_profiles.sh
```

`prestack_flux_profiles_2026-10-18_02/`에는 대상별 프로필, 지원 절단별 i/i2
차감 프로필, 공분산, CAMIRA 멤버 위치 수밀도와 탐색 상관이 저장된다. 자동
마스크·BCG 검토, extinction 및 passband 교정, 상류 sky/PSF 전달은 승인되지 않은
상태이므로 산출물은 `science_measurement=false`인 탐색 결과다.
정량·정성 품질 판정은 [`../stack_quality_assessment_2026-10-18.md`](../stack_quality_assessment_2026-10-18.md)에 기록했다.

첫 제출의 `_01` 경로는 800개 대상 프로필을 읽은 뒤 요약 단계의 대조장 필터
열 이름 오류로 실패했다. 해당 부분 산출물은 실패 기록으로 보존하며 분석에
사용하지 않는다. 수정된 재실행은 새 `_02` 경로를 사용한다.

## 2026-10-18 위성 위치 밀도 상관의 공간 블록 재표본화

`analyze_full_sample_satellite_association.py`는 전체 프로필의 기존 137개 지원쌍을
필터별로 나누고, 대상·대조장 중 어느 쪽이든 5도 하늘 셀을 공유하는 쌍을 같은
블록으로 묶어 5,000회 재표본화한다. 반경별 CAMIRA 멤버 위치 밀도와 차감광의
Spearman 상관 및 richness·적색편이를 제거한 Pearson 상관의 95% 구간을 기록한다.
양의 재표본 비율은 상관의 공간 안정도 요약이지 영가설 검정 p 값이 아니다.

```bash
python3 analysis/icl_wgl_membership/hsc_debug/analyze_full_sample_satellite_association.py \
  --work /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_association_full_2026-10-18_01
```

출력은 `block_bootstrap_metrics.tsv`와 `summary.json`이다. 입력 스택과 지원쌍의
기록 해시를 확인하며, 멤버 광도·질량·확률 및 위성 중심 외곽광은 이 단계에도
없으므로 `science_measurement=false`를 유지한다.

## 2026-10-18 약한렌즈 shape-catalog 겹침 사전계산

`prepare_wl_shape_preflight.py`는 800개 은하단의 3 Mpc 원뿔에 대해 HSC S19A
shape source 개수만 확인할 CAS SQL 배치를 만든다. 실행 시 네트워크 조회나 자료
다운로드를 하지 않으며, 준비된 쿼리도 아직 제출하지 않았다.

```bash
python3 analysis/icl_wgl_membership/hsc_debug/prepare_wl_shape_preflight.py \
  --manifest /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/sample_manifest.json \
  --work /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_shape_preflight_2026-10-18_01
```

32개 SQL 배치와 `preflight_manifest.json`이 생성된다. 렌즈별 source 행 수와 예상
전송량을 확인한 뒤 자료 전송 여부를 정하고, 사전계산 단계에서는 shape 행을
가져오지 않는다.

검증 및 count-only 실행기는 `submit_wl_shape_counts.py`다. `--plan`은 SQL 해시와
800개 렌즈 ID를 검증하며 네트워크를 사용하지 않는다. `--submit`은 `lageunha`의
대화형 터미널에서 한 번 숨김 비밀번호를 입력받아 32개 배치를 순차 실행하고,
렌즈당 정수 count 한 개씩만 저장한다. shape·좌표·photo-z 행은 반환하지 않는다.

```bash
python3 analysis/icl_wgl_membership/hsc_debug/submit_wl_shape_counts.py --plan
python3 analysis/icl_wgl_membership/hsc_debug/submit_wl_shape_counts.py --submit
```

기존 제출 결과는 `<preflight>/cas_counts/`에서 재개·해시 검증한다. count 합계는
중복 원천 은하가 여러 렌즈 원뿔에 들면 여러 lens-source 후보로 세는 기획용 수이며,
선택 후 실제 source 수나 다운로드 바이트가 아니다. shape source 행의 반출은 별도
용량 산정과 사용자 승인을 거친다.

2026-10-18에 사용자 승인으로 count-only 사전 조회를 수행했다. 32개 배치에서
렌즈별 count 800행을 받아 고유 렌즈 ID, 행 수, 배치 SHA-256을 검증했다. 렌즈별
count 합은 3,069,544이고 0개 count 렌즈는 168개다. 이 합은 겹치는 원천 은하를
렌즈마다 중복 집계한 값이며, source 품질·photo-z 선택 전 후보량이다. 저장 위치는
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_shape_preflight_2026-10-18_01/cas_counts/`.
shape·좌표·photo-z 원자료는 다운로드하지 않았으며 다음 원자료 조회에는 별도 승인과
전송량 산정이 필요하다.

## S19A shape pair stack adapter (2026-10-18)

`s19a_shape_stack.py`는 다운로드기가 아니다. 승인된 shape source 행과 보정 redshift
입력을 받은 뒤 렌즈–source 쌍의 접선·교차 전단 및 약한 전단 근사 `DeltaSigma`
프로필을 계산하는 함수형 adapter다. full per-object photo-z PDF 방식과 HSC Y3
tomographic `n(z)` 방식(`stack_s19a_tomographic_pairs`)을 모두 지원한다. 후자는
계획된 shape-row의 `source_zbin`을 직접 사용하고 full-PDF 대용량 tar를 요구하지 않는다.
`hsc_y3_zbin==0`은 제외하고 1..4를 보정 `n(z)` 행에 연결한다. `object_id`는 정수 또는 문자열로
정확히 보존해야 하며 float로 읽으면 거부한다. 같은 source가 여러 렌즈 원뿔에
들어가는 것은 허용하지만 같은 `(lens_id, object_id)` 쌍의 중복은 거부한다.

입력 `pairs` 열의 매핑은 다음과 같다.

| adapter 열 | HSC/분석 입력 |
|---|---|
| `e1`, `e2` | `i_hsmshaperegauss_e1/e2` |
| `sigma_e` | `i_hsmshaperegauss_derived_sigma_e` |
| `weight` | `i_hsmshaperegauss_derived_weight` |
| `rms_e` | `i_hsmshaperegauss_derived_rms_e` |
| `m`, `c1`, `c2` | `i_hsmshaperegauss_derived_shear_bias_m/c1/c2` |
| `b_mode_pass` | HSC `b_mode_mask`를 검증된 값 방향으로 decode한 통과 boolean |
| `mag_a10`, `resolution` | `i_apertureflux_10_mag`, `i_hsmshaperegauss_resolution` |
| `psf_shape11/22/12` | `i_sdssshape_psf_shape11/22/12` |
| 좌표·렌즈 적색이동 | 800개 렌즈 목록 및 `meas`의 천구좌표 |
| source redshift | 기본: `source_zbin` 1..4를 HSC Y3 SACC의 fiducial `n(z)` 곡선에 대응; 대체: `object_id`별 tract full PDF |

호출 함수는 `pdf_lensing_efficiency`로 각 lens-source 쌍의 PDF 적분
`<D_ls/D_s>`를 구하고, 각 반경 빈에서
`w_pair = shape_weight * SigmaCritInv**2`로
`R = 1 - sum(w_pair * e_RMS**2) / sum(w_pair)`를 다시 계산한다. per-object `m/c`와
HSC의 경계 밀도 선택편향 처방을 적용한다. 선택편향에는 반경 빈별
`shape_weight * SigmaCritInv**2` 가중치, 25.5 mag·`R2=0.3` 경계의 가중 밀도,
가중 평균 PSF ellipticity를 사용해 `m_sel` 및 `c_sel=a_sel*<e_PSF>`를 계산한다.
계수와 오차식은 HSC가 연결한
[hsc-y3-calib 보정 코드](https://github.com/mr-superonion/hsc-y3-calib)의
`get_sel_bias` 처방을 근거로 구현했으며, 저장소의 코드를 복사하거나 포함하지 않았다.
입력은 HSC shape catalog 행이어야 한다. adapter는 HSC Y3 논문의 선택 조건
`|e|<2`, `0<=sigma_e<=0.4`, `resolution>=0.3`, `mag_a10<=25.5`를 스택 전에
적용하고 각 shape/측광 선택 탈락 건수를 기록한다. 이전 어댑터는 수학적 ellipticity
상한 1을 강제했지만, HSC 논문이 실제 catalog 선택으로 정한 상한은 2다. B-mode mask도
기본 적용하며, 명시적 민감도 실행에서만 끌 수 있다.

HSC PDR3 문서는 `b_mode_mask`가 큰 B-mode 영역을 제거하는 데 쓰인다고 설명하지만,
저장값의 0/1 극성은 명시하지 않는다. Astrafek은 HSC Y3 PSF-systematics 논문에
실린 GAMA09H 제외 사각형과 승인된 shape-row 좌표를 교차해 이 릴리스에서 raw 0이
제외 영역임을 확인했다. 해석 함수는 여전히 `masked_value`를 필수로 받아 자동 추측을
막는다. 검증과 이후 stack은 아래 B-mode 검증 절의 보고서를 참조한다.

공식 HSC Y3 SACC 보정 파일을 로컬에 받은 뒤 `load_hsc_y3_nz_sacc(path)`로
`wl_0`–`wl_3` tracer를 읽어 공통 `z` 격자와 네 개 `n(z)` 배열을 만들 수 있다.
이 함수는 로컬 파일만 읽고 네트워크에 접속하지 않는다. tracer 곡선은 fiducial
`n(z)` 입력이며 posterior 표본이나 그 불확실성 전파를 대신하지 않는다. 실제 파일이
확보되면 다음처럼 adapter에 전달한다.

```python
from s19a_shape_stack import load_hsc_y3_nz_sacc, stack_s19a_tomographic_pairs

z_nz, nz_by_bin = load_hsc_y3_nz_sacc("/path/to/hsc_y3_fourier_space_data_vector.sacc")
result = stack_s19a_tomographic_pairs(pairs, z_nz, nz_by_bin, radial_edges_mpc, cosmology)
```

현재 결과는 weak-shear 근사 프로필과 독립 쌍으로 가정한 shape-noise 오차, cross-shear
null만 낸다. full P(z) 방식은 PDF와 기하 효율을 적분하고, tomographic 방식은 HSC Y3
보정 `n(z)`를 적분한다. `load_hsc_y3_nz_sacc`가 공식 SACC의 네 `wl_0`–`wl_3` tracer를
검증해 연결하며, 이 경로는 계획된 source-row export와 호환된다. fiducial `n(z)`에 대한
posterior/systematic 불확실성 전파는 아직 없으며, 개별 source PDF 방식은 HSC가 요구하는 photo-z bias/dilution
보정도 추가해야 한다. reduced shear,
random-lens 차감, 겹치는 source/하늘 블록 공분산, 선택 보정 계수의 오차 전파도
최종 공분산에 결합하지 않았다. 따라서 출력은 항상 `science_measurement=false`이고
최종 질량·과학 결과로 쓰지 않는다. 합성 주입 회수와 잘못된 ID/PDF/mask/selection 입력
거부 시험은 `test_s19a_shape_stack.py`에 있다. HSC 열·tomographic-bin·PDF 연결 규약은
[공식 S19A 안내](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s19a-shape-catalog-pdr3/),
응답도·per-object bias·selection correction은
[HSC 3-year shape catalog 방법론](https://arxiv.org/abs/2107.00136)을 따른다.

실제 row export와 SACC `n(z)` 파일이 준비되면 `run_wl_shape_stack.py`가 이 adapter를
연결한다. 후처리기는 네트워크를 쓰지 않으며 32개 row-batch 해시, 렌즈별 행 수,
800개 렌즈, count-only 계획 연결을 다시 검증한다. HSC `b_mode_mask`에서 제외할 저장값을
`--b-mode-masked-value`로 명시해야 하며, 현재 승인된 S19A 내보내기에서는 검증된 값 `0`을
사용한다. 다른 릴리스나 내보내기에 자동 일반화하지 않는다. 예시:

```bash
python3 analysis/icl_wgl_membership/hsc_debug/run_wl_shape_stack.py \
  --rows /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_shape_rows_export_2026-10-18_01 \
  --nz-sacc /path/to/hsc_y3_fourier_space_data_vector.sacc \
  --b-mode-masked-value 0 \
  --radial-edges-mpc 0.1 0.2 0.4 0.7 1.0 1.5 2.0 3.0 \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_cluster_stack_01.json
```

예시는 입력 경로와 마스크 극성을 검증한 뒤에만 실행한다. 산출물은 cluster-only 탐색
프로필이며 random-lens 차감, `n(z)` 오차 전파, paired spatial jackknife 공분산이 빠진
상태로 기록된다. source-row export 자체는 별도 명시 승인이 있을 때만 실행한다.

별도 후처리 `subtract_random_lens_profile`은 cluster와 random stack의 접선·교차
`DeltaSigma` 점 추정량을 bin별로 차감한다. `paired_jackknife_covariance`는 양쪽에
동일하게 적용한 spatial delete-one-block replicate의 차이에서 공분산을 계산한다.
이들은 실제 random source 쌍을 만들거나 각 jackknife replicate를 실행하지 않는다.
공통 sky block 및 source overlap을 반영한 실제 공분산이 없으면 두 stack의 개별
shape-noise 오차를 독립으로 합치지 않는다. 테스트에는 신호 회수·n(z) bin 처리·SACC
로딩/검증·선택 컷·무작위 차감 및 공분산 산식이 포함되며 현재 22개가 통과한다.

`paired_tomographic_spatial_jackknife` now runs that missing pair-level step: it
assigns source RA/Dec to the same equal-area RA/sin(Dec) cells in lens and random
pair tables, deletes each cell from both, and recomputes responsivity and edge
selection calibration for every replicate. It returns the full-sample
random-subtracted profile plus the joint tangential/cross covariance. This
calculation is not yet run because the random-source pairs have not been
exported. Its grid is an exploratory fallback and does not guarantee contiguous
independent HSC regions; region-based jackknife remains a science gate.

`prepare_wl_random_shape_count_plan.py` prepares (without submitting) a count-only
S19A support query for a finalized `astrafek-wl-random-centers-v1` manifest. It
retains the association between each unique random-center ID and its matched
cluster ID, enforces a common target redshift per matched cluster, and counts
distinct shape objects within that target's proper 3 Mpc aperture. The result
contains no source shape rows and no images. Its SQL/count plan must be reviewed
for a separate explicit approval before execution; the eventual source-row
export will require another approval after observed counts and transfer size
are known. `submit_wl_random_shape_counts.py --check` verifies the manifest and
all exact SQL text without network access. Its `--submit` mode also requires
`--approved-random-shape-count-query`; this guarded submitter has not been run.

## CAMIRA member-mass and weight sensitivity (2026-10-18)

`analyze_member_weighted_satellite_profiles.py` joins existing per-annulus satellite
profiles to the 800-cluster CAMIRA member file using saved member coordinates. It
verifies all 800 cluster identities and 36,513 member-position sets, then audits
5,832 profile rows in 543 clusters. It compares equal-cluster unweighted profiles,
CAMIRA-weighted profiles, and candidate stellar-mass bins. The source file is
headerless. The hash-matched acquisition summary identifies column 8 as
`weight`; column 7 is the undocumented `field7`, whose values are only treated
as candidate `log10(M*/Msun)`. The output therefore remains descriptive and
keeps `science_measurement=false`.

```bash
python3 analysis/icl_wgl_membership/hsc_debug/analyze_member_weighted_satellite_profiles.py \
  --work /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_member_weight_sensitivity_2026-10-18_04
```

See [`member_weighted_profile_audit_2026-10-18.md`](member_weighted_profile_audit_2026-10-18.md)
for results and limitations. Regression tests are in
`test_analyze_member_weighted_satellite_profiles.py`.

`analyze_weighted_satellite_icl_association.py` then compares the saved
50–75 kpc satellite profiles with the K2.5/40% 200–300 kpc ICL sample, using
equal-satellite versus CAMIRA-weighted means and a paired spatial-block
bootstrap. It is a post hoc sensitivity analysis and keeps
`science_measurement=false`. See
[`camira_weighted_icl_association_2026-10-18.md`](camira_weighted_icl_association_2026-10-18.md).

## ICL annulus sample-selection audit (2026-10-18)

`audit_icl_sample_selection.py` compares the 137-cluster K2.5/40% all-annulus sample
and the 224-cluster 200–300 kpc annulus-only sample with the 800-cluster parent,
within each filter. It reports descriptive z/richness shifts, central masking,
center-candidate fractions, and 5-degree footprint occupancy. It does not correct
selection or authorize a science measurement. See
[`icl_sample_selection_audit_2026-10-18.md`](icl_sample_selection_audit_2026-10-18.md)
and the hash-recorded output at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_sample_selection_audit_2026-10-18_02/`.

`audit_icl_center_quality.py` then separates the 224 outer-annulus pairs by the
automatic unique-center and K2.5 central-mask flags, using connected spatial-block
bootstrap intervals. The result is filter-dependent and remains exploratory; see
[`icl_center_quality_sensitivity_2026-10-18.md`](icl_center_quality_sensitivity_2026-10-18.md)
and `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_center_quality_sensitivity_2026-10-18_02/`.

`audit_satellite_icl_pixel_overlap.py` refines the circular-annulus geometry check by measuring
the actual shared valid pixels in the saved K2.5 ICL and SEP/Kron satellite masks. It can run on
one or more matched clusters with `--cluster-id`; the complete supported sample is 63 clusters
and 162 satellite profiles. It reads existing images only and keeps
`science_measurement=false`. See
[`satellite_icl_pixel_overlap_2026-10-18.md`](satellite_icl_pixel_overlap_2026-10-18.md) and
the hash-recorded output at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_icl_pixel_overlap_2026-10-18_02/`.

`audit_satellite_icl_shared_pixel_covariance.py` reconstructs the accepted 50–75 kpc satellite
target and azimuthal-control annuli, then contracts their signed equal-satellite weights with
the ICL target annulus and the HSC diagonal variance map. It verifies saved target/control
support counts and reports only a nominal independent-pixel noise covariance; coadd-correlated
noise and the ICL random-control term remain absent. It does not correct the measured relation.
See [`satellite_icl_shared_pixel_covariance_2026-10-18.md`](satellite_icl_shared_pixel_covariance_2026-10-18.md)
and `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_icl_shared_pixel_covariance_2026-10-18_01/`.

`audit_satellite_icl_leave_block_out.py` recalculates the redshift/richness-adjusted partial
correlation after omitting each connected 5-degree target/control sky block in turn. This
measures block influence only; it does not make a significance claim or change the science gate.
See [`satellite_icl_lobo_assessment_2026-10-19.md`](satellite_icl_lobo_assessment_2026-10-19.md)
and `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_icl_lobo_2026-10-19_03/`.
The deterministic leave-block-out table was rerun with the same input and verified byte-for-byte
at `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_icl_lobo_2026-10-19_04/`.

`audit_satellite_icl_spatial_blocks.py` restores target/control coordinates for the 94 rows in the
per-annulus association using the hash-verified ICL profile and matched-control extraction records.
It compares 5,000 ordinary pair-bootstrap draws with 5,000 connected 5-degree target/control block
draws. For the primary i2 all-supported subset (68 pairs, 33 blocks), the adjusted partial-Pearson
estimate is 0.293 with spatial-block 95% interval [0.042, 0.497]; the ordinary-pair interval
includes zero, and the unadjusted Spearman spatial interval also includes zero. This remains a
selected-sample exploratory result, not a detection. See
[`satellite_icl_spatial_block_assessment_2026-10-19.md`](satellite_icl_spatial_block_assessment_2026-10-19.md)
and `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_icl_spatial_blocks_2026-10-20_02/`.

`plan_icl_spatial_block_power.py` converts that exploratory effect size into a
design approximation using the 33 connected spatial blocks as independent
units, two adjustment variables, and the 30-test family. For the i2 partial
correlation scenario `r=0.293`, estimated power is 0.359 at nominal α=0.05 and
0.061 at Bonferroni α=0.05/30. The Fisher-z approximation requires 180
independent blocks for 80% power at the adjusted threshold. These are not 180
clusters: additional clusters inside the same connected sky blocks do not meet
the independent-block requirement. This is a planning bound, not a detection
or a replacement for the connected-block maxT analysis. The reproducible
report and JSON are at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_association_power_plan_2026-10-19_01/`.

`analyze_extinction_sensitivity_icl_satellite.py` applies center-based SFD extinction corrections
separately to the saved ICL targets and matched controls, and compares them with a target-only
HSC PDR3 forced-source `a_i` sensitivity. It verifies the existing profile/catalogue hashes and
uses 5,000 connected 5-degree block-bootstrap draws for the 200–300 kpc ICL means and the
200–300 kpc ICL / 50–75 kpc satellite association. Results are exploratory because extinction is
not corrected per pixel and other sky/PSF/mask/injection gates remain open. See
[`icl_extinction_sensitivity_2026-10-19.md`](icl_extinction_sensitivity_2026-10-19.md) and
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_extinction_sensitivity_2026-10-19_03/`.

`measure_pixel_extinction_profiles.py` checks the center-based correction against separate
pixel-level SFD corrections for target and matched-control images; local sky-plane subtraction is
reported only as a sensitivity test. The completed 224-pair profile and 74-pair
satellite-association audit remains exploratory. See
[`icl_pixel_extinction_assessment_2026-10-19.md`](icl_pixel_extinction_assessment_2026-10-19.md)
and `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_pixel_extinction_2026-10-19_03/`.

`audit_icl_satellite_sample_retention.py` now audits why 398 candidates that passed the
target-side gate failed the complete seven-annulus control gate. For 396, eight azimuthal
positions were generated; the attrition occurs when a control must pass the 60% coverage floor
in every annulus. Of these failures, 102 still have at least three controls in the primary
50–75 kpc annulus and match exactly to the existing independent-annulus profiles (25 i, 77 i2).
The complete-profile selection remains the reference; no threshold was relaxed. See
[`icl_satellite_control_gate_assessment_2026-10-19.md`](icl_satellite_control_gate_assessment_2026-10-19.md)
and `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_sample_retention_2026-10-19_02/`.

The approved HSC weak-lensing source-row export is complete for 800 lenses and 3,069,544
lens–source rows. It contains 22 columns, excludes full photo-z PDF rows, and remains
`science_measurement=false` pending random-lens subtraction and the remaining selection and
covariance gates. Its row summary is at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_shape_rows_export_2026-10-18_03/row_summary.json`.
The streaming quality audit finds 88,460 rows with `1<|e|<2` and none with `|e|>=2`,
so the adapter was corrected to follow the published `|e|<2` and `sigma_e<=0.4`
selection instead of rejecting the whole input at `|e|>=1`. The four published shape
cuts plus tomographic assignment leave 2,234,510 candidate pairs before B-mode and
radial selection. The raw B-mode bit was later decoded by coordinate matching to the
published HSC Y3 excluded footprint; see the validation record below. See
[`wl_shape_row_quality_assessment_2026-10-19.md`](wl_shape_row_quality_assessment_2026-10-19.md)
and `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_shape_row_quality_2026-10-19_03/`.

## Raw tangential-shape diagnostic (2026-10-19)

`audit_wl_raw_shape_profiles.py` calculates shape-weighted mean raw `e_t` and
cross `e_x` from the approved row export after the published HSC Y3 cuts and
tomographic bins 1–4. It uses 0.1–3 Mpc proper-radius bins and reports three
mask cases: no B-mode exclusion, exclude raw zero, and exclude raw one. The
last two were sensitivity cases before the coordinate-based polarity validation below.
The loader validates the raw bit without assigning a polarity.

The run is descriptive only: it has no shear response/bias calibration, source
`n(z)` dilution correction, random-lens subtraction, or source-overlap/spatial
covariance. It reports neither calibrated shear nor `DeltaSigma`, and keeps
`science_measurement=false`. See
[`wl_raw_shape_profile_assessment_2026-10-19.md`](wl_raw_shape_profile_assessment_2026-10-19.md)
and the hash-recorded output at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_raw_shape_profiles_2026-10-19_02/`.

Reproduce on `lageunha` without network access:

```bash
python3 analysis/icl_wgl_membership/hsc_debug/audit_wl_raw_shape_profiles.py \
  --rows /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_shape_rows_export_2026-10-18_03 \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_raw_shape_profiles_NEW
```

## Paired ICL sky-estimator contrast (2026-10-19)

`audit_icl_sky_estimator_delta.py` compares pixel-SFD and plane-first ICL
estimators on the same 74 matched pairs. It keeps each filter separate,
pre-registers `i2` as primary, residualizes on redshift and richness, and
uses identical connected target/control 5-degree block draws for both
estimators. In `i2`, the partial-correlation point estimate changes from
0.239 to 0.124; the paired 95% block-bootstrap interval for the difference
is −0.329 to 0.063 and includes zero. This is an exploratory sensitivity
result; the ICL science gate remains closed. Reproduction code, tests, and
interpretation are in
[`icl_sky_estimator_delta_assessment_2026-10-19.md`](icl_sky_estimator_delta_assessment_2026-10-19.md)
and [`test_audit_icl_sky_estimator_delta.py`](test_audit_icl_sky_estimator_delta.py).
The hash-recorded run is at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_sky_estimator_delta_2026-10-19_01/`.

## Paired-field plane-estimator injection (2026-10-19)

`audit_icl_plane_injection.py` injects four diffuse-light profiles into the
existing 74 matched target/control fields and reruns the downstream local
plane estimator on each image. This evaluates plane-subtraction transfer on
delivered coadd/bg data; it does not test HSC's upstream sky subtraction. The
run and limits are recorded in
[`icl_plane_injection_assessment_2026-10-19.md`](icl_plane_injection_assessment_2026-10-19.md).
The corrected Slurm run completed on 74 pairs with four shapes and 5,000
connected-block resamples. The 100 kpc profiles recovered 99.1–99.5% after
plane subtraction; the 250 kpc profiles recovered 86.0–90.8%. Cluster-control
recovery differences were consistent with zero for all shapes. The first
attempt's aggregation failure remains preserved in `_01`; the regression
test now checks the fixed pair count per shape. Verified results are at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_plane_injection_2026-10-19_02/`.

## Weak-lensing calibration input gate

The approved shape-row export and the official HSC Y3 four-bin SACC file are
available. The SACC retrieval, checksum, duplicated tracer tables, and
`source_zbin` count comparison are documented in
[`hsc_y3_nz_sacc_validation_2026-10-19.md`](hsc_y3_nz_sacc_validation_2026-10-19.md).
The official HSC page describes the four `wl_0`–`wl_3` n(z) tracers and calls
for dilution and selection-bias corrections. This input gate now passes. The B-mode
polarity has since been resolved and a calibrated cluster-only stack completed, but it
is not a science measurement: random-lens subtraction, n(z) uncertainty, source overlap,
and spatial covariance still need resolution.

### B-mode polarity sensitivity stack

The earlier two-hypothesis calibrated stack used only the previously approved row export
and SACC input; it did not query HSC or download anything. Its competing bit assignments
were resolved by matching exported source coordinates to the excluded footprint published
in the HSC Y3 systematics paper, not by choosing the more physical-looking shear profile.

```bash
sbatch analysis/icl_wgl_membership/hsc_debug/run_wl_bmode_polarity_sensitivity.sbatch
```

The script refuses to overwrite its output path. The historical sensitivity run and its
interpretation are recorded in
[`wl_bmode_polarity_sensitivity_2026-10-19.md`](wl_bmode_polarity_sensitivity_2026-10-19.md).

The coordinate validation exact-matches raw `b_mode_mask=0` to the published GAMA09H
rectangle (RA 132.5–140°, Dec 1.6–5°) and raw 1 to positions outside it after the
published shape and tomographic cuts. The audit records 77,021 inside-footprint pairs,
zero outside, and does not use the observed shear to infer bit meaning. Using raw 0 as
the excluded value, the calibrated cluster-only stack retained 2,153,004 lens-source
pairs across 609 unique lenses. Tangential `DeltaSigma` declines from 161.5 ± 11.5 to
15.3 ± 0.93 M_sun pc^-2 over 0.1–3 Mpc; quoted errors are independent-pair shape noise
only. The 1.5–2.0 Mpc cross component is 3.1 ± 1.6 and needs null-test follow-up. These
are preliminary cluster-only values: random subtraction, source overlap, photo-z
uncertainty, and spatial covariance are absent, so `science_measurement=false`. Full
checksums, counts, and the profile table are in
[`wl_bmode_coordinate_validation_2026-10-19.md`](wl_bmode_coordinate_validation_2026-10-19.md);
the exact output is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_cluster_stack_bmode0_2026-10-19_01/`.

### Random-control availability query (CAS count-only query complete)

HSC documents `pdr3_wide.random` as the public random-point table, with a
nominal density of 100 points/arcmin², an `adjust_density` thinning column,
non-overlap flags, and local i-band mask/coverage fields. The Astrafek planner
[`prepare_wl_random_control_counts.py`](prepare_wl_random_control_counts.py)
uses only the verified 800-lens preflight and emits 32 SQL batches. Each batch
returns one `lens_id` and one integer candidate count per lens. The query
counts thinned (`adjust_density < 0.01`) random points in the local 30′ annulus
outside twice the lens's angular 3 Mpc aperture, requiring primary points and
unmasked i-band pixels. HSC defines `isprimary=True` as objects in the inner
patch and inner tract with no children in its [PDR3 FAQ](https://hsc-release.mtk.nao.ac.jp/doc/index.php/faq__pdr3/).
The query returns no coordinates, shape rows, or photo-z values. HSC specifies
the random table, density, and mask fields in its [random-points documentation](https://hsc-release.mtk.nao.ac.jp/doc/index.php/random-points__pdr3/).

The first prepared plan at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_point_count_plan_2026-10-19_01/`
was rejected by CAS before a job was created because `patch_inner` is not a
column in the PDR3 random table. The corrected plan is at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_point_count_plan_2026-10-19_02/`;
it uses the documented `isprimary` selection for the same non-overlap intent.
On 2026-10-19, all 32 CAS batches completed and passed validation: 800 unique
lenses have one integer count each, totaling 780,191 random-point candidates;
all 800 lenses have at least one candidate. Only the two aggregate columns
(`lens_id`, `n_random_candidates`) were returned. No random coordinates, shape
rows, source rows, photo-z values, or images were returned, and this was not a
science measurement. A sufficient candidate count establishes availability
only. Any later coordinate export still needs all-CAMIRA exclusion,
full-aperture footprint checks, HSC-I depth/PSF matching, uniqueness and overlap
checks, and spatial covariance validation before random centers can serve as
random lenses. The count summary and per-batch checksums are stored outside Git
under `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_point_count_plan_2026-10-19_02/cas_counts/`.

#### Approved random-coordinate export and local geometry gate

After the count-only result, the user approved exporting only the candidate
coordinates needed for local random-lens screening. The immutable plan is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_point_coordinate_plan_2026-10-19_01/`.
It asks for 780,191 candidate associations in 32 CAS batches, with exactly
`lens_id, ra_deg, dec_deg`; it requests no object/source rows, shape fields,
photo-z values, or images. On 2026-10-19 the live export checkpoint had 2/32
batches complete (33,514 rows) and batch 002 submitted as CAS job 1062312. The
checkpoint and batch hashes, rather than this progress note, are authoritative
for completion.

`screen_wl_random_control_coordinates.py` is the offline next step after that
checkpoint is complete. It verifies each coordinate CSV against the approved
per-lens counts and checksums, then audits overlap between each candidate's
3-proper-Mpc aperture and every supplied CAMIRA cluster's 3-proper-Mpc
aperture. It also flags exact repeated sky coordinates across lens-associated
candidate rows. Its output is only a geometry-eligible pool: full-aperture HSC
coverage, i-band mask completeness, HSC-I depth/PSF/sky matching, random-center
selection, and spatial-block covariance remain separate gates. It is not a
science measurement. The coordinate export was active when recorded; the
offline screen has not been run. Any later HSC shape-row or image export needs
its own explicit user approval.

The 780k-coordinate geometry screen is available as the Slurm job
`run_wl_random_coordinate_geometry.sbatch`; its analysis code is CPU-only.
Submit it only after the approved
coordinate export has a terminal `complete` checkpoint; the runner checks that
state and refuses a partial result. It makes no HSC query.

`prepare_wl_random_centers.py` is the local step after the geometry screen. It
verifies the coordinate plan/checkpoint hashes, the geometry output hashes and
the per-lens counts, then selects up to an explicitly supplied number of
centers per matched lens by a recorded SHA-256 ranking seed. It writes the
`astrafek-wl-random-centers-v1` manifest consumed by
`prepare_wl_random_shape_count_plan.py`. The required `--per-lens` argument has
no default, so the provisional sample size is an explicit analysis choice. The
manifest remains `science_measurement=false`: full-aperture coverage, depth,
PSF, sky matching, random-center overlap and spatial covariance are still
unverified. Preparing the later S19A count-only SQL is offline; submitting it
requires separate user approval, and exporting any source-shape rows requires
another approval after the observed counts and transfer size are known. Tests
for deterministic balanced selection are in
`test_prepare_wl_random_centers.py`.

`plan_wl_random_center_scale.py` compares candidate random-to-lens ratios using
the measured 800-lens shape-row count and CSV bytes; it performs no network
request. The current export has 3,069,544 pairs in 909,833,598 bytes (296.4
bytes per pair). An 8:1 scenario would mean 6,400 random centers, 256 count-only
SQL batches, and an estimated 24,556,352 source-pair rows (6.78 GiB), if
random fields match the lens export's mean density and row width. Under equal,
independent per-center shape variance, the random stack would contribute 1/8 of
the lens-stack variance, or 11.1% of their combined variance. This is a
planning candidate only: spatial overlap, depth and the eventual per-center
CAS counts can change the estimate. The reproducible scenario table is saved
outside Git at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_center_scale_2026-10-19_01/scale_plan.json`.
Exact random centers and count-only SQL will be prepared after the coordinate
geometry screen; no query has been submitted for this sample-size plan.

## ICL–satellite spatial-block scale sensitivity (2026-10-19)

`audit_satellite_icl_block_scale_sensitivity.py` reruns the same verified
geometric-strata sample with 2.5°, 5°, and 10° connected target/control RA/Dec
cells. The 5° output exactly reproduces the archived 5,000-draw bootstrap in 60
numeric comparisons. For i2, the full-sample partial r remains 0.234, but the
95% spatial-block interval crosses zero at 2.5° and not at 5° or 10°. Thus the
interval depends on the assumed correlation scale; this remains exploratory,
has no multiplicity correction, and is not a detection. Details and hashes are
in [`icl_satellite_block_scale_assessment_2026-10-19.md`](icl_satellite_block_scale_assessment_2026-10-19.md).
The output is outside Git at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_icl_block_scale_sensitivity_2026-10-19_04/`.

The follow-up [`audit_satellite_icl_perbin_block_scale.py`](audit_satellite_icl_perbin_block_scale.py)
repeats the same test on the primary 94-row per-annulus sample (68 i2 pairs),
checking its 5° result against the archived bootstrap. The i2 partial r is
0.293 at each width, but its 95% block interval changes from [−0.003, 0.504]
at 2.5° to [0.042, 0.497] at 5° and [0.057, 0.491] at 10°. This remains an
exploratory, width-sensitive interval, not a detection; the report is
[`icl_satellite_perbin_block_scale_assessment_2026-10-19.md`](icl_satellite_perbin_block_scale_assessment_2026-10-19.md).
