# HSC 디버그 절차

이 디렉터리의 스크립트는 ICL 스케일 스택을 실제 목록으로 준비하는 절차이다. 논문에 함께 제공한다. 실행 결과는 과학 측정이 아니다. `science_measurement` 는 거짓으로 남는다.

두 연구 갈래의 현재 표본, ICL–위성 연관성 검증, 표적 전용 weak-lensing 결과와 남은 과학 관문은 [통합 연구 평가](icl_wgl_integrated_research_assessment_2026-10-19.md)를 참조한다.

2026-10-20에는 HSC patch-QA 트랙 연결요소를 하나씩 제외하는 약한렌즈
영향도 진단을 추가했다. 재현 명령은
[`run_wl_hsc_tract_region_influence.sbatch`](run_wl_hsc_tract_region_influence.sbatch),
결과와 해시는 [tract-region 평가](wl_hsc_tract_region_influence_assessment_2026-10-20.md)에
있다. 세 연결 영역만으로 공분산을 추정하지 않으며, 무작위 중심 차감도
대신하지 않는다.

같은 연결요소 민감도를 783개 matched target-control ICL 프로파일에 적용한
오프라인 감사는 [`audit_icl_tract_region_influence.py`](audit_icl_tract_region_influence.py)
로 재현한다. 쌍별 프로파일 재집계가 저장 스택 평균과 일치하는지 확인하며,
세 영역의 제외 결과는 공분산이 아닌 영향도로 기록한다. 상세 결과는
[ICL tract-region 평가](icl_hsc_tract_region_influence_assessment_2026-10-20.md)에 있다.

위성–ICL 연관성 표본 94개에도 양 끝점 영역 삭제를 적용했다. 재현기는
[`audit_satellite_icl_tract_region_influence.py`](audit_satellite_icl_tract_region_influence.py)이며,
i2 partial Pearson이 전체 `0.293`에서 한 연결영역 제거 후 `0.048`로 변하고
Spearman 부호도 바뀌는 영향도를 기록한다. 결과와 한계는
[연관성 tract-region 평가](satellite_icl_tract_region_influence_assessment_2026-10-20.md)에
있다.

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

2026-10-20에는 이미 저장된 PSF 배치를 재요청하지 않고 승인 계획 기준으로
744개 파일을 다시 감사했다. Astrafek `PositionPSFBank`가 57개 cluster/filter
프레임과 688개 위성 위치를 연결하고, 시험 중심 PSF를 포함한 745개 고유 FITS의
체크섬과 영상 배열을 모두 읽는 통합 검증도 수행했다. 재현 스크립트는
[`verify_position_psf_bank.py`](verify_position_psf_bank.py)이며, 결과는
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_position_psf_bank_validation_2026-10-20_01/`
에 있다. 요약 SHA-256은
`384b859fd59169446515027b80990c9482f64bf7e51d108f8008014c5d024207`, 스크립트
SHA-256은 `59ccba1aedc76a7cd43110f969c17cdd3f15b51307d3602428b5b3cf7c128302`이다.
이 검증은 PSF 코어 위치·파일 연결만 확인하며 확장 날개나 과학 측정의 검증은 아니다.

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

744개 승인 배치와 시험 PSF 1개의 중심부 크기 분포를 해시 재검증한 뒤 요약했다. 전체 745개 중 양수 성분 2차 모멘트의 가우스 등가 장축 FWHM 중앙값은 1.084 arcsec이고, 5–95백분위는 0.965–1.198 arcsec다. 프레임 중심 57개와 위성 위치 688개의 중앙값은 각각 1.087, 1.084 arcsec다. 이들은 7 arcsec 안팎 모델 코어의 설명 지표이며 target-specific 확장 PSF 날개나 seeing 보정값이 아니다. 음수 보간 픽셀을 보존하고 `science_measurement=false`를 유지한다. 재현 코드·테스트·해시는 [PSF 코어 분포 평가](psf_core_variation_assessment_2026-10-20.md)에 있다.

같은 57개 주입 프레임에서 ICL 200–300 kpc target-minus-control 신호와 중심 PSF 코어 폭의 양의 공변동을 점검했다. redshift·richness·저장 mask common coverage를 조정한 부분상관은 i `r=0.500`(13개), i2 `r=0.247`(44개)였다. 블록 재표집은 i에서 매우 넓고 i2에서 0 부근을 포함해 검출로 해석하지 않는다. 저장 annulus 마스크와 재검출 마스크의 지원률 차이도 기록했다. 상세 방법과 SHA-256은 [ICL–PSF 코어 공변동 평가](icl_psf_core_dependence_assessment_2026-10-20.md)에 있다.

2026-10-20에는 동일한 57개 프레임에서 100 kpc와 500 kpc 확산 성분을 각각 위치별 모델 PSF로 합성하고 SEP/K2.5 마스크를 시나리오별로 재생성했다. 평면 보정/원시 ICL 응답비의 중앙값은 i에서 `1.002 → 0.644`, i2에서 `0.997 → 0.633`으로 내려갔다. 이 q16–q84는 프레임 간 산포이며 신뢰구간이 아니다. 500 kpc i2 한 프레임은 `icl_only`에서 위성 측정 지원을 잃었지만 ICL 응답은 보존했다. 요약·감사·제한은 [PSF 마스크 재검출 규모 비교](icl_psf_redetection_scale_comparison_assessment_2026-10-20.md), 재현 집계 코드는 [`summarize_joint_icl_satellite_psf_scale_comparison.py`](summarize_joint_icl_satellite_psf_scale_comparison.py)에 있다. 결과는 주입 전달 민감도이며 `science_measurement=false`다.

`estimate_hsc_coadd_noise_correlation.py` compares short-range pixel residual
covariance in those 57 target coadds with their 57 already-downloaded matched
control coadds. It applies the saved K2.5 target source masks and regenerates
the same SEP/Kron-2.5 masks for matched controls. It samples masked 64×64
patches, removes a clipped local plane, normalizes by the delivered variance,
and records lag covariance plus an approximate large-top-hat variance factor.
It makes no data request and does not supply a source-photon-noise model;
unresolved sources and residual diffuse structure can still contribute. Run it through
[`run_icl_coadd_noise_correlation.sbatch`](run_icl_coadd_noise_correlation.sbatch).
This technical factor cannot replace the exact ICL/satellite estimator
covariance or total science error.
The matched-mask run found a 9×9 lag-grid variance factor near 1.59 in targets
and 1.54 in controls, with paired differences spanning zero. See
[`icl_coadd_residual_covariance_assessment_2026-10-19.md`](icl_coadd_residual_covariance_assessment_2026-10-19.md)
for the independent checksum and numerical audit; no correction was applied
to the science profile.

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

The cluster-only spatial jackknife has now run for the existing 800-lens stack.
Job 428563 used 52 occupied 5-degree equal-area cells and 2,153,004 selected
pairs. The tangential `DeltaSigma` profile declines from 1.61×10^14 at 0.1–0.2
Mpc to 1.53×10^13 Msun Mpc^-2 at 2–3 Mpc; the joint cross-component diagnostic
is `chi2=9.30` for 7 modes (`p=0.232`). This remains exploratory because the
covariance is cluster-only, random subtraction and matched-field variance are
absent, and the grid is not a contiguous HSC-region partition. See
[`wl_cluster_spatial_jackknife_assessment_2026-10-19.md`](wl_cluster_spatial_jackknife_assessment_2026-10-19.md)
and its recorded machine-readable output hash.

`audit_wl_tomographic_consistency.py` extends that cluster-only audit by
stacking HSC Y3 source bins 1–4 separately and using the same delete-one 5°
sky-cell replicates across bins and radial ranges. It stores a joint
tomographic/radial tangential-cross covariance and tests, per radial bin,
whether calibrated `DeltaSigma` agrees across source bins using that covariance.
This is a calibration consistency diagnostic, not a mass detection; it still
needs random-lens subtraction, HSC `n(z)` uncertainty and contiguous-region
covariance. The Slurm wrapper is
[`run_wl_tomographic_consistency.sbatch`](run_wl_tomographic_consistency.sbatch).
It reads only the already approved cluster lens-source rows and local SACC
`n(z)` input; no new HSC request or export is involved.
The successful run was Slurm job 428830 (12m 41s, 1.83 GiB peak RSS), with
52 shared sky-cell deletions and 2,153,004 selected pairs. Source-bin
`DeltaSigma` differs across bins at 0.1–0.7 Mpc; the GLS consistency tests
remain significant after Holm correction across seven radial intervals. The
global cross-component test also flags a joint-bin residual, although no
single-bin cross test survives four-test Holm correction. These are
cluster-only calibration diagnostics, not detections. Full values, covariance,
limitations, and reproducibility hashes are in
[`wl_tomographic_consistency_assessment_2026-10-20.md`](wl_tomographic_consistency_assessment_2026-10-20.md);
the saved figure is outside Git beside `summary.json` as
`tomographic_profiles.png`.

`audit_wl_source_reuse_covariance.py` measures how reusing the same source
galaxy in multiple lens stacks changes delete-one covariance when blocks are
deleted by lens endpoint, source endpoint, or either endpoint jointly. The
independent `validate_wl_source_reuse_covariance.py` checks each saved matrix's
symmetry, positive semidefiniteness to `1e-12` of the largest matrix-entry
scale, rank and profile-error diagonal; it also
recomputes the four-bin GLS consistency tests and Holm adjustments and checks
the source-reuse counts. Run it only after the producer writes `summary.json`,
using a fresh output path:

```bash
python3 analysis/icl_wgl_membership/hsc_debug/validate_wl_source_reuse_covariance.py \
  --summary /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_source_reuse_covariance_2026-10-20_01/summary.json \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_source_reuse_covariance_validation_2026-10-20_02/validation.json
```

The validator audits internal numerical consistency of the producer output; it
does not reload the catalogue rows or replace random-lens subtraction and
connected-survey-region covariance.

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

After a completed count-only export, `prepare_wl_random_shape_rows.py` builds an
offline, hash-bound source-row plan with the exact observed per-center row
totals. Its 23 explicit columns contain calibrated shape, source coordinates,
tomographic-bin and B-mode fields; it requests neither images nor full photo-z
PDFs. `submit_wl_random_shape_rows.py --check --plan-dir ...` audits the plan
and count evidence without network access. Submission requires a distinct
`--approved-random-shape-row-export` flag and remains a separate approval from
the count-only query. `run_wl_random_shape_stack.py` accepts only complete,
hash-matched row exports, checks each center-to-lens association and count, and
runs the random-center S19A stack. With the verified cluster row export, it
recomputes the stored target profile and forms paired delete-one 5-degree
equal-area sky-cell replicates, removing a pair whenever either endpoint is in
the deleted cell. This carries shared source positions through the same block
in both samples. The covariance is still conditional on this approximate grid;
block-width and connected HSC-region sensitivity, plus photo-z uncertainty,
remain to be assessed, so the output keeps `science_measurement=false`. None of
these shape-row requests or stacks has been run. Focused offline tests cover
count binding, SQL scope, submission gates, pair mapping, endpoint deletion and
existing stack math. The target/random provenance check passes against the
existing B-mode-zero cluster-only stack summary and compares its stored
`masked_value` and `H0_km_s_mpc` fields before subtraction.
The compute wrapper is `run_wl_random_shape_stack.sbatch`; it requests bounded
Slurm CPU and memory resources and only consumes completed local exports. It
does not contact CAS or request data.

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

## Full-sample 2-D ICL target/control stack and support sensitivity (2026-10-19)

The Astrafek 2-D stack uses the 800-cluster parent and 798 matched-control images;
783 same-filter pairs remain (302 i, 481 i2). An independent-support sensitivity
was completed for SEP/Kron mask factors K=2.5, 3.5, and 5.0, with all 40 filter/radial
aggregates in each run independently verified against per-pair 25 kpc profiles.
Same-pair changes versus common support vary with radius and reach 32.9% at K=5,
i2, 175–200 kpc. This estimator comparison is exploratory, leaves
`science_measurement=false`, and does not resolve target/control environment bias
or ICL–satellite association significance. See
[`icl_2d_target_control_stack_README.md`](icl_2d_target_control_stack_README.md)
and [`icl_2d_independent_support_assessment_2026-10-19.md`](icl_2d_independent_support_assessment_2026-10-19.md).

The 2026-10-20 outer-annulus test subtracts each matched pair's own 400–500 kpc
mean from its inner profile, requiring at least 5% common-pixel support in both
regions. In this support-selected sensitivity the i-band outer level is
`1.66–1.78` nJy/arcsec² across the three inner annuli, while i2 is near zero;
the resulting inner excesses remain positive in the connected 5° block
intervals. This does not show whether the i-band offset is residual sky or real
extended light, and the outer annulus itself may contain ICL. It is not a
corrected ICL profile or a detection. Reproduce with
[`audit_icl_outer_annulus_baseline.py`](audit_icl_outer_annulus_baseline.py);
the 5,000-replicate output and paired rows are in
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_outer_annulus_baseline_2026-10-20_02/`.
An offline regional influence audit found that deleting any one of the three
approximate occupied tract-center regions leaves the `i` outer baseline
positive (1.08–2.20 nJy/arcsec² across the inner annuli); `i2` remains near zero
overall and its leave-one-region-out values cross zero. Region labels use target
coordinates only. This sensitivity check does not separate residual sky from
real outer ICL and does not estimate covariance. Details, row counts, and the
verified result hash are in
[`icl_outer_baseline_region_influence_assessment_2026-10-20.md`](icl_outer_baseline_region_influence_assessment_2026-10-20.md),
reproduced with [`audit_icl_outer_baseline_regions.py`](audit_icl_outer_baseline_regions.py).
An additional pair-level coherence analysis finds that the same 400–500 kpc
offset covaries with inner target-minus-control light after adjustment for
redshift, richness, and area support, most clearly at 200–300 kpc (partial
`r=0.455` in i and `0.421` in i2). Connected-block intervals are descriptive
across six exploratory comparisons. The relation cannot distinguish residual
sky from real extended ICL; see
[`icl_outer_baseline_coherence_assessment_2026-10-20.md`](icl_outer_baseline_coherence_assessment_2026-10-20.md)
and reproduce it with
[`analyze_icl_outer_baseline_coherence.py`](analyze_icl_outer_baseline_coherence.py).

Keep the `i` and `i2` stacks separate. HSC PDR3 documents those as similar but
non-identical response curves; the combined coadds use object-specific
effective-response corrections, and the uncorrected photometric difference is
about 1% on average but can be larger for very red objects. Since diffuse ICL
has no per-pixel color information for that correction, pooling the two stacks
would add an uncontrolled color term. The current pipeline verifies identical
filter labels within every target/control pair and reports each filter apart
([HSC PDR3 processing](https://hsc-release.mtk.nao.ac.jp/doc/index.php/processing__pdr3/),
[filter-fraction note](https://hsc-release.mtk.nao.ac.jp/doc/index.php/filterfrac-2/)).

On 2026-10-20, `verify_icl_target_control_2d_stack.py` independently rechecked
the common-support K=2.5 stack and K=3.5/K=5.0 stacks: all three retain 783
same-filter pairs, 40 radial rows, and 201×201 maps whose target-minus-control
identity and pair counts agree. The audit resolves the omitted K=2.5 numeric
factor from the summary definition and the mask in its archived run snapshot,
then resolves K=3.5/K=5.0
runner versions from Git history by exact SHA-256. Verification records are
saved alongside each stack as `independent_verification_2026-10-20.json`.
Their SHA-256 values are K=2.5 `0da3238fa4105118667d5e67392e09c07dc1b59db68adc4344aada464447d740`,
K=3.5 `113645d19e2ca56a9f091d6cca38b5c48cf9bc15c534fdc319fb172d9f229997`,
and K=5.0 `34e0019cfeee01a4aed239c354e5f9035956798ea1cd9fdb9a52cdeb29a2a906`.
The verifier and regression tests are
[`verify_icl_target_control_2d_stack.py`](verify_icl_target_control_2d_stack.py)
and [`test_verify_icl_target_control_2d_stack.py`](test_verify_icl_target_control_2d_stack.py).
These checks confirm product integrity and provenance; all three products
remain exploratory with `science_measurement=false`.

The fixed 94-row ICL–satellite association sample has also been compared under
three independent-support mask factors using paired connected-sky-block
resampling. The full supported-pair i2 estimate remains near partial r=0.29,
while i remains near zero; the paired intervals for the estimator changes
include zero. This is a fixed-sample estimator sensitivity, not a population
reanalysis or a science measurement. Details and limitations are in
[`satellite_icl_independent_support_association_2026-10-20.md`](satellite_icl_independent_support_association_2026-10-20.md).

`audit_satellite_icl_leave_one_block_out.py` omits each connected target/control
sky block in turn and refits the redshift/richness partial correlation. The i2
estimate stays positive at 2.5°, 5°, and 10° widths, though its minimum falls
from 0.293 to 0.166–0.177; the i estimate changes sign at 2.5° and 10°. These
are influence ranges, not confidence intervals or detection tests. Full results
and cluster IDs are in
[`satellite_icl_leave_one_block_out_assessment_2026-10-20.md`](satellite_icl_leave_one_block_out_assessment_2026-10-20.md).

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
support counts and reports the nominal independent-pixel covariance. With
`--lag-covariance-dir`, it also contracts the exact signed weights through the
checksum-verified 9×9 per-frame short-range lag kernel measured from matched
target/control coadds. This is a finite-lag sensitivity, not a full coadd model;
the external ICL random-control term remains absent and no measured relation is
corrected.
See [`satellite_icl_shared_pixel_covariance_2026-10-18.md`](satellite_icl_shared_pixel_covariance_2026-10-18.md)
and `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_icl_shared_pixel_covariance_2026-10-18_01/`.
The lag-kernel sensitivity runner is
[`run_satellite_icl_lag_covariance.sbatch`](run_satellite_icl_lag_covariance.sbatch).

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

On 2026-10-20 the retention audit was rerun against the current 800-frame mask table,
corrected full-profile diagnostics, per-annulus profiles, and the spatial-block association
table; it verifies their SHA-256 links. The ICL–satellite statistics are now read from the
association result rows rather than copied as constants. The saved audit is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_sample_retention_2026-10-20_02/summary.json`
(SHA-256 `7952f53bd29289dc182b24abdd7b0356b9a09e40c4e15055a0ea1f6aa1ca354f`).
It confirms that independent 50–75 kpc satellite stacks have means
`2.51 [-0.30, 5.50]` nJy/arcsec² (i; 70 satellites/42 clusters) and
`1.19 [-0.24, 2.49]` (i2; 308/140), both with cluster-bootstrap intervals crossing zero.
For the 94 joined ICL–satellite pairs, i2 partial Pearson `r=0.293` has a 5° block interval
`[0.042, 0.497]`, while the ordinary-pair interval `[-0.007, 0.526]` and Spearman block
interval `[-0.088, 0.301]` include zero. The i partial Pearson estimate is `0.079` with a
block interval `[-0.186, 0.479]`. This estimator-dependent evidence remains exploratory.
The audit regression tests cover dynamic statistic extraction and reject inconsistent input
rows in [`test_audit_icl_satellite_sample_retention.py`](test_audit_icl_satellite_sample_retention.py).

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

The 2026-10-20 extension adds 500 kpc exponential stress profiles while
retaining the 100/250 kpc cases. After plane subtraction, the 500 kpc circular
profiles retain about 63% of injected signal and q=0.65 profiles retain
70–71%, similarly in target and matched controls; all paired differences
include zero. This quantifies downstream attenuation for these assumed
profiles, not a correction to observed ICL. Job 428956 completed on all 74
pairs. The full table, output hashes, and limitations are recorded in
[`icl_plane_injection_extended_assessment_2026-10-20.md`](icl_plane_injection_extended_assessment_2026-10-20.md);
the result is at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_plane_injection_extended_2026-10-20_01/`.

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

The tomographic disagreement investigation identified a newer DESI DR1/DR2
clustering-redshift calibration performed directly on the HSC PDR3 S19A shape
catalog. Its authors publish calibrated distributions under
[`desi-y3-hsc/results/distributions`](https://github.com/JeanCHDJdev/desi-y3-hsc/tree/main/results/distributions).
It includes high-redshift tracers absent from the older HSC Y3 calibration and
reports upward shifts in bins 3 and 4 relative to Rau+2022. See
[`wl_pdr3_desi_nz_calibration_update_2026-10-20.md`](wl_pdr3_desi_nz_calibration_update_2026-10-20.md).
Those result files have not been imported; obtain approval before downloading
the four `_v_1p1` distribution/metadata files (549,818 bytes total) identified
in the calibration note, then compare both published scale cuts to the SACC
baseline with the same random-subtracted stack and spatial jackknife.
While those files await approval, `audit_wl_published_nz_shift_sensitivity.py`
tests a no-download approximation: it translates the local SACC source-bin 3
and 4 curves toward higher redshift by the paper's central mean shifts (+0.043
and +0.050), renormalizes them, and repeats the same joint jackknife. The
[paper's v3](https://arxiv.org/abs/2511.18133v3) reports these shifts relative
to Rau+2022; this constant-translation sensitivity is not the authors'
calibrated spline distributions or their uncertainty.
The same Astrafek audit now also samples the published marginal 68% shift
interval endpoints: higher-redshift translations of 0.021–0.063 for bin 3 and
0.038–0.062 for bin 4, using the central point and all four rectangle corners.
The grid is a sensitivity envelope, not a joint posterior; the calibrated
distribution download remains separately approval-gated. The reproducible
Slurm runner writes a new `_02` result directory and preserves `_01`. After the
run completes, [`validate_wl_published_nz_shift_grid.py`](validate_wl_published_nz_shift_grid.py)
checks all five shift inputs, the saved SACC and baseline hashes, each 56×56
covariance's rank and symmetry, and the recomputed radial significance envelope.
The 2026-10-20 grid completed and passed that independent validator: every
scenario retained 52 jackknife cells and a rank-51 covariance. The source-bin
consistency discrepancy persisted in all five cases at 0.1–0.7 Mpc; the
seven-radius Holm-adjusted p-value envelope is 7.6×10⁻⁵–1.04×10⁻³ there.
This supports prioritizing the published calibrated distributions, but the
constant translations are not a posterior propagation and the cluster-only
profiles still omit random subtraction. The run and validation checksums are
recorded in [`wl_published_nz_shift_grid_assessment_2026-10-20.md`](wl_published_nz_shift_grid_assessment_2026-10-20.md).

The current S19A responsivity, object-level shear corrections, and HSC
edge-selection coefficients were also checked against the HSC Y3 paper and its
linked calibration utility. The source-level formula review is recorded in
[`wl_s19a_shape_calibration_formula_review_2026-10-20.md`](wl_s19a_shape_calibration_formula_review_2026-10-20.md).

To test whether the inner source-bin disagreement is concentrated in one lens-
redshift range, `audit_wl_tomographic_lens_redshift.py` splits lenses into equal-
count lower/upper halves sorted only by `(lens_z, lens_id)`. It recalculates all
four source-bin stacks in one shared spatial jackknife and Holm-adjusts the 14
radial consistency tests. The result plotter
`plot_wl_tomographic_lens_redshift.py` displays both lens cohorts' tangential and
cross profiles. The initial Slurm job 428849 completed using only the approved
shape rows and local SACC input; it represented 632 lenses, 316 per half, with
609 lenses surviving the published shape/radius selections. Both halves show
source-bin disagreement at 0.2–0.4 Mpc after 14-test Holm correction. The lower-
redshift half also shows corrected disagreement at 0.1–0.2 Mpc; its 0.4–0.7 Mpc
result does not pass correction. That initial summary did not persist the full
112×112 covariance, so it cannot support a direct test of whether the relative
bin offsets differ between lens-redshift halves.

The covariance-complete rerun and
`audit_wl_tomographic_redshift_interaction.py` now test the difference-in-
differences between source-bin offsets with the saved covariance. There is no
family-wise-significant change of the relative source-bin offsets across the
lens-redshift halves (global p=0.106); the 0.2–0.4 Mpc within-cohort mismatch
persists in both halves. Full statistics, limits, and figures are documented in
[`wl_tomographic_lens_redshift_assessment_2026-10-20.md`](wl_tomographic_lens_redshift_assessment_2026-10-20.md).
This remains a cluster-only diagnostic with no random-lens subtraction and
`science_measurement=false`.

`audit_wl_source_reuse_covariance.py` uses the same approved pair rows to
compare lens-endpoint, source-endpoint, and joint-endpoint spatial jackknife
covariances. It also counts reused HSC source IDs and verifies that one source
ID maps to one 5° source cell, so all repeated lens pairs are deleted together
in the source-aware variants. It makes no HSC request and keeps
`science_measurement=false`; output is written to a new directory by
[`run_wl_source_reuse_covariance.sbatch`](run_wl_source_reuse_covariance.sbatch).

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

`screen_wl_random_control_coordinates.py` is the offline full-sample step after
that checkpoint is complete. It verifies each coordinate CSV against the
approved per-lens counts and checksums, then audits overlap between each
candidate's 3-proper-Mpc aperture and every supplied CAMIRA cluster's
3-proper-Mpc aperture. It also flags exact repeated sky coordinates across
lens-associated candidate rows. Its output is only a geometry-eligible pool:
full-aperture HSC coverage, i-band mask completeness, HSC-I depth/PSF/sky
matching, random-center selection, and spatial-block covariance remain
separate gates. It is not a science measurement.

While the export is incomplete, `screen_wl_random_completed_batches.py` can
run an explicitly partial geometry pilot on only hash-verified completed
batches. Slurm job 428780 screened the first 50 lenses and 33,514 coordinates:
2,129 passed the geometry rules, only 31 lenses had any such candidate, and
the median was 5.5 candidates per lens. This early subset is not representative
of all 800 targets; it does show why per-lens support must be checked before
random-center selection. See
[`wl_random_coordinate_geometry_pilot_2026-10-20.md`](wl_random_coordinate_geometry_pilot_2026-10-20.md).
The pilot note records batch 002 as submitted while its authenticated
submitter was alive. A later SSH recheck on `lageunha` confirmed that
submitter PID 1049265 remains live; its older in-memory code does not persist
successful CAS status polls. The checkpoint is still at 2/32 complete with
batch 002 submitted. Its ID 1062312 is an HSC CAS
job ID, not a Slurm job ID; do not query it with `squeue` or submit it again.
The 2026-10-20 04:45 KST recheck found PID 1049265 alive after about 17 h 37 m.
The checkpoint remains at 2/32 batches complete (33,514 rows), with batch 002
submitted. Its file modification time is still 2026-10-19 11:07:53 KST, five
seconds after that process started. The live process has the earlier submitter
code loaded and does not persist successful CAS status polls, so its current
status is unknown. Keep the live submitter and do not resubmit batch 002. Any
later HSC shape-row or image export needs its own explicit user approval.

The Astrafek audit [`audit_wl_random_coordinate_checkpoint.py`](audit_wl_random_coordinate_checkpoint.py)
now verifies every completed coordinate CSV against the immutable count plan,
per-lens candidate totals, row schema, and checkpoint SHA-256 without contacting
CAS or inferring process liveness. Its 2026-10-20 run independently verified
2/32 batches (33,514/780,191 coordinate rows); both saved files matched their
recorded byte lengths and hashes, and the remaining 30 batches were still
incomplete. The report is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_coordinate_checkpoint_audit_2026-10-20_02/`.
The remote process was separately confirmed alive as PID 1049265 during this
check. No CAS request was made, and the output retains
`science_measurement=false`.

On 2026-10-20, the offline tests for coordinate submission, completed-batch
screening, full geometry screening, random-center selection, and count-plan
preparation passed together (33 tests). They exercise planning and
verification paths; they do not resolve the live CAS job status.

The submitter now clears stale top-level errors on entry and stores each
successful batch status, poll time, and poll count atomically in the checkpoint.
This code change applies only to a future invocation; the already-live PID
1049265 continues running its previously loaded code and is left untouched.

2026-10-20 07:06 KST에 `lageunha`에서 PID 1049265를 다시 확인했다. 프로세스는
여전히 살아 있었고, 반복 확인에도 같은 Astrafek 체크포인트 감사는 2/32 배치
(33,514/780,191행)만 확인했다. 감사 요약 SHA-256은
`d5072ee7cecd4ec4798892dc6d3be4acf1d916738e0ff65564666234a3ec8212`, 체크포인트
SHA-256은 `84cccc6443ec8d2fc8bee1611ab6b9b57a3aabdc90978fcdfc558151d2bb92b2`다.
이번 감사 디렉터리는
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_coordinate_checkpoint_audit_2026-10-20_07/`다.
이는 저장 체크포인트만 검증하며 살아 있는 프로세스의 CAS 작업 완료 여부는 보여주지
않는다. 프로세스와 기존 CAS 요청은 그대로 두었고 재제출하지 않았다.

### Reuse test for the existing ICL matched-control centers

The 798 already downloaded ICL matched-control positions were tested as a
possible shortcut for weak-lensing random subtraction. Astrafek's
[`audit_wl_matched_control_centers.py`](audit_wl_matched_control_centers.py)
rechecks the assigned target redshifts and audits each 3-proper-Mpc aperture
against all 7,319 entries in the supplied S20A CAMIRA catalogue. Among the
798 unique controls, 768 overlap a CAMIRA 3-Mpc aperture; only 30 remain free
of both CAMIRA overlaps and mutual overlap with another clean control aperture
(3.8%). The original ICL-control selection therefore cannot represent the
full 800-lens weak-lensing random sample. The 30 clean positions are only
geometry candidates, not approved WL random centers; their shape-source
coverage, selection, and spatial covariance remain unverified.

This audit used no HSC query or source-shape rows. It ran on `lageunha` against
the existing manifests and the hash-verified full CAMIRA catalogue. Results:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_icl_matched_control_geometry_2026-10-20_04/`.
The 798-row table SHA-256 is
`74dd333e7af2b35c4da305d747f585e5ab4ae2e9cde19753550a2f17df1db0dc`; the
summary SHA-256 is
`12eec246699c24a6f194e271758bfdc09d85660b1eeedcc73a54c7aec1109cc9`.
Reproduce on `lageunha` with a new output directory:

```bash
python3 analysis/icl_wgl_membership/hsc_debug/audit_wl_matched_control_centers.py \
  --controls /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/control_download_01/download_manifest.json \
  --targets /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/download_01/download_manifest.json \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_icl_matched_control_geometry_NEW
```

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

## Finite-lag ICL–satellite pixel covariance (2026-10-19)

`audit_satellite_icl_shared_pixel_covariance.py --lag-covariance-dir ...`
contracts the actual signed ICL 200–300 kpc and satellite 50–75 kpc estimator
weights with a matched 9×9 empirical lag kernel. Across 57 supported frames
(13 i, 44 i2), this increases each estimator's median target-side noise sigma
by about 25% relative to the independent-pixel diagonal calculation; the
normalized ICL–satellite pixel-noise correlation remains close to zero and
close to its diagonal estimate. This is a noise sensitivity, not a correction
to the observed association or final science error: the ICL random-control
term, full position-dependent covariance, source Poisson noise, PSF, and sky
transfer remain outside this calculation. The row-level output, hashes, and
interpretation are in
[`satellite_icl_lag_covariance_assessment_2026-10-19.md`](satellite_icl_lag_covariance_assessment_2026-10-19.md).

`audit_satellite_icl_target_control_covariance.py` adds the matched-control
image's 200–300 kpc aperture variance to that target-side calculation. All 57
control apertures pass the same 40% coverage floor and are disjoint from their
target annuli. The total ICL sigma inflation remains about 24%; the normalized
noise-only ICL–satellite covariance stays small. A subsequent geometry audit
validated target and control image pixel scales at 0.168 arcsec/pixel. The
nearest target/control annulus boundaries are at least 297.497 arcsec apart,
versus 0.9504 arcsec maximum separation in the measured 9×9 lag kernel. Thus
the finite-lag model cannot couple pixel noise between these apertures;
large-scale sky and calibration covariance remains unmeasured. Other major
error terms remain outside the model. Results, row/table hashes, and limitations are in
[`satellite_icl_target_control_covariance_assessment_2026-10-19.md`](satellite_icl_target_control_covariance_assessment_2026-10-19.md).
The saved row table and summary can be independently checked with
[`verify_satellite_icl_target_control_covariance.py`](verify_satellite_icl_target_control_covariance.py).
The target-WCS-verified geometry audit is saved as
`finite_lag_separation_target_wcs_audit.json` beside the covariance table
(SHA-256 `79937c8aca5e90e3cc9cdc99ad64b5b18668e67f1683c1a0650649a081a183eb`).

## ICL–위성 외곽광 연관성의 외곽 기준 민감도 (2026-10-20)

기존 동일 cluster/filter/control 쌍에서 50–75 kpc 위성 외곽광과 ICL 200–300 kpc
측정의 관계를, 쌍별 400–500 kpc 외곽 고리를 빼기 전후로 비교했다. 연결된 5도
target/control 블록 5,000회 재표집 후 i 77쌍, i2 183쌍에서 raw 및 corrected 관계와
그 변화 구간이 모두 0을 포함했다. outer annulus의 공통 면적 중앙값은 약 17%이며
20% 미만인 쌍은 i 47/77, i2 111/183이다. 따라서 차감 결과는 sky 보정 또는 검출로
볼 수 없다. 코드는 [`analyze_satellite_outer_icl_baseline_association.py`](analyze_satellite_outer_icl_baseline_association.py),
입력·출력 해시와 정량 해석은 [민감도 평가](satellite_outer_icl_baseline_association_assessment_2026-10-20.md)에
있다. 위성 수와 ICL 면적 지지를 추가 조정한 분석도 같은 결론을 보였다. 결과는
`science_measurement=false`다.

2026-10-20에는 같은 77개 i·183개 i2 은하단 공통표본에서 위성 5–75 kpc의 7개
반경과 ICL의 관계를 추가 비교했다. 원시 상관과 raw-minus-outer 상관 차이의 95%
공간블록 구간은 14개 필터·반경 조합 모두 0을 포함했다. i의 5–10 kpc 구간에서만
외곽 차감 후 상관의 미보정 구간이 0보다 근소하게 높았으나 사후 다중비교를 보정하지
않았고, paired 차이 구간은 0을 포함한다. 이를 검출로 해석하지 않는다. 자세한 반경별
표와 공통표본 검증은 [고정 공통표본 평가](satellite_radial_icl_baseline_association_assessment_2026-10-20.md)에
기록했다.

이어서 위성 7개 반경별 raw·외곽차감 ICL 상관의 은하단 하나/연결 블록 하나 삭제
영향도를 계산했다. 외곽 차감 i-band 안쪽 상관은 삭제 후에도 양수였지만, 이는
신뢰구간이 아니며 기존 bootstrap 변화 구간은 0을 포함한다. 표와 해시는
[영향도 평가](satellite_radial_icl_association_influence_assessment_2026-10-20.md)에
있다. 결과는 여전히 `science_measurement=false`다.

`propagate_satellite_icl_measurement_noise.py`는 공분산이 있는 동일 57쌍(i 13,
i2 44)의 ICL·위성 관측값에 diagonal 및 finite-lag 측정 잡음을 각각 전파한다.
이는 공간 블록 재표본화와 합치지 않는 조건부 잡음 민감도다. i2의 44쌍 부분표본은
전체 68쌍과 다른 `r`을 보이므로 전체 관계로 확대하지 않는다. 산출물·한계와
체크섬은 [`icl_satellite_measurement_noise_mc_assessment_2026-10-19.md`](icl_satellite_measurement_noise_mc_assessment_2026-10-19.md)에 기록한다.

## ICL–위성 HSC 영역 고정효과 민감도 (2026-10-20)

같은 94쌍에서 redshift·richness를 보정한 뒤 표적/대조 HSC 연결영역 지시변수를 추가했다. i2 부분상관은 `r=0.29264`에서 `0.29105`로 거의 유지됐고, i는 `0.07898`에서 `0.13860`으로 변했으나 26쌍뿐이다. 표본의 target/control endpoint가 모두 같은 근사 영역에 배정되어 두 영역 효과 열은 중복되므로 하나만 적합했다. 이는 영역 간 평균 오프셋에 대한 민감도이며 검출이나 공분산 추정이 아니다. 방법과 한계는 [`satellite_icl_region_fixed_effects_assessment_2026-10-20.md`](satellite_icl_region_fixed_effects_assessment_2026-10-20.md)에 있다.

각 필터의 cluster 단위 leave-one-out에서 전체 지원표본 i2 `r`은 68회 삭제 모두 양수(`0.151–0.338`), i는 26회 중 한 번 부호가 바뀌었다. 더 엄격한 중심 적격 표본은 i2 15쌍(`0.254–0.433`), i 8쌍(`−0.787`–`−0.182`)이지만, 기존 공간 블록 구간은 각각 `[-0.059, 0.812]`, `[-1, 1]`로 넓고 0을 포함한다. 이는 선택별 영향도 민감도이지 검출이 아니다. 결과와 요약 해시는 [`satellite_icl_cluster_influence_assessment_2026-10-20.md`](satellite_icl_cluster_influence_assessment_2026-10-20.md)에 기록했다.

동일한 영향도 계산에 redshift·richness 조정 순위상관을 추가했다. 전체 지원 i2 순위 부분상관은 `0.139`(Pearson `0.293`), i는 `0.248`(Pearson `0.079`)다. i 전체 지원 Pearson은 한 은하단 제외 시 부호가 바뀌지만 순위 부분상관은 유지됐다. 이 차이는 관계의 크기·형태가 상관 정의에 민감함을 보여주는 기술적 민감도일 뿐 검출 검정은 아니다.

94개 결합 행은 200–300 kpc ICL 지원 224개와 독립 50–75 kpc 위성 프로파일 182개의 정확한 교집합이다(i 26, i2 68). 입력 관문별로 제외된 행은 i에서 ICL 전용 45/위성 전용 16, i2에서 85/72였다. 결합률은 i 26/71(36.6%), i2 68/153(44.4%)이며, 표적 면적지원 차이는 i `+0.60 SD`, i2 `+1.05 SD`다. i2에서 redshift `+0.32 SD`, richness `−0.33 SD` 차이도 나타나 현재 상관을 전체 800개 모집단에 일반화할 수 없다. [표본 선택 감사](icl_satellite_association_selection_assessment_2026-10-20.md)는 선택 가중치 보정이나 검출을 주장하지 않는다.

면적지원 불균형의 영향을 점검하기 위해 동일 94행에서 공통지원률 또는 표적·대조 지원률을 추가 공변량으로 둔 부분상관도 계산했다. i2 기준 `r=0.293`은 각각 `0.320`, `0.305`였고, 같은 5° 연결 공간블록의 짝지은 조정 전후 차이 분포는 모두 0을 포함했다. i는 26쌍뿐이며 기준 `0.079`가 `0.086`, `0.129`로 변했지만 차이 범위가 넓다. 이는 선택 편향 교정이나 검출이 아닌 지원 민감도다. 재현 코드·시험·원자료 해시는 [면적지원 조정 민감도](satellite_icl_coverage_adjustment_assessment_2026-10-20.md)에 기록했다.

2026-10-20 전체 `hsc_debug` 테스트 스위트는 282개 테스트와 6개 하위 테스트를 모두 통과했다. 이는 재현 도구의 코드 검증이며 ICL·렌즈 신호의 과학 승인이나 랜덤 차감 관문을 대신하지 않는다.

## ICL–위성 독립 하늘표본 확장 가능성 (2026-10-20)

`assess_icl_catalog_expansion_feasibility.py`는 기존 공간블록 power 계획과 CAMIRA 오프라인 후보 풀을 연결한다. 탐색 효과 시나리오 `r=0.293`, 30개 비교 보정 아래 80% 검정력에는 독립 연결 블록 180개가 필요하지만 현재는 33개다. 기존 S20A 카탈로그의 가장 넓은 후보 시나리오는 2,596개를 더할 수 있으나 새 target 5도 셀은 6개, 전체는 84개다. 현재 pair-block 구현은 target 셀을 공유하는 행을 연결하므로 84개는 연결 블록 수의 상한이고, 계획치보다 최소 96개 적다. 이 상한은 usable 또는 독립 블록을 보장하지 않는다. HSC 공식 PDR3 CAMIRA 목록도 PDR3(S20A) 기반 Wide/Deep/UD이며, 현재 입력은 이미 S20A Wide 별마스크 목록이다. 따라서 Wide 재조회는 새로운 독립 자료가 아니고, Deep/UD의 명목 면적(27+3.5 deg²)은 표본 보조에는 쓸 수 있어도 블록 부족을 채우지 못한다. 최신 실행 기록은 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_catalog_expansion_feasibility_2026-10-20_04/`에 있다. 재현 방법과 공식 HSC 출처는 [`assess_icl_catalog_expansion_feasibility.py`](assess_icl_catalog_expansion_feasibility.py) 및 산출물 `assessment.md`에 기록했다. 모든 결과는 `science_measurement=false`다.

별도 S23B 최종연도 CAMIRA Wide 별마스크 목록은 논문상 12,162개 후보와 1,220 deg²를 포함한다. 저장소 API 메타데이터로 정확한 자산을 `camira_s23b_wide_sm.tbl` (827,339 bytes; Git blob SHA-1 `066c2c8beca8f3a65186db7511fddca52d298908`)로 확인했지만, 카탈로그 바이트는 받지 않았고 해당 단일 파일의 사용자 승인을 기다린다. 승인 뒤에만 S20A 중복, 기존 고정 선택, target/control 연결 블록 및 현재 PDR3 영상 지원을 교차검사한다. S23B 명목 목록이나 면적은 PDR3 영상 지원 또는 새 독립 블록을 보장하지 않는다. 범위와 사후 감사 단계는 [S23B CAMIRA 수신 계획](s23b_camira_catalog_acquisition_plan_2026-10-20.md)에 기록했다.
