# HSC 디버그 절차

이 디렉터리의 스크립트는 ICL 스케일 스택을 실제 목록으로 준비하는 절차이다. 논문에 함께 제공한다. 실행 결과는 과학 측정이 아니다. `science_measurement` 는 거짓으로 남는다.

두 연구 갈래의 현재 표본, ICL–위성 연관성 검증, 표적 전용 weak-lensing 결과와 남은 과학 관문은 [통합 연구 평가](icl_wgl_integrated_research_assessment_2026-10-19.md)를 참조한다.

2026-11-06에는 상대적인 200–300 대 300–500 kpc ICL 프로파일 대비와
검출된 CAMIRA 멤버의 50–75 kpc 위성 중심 잔차 프로파일을 적색편이·풍부도
보정 후 비교했다. i 표본 20개/10개 블록의 부분상관은 약 −0.13이고,
i2 표본 47개/30개 블록은 약 +0.39–0.42이나 네 95% 블록 부트스트랩
구간이 모두 0을 포함한다. 검출 위성 표본은 완전한 위성 목록이 아니며,
측정오차·source-Poisson·coadd 공분산·하늘 전달·PSF 날개 불확실성을 전파하지
않았다. 세부 결과: [ICL–위성 외곽 프로파일 평가](icl_satellite_outer_profile_relation_assessment_2026-11-06.md).

2026-11-06에는 94개 지원쌍에서 ICL 신호가 위성 외곽광을 공간 블록 밖으로
예측하는지 leave-one-block-out 교차검증했다. 풍부도·적색편이 모형에 ICL을
추가했을 때 `i` 오차는 6.0% 커지고 `i2`는 3.5% 줄었으나, 두 블록 부트스트랩
구간 모두 0을 포함한다. 따라서 현재 표본에서는 새 공간 블록 예측 개선이
확인되지 않았다. 재현 스크립트와 결과는
[공간 교차검증 평가](icl_satellite_spatial_cross_validation_assessment_2026-11-06.md)에 있다.

2026-10-21에는 S20A·S23B CAMIRA 목록을 같은 `z_cl`·BCG 분광 기준으로 합친 오프라인 표본 상한을 계산했다. 6′ 최소 간격 뒤 1,620개 후보(22개 새 5° 셀)가 남지만 현재 joint-support yield를 적용하면 예상 i2 공통지원은 총 약 206행으로, 다중비교 표본력 시나리오의 323–804행보다 작다. 영상·마스크·PSF 지원은 미검증이다. 재현기는 [`audit_icl_cross_release_expansion.py`](audit_icl_cross_release_expansion.py), 결과와 해시는 [cross-release 표본 확장 평가](icl_cross_release_expansion_assessment_2026-10-21.md)에 있다.

2026-10-22에는 S23B 후보 1,333개를 현재 parent와 ICL–위성 결합표에 좌표 기반으로 연결해 행별 지원 장부를 만들었다. i2 표본 33개 target 셀에서 최대 67개 새 후보 셀을 더할 수 있어 총 100개 셀이 기하 상한이지만, 보수 180-block 설계 목표와 영상·대조장·PSF·위성 완전도 검증을 충족하지 않는다. 재현기와 결과 해시는 [S23B 공간 후보 장부](icl_s23b_spatial_candidate_ledger_2026-10-22.md)에 있다.
2026-11-06 tract graph 감사에서는 1,333개 후보 중 1,330개가 기존 표본이 차지한 세 HSC 연결 성분에 속하고, discovery 표본이 비어 있는 네 번째 성분에는 3개만 속했다. 따라서 S23B pool은 22개의 새 5° 셀을 제시하지만 실질적인 독립 HSC field support를 거의 늘리지 않는다. 재현기는 [`audit_s23b_candidate_tract_region_support.py`](audit_s23b_candidate_tract_region_support.py), 결과는 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_tract_region_support_2026-11-06_01/summary.json`이다. 좌표 배정은 근사이며 후보 영상은 아직 없다.

이 1,333개 후보의 PDR3 HSC-I 영상 패치 가능성을 검사하기 위해 중심과 proper 1-Mpc 원주 8점을 대상으로 27개 batch metadata-only SQL을 준비했다(총 11,997 probe, CAS 미제출). 이 별도 S23B query 범위는 승인 대기이며 기존 2,596개 S20A 계획과 합쳐서 취급하지 않는다. [S23B coverage metadata 계획](icl_s23b_coverage_metadata_plan_2026-10-22.md)

공통 제출기 `submit_icl_expansion_coverage_metadata.py --check --plan <plan-dir>`는 양쪽 계획을 오프라인 검증한다. S23B `--submit`은 승인한 JSON의 정확한 SHA-256이 별도로 주어지지 않으면 거부한다. metadata 반환 후 기존 9점 화면기는 S23B의 `S23B:<catalog_id>` 키도 허용한다.

별도 검정력 관문은 기존 탐색 표본을 독립 재현에 중복 사용하지 않고 S23B가 제공할 수 있는 새 target 셀 상한을 계산한다. 최대 67개 S23B 셀은 r=0.293 가정의 단일 고정 endpoint 80% 목표(92 block)에도 약 25개 부족하다. [S23B 독립 재현 표본력 평가](icl_s23b_replication_power_gate_2026-10-22.md)

미제출 S20A broad-redshift 후보 2,596개는 discovery 밖 51셀을 차지하지만 S23B와의 합집합에 단 2셀만 더해 총 69개 새 target 셀이다. 따라서 기존 S20A metadata 계획을 실행해도 독립 footprint 부족을 해소하지 못한다. [cross-release 표본력 상한](icl_cross_release_replication_capacity_2026-10-22.md)

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

`prepare_wl_random_shape_count_plan.py` and the guarded count submitter build
and execute the count-only support query for an
`astrafek-wl-random-centers-v1` manifest. The user separately approved the
subsequent shape-row export plan
`fc883444d2b23663d5f6034d773b76db948add29f1cda0c23384d17cfed2996a`.
All 124 batches completed: 3,100 random centers tied to 620 matched cluster
IDs, 9,010,324 center-source rows, and 2,962,522,918 CSV bytes. The exported
23 columns include calibrated shape, source coordinates, tomographic-bin and
B-mode fields; images and full photo-z PDFs were not requested.

`run_wl_random_shape_stack.py` verifies the complete hash-bound export and
center associations, reproduces the saved 800-lens target profile, and then
forms paired cluster-minus-random profiles. The count preflight found only 435
of the 620 random-matched lenses had both target-side and random-side source
support. The estimator and paired delete-one 5-degree equal-area sky-cell
covariance therefore use the same 435 lenses and their 2,175 associated
centers (2,046 with source rows, 129 count-zero). The final output is
exploratory: `science_measurement=false` because the random centers' full
aperture depth/mask/PSF/sky matching, connected HSC-region covariance,
block-width sensitivity, photo-z calibration uncertainty and reduced-shear
correction are not established.

The block-width comparison is documented in [`wl_random_shape_jackknife_sensitivity_assessment_2026-11-05.md`](wl_random_shape_jackknife_sensitivity_assessment_2026-11-05.md). The 2.5-degree run uses 154 cells versus 51 at 5 degrees; the profile central values match, while diagonal errors shift by −18% to +7%. This is a partial covariance sensitivity only; disconnected-region and random-aperture matching gates remain open. The full-estimator-per-block implementation took about 90 minutes for 2.5 degrees.

`analyze_wl_random_tract_region_influence.py` adds a leave-one-connected-tract-region-out check to the paired lens/random stack. The run uses all verified target and random source rows, recreates the saved seven-bin profile exactly, and assigns both lens/random centers and sources to approximate connected components of the HSC patch-QA tract-center graph. Only three components overlap the stack. Removing the dominant component removes 261/435 lenses and shifts ΔΣ by −28.6% at 0.1–0.2 Mpc and −26.8% at 0.4–0.7 Mpc. The coarse region graph is only an influence diagnostic and cannot provide a stable covariance. See [`wl_random_tract_region_influence_assessment_2026-11-06.md`](wl_random_tract_region_influence_assessment_2026-11-06.md). The science flag remains false.

`audit_wl_random_tract_partition_sensitivity.py` repeats the tract-center graph at 1.6°, 2.0°, 2.4°, and 3.0° adjacency. All thresholds preserve the same four components, with identical assignments for the 435 lenses and 2,175 random centers. The audit does not validate nearest-center boundaries or polygon geometry. Its output is `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_tract_partition_sensitivity_2026-11-06_01/summary.json`.

The completed stack summary is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_shape_stack_lens_reuse_2026-11-04_01/summary.json`
(SHA-256 `83c5484454033b2652adcba5bd0699f28082fa6983eefa960613da75491bafc8`).
The scientific assessment and radial-profile figure are documented in
[`wl_random_shape_stack_assessment_2026-11-04.md`](wl_random_shape_stack_assessment_2026-11-04.md).
Recreate the figure from the saved summary with:

```bash
python3 analysis/icl_wgl_membership/hsc_debug/plot_wl_random_shape_stack.py \
  --summary /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_shape_stack_lens_reuse_2026-11-04_01/summary.json \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_shape_stack_lens_reuse_2026-11-04_01/random_subtracted_profile.png
```

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
the actual shared valid pixels in the saved K2.5 ICL and SEP/Kron satellite masks. Its support
input is the per-annulus ICL table, not the smaller common-stack table. With this scope it covers
94 cluster/filter pairs and 235 satellite profiles; the fixed 74-pair association sample is fully
contained. It can run on one or more matched clusters with `--cluster-id`, reads existing images
only, and keeps `science_measurement=false`. See
[`satellite_icl_pixel_overlap_2026-10-18.md`](satellite_icl_pixel_overlap_2026-10-18.md) and
the corrected full-support and shared-pixel sensitivity assessment
[`icl_satellite_full_support_overlap_assessment_2026-11-03.md`](icl_satellite_full_support_overlap_assessment_2026-11-03.md).
The 94-pair output is at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_icl_pixel_overlap_2026-11-03_01/`.

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

`assess_icl_power_family_sensitivity.py` extends that design approximation
across 1, 6, or 30 hypothetical tests and `r=0.2, 0.293, 0.4`, comparing
requirements with the current 33 connected blocks and the audited candidate
pool upper bound of 84. For `r=0.293`, a single prospective primary test needs
92 independent blocks and has approximate power 0.765 even at the 84-block
upper bound; 30 tests need 180 blocks. At `r=0.2`, one test needs 196 blocks.
The single-test scenario cannot be applied retroactively to the current
exploratory analysis: any primary endpoint must be fixed before new data are
examined. The output and input hashes are in
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_power_family_sensitivity_2026-10-20_01/`.
This Fisher-z design sensitivity does not estimate usable independent blocks or
establish a detection; `science_measurement=false` remains in force.
The prospective follow-up design candidate in
[`icl_satellite_confirmatory_replication_design_candidate_2026-10-20.md`](icl_satellite_confirmatory_replication_design_candidate_2026-10-20.md)
sets the current i2 200–300 kpc ICL / 50–75 kpc satellite observable pair as a
reviewable candidate for an independent sample, and gives 196 usable blocks as
the conservative `r=0.2`, one-test Fisher-z planning target. It is not approved
or preregistered and cannot retroactively validate the current exploratory sample.

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
The user approved exactly the four `_v_1p1` distribution/metadata files
(549,818 bytes total). They are stored outside Git under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/desi_y3_hsc_nz_abc70af2_2026-10-20/`,
with pinned commit, source URLs, byte counts, SHA-256 values, and schema recorded
in `manifest.json`. Only these four files were fetched. Their NPZ schema contains
mean/median and pointwise lower/upper curves for four bins, including
magnification-corrected variants, but no posterior samples. The completed
comparison reuses the verified SACC baseline and applies a matched 52-cell
spatial jackknife to the four DESI mean curves. Median `DeltaSigma` shifts
range from +12.5–14.8% in source bin 1 and −1–3.5% in bins 2–4. The 0.1–0.7 Mpc
tomographic inconsistency becomes less extreme but remains; random-lens
subtraction is unavailable. Detailed values, caveats, and hashes are in
[`wl_desi_nz_distribution_sensitivity_assessment_2026-10-20.md`](wl_desi_nz_distribution_sensitivity_assessment_2026-10-20.md),
with independent validator output in the GPFS result directory.
`audit_wl_published_nz_shift_sensitivity.py` remains a no-download approximation:
it translates the local SACC source-bin 3
and 4 curves toward higher redshift by the paper's central mean shifts (+0.043
and +0.050), renormalizes them, and repeats the same joint jackknife. The
[paper's v3](https://arxiv.org/abs/2511.18133v3) reports these shifts relative
to Rau+2022; this constant-translation sensitivity is not the authors'
calibrated spline distributions or their uncertainty.
The same Astrafek audit now also samples the published marginal 68% shift
interval endpoints: higher-redshift translations of 0.021–0.063 for bin 3 and
0.038–0.062 for bin 4, using the central point and all four rectangle corners.
The grid is a sensitivity envelope, not a joint posterior. The author-published
distribution comparison evaluates both scale cuts, the mean and
magnification-corrected means with spatial jackknife, and the pointwise lower /
upper curves as profile sensitivities. It does not treat those marginal curves
as posterior samples. The reproducible
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


The approved random-lens support gate is now running from the existing
3,100-center, 620-lens plan. It queries only integer distinct-shape-source counts
inside 3 proper Mpc, in 124 batches of 25 centers. Checkpoint and verified CSVs
are under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_shape_count_results_lens_reuse_2026-11-03_01/`.
Use `count_status.json` as the live state; completed batch IDs and checksums are
retained for safe resumption. This approval does not include exporting random
source rows. The plan and query guard are
[`prepare_wl_random_shape_count_plan.py`](prepare_wl_random_shape_count_plan.py)
and [`submit_wl_random_shape_counts.py`](submit_wl_random_shape_counts.py).

The offline support audit is
[`audit_wl_random_shape_count_support.py`](audit_wl_random_shape_count_support.py).
At the 2026-11-03 interim snapshot it verified 68/124 completed batches
(1,700 centers); 600 centers had zero shape sources, and 107/340 fully covered
matched lenses had zero support across all five centers. These first manifest-
ordered lenses are not representative of the full sample, so the support rates
are provisional. The audit and limits are documented in
[`wl_random_shape_count_support_assessment_2026-11-03.md`](wl_random_shape_count_support_assessment_2026-11-03.md).
The same audit can compare completed random-center counts to the already
verified lens-position count preflight with `--lens-preflight-dir`. At the
104-batch snapshot, 118/134 matched lenses with zero lens-position counts also
had zero counts at all five random centers; among 386 count-positive lenses,
309 had positive counts at every center and 34 had none. Their per-lens random-
mean/lens-position count ratio has median 0.850 (q16–q84: 0.302–1.115). This
remains a partial,
pre-quality-cut support check, not a final retention result. The count-only
query remains separate from any random source-row export.

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
verifies every completed coordinate CSV against the immutable count plan,
per-lens candidate totals, row schema, and checkpoint SHA-256 without contacting
CAS or inferring process liveness. The 2026-10-20 audit verified 2/32 batches;
the 2026-10-22 recheck again verified the same 33,514/780,191 coordinate rows.
The former authenticated submitter PID 1049265 is now absent on `lageunha`,
but HSC CAS job 1062312's last persisted status is still `running`, so its
current state is unknown. Do not resubmit that job or start remaining batches
until its authenticated CAS status is checked. The dated audit and exact hashes
are in [`wl_random_coordinate_recheck_2026-10-22.md`](wl_random_coordinate_recheck_2026-10-22.md).
No CAS request was made, and `science_measurement=false` remains in force.

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

2026-10-20 09:08 KST 재확인에서도 PID 1049265는 살아 있었다(경과 21시간 59분,
`hrtimer_nanosleep`). 체크포인트 시각은 여전히 2026-10-19 11:07:53 KST였고,
2/32 배치(33,514/780,191행) 완료, batch 002 CAS 작업 `1062312` 제출 상태로
변화가 없었다. 최상위 `error_type=HTTPError`가 저장되어 있으나 발생 시각과 대상
요청이 없어 특정 CAS 작업의 종료 오류라고 단정할 수 없다. CAS 상태 poll 시각은
기록되어 있지 않다. 동일 프로세스와 제출 건을 유지했고 중복 제출하지 않았다.

09:13 KST에 같은 승인 좌표 계획의 오프라인 checkpoint audit를 다시 실행했다.
32개 중 기존 2개(33,514/780,191행)만 행 수·계획별 lens count·SHA-256이 일치했다.
체크포인트 SHA-256은 앞선 감사와 동일한
`84cccc6443ec8d2fc8bee1611ab6b9b57a3aabdc90978fcdfc558151d2bb92b2`이고,
이번 실행은 CAS에 연결하지 않았다. 보고서는
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_coordinate_checkpoint_audit_2026-10-20_10/`에
있으며 JSON SHA-256은 `d5072ee7cecd4ec4798892dc6d3be4acf1d916738e0ff65564666234a3ec8212`다.

2026-10-20 09:18 KST에 SSH로 다시 확인했을 때 같은 PID 1049265가 살아 있었고,
프로세스는 60초 대기 호출(`hrtimer_nanosleep`) 중이었다. 체크포인트 수정 시각은
여전히 2026-10-19 11:07:53 KST이며 batch 000·001만 완료(33,514/780,191행),
batch 002는 CAS 작업 `1062312` 제출 상태였다. 최상위 `error_type=HTTPError`도
남아 있지만 요청 시각과 응답이 없어 CAS 작업 상태를 판정할 수 없다. 이 확인은
프로세스와 기존 로컬 체크포인트만 읽었으며 CAS에 재조회·재제출하지 않았다.
프로세스가 살아 있어 종료하거나 같은 배치를 다시 제출하지 않는다.

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

The 2026-11-03 63-pair extension is a historical subset because it used the
137-row common-stack table to select per-annulus support. The corrected full-support rerun covers
94 pairs and all 74 fixed association pairs. Actual masks show shared valid ICL/satellite pixels
in 140/206 association profiles. Excluding the union of the shared satellite target-shell pixels
reduces the i2 redshift/richness-adjusted partial correlation from 0.217 to 0.176 in the common
coverage sample; the paired 95% block-bootstrap change is [−0.065, −0.018]. This is a strong
measurement-overlap sensitivity, not a causal correction or science detection; details and hashes
are in [`icl_satellite_full_support_overlap_assessment_2026-11-03.md`](icl_satellite_full_support_overlap_assessment_2026-11-03.md).

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

## 초기 20개 Gaia 별 우표의 프레임 가중 PSF 파일럿 (2026-10-20)

초기 pilot은 20개 별 우표와 15개 HSC 프레임을 프레임 동등 가중 방사 프로필로 쌓았다. 네 설정 모두 48–96픽셀 외곽 평균이 음수였으나 이 결과는 뒤에 실행한 333개 stamp/227개 프레임의 확장 평가로 대체됐다. 확장 표본에서는 부호가 local-sky/source-mask 설정에 따라 달랐다. 자세한 확장 결과는 [프레임 가중 PSF 평가](empirical_psf_frame_weighted_stack_assessment_2026-10-20.md)에 있다. 두 결과 모두 `science_measurement=false`, `psf_approved=false`다.

## ICL–위성 독립 하늘표본 확장 가능성 (2026-10-20)

`assess_icl_catalog_expansion_feasibility.py`는 기존 공간블록 power 계획과 CAMIRA 오프라인 후보 풀을 연결한다. 탐색 효과 시나리오 `r=0.293`, 30개 비교 보정 아래 80% 검정력에는 독립 연결 블록 180개가 필요하지만 현재는 33개다. 기존 S20A 카탈로그의 가장 넓은 후보 시나리오는 2,596개를 더할 수 있으나 새 target 5도 셀은 6개, 전체는 84개다. 현재 pair-block 구현은 target 셀을 공유하는 행을 연결하므로 84개는 연결 블록 수의 상한이고, 계획치보다 최소 96개 적다. 이 상한은 usable 또는 독립 블록을 보장하지 않는다. HSC 공식 PDR3 CAMIRA 목록도 PDR3(S20A) 기반 Wide/Deep/UD이며, 현재 입력은 이미 S20A Wide 별마스크 목록이다. 따라서 Wide 재조회는 새로운 독립 자료가 아니고, Deep/UD의 명목 면적(27+3.5 deg²)은 표본 보조에는 쓸 수 있어도 블록 부족을 채우지 못한다. 최신 실행 기록은 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_catalog_expansion_feasibility_2026-10-20_04/`에 있다. 재현 방법과 공식 HSC 출처는 [`assess_icl_catalog_expansion_feasibility.py`](assess_icl_catalog_expansion_feasibility.py) 및 산출물 `assessment.md`에 기록했다. 모든 결과는 `science_measurement=false`다.

S23B 최종연도 CAMIRA Wide 별마스크 표는 승인된 단일 파일만 수신했다. 크기는 827,339 bytes이고 Git blob SHA-1 `066c2c8beca8f3a65186db7511fddca52d298908`, SHA-256 `5af057028feab7976cd3ed389a7246cbe73add8bbb96891101f0f692b1de4d71`로 검증했다. 12,162행 중 4,492행이 S20A 카탈로그와 1초각 이내에서 겹쳤고, 기존 800개 영상 표본의 위치와 1초각 이내로 겹친 607행 중 559행이 현재 적색편이/BCG 분광 컷도 통과했다. 같은 컷과 S20A 중복 제거 뒤에는 1,440개 후보, 기존 표본과 신규 후보 사이 6분각 간격을 적용한 richness 우선 탐욕 선택에서는 1,333개 후보가 남았고, 새 5도 좌표 셀은 22개(기존과 합계 100개)였다. 단, 이미 검증된 PDR3 영상 manifest는 기존 표본 800개만 포함하며 새 위치의 coadd·mask·PSF·깊이·전체 aperture 및 matched-control 지원은 확인되지 않았다. 따라서 이 값은 목록 수준 상한이지 PDR3 영상 가능성 또는 독립 공분산 블록 수가 아니다. 전체 절차·해시는 [S23B CAMIRA 감사 결과](s23b_camira_catalog_audit_2026-10-20.md)에, 재현 스크립트는 [`acquire_s23b_camira_catalog.py`](acquire_s23b_camira_catalog.py) 및 [`audit_s23b_camira_catalog.py`](audit_s23b_camira_catalog.py)에 기록했다. 원본과 요약은 `/gpfs/kjhan/Astrafek/data/catalogue_s23b_camira_2026-10-20_01/`, 감사 표는 `/gpfs/kjhan/Astrafek/data/catalogue_s23b_camira_audit_2026-10-20_02/`에 있으며 모두 `science_measurement=false`다.

2026-10-20 전체 Gaia 별 우표 프레임 가중 스택은 333개 stamp/227개 HSC-I 프레임에서 완료했다. 64–96픽셀 고리의 부호는 국소 배경 및 source-mask 선택에 따라 달랐고 양의 확장 PSF 날개는 검증되지 않았다. 원래 20-stamp 시험 결과는 보존하되, 전체 표본 평가를 현재 해석에 사용한다. 수치와 해시는 [프레임 가중 PSF 평가](empirical_psf_frame_weighted_stack_assessment_2026-10-20.md)의 확장 표본 절에 있다. `science_measurement=false`, `psf_approved=false`.

같은 333개 별 우표에서 319개는 매칭 control의 같은 하위 필터 `i`/`i2`와 401×401 전체 footprint를 지원한다(217개 표적 프레임). 14개 후보는 하위 필터가 달라 제외했다. 이 WCS/header preflight는 영상 픽셀을 읽지 않았다. 다음 짝지은 control-profile 진단의 입력 지원 결과이며 PSF 날개나 ICL sky 보정은 아니다. 상세 방법·산출물·SHA-256은 [matched-control 지원 평가](empirical_psf_matched_control_support_assessment_2026-10-20.md)에 있다.

2026-10-20 matched-control profile 진단은 같은 하위 필터의 319개 표적 별 우표(217 프레임)와 대응 대조장 위치를 세 가지 mask/coverage 설정으로 비교했다. 중심-source 제한을 외곽 링에 적용하지 않아도 가장 큰 outer-bin 표본은 i2 64–96픽셀의 12개 프레임뿐이며, 수치는 sky annulus 설정에 민감하다. 대조점은 별이 없는 placebo이고 target에는 ICL 곡률이 남으므로 이를 PSF wing으로 분리할 수 없다. 요약과 세 실행의 해시는 [matched-control profile 평가](empirical_psf_matched_control_profiles_assessment_2026-10-20.md)에 기록했다. `psf_approved=false`, `science_measurement=false`.

Gaia↔HSC 후보 연결의 실제 관측 epoch를 확인하기 위해 기존 표적 coadd 227개 전체의 FITS header를 감사했다. 표준 date/MJD 카드가 0/227이어서 실제 epoch의 Gaia cross-match를 이 파일들만으로 확인할 수는 없다. PSF 우표 중심 자체는 Gaia 위치가 아니라 HSC 검출의 `X_IMAGE/Y_IMAGE`를 사용한다. 픽셀과 HSC 서비스는 접근하지 않았다. 재현 스크립트·해시는 [HSC coadd epoch 감사](empirical_psf_hsc_coadd_epoch_audit_2026-10-20.md)에 있다.

같은 기존 333개 별 우표의 품질 감사에서 중심 12픽셀 안의 포화 비트는 0개, `BRIGHT_OBJECT` 겹침은 3개였다. 중심 마스크 개방률 99% 기준을 통과한 표본은 Kron 2.5에서 290개/201프레임, Kron 3.5에서 242개/179프레임이다. 그 중심 적격 표본도 48–96픽셀 고리의 개방률 중앙값이 각각 0.69–0.72 및 0.58–0.60이며, 다수에 검출 이웃이 있어 날개 스택 해석에 혼잡도·마스크 검증이 필요하다. 상세 수치·해시는 [별 우표 품질 감사](empirical_psf_star_stamp_quality_assessment_2026-10-20.md)에 기록했다. `psf_approved=false`, `science_measurement=false`.

약한렌즈 random 중심 승인 count-only 결과 32개 배치를 모두 다시 검증해 800개 렌즈별 원시 후보 수를 요약했다. 각 렌즈는 계획한 200개보다 많은 364–1,443개(중앙값 1,000개)를 보유하지만, 이는 전체 CAMIRA aperture 중첩과 영상/형상 품질을 걸러내기 전 숫자다. 첫 50개 부분 기하 검사에서 6.35%만 통과해, 전체 geometry screen이 다음 병목이다. 고정된 표적 표본에서 후보 수와 렌즈 적색편이의 Spearman 상관은 0.595이며 유효 annulus 면적으로 나눈 후보 밀도에서는 0.080으로 줄었다. 입력 및 per-lens 산출물 해시는 [random 중심 후보 수 감사](wl_random_coordinate_candidate_count_assessment_2026-10-20.md)에 기록했다. 새 CAS 질의는 없고 `science_measurement=false`다.

333개 Gaia 별 우표의 물리적 반경을 실제 227개 cluster frame의 적색편이에서 재계산했다. 현재 경험 PSF 프로필은 최대 96픽셀, 중앙 물리반경 약 88 kpc로 ICL 측정 200–300 kpc에 도달하지 않는다. 원본 401×401 우표의 inscribed circle도 200 kpc에 닿는 프레임은 52/227, 300 kpc는 0/227이다. 따라서 현재 경험 별 스택은 satellite 50–75 kpc 척도와 기하학적 겹침은 있지만 ICL 고리의 PSF wing을 직접 보정할 수 없다. 재현 절차와 체크섬은 [PSF 물리 반경 coverage 감사](empirical_psf_physical_radius_coverage_assessment_2026-10-20.md)에 있다. `psf_approved=false`, `science_measurement=false`.

공개 HSC PDR3 확장 PSF 연구는 현재 별 우표의 반경 한계를 보완할 수 있는 `coadd/bg`용 5.6분각 모델을 제공한다. HSC 공식 처리문서는 비슷하지만 동일하지 않은 i/i2 노출이 공합성된다고 설명한다. Zenodo v4.0.0의 대칭·비대칭 HSC-I FITS 두 파일(합계 19,618,560 bytes)을 사용자 승인 범위로 수신해 크기·MD5·SHA-256을 검증했다. 두 4001×4001 모델은 총 광량을 거의 1로 정규화하며, z=0.30–0.55에서 200–300 kpc 고리에 점광원 광량의 0.222–0.242%가 들어간다. 독립 검증기가 원본 FITS와 aperture 수치를 재계산했다. 중심 43×43픽셀은 기존 744개 위치별 HSC PSF 코어 중앙값보다 FWHM이 21.5% 넓다. 1.5–3.3″ 중첩고리에서 public wing 대 위치별 core scale 중앙값은 0.739(5–95백분위 0.467–1.045), profile 잔차 RMS 중앙값은 12.6%였다. 따라서 hybrid 감도분석의 입력은 얻었지만 위치별 외곽 보정은 아직 승인되지 않았다. 6개 합성시험의 [`audit_public_extended_psf_reference.py`](audit_public_extended_psf_reference.py), 독립 검증기 [`validate_public_extended_psf_reference_audit.py`](validate_public_extended_psf_reference_audit.py), 코어 비교기 [`compare_public_extended_psf_to_picker_cores.py`](compare_public_extended_psf_to_picker_cores.py), 중첩 감사기 [`audit_public_psf_picker_core_overlap.py`](audit_public_psf_picker_core_overlap.py), 수치·한계는 [공개 확장 PSF 적합성 평가](extended_psf_public_reference_assessment_2026-10-21.md)에 기록했다. `psf_approved=false`, `science_measurement=false`다.

## Random coordinate radial-yield audit (2026-10-20)

첫 두 승인 좌표 batch의 50개 렌즈에서 inner exclusion 경계부터 30′까지 실제 query 폭은 7.14–8.76′(중앙값 8.10′)였다. 33,514 후보 중 2,129개(6.35%)가 CAMIRA aperture overlap·중복 제거 뒤 남았고 렌즈별 지원은 편중됐다. 이 결과는 첫 50개 비대표 pilot이며 random-center 승인이나 shape 반출 결과가 아니다. 보충 30–60′ 영역의 이상적 구면 면적은 현재 annulus의 중앙값 6.42배지만 HSC footprint/depth는 반영하지 않았다. 재현 스크립트와 산출물 hash는 [radial-yield 평가](wl_random_coordinate_radial_yield_assessment_2026-10-20.md)에 있다. 준비된 count-only 쿼리는 별도 승인 전 제출하지 않는다.

## ICL–위성 표본수 계획 민감도 (2026-10-21)

현재 800개 parent 은하단 중 joint ICL 200–300 kpc / 위성 50–75 kpc 지원은 i 26개, i2 68개다. 부분상관 가정치 0.10/0.20/0.293, redshift·richness 두 공변량, 168개 가족 Bonferroni 보정에서 80% 검정력에 필요한 독립 유효 행은 각각 2,849/536/323이다. 현재 joint-support 수율을 유지한다고 단순 외삽하면 i2 parent 규모는 약 33,518/6,306/3,800개이며, 설계효과 1.5에서는 약 50,283/9,459/5,706개다. `r=0.293`은 선택된 탐색 endpoint 점추정치라 보수적이지 않을 수 있다. 이 계산은 영상 지원·PSF/sky 체계오차를 해결하지 않으므로 추가 HSC 표본이 현재 측정을 자동으로 유효하게 만들지 않는다. 재현기·가정·제한은 [표본수 계획 평가](icl_sample_size_planning_assessment_2026-10-21.md)에 있다.

같은 S23B 카탈로그 오프라인 감사의 신규 후보 상한 1,333개를 현재 표본에 더해도 총 2,133개다. 기존 i2 joint-support 수율을 적용한 기대 지원은 총 약 181행이어서 가족 보정 `r=0.293` 시나리오의 323행 목표에도 못 미친다. 후보는 PDR3 영상·매칭 대조 지원이 확인되지 않아 실제 usable 표본으로 세지 않았다. 계산 스크립트와 caveat는 [표본수 평가](icl_sample_size_planning_assessment_2026-10-21.md#현재-공개-후보-풀과의-비교-2026-10-21), 산출물은 `icl_s23b_expansion_feasibility_2026-10-21_02/summary.json`에 있다.

로컬 S20A CAMIRA 부모 카탈로그에도 cluster photo-z만 쓰는 넓은 `0.20–0.70` 선택에서 2,596개 후보가 있지만, 5도 좌표 셀은 6개만 늘어난다. 현 i2 지원률을 적용한 총 기대 지원은 약 289행으로, 사후 효과 시나리오의 323행보다 작다. redshift 선택과 PDR3 이미지 지원이 달라 실제 확장 표본으로 인정하지 않았다. 상세는 [확장 표본수 평가](icl_sample_size_planning_assessment_2026-10-21.md#현재-공개-후보-풀과의-비교-2026-10-21)에 적었다.

위 후보 2,596개의 local CAMIRA member 좌표를 130,148개 검증하고, 중심 및 1-Mpc 원주 8개 위치에 대한 PDR3 HSC-I mosaic/patch-QA metadata SQL 계획을 52개 배치로 생성했다. 총 23,364 probe positions이며 **CAS 제출은 하지 않았다**. 새 영상·source row·shape·photo-z는 요청하지 않는다. 기존 승인 좌표 CAS 프로세스가 끝나지 않아 겹치는 요청은 보류한다. 세부 열·입력 및 산출물 해시·한계는 [coverage query plan](icl_expansion_coverage_query_plan_2026-10-21.md)에 기록했다.

가족 전체를 반경별로 검토하도록 [`plot_icl_satellite_family_matrix.py`](plot_icl_satellite_family_matrix.py)가 168개 hash-pinned 결과를 4패널로 그린다. 각 셀에 부분상관과 cluster 수를 표시하고 색상은 `−0.5…+0.5`로 고정했다. 168개 전체에서 Holm 또는 max-T `p<0.05`인 셀은 0개다. PNG/PDF와 입력·출력 해시는 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_radial_family_plot_2026-10-21_01/plot_manifest.json`에 있다. 이 그림도 탐색 진단이며 `science_measurement=false`를 유지한다.

## Weak-lensing richness split: redshift-matched sensitivity (2026-11-03)

The 400/400 parent-sample richness split has median lens redshifts 0.442 and
0.427. A maximum-cardinality 1:1 match within `|delta z|<=0.02`, restricted to
matched lenses with source-row support on both sides, retains 251 lens pairs.
The upper-minus-lower `DeltaSigma` difference from 0.2–3 Mpc remains similar
to the unmatched split, while the innermost-bin difference is broader. This
only checks lens-redshift imbalance: 141/392 candidate matches lacked complete
source-row support, and random-lens subtraction, matched-field variance, and
full calibration systematics remain open. Keep `science_measurement=false`.
See [`wl_redshift_matched_richness_assessment_2026-11-03.md`](wl_redshift_matched_richness_assessment_2026-11-03.md)
and reproduce with [`audit_wl_redshift_matched_richness_split.py`](audit_wl_redshift_matched_richness_split.py).

## ICL–satellite shell-excluded injection response (2026-11-03)

The fixed 74-pair association sample intersects the approved compact-PSF joint
support in 48 frames (12 `i`, 36 `i2`). On these exact frames, the supported
satellite 50–75 kpc target-pixel union was excluded from the ICL 200–300 kpc
measurement while running diffuse-only, satellite-only, and combined PSF-core
injections. Diffuse-to-ICL response is stable; paired frame-level changes in
satellite-to-ICL response span zero and are broad. These synthetic responses do
not justify a scalar correction. Extended PSF wings, source Poisson noise,
coadd covariance, upstream sky transfer, mask re-selection, matched-control
variance, and full spatial covariance remain open. Keep
`science_measurement=false`; see
[`icl_satellite_shell_exclusion_transfer_assessment_2026-11-03.md`](icl_satellite_shell_exclusion_transfer_assessment_2026-11-03.md).

The same 48 exact association/PSF frames were then rerun with the public HSC-I
extended wing spliced onto the local PSF core (2.4 arcsec splice, 0.6 arcsec
blend, 0.65 wing amplitude). In the matched fixed-mask comparison, diffuse
response shifts are small and satellite-envelope response shifts have broad
frame scatter with q16–q84 spanning zero; no scalar response correction is
supported. The model point-source wing flux left in the 200–300 kpc annulus
falls slightly when the supported satellite shells are excluded, while median
coverage also falls by 2.46 percentage points in `i` and 3.70 points in `i2`.
The `i2` wing is an HSC-I proxy, not an i2 calibration. The independent
checksum/sample audit is
[`summarize_icl_hybrid_shell_exclusion.py`](summarize_icl_hybrid_shell_exclusion.py),
with full numbers and limitations in
[`icl_satellite_shell_exclusion_hybrid_psf_assessment_2026-11-03.md`](icl_satellite_shell_exclusion_hybrid_psf_assessment_2026-11-03.md).
Keep `science_measurement=false`, `psf_approved=false`, and
`observational_validation_passed=false`.

## Joint sky and coadd-covariance sensitivity (2026-11-03)

For the same fixed 74-pair association set, an offline four-model comparison
adds the paired 600–900 kpc outer-residual proxy and/or a local target-minus-
control coadd sigma-inflation proxy to redshift and richness. In `i2` the partial
correlation changes from 0.239 to 0.223 with both proxies; in `i` it changes
from 0.112 to 0.098. The paired spatial-block bootstrap intervals are broad,
especially for `i`; they are not detection intervals. The outer annulus is not
pure sky and the local lag statistic is not radial-estimator covariance. No
residual was subtracted and no covariance correction was applied. Reproduce
with [`analyze_icl_satellite_joint_sky_covariance.py`](analyze_icl_satellite_joint_sky_covariance.py);
see [`icl_satellite_joint_sky_covariance_assessment_2026-11-03.md`](icl_satellite_joint_sky_covariance_assessment_2026-11-03.md).
Keep `science_measurement=false`.

## S23B spatial-block uncertainty for the target/control ICL stack (2026-11-07)

Run the saved pair-profile block bootstrap with:

```bash
python3 analysis/icl_wgl_membership/hsc_debug/analyze_s23b_icl_spatial_blocks.py \
  --replicas 5000 \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_spatial_blocks_2026-11-07_02
```

The runner pins the S23B target/control manifests and exploratory stack,
validates the pair-to-control join, and refuses to overwrite an output
directory. It reads saved annular profiles only, then writes the spatial-block
table, leave-one-block-out contrasts, block-bootstrap intervals, figure, and
assessment. Six layouts combine 2.5/5/10-degree RA-Dec widths with zero/half
cell phase shifts. Treat them as spatial dependence sensitivities, not as a
survey-region covariance. Result SHA-256:
`423bb5e3ad52f3b349dedc74b00a2c05ebb4589e20402abb68b4e99d17199afc`.
The S23B target/control stack remains exploratory and
`science_measurement=false`.

### 2026-11-07 S23B target-minus-random-control ICL path

The user approved exactly 1,041 HSC PDR3 `pdr3_wide` HSC-I `coadd/bg`
IMAGE/MASK/VARIANCE control cutouts, pinned to
`s23b_icl_control_image_plan_2026-11-07_02/control_manifest.json`
(SHA-256 `ae9eb406cfc179d9980031bb67f85f8ee71c4580ba59e97713278823db5c06df`).
The resumable download is in
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_control_image_plan_2026-11-07_01/control_download_01/`.
The 628th response had the valid, header-declared `MP_NOT_DEBLENDED` /
`MP_UNMASKEDNAN` bit positions swapped relative to the other files. HSC
documentation says to read the mask header because bit/event correspondence
can vary. `download_pilot.validate_fits` now checks the known plane names and
unique bit assignment; mask generation and control masking use each FITS
header's named mapping. Do not replace or expand the approved control list.

After transfer, perform these steps in order:

1. Run `audit_s23b_icl_control_images.py` for complete file, HDU, filter, and
   WCS integrity.
2. Run `prepare_full_cluster_mask_support.py --download-dir ...` with the
   approved v2 manifest, its pinned SHA, and that audit summary, writing masks
   under `s23b_icl_mask_support_2026-11-07_01/control/`.
3. Run `audit_s23b_icl_pair_common_support.py` to report target/control/common
   radial mask coverage on the proper-kpc grid. It reads FITS headers and mask
   planes, not science-image pixel values or flux.
4. Review filter counts and support yield. Only then run
   `stack_s23b_icl_target_control_2d.py`, which requires the completed
   mask-support audit and control masks and keeps `i` and `i2` separate.

The pre-stack and flux-stack scripts exclude the four target images already
identified by the zero-usable-pixel audit, and the flux stack also excludes
any pair with zero common mask support inside 1,000 proper kpc. No minimum
coverage cut has been selected. Outputs stay `science_measurement=false`; no
extinction, sky-transfer, PSF-wing, or coadd-covariance correction is implied.

### 2026-11-07 S23B image-transfer continuation

The approved target-image transfer remains active in session 5633. At the
latest checkpoint it has validated 850/1,043 files (47,507,889,600 bytes;
685 `i2`, 165 `i`). The separately approved control transfer uses the corrected
v2 manifest SHA-256
`ae9eb406cfc179d9980031bb67f85f8ee71c4580ba59e97713278823db5c06df`; session
56526 has validated 86/1,041 files (5,233,340,160 bytes), including the
replacement for the earlier HTTP 404. Both transfers retain the exact approved
scopes and `science_measurement=false`.

After the control transfer completes,
[`audit_s23b_icl_control_images.py`](audit_s23b_icl_control_images.py) will
verify every checksum, FITS/HDU layout, HSC-I filter, WCS footprint, and
one-to-one target mapping. It is an integrity audit only. The target and
control image audits must both pass before any S23B stack; pixel-level usable
coverage, residual sky, BCG centering, member positions, noise covariance, and
extended PSF behavior remain separate scientific gates.
The S23B methods paper says the algorithm selects a BCG and recomputes richness
and redshift around it; the released cluster table lists RA/Dec without a
separate BCG coordinate pair. Therefore the v1 manifest's “not a BCG position”
label is unverified, not a confirmed offset. Its approved hash and coordinates
remain untouched; the manifest generator's wording is corrected for future
plans. See https://arxiv.org/abs/2512.13954 and
https://github.com/oguri/cluster_catalogs/tree/main/hsc_s23b_camira.

## ICL aperture-specific finite-lag covariance (2026-11-03)

[`analyze_icl_annulus_aperture_covariance.py`](analyze_icl_annulus_aperture_covariance.py)
contracts the independently audited per-frame 9×9 coadd lag kernels with the
actual valid-pixel weights in the saved 200–300 kpc target and control annuli.
All 224 supported pairs (71 `i`, 153 `i2`) and 448 frames completed; the
reconstructed means and coverage were checked against the saved profiles. The
finite-lag pixel-noise stack sigma is 0.07869 nJy arcsec⁻² in `i` and 0.07356
in `i2`, respectively 1.254 and 1.225 times the diagonal-only estimate. This
is about 7.56% and 6.94% of the corresponding spatial-block-bootstrap interval
half-widths. It is an approximate pixel-noise floor, not the total estimator
covariance or a profile correction. Stationarity within each aperture is
assumed; long-range lags, sky transfer, source Poisson noise, PSF wings, mask
uncertainty, and cross-frame covariance remain omitted. Full results and
hashes are in
[`icl_aperture_covariance_200_300_assessment_2026-11-03.md`](icl_aperture_covariance_200_300_assessment_2026-11-03.md).
Keep `science_measurement=false`.

## ICL synthetic-source Poisson-noise input gate (2026-11-04)

[`audit_icl_poisson_noise_inputs.py`](audit_icl_poisson_noise_inputs.py)
header-audited all 800 target and 798 matched-control coadds. Every file has
the recorded IMAGE/MASK/VARIANCE planes and byte size, but none of the 1,598
headers carries an explicit gain or exposure key. The coadd variance plane and
`FLUXMAG0` alone do not determine the effective gain for an injected source's
Poisson variance. Do not label the current ICL injection as noise-aware; see
[`icl_coadd_poisson_noise_assessment_2026-11-04.md`](icl_coadd_poisson_noise_assessment_2026-11-04.md).
No pixels were read and no new HSC data were requested. Keep
`science_measurement=false`.

## ICL--satellite evidence synthesis (2026-11-04)

[`build_icl_association_evidence_report.py`](build_icl_association_evidence_report.py)
combines saved JSON diagnostics, hashes every input summary, and writes a
reviewable report without rereading images or querying HSC. Reproduce the
current version with:

```bash
python3 analysis/icl_wgl_membership/hsc_debug/build_icl_association_evidence_report.py \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_association_evidence_report_2026-11-04_02/report.md
```

The fixed 74-pair association sample has a positive exploratory `i2` partial
correlation (raw `r=0.234`, finite-lag pixel-noise diagnostic `r=0.249`), while
removing supported satellite-shell pixels from the ICL aperture lowers it to
`r=0.193`; the paired change interval is `[-0.064,-0.021]` across 28 approximate
spatial blocks. The `i` sample has 21 pairs and broad intervals. An outer-sky
residual covariate changes the `i2` estimate by `+0.013` with a paired interval
that includes zero. This is suggestive but sensitive to pixel overlap and
limited by PSF, sky-transfer, completeness, and covariance systematics. It is
not a confirmed science detection; all gates remain false. The generated
hash-provenance report is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_association_evidence_report_2026-11-04_02/report.md`.

The revised synthesis connects this 74-pair diagnostic subset to the separate
support-selection audit (94 joint-support rows: 26 in `i`, 68 in `i2`) and the
sample-size plan. Under an assumed partial `r=0.20`, 80% power, and conservative
168-hypothesis Bonferroni control, the plan needs 536–804 usable `i2` rows
(design effect 1.0–1.5), corresponding to about 6,306–9,459 parent clusters at
the observed support yield. The offline S23B catalog-only upper bound projects
about 181 supported `i2` rows if that yield transfers, and its image/control
coverage is unverified. These are design sensitivities, not detection claims;
the current 800-cluster sample is not adequately powered for that scenario.

## HSC coadd contributor gain metadata route (2026-11-04)

The official PDR3 schema browser exposes `pdr3_wide.mosaicframe`, whose rows
map a `(skymap_id, band)` coadd to contributing CCD `frame_num`, visit, and
CCD. The PDR3 [`mosaic` schema](https://hsc-release.mtk.nao.ac.jp/schema/schemata/pdr3.pdr3_wide.mosaic.json)
exposes `areacube`, used by v3 to match a requested sky coordinate. The related
`pdr3_wide.frame` schema includes `exptime` and `gain1` through
`gain4`. See the [mosaicframe schema](https://hsc-release.mtk.nao.ac.jp/schema/schemata/pdr3.pdr3_wide.mosaicframe.json),
[frame schema](https://hsc-release.mtk.nao.ac.jp/schema/schemata/pdr3.pdr3_wide.frame.json),
and [PDR3 database notes](https://hsc-release.mtk.nao.ac.jp/doc/index.php/database__pdr3/).
This is a route to investigate per-visit inputs missing from the coadd headers;
it does not by itself recover pixel-specific coadd weights or establish an
effective source-Poisson variance.

[`prepare_icl_frame_metadata_counts.py`](prepare_icl_frame_metadata_counts.py)
builds a local count-only SQL plan from existing target/control manifests. It
maps the 800 target and 798 control image centers to their nearest containing
HSC-I patch and counts its contributing frame/visit/CCD rows. The plan has 16
batches and returns one compact count row per center; it contains no image,
shape, or photo-z rows. The user approved plan v1 on 2026-11-05. CAS rejected it
with HTTP 406 before a job ID was returned because `public.skymap.patch_area`
was unavailable. The user then approved plan v2 using the unqualified
`patch_area, wcs` form shown in the PDR3 stored-functions guide. CAS again
rejected that SQL before job creation, this time reporting `column "patch_area"
does not exist`. The HTTP response body was inspected with credentials redacted;
both checkpoints confirm `query_submitted=false`. The published PDR3 `mosaic`
schema exposes `areacube`, so a v3 plan uses its documented spatial-search
operator with `radec_to_coord()` and removes the `public.skymap` dependency.
Plan v3 is offline validated and awaits approval because its SQL hash differs.

```bash
python3 analysis/icl_wgl_membership/hsc_debug/prepare_icl_frame_metadata_counts.py
```

The prepared plan is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_metadata_count_plan_2026-11-04_01/plan.json`
(SHA-256 `7d8b04df45489518b24095910bf201a93687f9e317aa2bf350ef0747c1c97570`).
The corrected, offline-validated v2 is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_metadata_count_plan_2026-11-05_02/plan.json`
(SHA-256 `b92a3c4931de25cf560abf5f6b0cd39af86b4e6d8cd3001a19cbeaa681f4316a`,
SQL aggregate SHA-256 `403869868c336b7f2c3ff8ce02430065b55345289396445e6feb4d4f8f0b4645`,
174,102 SQL bytes). It keeps the same 1,598 centers, 16 batches, 13 output
columns, and count-only scope. The approved attempt failed before job creation;
columns, and count-only scope. The exact error was `column "patch_area" does
not exist`; no job was created.
The revised, offline-validated v3 is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_metadata_count_plan_2026-11-05_03/plan.json`
(SHA-256 `3dac4db631d0dcbee5c5f4509e639f2bd90ca2af99fa7cf5f66c2a9f60bd9014`,
SQL aggregate SHA-256 `924c478ec02b5e703027f0446b834c381cf479f7e4ac2ab7dcf27b3af9d999d6`,
173,398 SQL bytes). It retains the 1,598 centers, 16 batches, 13 output
columns, and count-only scope; its new SQL awaits explicit approval.
After the count result is reviewed, any query returning actual gain, exposure,
and calibration values will need its own bounded plan and approval. Keep the
shot-noise and science gates false until weights, variance semantics, and the
injected-source transfer are validated.

## Supported 94-pair ICL–satellite association scope (2026-11-05)

[`analyze_icl_satellite_supported_94.py`](analyze_icl_satellite_supported_94.py)
recomputes the raw ICL–satellite partial correlation for the existing 94-pair
joint-support sample and compares it with the nested 74-pair diagnostic sample.
The 94-pair result is `r=0.079` for `i` (26 pairs) and `r=0.293` for `i2`
(68 pairs). The paired 94-minus-74 block-bootstrap intervals include zero in
both filters. This is a selected-support sensitivity, not a correction for
selection or a science detection. The analysis reads saved rows and coordinate
manifests only, makes no HSC request, and keeps `science_measurement=false`.
The hash-bound assessment is
[`icl_satellite_supported_94_assessment_2026-11-05.md`](icl_satellite_supported_94_assessment_2026-11-05.md).
The regenerated nine-input synthesis, including the paired 74-versus-94
sample comparison, is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_association_evidence_report_2026-11-05_02/report.md`.
The i2 point estimate changes from `r=0.234` to `r=0.293`, while the paired
95% block interval `[-0.188, 0.343]` includes zero. The result is a support
scope sensitivity, not evidence for a stronger association.

## Weak-lensing source-density boost gate (2026-11-05)

[`audit_wl_source_density_boost.py`](audit_wl_source_density_boost.py) compares
eligible source-pair densities in the 435-lens stack and 2,175 matched random
centers, using only the saved 5°/2.5° stack summaries. The count-based boost
candidate ranges from 1.765 at 0.1–0.2 Mpc to 1.148 at 2–3 Mpc. This may signal
non-negligible source dilution, but coverage/depth matching and spatial
uncertainty are not resolved. It is not applied to `DeltaSigma`, has no error
interval, and retains `science_measurement=false`. See the
[`source-density boost assessment`](wl_source_density_boost_assessment_2026-11-05.md).

The follow-up [`analyze_wl_source_density_boost_blocks.py`](analyze_wl_source_density_boost_blocks.py)
recomputes the boost by paired spatial block from the checksum-verified local
shape rows. In the 435-lens and 2,175-center sample, the 0.1–0.2 Mpc count ratio
is `1.765 ± 0.064`, while the corrected S19A response-adjusted estimator-weight
ratio using `SigmaCritInv` is `1.139 ± 0.048`. At 2–3 Mpc, the ratios are
`1.148 ± 0.030` and `1.126 ± 0.031`. The errors are delete-one block jackknife
standard errors, not confidence intervals. The script reconstructs the saved
`DeltaSigma` stack and records a 42-by-42 covariance across four boost metrics,
seven radial bins, and fourteen tangential/cross `DeltaSigma` parameters.
Boost–`DeltaSigma` correlation coefficients range from -0.346 to 0.325 at 2.5
degrees and from -0.393 to 0.465 at 5 degrees. The corrected ratios remain
diagnostics because full-aperture coverage matching and photo-z dilution
calibration are unresolved. No correction was applied and
`science_measurement=false` remains in force. The hash-bound results and limits
are in the [`joint boost/DeltaSigma assessment`](wl_source_density_boost_joint_covariance_assessment_2026-11-05.md).
The earlier dimensionless-`beta` ratios in the [preliminary assessment](wl_source_density_boost_spatial_assessment_2026-11-05.md)
are superseded.

The saved result is also shown in the Astrafek figure
[`plot_wl_boost_delta_sigma_joint.py`](plot_wl_boost_delta_sigma_joint.py) and
the [joint visualization note](wl_source_density_boost_joint_visualization_2026-11-05.md).
The figure places all four density ratios beside the random-subtracted
`DeltaSigma` and the complete 42-parameter jackknife correlation matrix. It
uses only the saved diagnostic JSON; no boost correction or fit is applied.

## 800-parent ICL–satellite support yield (2026-11-05)

[`audit_icl_satellite_support_yield.py`](audit_icl_satellite_support_yield.py)
joins the saved 200–300 kpc ICL support to satellite profiles at each radius.
The overlap decreases from 198 cluster/filter pairs at 5–10 kpc to 94 at
50–75 kpc; only 23 of the latter pass the existing center-eligibility flag.
Among 9,468 satellite candidates, 6,383 fail the target 60% annulus-coverage
floor and 1,306 have insufficient azimuthal controls. This audit does not
loosen either cut. Full results and provenance are in the
[support-yield assessment](icl_satellite_support_yield_assessment_2026-11-05.md).

The same 800 images were then remeasured at a 40% coverage floor as a
selection sensitivity. At 50–75 kpc the ICL/satellite overlap grows from 94
to 204 cluster/filter pairs, with center eligibility growing from 23 to 52.
The i2 partial correlation changes from `0.293` (68 pairs) to `0.249` (140
pairs); its relaxed-floor interval includes zero and the paired floor-change
interval includes zero. The larger sample is not adopted because the lower
floor admits partly masked apertures and additional controls. Details and
provenance are in the
[coverage-threshold assessment](icl_satellite_coverage_threshold_assessment_2026-11-05.md),
reproduced by
[`analyze_icl_satellite_coverage_threshold.py`](analyze_icl_satellite_coverage_threshold.py).
The support-versus-association summary figure is generated by
[`plot_icl_satellite_coverage_threshold.py`](plot_icl_satellite_coverage_threshold.py).
The center-eligible i2 subset at 40% has 30 pairs and `r=0.242` with a 95%
block interval `[-0.024, 0.684]`; at 60%, 15 pairs give `[-0.081, 0.785]`.
This central-mask-restricted comparison is also inconclusive.

## Outer-annulus coverage-loss components (2026-11-05)

[`audit_icl_satellite_outer_coverage_loss.py`](audit_icl_satellite_outer_coverage_loss.py)
reconstructs 50–75 kpc coverage for 2,691 candidates in the 224 frames with
ICL 200–300 kpc support. It separates image footprint, base quality mask, SEP
segmentation, neighboring-source masks, and fallback member circles. The
exact-pixel result has 1,171 candidates above the 40% target-coverage floor;
among 1,520 failures, 465 first cross the floor at the base mask, 81 at
segmentation, 973 at neighboring-source masks, and one at a fallback member
circle. No candidate fails first at the image boundary. These are sequential
mask-selection diagnostics, not flux corrections or a science detection.
Detailed hashes and interpretation are in the
[outer-annulus coverage-loss assessment](icl_satellite_outer_coverage_loss_assessment_2026-11-05.md).

The [mask-selection analysis](analyze_icl_satellite_mask_selection.py) then
compared coverage passage and retained-profile rates across candidate quartiles
in cluster-centric radius, SEP S/N, redshift, and richness, with cluster-ID
bootstrap intervals. The retained-profile fraction is 20.6–22.7 percentage
points higher in the largest-radius quartile in both filters. This indicates
non-random radial support and remains an uncorrected selection diagnostic;
details are in the
[mask-selection assessment](icl_satellite_mask_selection_assessment_2026-11-05.md).
An adjusted cluster-robust logistic diagnostic retains positive radius and
negative log-richness associations in both filters for coverage passage and
profile retention. The model converged but is not injection-calibrated and
must not yet be used as a science weight. It is reproducible with
[`analyze_icl_satellite_mask_selection_adjusted.py`](analyze_icl_satellite_mask_selection_adjusted.py);
the output is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_mask_selection_adjusted_2026-11-05_02/summary.json`.

The 94-pair supported ICL/satellite association was also split by the saved
central-mask eligibility flag. Only 8 `i` and 15 `i2` pairs pass that flag.
The `i2` partial-r point estimates are similar for eligible/ineligible pairs
(0.307/0.312), but both subgroup block intervals include zero; the `i` eligible
estimate is negative with a nearly full-range interval. This subgroup check is
too small to establish a mask effect or validate the current association. The
paired eligible-minus-ineligible block-bootstrap intervals include zero in
both filters (`i`: `[-1.665,0.632]`; `i2`: `[-0.516,0.693]`); only 86% of
paired `i` resamples meet the minimum subgroup-size rule. It does not relax
masks or correct fluxes. See
[`analyze_icl_satellite_center_mask_sensitivity.py`](analyze_icl_satellite_center_mask_sensitivity.py)
and its assessment
[`icl_satellite_center_mask_sensitivity_assessment_2026-11-05.md`](icl_satellite_center_mask_sensitivity_assessment_2026-11-05.md).
Machine result:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_center_mask_sensitivity_2026-11-05_02/summary.json`.

## Paired 40% and 60% satellite-support injection comparison (2026-11-05)

[`compare_icl_satellite_support_thresholds.py`](compare_icl_satellite_support_thresholds.py)
compares the completed 40% run and the earlier 60% run on identical image and
position-dependent PSF contents. It verified 57 paired frames, matching
re-detection counts/masks, and matching injection profiles. Median satellite
support increased from 2 to 5 per frame in `i` and 2 to 6 in `i2`; median
control profiles increased from 7 to 26 and 9 to 34. The paired satellite
unit-injection response change has cluster-bootstrap intervals containing zero
in both filters. This is a support/estimator sensitivity on a PSF-supported
subset, not a completeness correction or detection. See the
[detailed assessment](icl_satellite_support_threshold_comparison_2026-11-05.md)
and paired outputs under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_support_threshold_comparison_2026-11-05_02/`.

## Joint-support selection weighting sensitivity (2026-11-05)

[`analyze_icl_satellite_joint_support_ipw_sensitivity.py`](analyze_icl_satellite_joint_support_ipw_sensitivity.py)
models whether each of the 224 ICL-supported cluster/filter pairs has a saved
50–75 kpc satellite profile. The 94 supported profiles differ from unsupported
pairs in measured ICL coverage, richness, and redshift; the largest standardized
mean difference is 0.70 for i2 ICL coverage. Stabilized inverse-probability
weighting changes the z/richness-adjusted partial correlations from `0.079` to
`−0.011` in `i`, and from `0.293` to `0.227` in `i2`. Connected 5-degree block
bootstrap intervals are `[-0.566, 0.581]` and `[-0.034, 0.416]`; both include
zero. Effective sample sizes remain
22/26 and 60/68. The original 216 BFGS precision-loss messages were checked:
their score gradients met the stationarity threshold (maximum `1.40×10⁻⁶`), so
none were failed fits. Omitting central-mask state, adding quadratic coverage,
or allowing filter-specific logistic slopes gives point estimates in
`[-0.014, 0.021]` (`i`) and `[0.226, 0.231]` (`i2`). The method assumes profile
availability is conditionally independent of the unmeasured satellite signal,
which is unverified. Treat IPW as a sensitivity only, not a completeness
correction or a detection. The
[assessment](icl_satellite_joint_support_ipw_sensitivity_2026-11-05.md)
records the assumptions and optimizer audit. Output summary:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_joint_support_ipw_2026-11-05_06/summary.json`
(SHA-256 `8dbfae2169361c9fb6ce5ccfcce684b486162dbf4820e531112cdfc4f39fdbc7`;
`science_measurement=false`).

## HSC frame-count query submission update (2026-11-05)

The user approved the exact count-only plan v2 (SHA-256
`b92a3c4931de25cf560abf5f6b0cd39af86b4e6d8cd3001a19cbeaa681f4316a`). Its
first batch was attempted, but CAS returned HTTP 406 because the SQL's
unqualified `patch_area` column is not present in the referenced schema. No
job was created and `query_submitted=false`; no result rows or credentials were
saved. Do not retry this unchanged plan. Resolve the official schema/table
columns first; any changed SQL plan needs approval for its new hash. The query
scope remains count-only, with no images, shapes, or photo-z rows.

## PSF coverage plan for the full ICL-supported pair sample (2026-11-05)

[`prepare_icl_psf_picker_supported_pairs_plan.py`](prepare_icl_psf_picker_supported_pairs_plan.py)
reconstructs the member positions from the existing 800-image sample and
40%-coverage profile product. It finds 2,693 unique satellite injection
positions across the 224 cluster/filter pairs with ICL 200–300 kpc support.
Existing, audited PSFs cover 688 satellite positions in the 57-frame injection
subset and all 224 cluster centers. The resulting offline plan therefore
contains 2,005 additional satellite PSF coordinates, estimated at 40,420,800
bytes (38.55 MiB) from the actual existing PSF size distribution. It requests
only centered HSC-I `pdr3_wide` coadd model PSF cores: no images, catalog shape
rows, or photo-z data. The plan and coordinate files are under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_psf_picker_satellite_support224_plan_2026-11-05_01/`.
The user approved the exact plan hash
`1d7c0189721600fc31e8524a4258dc70a5c0f0d436df45b57b9a0c3dd0058b65`; all
2,005 requested PSFs were received and audited on 2026-11-06. The plan did not
change the cluster image sample. The resulting injection analysis remains a
selection/systematics test, not a science detection.

## Injection-amplitude transfer check (2026-11-05)

The PSF-supported 57-frame injection/re-detection was repeated with satellite
amplitudes 0.5 and 2 nJy arcsec⁻² around the existing 1 nJy arcsec⁻² run. The
median `i2` unit-response stayed near 0.56–0.59 and paired cluster-bootstrap
intervals for both amplitude changes include zero. Detection masks still grow
with injected brightness, so stable ensemble transfer does not imply
amplitude-independent selection. Reproduce and audit with
[`analyze_icl_satellite_injection_amplitude_grid.py`](analyze_icl_satellite_injection_amplitude_grid.py);
see the [assessment](icl_satellite_injection_amplitude_grid_assessment_2026-11-05.md)
and output under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_injection_amplitude_grid_2026-11-05_02/`.
The `i` subsample has only 13 frames; all outputs remain `science_measurement=false`.

## Global weak-lensing profile null audit (2026-11-05)

[`analyze_wl_random_shape_stack_detection.py`](analyze_wl_random_shape_stack_detection.py)
uses the complete 14-component paired lens/random endpoint-jackknife covariance
to test the 7-bin tangential ΔΣ profile against zero, including cuts that
drop the inner 0.1–0.2 and 0.4–0.7 Mpc bins. The current 435-lens common-support
profile is strongly inconsistent with zero under both saved covariance
partitions, while the 7-bin cross component is consistent with zero. The point
estimates are identical for 5° and 2.5° blocks, but nominal significance and
uncertainties are conditional on the partition. This is not yet a publication
detection: contiguous survey regions, random-point mask/depth matching, n(z)
calibration, reduced shear, and coadd covariance remain unverified. See the
[assessment](wl_random_shape_stack_detection_assessment_2026-11-05.md) and
machine output at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_shape_stack_detection_audit_2026-11-05_01/`.
The output keeps `science_measurement=false`.

## ICL–weak-lensing sample-overlap audit (2026-11-05)

[`audit_icl_wl_sample_overlap.py`](audit_icl_wl_sample_overlap.py) streamed and
hash-verified the existing 800-lens target shape export and all 124 random
count batches, then reproduced the 435-lens common random-subtracted WL set and
its stored lens-ID hash. Of the 224 ICL profile-supported filter pairs, the
common WL set contains 38 `i` and 70 `i2` pairs. For the narrower 94-pair
50–75 kpc satellite-profile subset, only 12 `i` and 27 `i2` pairs have both
target and random shape support. Treat the joint subset as exploratory; it is
too small for confirmatory filter-specific mass–satellite inference. An
ICL-supported WL stack on the larger 38/70 sets is a separate exploratory
measurement and cannot establish the satellite relation. No HSC query was
submitted. Details and limitations are in
[`icl_wl_sample_overlap_assessment_2026-11-05.md`](icl_wl_sample_overlap_assessment_2026-11-05.md);
machine output is under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_sample_overlap_audit_2026-11-05_03/`
(summary SHA-256 `98d94834b895f1275337f494a801aab7e021048ef502d4f7074635bb26b25e12`).
All outputs keep `science_measurement=false`.

## ICL-supported weak-lensing subset stacks (2026-11-05)

[`run_icl_supported_wl_subset_stack.py`](run_icl_supported_wl_subset_stack.py)
reuses the verified S19A tomographic estimator, matched random centers, and
paired endpoint jackknife to stack the `i` 38 / `i2` 70 ICL-profile-supported
common lenses separately. It also supports the narrower joint ICL plus
50–75 kpc satellite-profile subset (`i` 12, `i2` 27). The 5-degree and
2.5-degree ICL-only runs and their full-covariance null checks are under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_supported_wl_subset_stack_2026-11-05_01/`
and `_02p5/`. Tangential point estimates are unchanged; median diagonal error
ratios (2.5°/5°) are 1.05 (`i`) and 1.12 (`i2`). The `i` cross component
rejects zero at both widths (`p=8.0×10⁻⁴` and `3.4×10⁻¹¹`), so its tangential
profile is not interpretable as a detection. The `i2` cross component is
consistent with zero (`p=0.42`, `0.31`) and its tangential stack is promising
but still exploratory. The [block-width comparison](compare_icl_supported_wl_blockwidth.py)
keeps both covariance partitions and their assumptions. A separate joint
satellite/WL stack uses `i` 12 and `i2` 27. Its `i` cross null is also rejected
at both widths (`p=3.5×10⁻²²`, `2.5×10⁻²⁷⁴`), while `i2` remains consistent
with zero (`p=0.575`, `0.282`). The joint subset is diagnostic only. The full
interpretation, profile figure and covariance outputs are in
[`icl_supported_wl_subset_stack_assessment_2026-11-05.md`](icl_supported_wl_subset_stack_assessment_2026-11-05.md).
No subset product changes the science gate.

## HSC PDR3 frame-count query v3 (2026-11-05)

The exact approved v3 plan (`3dac4db631d0dcbee5c5f4509e639f2bd90ca2af99fa7cf5f66c2a9f60bd9014`)
completed 16 count-only batches for 800 target and 798 matched-control centers.
All 1,598 map to a patch. Target/control medians are 42/42 frames, 9/9 visits,
and 38/37 CCDs; the median paired difference is zero. These counts show no
gross central coverage offset, but do not measure pixel depth or shot noise.
The [Astrafek auditor](audit_icl_frame_metadata_count_results.py) verifies all
batch hashes and paired mappings. See the
[assessment](icl_frame_metadata_count_assessment_2026-11-05.md) and output at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_metadata_count_audit_2026-11-05_02/`
(summary SHA-256 `63e7424cb6abe3f89b7138bb99ea1072fbb367ab1e18b014974981334522810b`).
No images, shape rows, photo-z, gains, exposures, sky values, or credentials
were returned; `science_measurement=false`.

The combined position-PSF bank path was exercised after approval of the exact
2,005-position plan (`1d7c0189721600fc31e8524a4258dc70a5c0f0d436df45b57b9a0c3dd0058b65`).
`verify_supported_pairs_position_psf_bank.py` verified the original 744
injection PSFs, 1,542 target-center PSFs, and 2,005 additional satellite PSFs:
224 ICL-supported frames, 2,693 satellite roles, and 2,917 unique PSF files.
The expanded redetection run then completed 137 pairs; see
[`icl_supported_pairs_psf_redetection_assessment_2026-11-06.md`](icl_supported_pairs_psf_redetection_assessment_2026-11-06.md).

## PSF-injection sample representativeness (2026-11-05)

[`audit_joint_injection_psf_sample_representativeness.py`](audit_joint_injection_psf_sample_representativeness.py)
compares the existing PSF-convolved 57-pair injection subset with the 224-pair
ICL-supported parent using only saved summary/profile tables. The 57-pair
subset has higher 200–300 kpc ICL coverage (overall standardized mean
difference +0.80; `i` +1.06, `i2` +0.74) and substantially higher 50–75 kpc
satellite support at the 60% floor (54/57 versus 98/224). Its injection
transfer should not be generalized to the parent sample. The offline report
and checksummed table are in
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_joint_injection_psf_sample_representativeness_2026-11-05_01/`;
all outputs retain `science_measurement=false`.

The batch runner also exposes `--icl-supported-only` so the injection
selection can include every ICL-supported frame and preserve satellite-profile
support failures as measured outcomes. For the 224-pair run, use the
`SupportedPairsPositionPSFBank` with the original injection batch, audited
target-center batch, and approved satellite-position batch. The combined bank
passed verification on 2026-11-06: all 224 frames, 2,693 satellite roles, and
2,917 unique PSF files are covered. The official PSF Picker batch returned
2,005/2,005 files and passed the exact approved-plan and coordinate audit.

The post-download bank audit and batch runner also require the exact approved
satellite-plan SHA-256 through `--approved-satellite-plan-sha256` and
`--approved-supported-satellite-psf-plan-sha256`, respectively. Both reject a
plan that differs from the explicitly approved 2,005-position request. The
result records request-approval verification separately from `psf_approved`;
the latter remains false until scientific PSF validation passes.

The expanded redetection injection completed 137 pairs (44 `i`, 93 `i2`) after
the measured 40%-coverage selection from the 224-frame bank; 87 bank frames
did not pass that profile-support threshold. All 137 frame records and the
position-specific PSF hashes passed
[`audit_icl_supported_pairs_psf_redetection_batch.py`](audit_icl_supported_pairs_psf_redetection_batch.py).
At least one satellite support gap remains in 79/137 frames. Median smooth-ICL
recovery is about 0.226 in both filters, while satellite-profile recovery is
wide and support-limited (14 valid `i`, 44 valid `i2`). This is not a science
measurement or a calibration matrix: injected Poisson noise, upstream sky
transfer, and extended PSF wings are absent, and `mask_all` can remove BCG
light. See the
[assessment](icl_supported_pairs_psf_redetection_assessment_2026-11-06.md)
and its checksummed summaries; keep `science_measurement=false`.

## Random-center shape-row export integrity (2026-11-05)

[`audit_wl_random_shape_export_integrity.py`](audit_wl_random_shape_export_integrity.py)
reconciles the approved 124-batch source-shape export against its exact plan,
checkpoint, row summary, and the 5°/2.5° stack provenance. All 124 CSV hashes
and byte counts match, totaling 9,010,324 rows and 2,962,522,918 bytes. Both
stacks reference the same plan hash and all current per-batch hashes. A stale
`KeyboardInterrupt` string remains in the top-level checkpoint although its
status, all jobs, row summary, and stack references are complete; the audit
records but does not modify that historical field. This is a file-integrity
result only: random-aperture depth/mask/PSF/sky matching and covariance
footprint assumptions remain open, and `science_measurement=false`.
Assessment:
[`wl_random_shape_export_integrity_assessment_2026-11-05.md`](wl_random_shape_export_integrity_assessment_2026-11-05.md).
Machine summary:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_shape_export_integrity_2026-11-05_01/summary.json`.

## S23B support-yield and spatial-block ceiling (2026-11-05)

[`assess_s23b_support_yield_and_blocks.py`](assess_s23b_support_yield_and_blocks.py)
combines the existing S23B catalog-only candidate ceiling, spatial-cell ledger,
observed 60% and 40% outer-satellite support, and Fisher-z block-power plan.
Transferring the observed i2 support rates to all 1,333 candidates projects
about 181 pairs at 60% coverage or 373 at 40%. The 40% point estimate is
`r=0.249` with a 95% block interval that includes zero. The S23B catalog spans
100 target 5° cells, but only 22 lie outside the existing discovery footprint;
under the current block rule these cells cannot meet even the 92-block
single-endpoint target, much less 180 blocks for 30 comparisons. The counts are
planning upper bounds: candidate images, masks, depth, PSF, and matched controls
remain unverified. No query or download occurred; `science_measurement=false`.
See
[`icl_s23b_support_yield_and_blocks_assessment_2026-11-05.md`](icl_s23b_support_yield_and_blocks_assessment_2026-11-05.md)
and machine output
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_support_yield_and_blocks_2026-11-05_01/`.

The follow-up catalog-only redshift/BCG-spec cut comparison retained the
exact 1,333-ID baseline and evaluated five alternatives. Widening the redshift
range while retaining BCG spectroscopic consistency gives 1,983 or 2,339
greedy additions but only 25 or 26 new 5-degree cells. Removing the BCG-spec
condition gives up to 4,137 additions and 31 new cells, with additional lens
redshift calibration required. These cell counts are not independent blocks
or proof of PDR3 image/control support. Keep the current selection as the
replication baseline; evaluate any wider sample as a prespecified,
redshift-stratified extension after coverage and calibration audits. No HSC
query or download occurred. Runner and full scenario table:
[`audit_s23b_camira_selection_scenarios.py`](audit_s23b_camira_selection_scenarios.py),
[`icl_s23b_camira_selection_scenarios_2026-11-05.md`](icl_s23b_camira_selection_scenarios_2026-11-05.md).
Machine summary:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_camira_selection_scenarios_2026-11-05_02/summary.json`.

The random-subtracted weak-lensing stack was decomposed by HSC Y3 source
tomographic bin. Source bin 1 has an exact estimator-weight excess of 4.61
in 0.1–0.2 Mpc, declining to 1.26 at 2–3 Mpc; bins 2–4 remain near 0.97–1.16.
The difference from each higher bin survives Holm correction over 21
radial-bin contrasts for both 5° and 2.5° endpoint jackknives. Excluding bin 1
retains a positive tangential profile with a cross profile consistent with
zero, while shifting the inner tangential profile materially. This points to
a low-redshift source-selection/dilution issue, but does not identify member
contamination or justify a boost correction. The four-bin counts and profile
reproduce the saved stack; `science_measurement=false` remains. Reproduction
and limitations are in
[`analyze_wl_tomographic_weight_boost.py`](analyze_wl_tomographic_weight_boost.py)
and
[`wl_tomographic_weight_boost_assessment_2026-11-05.md`](wl_tomographic_weight_boost_assessment_2026-11-05.md).
Machine output:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_tomographic_weight_boost_2026-11-05_04/summary.json`.

## HSC frame-level noise metadata plan (2026-11-05)

The approved v3 frame-count query completed for all 1,598 target/control
centers. Its 16 batches show a median of 42 contributing frames per center in
both samples and no median paired offset, but counts cannot provide a
pixel-level Poisson variance. The separate [count audit](icl_frame_metadata_count_assessment_2026-11-05.md)
preserves the measured distributions and its `science_measurement=false`
limit.

[`prepare_icl_frame_metadata_rows.py`](prepare_icl_frame_metadata_rows.py)
uses those audited center-to-patch mappings and exact per-center frame counts
to prepare, offline, a possible row export from `pdr3_wide.mosaicframe` joined
to `pdr3_wide.frame`. The plan asks for contributing frame identity, original
filter, exposure time, four amplifier gains, `FLUXMAG0`, sky RMS, and MJD. It
contains 21 columns, 16 batches, and an exact estimate of 67,071 center-frame
rows (roughly 13–27 MB as CSV). It returns no image, shape, photo-z, or object
catalog rows. These are calibration metadata, so this distinct row export
requires its own explicit approval; the earlier count-only approval does not
authorize submission.

The [official mosaicframe schema](https://hsc-release.mtk.nao.ac.jp/schema/schemata/pdr3.pdr3_wide.mosaicframe.json)
defines `frame_num`, visit, and CCD membership for each coadd patch, while the
[official frame schema](https://hsc-release.mtk.nao.ac.jp/schema/schemata/pdr3.pdr3_wide.frame.json)
defines the requested exposure and calibration fields. Even if approved and
returned, these rows will not provide per-pixel coadd weights or the mapping
from image pixels to amplifier sections; they cannot alone establish an ICL
shot-noise model. Do not update injected-source variance or the science gate
from these metadata alone.

Prepared plan:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_metadata_row_plan_2026-11-05_01/plan.json`
(SHA-256 `49aa4b7601dd42243209d7b3b394a4dd3b4cbf00d87d28cbb9467553c29f4ce3`,
SQL aggregate SHA-256 `4270f3f0d290ac427b24694f7561f46cdee71cc8c3621a57c3fcf3f034583c76`).
The user explicitly approved this exact plan, and
[`submit_icl_frame_metadata_rows.py`](submit_icl_frame_metadata_rows.py)
submitted its 16 batches on `lageunha`. The first response exposed a validator
assumption: CAS returns original detector filters `i` and `i2`, not release
labels `HSC-I` and `HSC-I2`. The runner was corrected to accept only those two
original filter codes and resumed using the checkpointed first job; no extra
query scope was added.

All 67,071 rows for 1,598 requested centers were returned and verified, with
zero null values in exposure, amplifier gains, FLUXMAG0, sky RMS, or MJD. The
independent offline auditor
[`audit_icl_frame_metadata_row_results.py`](audit_icl_frame_metadata_row_results.py)
rechecked all 16 batch hashes, exact center/patch mappings, frame uniqueness,
and per-center membership counts. Among 798 matched target/control pairs,
median differences were zero for frame count, visit count, CCD count, total
exposure, median sky RMS, zeropoint, and all four median gains. This rules out
only a large systematic imbalance in these patch-level summaries; it does not
measure per-pixel depth. HSC returned original-frame filter rows: 20,341 `i`
and 46,730 `i2`. Images, shape, photo-z, object rows, and credentials were not
returned or saved. The export cannot supply per-pixel coadd weights or map
pixels to amplifier sections, so it does not establish a source-Poisson model;
keep `science_measurement=false`.

Result audit:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_metadata_row_audit_2026-11-06_01/summary.json`
(SHA-256 `0f361fcb496fbd7be7c1059f638414f6900535f73a9ac48696030480f400a237`). Raw batch checkpoint and rows:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_metadata_row_results_2026-11-06_01/`.

The audit shows matched target/control medians of zero for the number of
contributing frames, visits and CCDs, summed exposure time, median sky RMS,
median FLUXMAG0 and each amplifier gain. The central 68% ranges are small for
counts/exposure (frame rows −4 to +4; visits −1 to +1; CCDs −3 to +4;
exposure −800 to +800 s) but are not pixel-depth comparisons. The result
therefore rules out a gross imbalance in these patch-level summaries only.

## ICL frame-count offset sensitivity (2026-11-05)

[`analyze_icl_frame_count_icl_sensitivity.py`](analyze_icl_frame_count_icl_sensitivity.py)
joined the audited count-only HSC result to the saved 200–300 kpc ICL
target–control contrasts for 224 supported cluster/filter pairs. It tested
target-minus-control differences in frame, visit, and CCD counts against the
ICL contrast, adjusting for redshift, log richness, and common aperture
coverage; 5,000 connected 5-degree block resamples were run within each filter.
All six partial correlations are near zero and all intervals include zero.
This does not establish equal pixel depth, sky noise, or coadd weighting, and
does not alter the ICL estimate or science gate. The runner used only saved
tables and no image pixels or HSC requests. See the
[assessment](icl_frame_count_icl_sensitivity_assessment_2026-11-05.md) and
machine summary
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_count_sensitivity_2026-11-05_01/summary.json`.

## Random-subtracted weak-lensing NFW mass screen (2026-11-05)

`fit_wl_random_stack_nfw_mass.py` fits the existing 435-lens random-subtracted
profile with a single-halo NFW forward model using source bins 2–4 and proper
0.2–1.5 Mpc. Fixed-concentration fits (`c=3,4,5`) are repeated with the saved
5° and 2.5° endpoint-jackknife covariance. At `c=3`, the effective
`M200c` is `1.78e14` and `1.80e14 Msun`; at `c=5` it is `1.54e14` and
`1.44e14 Msun`. The roughly 20% concentration-grid shift and poor `c=5`
profile fit show material model dependence. These conditional intervals omit
concentration, reduced-shear, magnification, miscentering, two-halo, photo-z,
random-aperture matching, and covariance uncertainties. This is an exploratory
physical-scale screen, not a science-grade mass calibration; keep
`science_measurement=false`. Full method, caveats, and hashes are in the
[NFW screen assessment](wl_random_stack_nfw_mass_screen_assessment_2026-11-05.md).
No new HSC query or image read was used.

## Reduced-shear sensitivity (2026-11-05)

[`fit_wl_random_stack_reduced_shear_sensitivity.py`](fit_wl_random_stack_reduced_shear_sensitivity.py)
refits the same 435-lens random-subtracted profile using calibrated HSC Y3
`n(z)` moments through ninth power and the reduced-shear geometric series
through order 8. Its conservative omitted-tail bound is 0.012–0.022%. The
effective mass shifts down 1.7–5.5% relative to the weak-shear NFW screen;
fixed-concentration model dependence remains about 20%, and the `c=5` fits
remain poor. This addresses only one modeling term and does not calibrate
mass or the ICL–satellite relation. See the
[assessment](wl_random_stack_reduced_shear_sensitivity_assessment_2026-11-05.md)
and summary
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_stack_reduced_shear_sensitivity_2026-11-05_03/summary.json`.
No HSC query or image read was made; `science_measurement=false`.

## Random-subtracted richness contrast (2026-11-05)

[`analyze_wl_random_subtracted_richness_split.py`](analyze_wl_random_subtracted_richness_split.py)
retains the predeclared 400/400 parent richness split, filters to lenses with
estimator-usable target and matched-random support, redshift-matches 204 pairs,
and recomputes four stacks under a shared 50-block jackknife. The high-minus-low
random-subtracted `DeltaSigma` contrast is positive in every bin from 0.2 to
3 Mpc. The nominal seven-bin Wald diagnostic is not a detection significance;
selection, matched-aperture depth/mask/PSF, and contiguous survey covariance
remain unresolved. The earlier 214-pair run was superseded after discovering
10 lenses per group lacked target tomography support. See the
[assessment](wl_random_subtracted_richness_split_assessment_2026-11-05.md) and
machine summary
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_subtracted_richness_split_2026-11-05_02/summary.json`.
Keep `science_measurement=false`; no HSC request or image read was made.

## Two-mass NFW contrast screen (2026-11-05)

[`fit_wl_richness_split_nfw_mass_contrast.py`](fit_wl_richness_split_nfw_mass_contrast.py)
fits lower- and higher-richness effective `M200c` jointly to their
random-subtracted profile contrast with the saved shared covariance. For
fixed `c=3,4,5`, the higher/lower mass ratio is 1.21–1.89 using source bins
1–4 and 1.39–2.85 after excluding source bin 1. Profile intervals are broad
(upper limits up to 28); concentrations and source-bin choice move the masses
substantially. This is an exploratory hierarchy screen, not a richness-mass
calibration. See the
[assessment](wl_richness_split_nfw_mass_contrast_assessment_2026-11-05.md)
and two machine summaries listed there. No new HSC request or image read;
`science_measurement=false`.


## Integrated ICL–satellite evidence update (2026-11-05)

[`build_icl_association_evidence_report.py`](build_icl_association_evidence_report.py)
now consolidates the 94-pair support-selection weighting, PSF-injection
sample representativeness, and frame-count sensitivity alongside the earlier
overlap, sky, noise, and sample-size audits. The selection-weighted partial
correlation is `−0.011` in `i` and `0.227` in `i2`, versus raw values `0.079`
and `0.293`; the i2 57-pair PSF injection subset has 0.736 SD higher ICL
coverage than the parent support set. These diagnostics further weaken
population and transfer claims. Full values, caveats, input hashes, and
reproduction provenance are in
[`icl_association_evidence_status_2026-11-05_03.md`](icl_association_evidence_status_2026-11-05_03.md).
The report reads saved summaries only and keeps `science_measurement=false`.

## Reduced-shear richness contrast (2026-11-05)

[`fit_wl_richness_split_reduced_shear.py`](fit_wl_richness_split_reduced_shear.py)
refits the same 204 redshift-matched lenses per richness group and the same
random-subtracted profile/covariance using calibrated source `n(z)` moments
through power 17. The high/low effective mass ratios are 1.168–1.214 across
fixed `c=3,4,5` with source bins 1–4, and 1.150–1.198 after excluding source
bin 1. The weak-shear `c=3` ratio was 1.89 (or 2.85 without source bin 1), so
the forward-model choice materially affects the concentration sensitivity.
The ratios are diagnostic single-mass fits; the approximate covariance and
mass-distribution, photo-z, miscentering, magnification, and matched-random
systematics remain. Source-supported order-16 truncation bounds peak at 6.94%
per pair for `c=3` without bin 1; fitted ratios change by at most 0.3% between
series orders 8 and 16. See the
[assessment](wl_richness_split_reduced_shear_assessment_2026-11-05.md) and
its exact input/output hashes. No HSC query or image read; keep
`science_measurement=false`.

## Random-subtracted weak lensing by ICL strength (2026-11-05)

[`analyze_wl_random_subtracted_icl_split.py`](analyze_wl_random_subtracted_icl_split.py)
adds an exploratory ICL-to-weak-lensing analysis using the saved ICL annulus
support and exact HSC random-subtracted source rows. In the `i2` common sample,
70 lenses were split at the median 200–300 kpc target-minus-control ICL
contrast and 29 high/low pairs were matched in redshift and log richness.
Only source bins 2–4 were stacked. The random-subtracted high-minus-low
`DeltaSigma` profile is radially mixed; its conditional seven-bin tangential
Wald test changes from nominal `p=2.42e-4` with 26 approximate 5° blocks to
`p=0.099` with 54 approximate 2.5° blocks; cross-null values are `p=0.716` and
`p=0.633`. The point profile is identical, but its error bars change by
0.89–1.26. The nominal contrast does not pass this block-width sensitivity;
neither p-value is a calibrated detection significance. Image-side ICL
Poisson-noise, extended-PSF, sky-transfer and completeness corrections remain
open, as do validated contiguous-region covariance and lensing calibration.
The `i` common sample is only 38 lenses and is not included in this first run.
See [`wl_icl_split_assessment_2026-11-05.md`](wl_icl_split_assessment_2026-11-05.md)
for the full profile table, matched-sample details and provenance. The plotter
is [`plot_wl_random_subtracted_icl_split.py`](plot_wl_random_subtracted_icl_split.py).
The 5° and 2.5° summaries and figures are kept separately for reproducibility.
[`compare_wl_icl_split_blockwidth.py`](compare_wl_icl_split_blockwidth.py)
verifies that the matched pairs and profile points are identical and writes
the error-bar and nominal-null comparison.
This analysis reads saved files only, requests no HSC data, and keeps
`science_measurement=false`.

## ICL-supported satellite-amplitude response on 137 pairs (2026-11-06)

The earlier 0.5/1/2 nJy satellite-injection grid was limited to 57 frames.
[`analyze_icl_satellite_injection_amplitude_grid.py`](analyze_icl_satellite_injection_amplitude_grid.py)
now accepts a declared paired-frame count and excludes unavailable profile
responses from summaries and cluster bootstrap intervals. The grid was
repeated on the full 137-pair ICL-coverage-selected sample with the same
images, local PSFs, and baseline detections across amplitudes. All three
137-frame runs completed without frame errors and passed
[`audit_icl_supported_pairs_psf_redetection_batch.py`](audit_icl_supported_pairs_psf_redetection_batch.py).

Median satellite-profile response in `i2` is 0.612/0.409/0.411 at amplitudes
0.5/1/2, but paired cluster-bootstrap intervals for both amplitude changes
include zero; the intervals are broad and responses are available in only
43–44 frames. In `i`, only 12–15 frames contribute. Meanwhile median added
mask pixels grow from 151 to 871 in `i2` satellite-only injections and from
349 to 1,423 in `i`. At least one scenario-level satellite-support gap occurs
in 78/137, 79/137, and 82/137 frames at amplitudes 0.5, 1, and 2.

The results show amplitude-dependent mask growth but do not provide a stable
completeness correction or a science detection. Injected Poisson noise,
upstream sky transfer, extended-wing PSF calibration, and validated covariance
remain open. See the detailed
[assessment](icl_supported_pairs_satellite_amplitude_assessment_2026-11-06.md)
and paired checksummed output under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_injection_amplitude_grid_supported137_2026-11-06_01/`.
Keep `science_measurement=false`.

## Diffuse ICL injection amplitude response (2026-11-06)

[`analyze_icl_diffuse_injection_amplitude_grid.py`](analyze_icl_diffuse_injection_amplitude_grid.py)
compares the same 137 PSF-convolved, ICL-coverage-selected frames at diffuse
input amplitudes 0.5, 1.0 and 2.0 nJy arcsec⁻² at 100 kpc. The satellite
injection stays at 1.0 and SEP detections/masks are regenerated. It verifies
identical frame IDs, images, position-dependent PSFs and no-injection baseline
detections before paired 5° connected-sky-block bootstrap comparisons.

The median response per unit diffuse input is `i=0.240, 0.226, 0.224` and
`i2=0.225, 0.226, 0.226` at amplitudes 0.5, 1.0, 2.0. Paired response changes
relative to 1.0 are consistent with zero in both filters. The absolute response
is only about 0.23 of the injected light, and the response stability does not
correct for missing injected-source Poisson noise, coadd weights, upstream HSC
sky transfer, or extended PSF wings. Source masks and detections change in many
frames across all amplitudes, so this is a transfer-linearity diagnostic, not
a flux correction or astrophysical measurement; `science_measurement=false`.

Assessment and paired table are in
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_diffuse_injection_amplitude_grid_supported137_2026-11-06_01/`
(summary SHA-256 `538656e4b3c5141282ddd24c16ad8ce0d20454692228f38b437ad996b29fcd96`).
Inputs were the 0.5, 1.0, and 2.0 runs in their corresponding
`joint_icl_satellite_psf_redetection_icl_supported224*` directories; all three
completed 137/137 frames with zero errors.

## ICL–satellite support-sample block-width sensitivity (2026-11-05)

[`analyze_icl_satellite_supported_94_blockwidth.py`](analyze_icl_satellite_supported_94_blockwidth.py)
recomputes the raw ICL/satellite partial correlation for the prior 74-pair
diagnostic subset and the full 94-pair joint-support subset at 2.5°, 5°, and
10° connected sky-cell widths. The `i2` 94-pair estimate stays at `r=0.293`,
with bootstrap intervals `[0.002,0.504]`, `[0.048,0.498]`, and `[0.056,0.496]`.
The 2.5° interval is borderline, while the 74-pair estimate crosses zero at
2.5°. Paired 94-minus-74 intervals cross zero at all three widths in both
filters. These are selected-sample diagnostics under approximate cells, not a
population detection. See
[`icl_satellite_supported_94_blockwidth_assessment_2026-11-05.md`](icl_satellite_supported_94_blockwidth_assessment_2026-11-05.md)
for full values, limitations, and hashes. The plotter is
[`plot_icl_satellite_supported_94_blockwidth.py`](plot_icl_satellite_supported_94_blockwidth.py).
No HSC query or image read; keep `science_measurement=false`.

## ICL–satellite pixel-noise correction on 94 supported pairs (2026-11-06)

[`assess_icl_satellite_pixel_noise_bias_supported_94.py`](assess_icl_satellite_pixel_noise_bias_supported_94.py)
extends the saved finite-lag pixel-noise correction to all 94 existing
supported pairs (26 `i`, 68 `i2`). The independent support audit confirmed
exact 94-pair covariance coverage, and all 74 nested-pair raw measurements
match the earlier association table. Coordinate-connected block resampling
gives corrected `i`/`i2` partial correlations `0.081`/`0.306` at 2.5°, 5°,
and 10°; 95% intervals are `[-0.359,0.610]`, `[-0.195,0.511]`,
`[-0.501,0.524]` for `i`, and `[-0.004,0.523]`, `[0.050,0.519]`,
`[0.062,0.512]` for `i2`. The finite-lag correction is small, while the
smallest-width `i2` interval is borderline with zero. These selected-sample
intervals are not calibrated detections; support selection, sky,
source-Poisson, PSF-wing, calibration, and astrophysical covariance remain
unresolved. See
[`icl_satellite_pixel_noise_bias_supported_94_assessment_2026-11-06.md`](icl_satellite_pixel_noise_bias_supported_94_assessment_2026-11-06.md).
No HSC query, download, or image read was done; `science_measurement=false`.

## HSC warp-product scope assessment (2026-11-06)

The authenticated PDR3 DAS Cutout manual documents `type=warp`, which returns
visit-level warped CCD cutouts in a tar archive with selectable image, mask,
and variance layers. Deduplicating the approved 67,071-row frame metadata gives
12,839 patch/filter/visit tuples and 56,731 patch/filter/visit/CCD tuples.
These are not transfer counts; one spatial request can return multiple warps,
and overlapping target/control requests can repeat data. The full transfer is
not planned until a pilot establishes actual output sizes and archive layout.

A proposed one-position pilot (target:6966, 2′×2′, HSC-I, nine metadata visits,
43 frame rows, 40 CCDs, all three layers) is estimated at 55,212,300 bytes
(52.65 MiB) under full-square 32-bit planes. The exact offline plan is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_warp_pilot_plan_2026-11-06_02/plan.json`
(SHA-256 `a436df245a33bda09cd5599e564b2b9302e5d0d4c006b59fcf6fab1c2c33a65e`)
and caps a possible response at 70 MiB. `prepare_icl_warp_pilot_plan.py`
builds it from the already approved metadata; `download_icl_warp_pilot.py` is
plan-only by default and requires separate explicit approval of this exact
plan hash before a network request. No HSC query or download was made. See the
[scope assessment](icl_warp_product_scope_assessment_2026-11-06.md). Warp
variance plus frame metadata still needs calibration to a source-Poisson model;
keep `science_measurement=false`.

The offline FITS audit records per-HDU WCS alignment, requested-center
coverage, pixel scale, finite image fraction, mask occupancy, and finite and
positive variance fractions. A visit-complete archive with failed QA is
reported as `complete_with_qa_warnings`. Updated downloader SHA-256:
`766dbedb2b76951307b6c10a9d4be309f96ba1f509d6e09073c0716a6c71e291`.

## Radial coadd-variance balance in the supported ICL–satellite sample (2026-11-06)

[`analyze_icl_variance_radial_balance.py`](analyze_icl_variance_radial_balance.py)
checked the delivered variance plane in 94 existing matched target/control
pairs (26 `i`, 68 `i2`) across 0–900 proper kpc. The median paired variance
ratios remain near one. Three of 20 unadjusted radial intervals exclude equal
variance: `i` at 0–50 kpc (+5.8%) and `i2` at 50–75 and 75–100 kpc (−2.6%
and −4.8%). All ICL-versus-variance partial-correlation intervals include
zero. These pointwise, selected-sample diagnostics do not support a
multiplicity-aware scientific claim and do not supply effective gain or
source-Poisson variance. The 600–900 kpc median valid-pixel fraction is about
0.53. Detailed limits and hashes:
[`icl_variance_radial_balance_supported94_assessment_2026-11-06.md`](icl_variance_radial_balance_supported94_assessment_2026-11-06.md).
Results are under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_variance_radial_balance_supported94_2026-11-06_01/`.
No new HSC data request; `science_measurement=false`.

## i-band ICL-strength weak-lensing follow-up (2026-11-06)

The saved ICL/WL overlap has 37 usable `i` lenses. Applying the same median
ICL split and redshift/richness matching as the prior `i2` analysis retains
only 10 high/low pairs. The seven-bin high-minus-low random-subtracted profile
is radially mixed. Repeating it at 2.5°, 5°, and 10° block widths preserves the
same matched pairs and point profile but changes the conditional tangential
Wald diagnostic from χ²=15.79 to 161.92 to 101.28. The cross-component null is
also large at every width (χ²=24.70, 90.97, 105.50). These are not calibrated
p-values: there are only 9–24 approximate spatial blocks for 10 pairs, and the
covariance is highly width-dependent. This `i` sample does not strengthen the
ICL–lensing association. See
[`wl_icl_i_band_blockwidth_assessment_2026-11-06.md`](wl_icl_i_band_blockwidth_assessment_2026-11-06.md)
for the profile, errors, exact hashes, and limitations. No new HSC data or
image pixels were read; `science_measurement=false`.

## Tract-component balance for the ICL/WL split (2026-11-06)

The Astrafek ID audit maps the existing i and i2 ICL/WL common samples and
matched high/low pairs to the parent's approximate nearest-tract components.
The i common sample has 38 lenses in components 0/1 as 18/20; its 10 matched
pairs split 4/3 and 6/7 across those components, with five pairs crossing
components. The i2 common sample has 70 lenses in components 0/1/2 as 48/18/4;
its 29 matched pairs have low/high counts 22/21, 4/7, 3/1 and 12 cross-component
pairs. This rules out a gross low/high count imbalance for i, but does not
explain the strong cross-null: source and random-center coordinates were not
mapped, and these components approximate rather than define HSC fields.
Runner: [`audit_wl_icl_split_tract_region_balance.py`](audit_wl_icl_split_tract_region_balance.py).
Summary hashes are `8cba95a740c5d332b8c313c5687e58289652042497a5cc977061032203657130`
(`i`) and `adb313d390c9def350255dd8719d46ea9c536b1864f142b73283455421c0f6d1`
(`i2`). No new HSC data or image pixels; `science_measurement=false`.

## i-band cross-null by source and random-aperture tract component (2026-11-06)

[`audit_wl_icl_split_source_region_cross_null.py`](audit_wl_icl_split_source_region_cross_null.py)
reapplies the exact saved estimator cuts and recomputes target/random profiles
by the approximate component containing each lens or random center. It also
maps every source coordinate. The selection reproduces 53,399 target and
241,258 random used pairs; all used source-endpoint pairs remain within one
component, so cross-component source leakage is not the explanation. The
high-minus-low random-subtracted cross profiles are oscillatory in both
components, with the first two radial bins having the same signs and large
amplitudes. No region covariance is estimated, and only 3–7 lenses/group/
component are available, so the decomposition localizes the pattern but does
not diagnose its physical/instrumental cause. Summary SHA-256
`fb90c982da784dff60b270efdfe986d0fbe5b8f3a656e15dd5bc14629810862c`;
`science_measurement=false`.

## Target/random cross terms and radial-cut sensitivity (2026-11-06)

The saved endpoint profiles show that the innermost high-minus-low cross
contrast is mainly in the target stacks; random-center differences partly
modify it and reinforce it in component 1. A post-hoc radial-cut audit finds
that for `R≥1 Mpc` the three-bin cross null has conditional nominal p values
0.307/0.303/0.107 for 2.5°/5°/10° blocks, while the tangential high/low null
has p=0.0093/6.6×10⁻⁵/5.3×10⁻⁵. This flags a possible outer tangential contrast
after the inner cross-mode structure is removed. It is not a detection: the
cut followed inspection of the cross profile; the sample has 10 matched lens
pairs and only 9–24 approximate blocks. The larger i2 split does not reproduce
an outer tangential contrast at `R≥1 Mpc` (conditional nominal p=0.172/0.114
for 2.5°/5° blocks); its cross null remains consistent with zero. This i-band
outer signal is filter/sample-specific. See
[`wl_icl_i_band_blockwidth_assessment_2026-11-06.md`](wl_icl_i_band_blockwidth_assessment_2026-11-06.md)
for full endpoint terms, all radial cuts, exact hashes and limits. No new HSC
data; `science_measurement=false`.

## Continuous ICL predictor in the i2 weak-lensing stack (2026-11-06)

The continuous-lens analysis uses the exact 70 i2 lenses in the audited
ICL/WL common support. For each lens it subtracts the pooled profile of five
matched random centers, then regresses each radial tangential/cross component
on continuous target-minus-control ICL surface brightness, redshift and log
richness. The only nominal tangential-bin result is a negative slope in
2.0–3.0 proper Mpc (`p=0.021`); it fails the seven-bin multiplicity correction
(minimum Bonferroni/BH-adjusted value `0.149`). All cross-component slopes are
consistent with zero, and the outer result is not reproduced by the prior i2
median split. Radial support varies from 58 to 70 lenses, and the 27-cell
delete-one estimate uses approximate equal-area cells rather than validated
HSC regions. Treat this as an exploratory conditional association, not an
ICL–halo-mass detection. Runner:
[`analyze_wl_continuous_icl_relation.py`](analyze_wl_continuous_icl_relation.py).
Full quantitative results, caveats and hashes are in
[`wl_continuous_icl_relation_assessment_2026-11-06.md`](wl_continuous_icl_relation_assessment_2026-11-06.md).
Saved per-lens profiles and regression output are under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_continuous_icl_relation_i2_2026-11-06_02/`.
The coefficient figure is reproduced by
[`plot_wl_continuous_icl_relation.py`](plot_wl_continuous_icl_relation.py).
PNG and vector PDF files are in the result directory as
`wl_continuous_icl_relation_i2_v2.png` and
`wl_continuous_icl_relation_i2_v2.pdf`.
No new HSC data request; `science_measurement=false`.

## i2 ICL–satellite–weak-lensing triple overlap (2026-11-06)

[`analyze_icl_satellite_wl_triple_overlap.py`](analyze_icl_satellite_wl_triple_overlap.py)
joins the exact 27 i2 lenses listed by the ICL/WL overlap audit to their saved
satellite profiles and continuous, random-subtracted WL profiles. It checks
the pinned source hashes and lens-by-lens agreement of ICL, redshift, and
richness. The direct ICL–satellite partial correlation controlling redshift
and log richness is `r=0.067` (5,000-replicate connected 5° block-bootstrap
95% interval `[-0.144, 0.293]`). All 14 pointwise intervals for conditional
WL versus ICL include zero. Two of 14 conditional WL versus satellite
intervals exclude zero: cross at 0.4–0.7 Mpc and tangential at 2–3 Mpc. These
are unadjusted findings among 29 reported relationships in a selected sample
of 27 lenses; do not treat them as detections. Recomputed fits use 18–20
approximate connected blocks and no block contains more than two lenses.
Full interpretation and SHA-256 values are in
[`icl_satellite_wl_triple_overlap_assessment_2026-11-06.md`](icl_satellite_wl_triple_overlap_assessment_2026-11-06.md).
Outputs are under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_wl_triple_overlap_2026-11-06_01/`.
No HSC query or download; `science_measurement=false`.

## Full-sample radial coadd-variance balance (2026-11-06)

The supported-pair variance diagnostic was extended to all 783 same-filter
target/control coadds (302 `i`, 481 `i2`). A joint 10,000-resample max-t
bootstrap adjusted the 20 radial contrasts across filters. Three `i2` bins
retain simultaneous intervals below equal variance: 75–100 kpc (ratio 0.961,
interval 0.924–0.999), 100–150 kpc (0.955, 0.915–0.996), and 200–300 kpc
(0.961, 0.934–0.988). No `i` interval excludes equality. The joint partition
has 38 approximate blocks, with one block containing 107 pairs; these results
remain a selected-sample variance-plane diagnostic, not a total-noise or ICL
signal correction. Details and reproducible hashes are in
[`icl_variance_radial_balance_full783_assessment_2026-11-06.md`](icl_variance_radial_balance_full783_assessment_2026-11-06.md).
Pair tables and simultaneous intervals are under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_variance_radial_balance_full783_2026-11-06_01/`
and
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_variance_radial_balance_familywise_2026-11-06_02/`.
The initial `_01` familywise output is marked invalid in its summary because
it pooled filter rows; `_02` contains the corrected filter-separated result.
No new HSC data; `science_measurement=false`.

## Variance-adjusted ICL–satellite association sensitivity (2026-11-06)

[`analyze_icl_variance_adjusted_satellite_assoc.py`](analyze_icl_variance_adjusted_satellite_assoc.py)
joins the 94 saved ICL–satellite supported pairs (26 `i`, 68 `i2`) to the
existing full-sample delivered-variance table at 200–300 kpc. In paired
5,000-resample connected 5-degree block comparisons, the `i2` partial
correlation controlling redshift and richness is `r=0.293` (95% interval
0.043–0.498); adding log target/control variance ratio gives `r=0.285`
(0.031–0.499), a paired change of −0.008 (−0.078–0.034). Adding a coarse
mean-log-variance depth proxy as well gives `r=0.260` (−0.018–0.493), change
−0.032 (−0.148–0.027). All paired changes include zero. The `i` sample is
more weakly constrained (N=26) and all its model intervals include zero.
Thus the delivered-variance imbalance does not explain the supported-sample
association in this sensitivity analysis. These selected-sample correlations
are not causal or calibrated measurements; the variance plane omits source
Poisson/effective-gain terms and does not measure ICL-estimator covariance.
Full assessment and hashes are in
[`icl_variance_adjusted_satellite_assoc_assessment_2026-11-06.md`](icl_variance_adjusted_satellite_assoc_assessment_2026-11-06.md).
Output:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_variance_adjusted_satellite_assoc_2026-11-06_01/`.
No HSC query/download; `science_measurement=false`.

## Joint support-selection and pixel-noise sensitivity (2026-11-06)

[`analyze_icl_satellite_ipw_pixelnoise_joint.py`](analyze_icl_satellite_ipw_pixelnoise_joint.py)
combines stabilized IPW for the 94/224 satellite-profile support subset with
the saved finite-lag ICL/satellite pixel covariance. The baseline noise-
corrected `i2` partial correlation is `r=0.238` (2,000 connected-block
95% interval −0.023–0.446), versus 0.306 without the support IPW. Adding
variance balance to the support model gives `r=0.278` (−0.018–0.495); adding
a coarse mean-depth proxy gives `r=0.201` (−0.054–0.410). Paired changes from
the baseline include zero. `i` is imprecise, and the fuller propensity model
reduces its effective N to 13.5/26. This does not pass the science gate:
missing-at-random selection and total-noise calibration remain unresolved.
The exact covariance scope audit covers all 94 selected pairs. Methods and
hashes are in
[`icl_satellite_ipw_pixelnoise_joint_assessment_2026-11-06.md`](icl_satellite_ipw_pixelnoise_joint_assessment_2026-11-06.md).
Output:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_ipw_pixelnoise_joint_2026-11-06_01/`.
No HSC query/download; `science_measurement=false`.

## Variance terms in the ICL–satellite support model (2026-11-06)

[`analyze_icl_satellite_ipw_variance_sensitivity.py`](analyze_icl_satellite_ipw_variance_sensitivity.py)
extends the 224-to-94 satellite-profile support analysis by adding the
200–300 kpc target/control delivered-variance ratio and a coarse mean-log-
variance proxy to the propensity model. The `i2` IPW partial correlation is
0.227 (95% paired spatial-block interval −0.035–0.419) at baseline, 0.266
(−0.026–0.470) with variance ratio, and 0.193 (−0.048–0.401) with both
variance terms. Paired changes from baseline are +0.039 (−0.069–0.162) and
−0.034 (−0.163–0.092), both including zero. The richer model reduces `i`
effective N to 13.5/26 and raises its maximum stabilized weight to 5.08.
This does not establish an IPW correction: the missing-at-random assumption
is unverified, and delivered variance is not total noise. Full results and
hashes are in
[`icl_satellite_ipw_variance_sensitivity_assessment_2026-11-06.md`](icl_satellite_ipw_variance_sensitivity_assessment_2026-11-06.md).
Output:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_ipw_variance_sensitivity_2026-11-06_01/`.
No HSC query/download; `science_measurement=false`.

## ICL-strength split stacked satellite profiles (2026-11-06)

[`analyze_icl_strength_split_satellite_profiles.py`](analyze_icl_strength_split_satellite_profiles.py)
stacks the existing cluster-level satellite excess profiles in the fixed
94-pair ICL/satellite support sample (26 `i`, 68 `i2`), splitting each filter
at its median 200–300 kpc ICL signal. It estimates the redshift- and
richness-adjusted high-minus-low contrast at seven radii per filter. The
10,000-replicate restricted Webb wild-cluster bootstrap shares spatial-block
weights across both filters; Holm and max-T adjust all 14 contrasts. All
pointwise CR1 95% intervals include zero, every Holm-adjusted p-value is 1.0,
and the minimum max-T p-value is 0.520. The i2 contrasts are positive, but not
resolved from zero after clustered uncertainty; the i estimates are noisy and
change sign with radius. This post-hoc selected-sample split is not evidence
of a physical ICL–satellite relation. The figure is produced by
[`plot_icl_strength_split_satellite_profiles.py`](plot_icl_strength_split_satellite_profiles.py).
Methods, hashes, and limitations are in
[`icl_strength_split_satellite_profiles_assessment_2026-11-06.md`](icl_strength_split_satellite_profiles_assessment_2026-11-06.md).
Machine-readable outputs are in
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_icl_strength_stacked_profiles_2026-11-06_01/`.
No new HSC query/download; `science_measurement=false`.

## Continuous ICL–satellite profile relation (2026-11-06)

[`analyze_icl_satellite_radial_relations.py`](analyze_icl_satellite_radial_relations.py)
compares both the continuous one-standard-deviation ICL slope and the median
split contrast across seven radial bins, adjusting for redshift and log
richness in the same 94-pair sample. Holm and max-T corrections cover all 28
filter/radius/method comparisons with shared 5-degree block weights. The
smallest raw bootstrap p is 0.0369 (i2, 50–75 kpc), but its max-T p is 0.770;
the minimum max-T p in the family is 0.190 (continuous i, 5–10 kpc). Holm p is
1.0 for all 28. The three negative pointwise CR1 intervals in the inner i
bins are exploratory and do not pass the shared-block tests. No relation is
established; `science_measurement=false`. Figure and full methods are in
[`icl_satellite_radial_relations_assessment_2026-11-06.md`](icl_satellite_radial_relations_assessment_2026-11-06.md)
and `plot_icl_satellite_radial_relations.py`. Outputs:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_radial_relations_joint_2026-11-06_01/`.
No new HSC request/download.


## 2026-11-06: paired weak-lensing boost-profile sensitivity

The runner `analyze_wl_source_density_boost_blocks.py` now also reports the
explicitly unvalidated sensitivity profile `B_W(R) * DeltaSigma_lens-random`,
recomputing the nonlinear product inside paired 2.5- and 5-degree spatial
jackknife replicates. Across seven radii, the exact response-adjusted S19A
estimator-weight boost is 1.125–1.163 and raises the tangential profile by
12.5–16.3%; the two block-width central profiles agree to stored precision.
This is not an applied dilution correction: magnification, depth, mask, PSF,
sky transfer, source selection, source-redshift uncertainty, reduced shear,
and regional covariance remain unresolved. The saved profile is unchanged,
`boost_correction_applied=false`, and `science_measurement=false`. Methods and
hashes are in
[`wl_boost_profile_sensitivity_assessment_2026-11-06.md`](wl_boost_profile_sensitivity_assessment_2026-11-06.md).
Output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_source_density_boost_spatial_2026-11-05_06/`.
No HSC query, download, or image read occurred.


## 2026-11-06: forced CAMIRA-position coverage sensitivity

The forced-position runner now supports an explicit annular-coverage floor;
the science-default remains 60%. On existing cutouts, unmasking other
unmatched CAMIRA positions did not change support at 60% (20 retained
position-bin measurements versus 20) and changed target coverage by at most
0.025 across the candidate set. Lowering the exploratory floor to 20% while
keeping at least three rotated controls increased 50–75 kpc forced positions
from 9 to 65 (25 `i`, 40 `i2`). The forced-only `i2` association is imprecise
(`r=0.147`, block interval −0.485–0.425); the combined 74-pair estimate is
`r=0.257`, with paired 20%-vs-60% same-sample change +0.003 (−0.070–0.126).
This is post-hoc, coverage-selected and not a detection; missing pixels may be
non-random and satellite identity is unconfirmed. Full assessment:
[`icl_satellite_forced_coverage_sensitivity_assessment_2026-11-06.md`](icl_satellite_forced_coverage_sensitivity_assessment_2026-11-06.md).
Runner: `measure_forced_satellite_profiles.py`; paired comparison:
`analyze_icl_satellite_forced_coverage_sensitivity.py`. Results are under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_forced_coverage_sensitivity_2026-11-06_02/`
and `satellite_forced_outer_profile_unmatched_camira_unmasked_cov20_2026-11-06_01/`.
No HSC query/download; `science_measurement=false`.


## 2026-11-06: forced-position radial-selection audit

The 20%-coverage forced stack is strongly selected in host-centric radius. Of
65 retained 50–75 kpc profiles, 41 (63%) lie at 300–500 kpc from the cluster
center and none lies inside 100 kpc. Parent candidate-position median radii
are 200 kpc (`i`) and 195 kpc (`i2`), versus 388 and 339 kpc among retained
profiles. Only 18/282 candidate-bearing `i` clusters and 32/439 `i2` clusters
retain any 50–75 kpc profile. Positivity fails inside 100 kpc, so this sample
cannot be IPW-corrected into the parent radial population. Keep forced and
detection-based profiles separate; `i2 r=0.257` remains a sensitivity result,
not a corrected ICL–satellite relation. Details:
[`icl_satellite_forced_position_selection_assessment_2026-11-06.md`](icl_satellite_forced_position_selection_assessment_2026-11-06.md).
Runner `audit_icl_satellite_forced_position_coverage_selection.py`; output
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_forced_position_coverage_selection_2026-11-06_01/`.


## 2026-11-06: contributing-member count sensitivity

On the fixed 94-pair sample, the `i2` ICL–satellite partial correlation is
0.293 after linear redshift and richness adjustment. Adding the log number of
members in each 50–75 kpc stack changes the estimate to 0.295. The paired
5-degree block-bootstrap change is 0.002 (95% interval −0.023–0.044; 9,989
valid draws out of 10,000 requested). For `i`,
the estimate changes from 0.079 to 0.088, with change 0.009
(−0.081–0.192). Member count has no resolved partial association with either
ICL or satellite excess in the two filters. The analysis reads existing local
tables and leaves `science_measurement=false`.

Runner: `analyze_icl_satellite_membercount_adjustment.py`.
Assessment: `icl_satellite_membercount_adjustment_assessment_2026-11-06.md`.
Results: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_membercount_adjustment_2026-11-06_02/`.
No HSC query, download, or image-pixel read occurred.


## 2026-11-06: long-lag covariance propagated through ICL apertures

The ±8-pixel residual covariance kernels were contracted with the exact
200–300 proper-kpc aperture masks for all 224 supported target/control pairs.
The independent-frame noise sigma rose by 0.8% in `i` and 1.3% in `i2`
relative to the earlier 9×9 kernel. The resulting estimates are 0.0793 and
0.0745 nJy arcsec⁻², about 7.6% and 7.0% of the saved spatial-block interval
half-widths. The stationary-kernel calculation omits source Poisson, sky
transfer, PSF wings, and covariance beyond eight pixels. No correction was
applied and `science_measurement=false` remains.

Runner: `analyze_icl_long_lag_aperture_covariance.py`.
Assessment: `icl_long_lag_aperture_covariance_assessment_2026-11-06.md`.
Results: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_aperture_covariance_200_300_2026-11-06_01/`.
No HSC query or image download occurred.


## 2026-11-06: random-point variance result auditor prepared

`audit_icl_random_pixel_variance_results.py` is ready to validate the exact
approved-plan checkpoint, every batch hash and row, center/patch identity,
annulus membership, mask flags, and target/control pairing. The auditor writes
1,598 center summaries and 798 paired summaries. The exact HSC plan remains
unsubmitted pending its separate user approval, so no result audit has run.
The pipeline keeps `science_measurement=false` and does not interpret random
point variance as a coadd weight or source-Poisson correction.


## 2026-11-06: host-centric satellite-position mix sensitivity

The fixed 94-pair ICL–satellite sample was re-evaluated after adding the median host-centric radius of detected members contributing to each 50–75 proper-kpc satellite stack. With redshift and linear richness controlled, the `i2` partial correlation changes from 0.293 to 0.322 (paired change +0.030, shared 5-degree block-bootstrap 95% interval −0.019–0.101). In `i`, it changes from 0.079 to 0.029 (change −0.050, interval −0.330–0.044). All position-association and model-change intervals include zero. This measured position-mix adjustment does not resolve the selected-sample association and is not a population/completeness correction. Keep `science_measurement=false`; source Poisson, PSF wings, sky transfer, and total image covariance remain open.

Runner: `analyze_icl_satellite_clustercentric_position_adjustment.py` (SHA-256 `65e8557516e6df8633a85e75ad519bf343f997f52a5e223424f7e5e4f9e55f86`). Assessment: `icl_satellite_clustercentric_position_adjustment_assessment_2026-11-06.md`. Results: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_clustercentric_position_adjustment_2026-11-06_01/`; summary SHA-256 `78317905c3d0b2b0ada48c8b39a27ae564caaeb0bd744f365c12fedc0a7d8971`. Existing local tables only; no HSC request, download, or pixel read.


## 2026-11-06: S20A coverage-query submission audit

The previously approved metadata-only S20A plan covers 2,596 candidates and 23,364 center/boundary probes in 52 batches. Its first CAS submission returned HTTP 406. After correcting the mosaic join syntax, the regenerated plan preserved the candidate-center and probe files byte-for-byte, passed offline validation, and also received HTTP 406 on its first batch. Both checkpoints show no CAS job ID and no results; no later batches ran. The remaining 406 cause is undiagnosed. Do not retry before checking CAS job state and resolving the server rejection. No images, source rows, shapes, or photo-z were requested. Details: `icl_expansion_coverage_submission_assessment_2026-11-06.md`.


## 2026-11-06: current ICL and weak-lensing evidence synthesis

Added `build_icl_wl_current_evidence_report.py`, which checks 13 saved component summaries for `science_measurement=false`, hashes every source summary, and writes a single ICL/satellite plus weak-lensing evidence report. It reports the 94-pair ICL result under raw, finite-lag noise, and support-weighted models, and documents the 435-lens lensing sample's strong regional influence. The current results remain exploratory; no profile was promoted. No HSC query, download, or pixel read occurred.

Output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_03/icl_wl_current_evidence_report.md` (SHA-256 `fb061edb26f6ea25666eefcc644b464bc4b54eb68239f649a920ca8145c4439f`). Summary SHA-256 `1b6fd47fea010fbb8622def7326e9ca62e467aff8e98dc1f713b9f241ea6644e`. Runner SHA-256 `bc0a9dd6b85bb87eb59eaa6eee5ded8437fe3b40442cca2e8c9db8abf8c6041d`.

## 2026-11-06: frame-metadata sensitivity of the ICL--satellite relation

`analyze_icl_frame_metadata_satellite_relation.py` joined the re-audited
67,071-row CCD-frame metadata export to the fixed 94-pair ICL/satellite
support. It added five coarse target/control frame-summary offsets one at a
time to the existing redshift- and linear-richness-adjusted partial
correlation, with 5,000 shared connected-block bootstrap draws and one max-T
family across ten filter/metadata comparisons. No adjusted change is resolved:
all max-T p-values are >=0.895, and every simultaneous interval includes
zero. The i2 estimate remains 0.286--0.309 around the 0.293 baseline; the
larger changes in the small i sample are uncertain. This excludes neither
subpixel depth differences nor Poisson, sky-transfer, PSF-wing or selection
effects, and does not promote `science_measurement`.

The authoritative result is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_metadata_satellite_relation_2026-11-06_03/summary.json`
(SHA-256 `cd2b17e5b9e3579ae580b0e7531b54f5ee6e9508477f77d1b3b78bb7c9a02643`);
the joined pair table SHA-256 is
`368137cae49ff728f97fc97562721dd3cbcaaa5765aa83d32c65d664ccec5950`.
Assessment: `icl_frame_metadata_satellite_relation_assessment_2026-11-06.md`.
The earlier `_01` output used log richness and is superseded. `_02` used the
correct estimator but omitted the runner digest; `_03` is the final
reproducible output. No additional CAS request or image read was made.

The consolidated evidence report was rebuilt to include this result:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_07/icl_wl_current_evidence_report.md`
(SHA-256 `dba2cb79a2a3f30bea8cb1ce132bdf49539d822871aa7c6d765fbe4cf32e4efb`);
summary SHA-256 `e31797d0ad936f43a115817f589fa34d0308f2c97cf4e7bba561192c221afb4b`.
Builder SHA-256 `98b51e03b360a3ff89a0c231a31e71eb99b27b833a0329a79b2e3fe0660ba5ee`.

## 2026-11-06: nonlinear redshift/richness model-form sensitivity

`analyze_icl_satellite_nonlinear_adjustment.py` adds centered/scaled
quadratic redshift and richness terms and their interaction to the established
linear adjustment, on the same 94 supported ICL/satellite pairs. Four models
per filter are compared to baseline with 5,000 connected-block resamples; a
single max-T family covers all eight changes. The `i2` association remains
near its baseline (r=0.293; adjusted values 0.294--0.305), and every paired
change interval includes zero. The `i` estimate is less stable but remains
imprecise on 26 pairs. This post-hoc sensitivity does not resolve confounding
or selection and does not change `science_measurement=false`.

Assessment: `icl_satellite_nonlinear_adjustment_assessment_2026-11-06.md`.
Result summary:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_nonlinear_adjustment_2026-11-06_01/summary.json`
(SHA-256 `4777d93629348a8849a74af8aa3217b4553ff27142d7802724b8c6091a5c332c`).

The consolidated evidence report now includes the nonlinear model-form result
as well as the frame-metadata association sensitivity:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_08/icl_wl_current_evidence_report.md`
(SHA-256 `812239a91e0071ae8abcaa7470a925d83db58ca3f2f2d0993a089f666a0b7cdf`);
summary SHA-256 `67075fdc8b032075dd898eb59efddde03c045d83f26eb5edf0a4461a3279827c`.
Builder SHA-256 `a2d5d5ab297ee403d2417c8a28d76757dc88429f0d4f422c5dd533bdf1c0a996`.

## 2026-11-06: random-only weak-lensing placebo

`analyze_wl_random_only_placebo.py` revalidated the already approved random
and cluster source-shape exports, selected the exact 435-lens common-support
sample and 2,175 matched random centers, reproduced the stored random-only
profile, and estimated a 5-degree spatial delete-one-block covariance over 50
occupied cells. The tangential random-only profile is consistent with zero
(covariance-conditional χ²=6.90, rank 7, p=0.439), as is the cross profile
(χ²=4.06, p=0.773). This is a placebo check for large random-stack residuals;
it does not validate the cluster signal or establish a field covariance.
`science_measurement=false` remains.

Assessment: `wl_random_only_placebo_assessment_2026-11-06.md`.
Output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_only_placebo_2026-11-06_01/summary.json`
(SHA-256 `9c2b6ad3f1525f18e26c8bf3ec947e40de49abc40271dd75313dad32b2b50fa2`);
covariance SHA-256 `c312159b7efc861b299f7b513348fdbac1bebbcb2033e99d19ea3f6793d4ed99`.
Runner SHA-256 `44fb39796180f452dc6283f443cb5b4c0de21cc6fd6fcdcb7bcf3034ac35606b`.

The comparison figure is generated by `plot_wl_random_only_placebo.py`
(SHA-256 `0baf2ccb1ff376cae335c29a637883f958eb2de3477db356568cfa13d3d77fd7`):
[PNG](/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_only_placebo_2026-11-06_01/random_subtracted_vs_random_only.png)
and
[PDF](/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_only_placebo_2026-11-06_01/random_subtracted_vs_random_only.pdf).

The consolidated report including this placebo is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_10/icl_wl_current_evidence_report.md`
(SHA-256 `5e56375a4a16bff9b257ac31f6664631a2bb0beef2f9246345654c8a9c0b9f2e`);
summary SHA-256 `4b5073163ce0fb226f185306e6c2d1c4be2546dd8a7b14b1725e6c5e41754c99`.
Builder SHA-256 `e07a66ce4e4e3c5e42b4869033c02f5921a901a58c7115a69c9e760859aa370b`.

## 2026-11-06: spatially held-out ICL–satellite prediction

The new leave-one-connected-block-out analysis finds no clear prediction gain
from adding ICL to redshift and richness: held-out MSE changes by +6.0% in `i`
and −3.5% in `i2`; paired block-bootstrap intervals include zero for both.
See [the assessment](icl_satellite_spatial_cross_validation_assessment_2026-11-06.md)
and [the Astrafek runner](analyze_icl_satellite_spatial_cross_validation.py).
Runner SHA-256 `a068c2ac0330a012964596e7bafc11a3c3d2be7aca293866a0e90b4024b8b510`;
result summary SHA-256 `7059c7b196467afb581f0edf99ce69a65b9e132979f7a3284be8182c6cc1ceb7`.
The cumulative evidence report was rebuilt with this diagnostic:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_11/icl_wl_current_evidence_report.md`
(SHA-256 `3023a095d419a75c8ded3e4a2ee44ae736341d257101acd00e812322be92dc3b`);
summary SHA-256 `87dde195464425c552a05f4159f89224536065c9815e503128ab0ebeb00cf98b`.
Updated builder SHA-256 `cee435d9c6f5a72c2fc51115814e5a5a198f87d28bc041cb388b11df1120d44e`.

## 2026-11-06: fixed-support radial ICL–satellite relation

Astrafek now also checks all seven satellite radial bins on the same cluster
set within each filter: 21 `i` pairs and 53 `i2` pairs. The `i2` 50–75 kpc
pointwise slope is positive, but its 14-test max-T p-value is 0.410; the
smallest adjusted p across all filters/radii is 0.237. No radius passes the
familywise check. See
[the assessment](icl_satellite_fixed_profile_assessment_2026-11-06.md) and
[the runner](analyze_icl_satellite_fixed_profile.py). The figure is
`icl_satellite_fixed_sample_radial_slopes.png` in
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_fixed_profile_2026-11-06_02/`.
Runner SHA-256 `b19d18de0b54e4918463b55232ac260e6b06e109e585475ed42511cbdcde187c`;
summary SHA-256 `08e60ec574a624f4c177c7fbbc66c9ee93a57a4072cbea1755acf23aaeda113d`.
Assessment SHA-256 `57bfb9a077790a18b746f2835f12a98d4b68237c5cebbfd73ce11db87983b1b8`.
The updated consolidated evidence report is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_12/icl_wl_current_evidence_report.md`
(SHA-256 `9dc5065c0d4ebd4ff9f4e04c4e07b6bab8e6022eb533a0d962a8af7c89e1a011`);
summary SHA-256 `b458ce7120ec94f1b3c7b4d6f1d66732da87d81ca7db89450812b473f5bd04be`.
Updated report builder SHA-256 `cfe710f5909e48a172fe19357f4734396734df824e291dc57dad4575baa5d7b9`.

The fixed-sample high/low ICL split uses the same complete-profile clusters at
all seven radii. The adjusted outer high-minus-low differences are +4.61
nJy arcsec⁻² (`i`) and +1.35 (`i2`); neither is resolved after 14-test
max-T correction, whose minimum p-value is 0.452. Details:
[assessment](icl_satellite_fixed_profile_split_assessment_2026-11-06.md) and
[runner](analyze_icl_satellite_fixed_profile_split.py). Figure:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_fixed_profile_split_2026-11-06_02/icl_satellite_fixed_sample_icl_split.png`.
Runner SHA-256 `79a54e88e03d1e2827701d276264a368e4da149a2d5b74afd4fb2aa7bd382d30`;
summary SHA-256 `49403d38eb267f61bee9812c7b2a24694a2fa45b5463fb2db60c13cbd4d1f909`.
The updated consolidated report is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_13/icl_wl_current_evidence_report.md`
(SHA-256 `3982044b82aa1636531007544f8355986eb7525761c1e036cd45a684f698555e`);
summary SHA-256 `269248dbced884a5373c643f1eff3bfd3e387e1cb41e0d015a3e467fd2c95e26`.
Updated report builder SHA-256 `8e7daef251292fc186944b5eb6466dd56d996b717a60dec9796afb1d304346fa`.

## Approved-scope warp noise validation

The single-position target-6966 warp pilot has been downloaded and audited. Its
assessment and result hashes are recorded in
`icl_warp_pilot_coadd_variance_assessment_2026-11-06.md`. To prepare the next
request, covering the 94 existing ICL–satellite supported target/control pairs
(26 `i`, 68 `i2`), run:

```bash
python3 prepare_icl_supported94_warp_noise_plan.py
python3 download_icl_supported94_warp_noise.py
```

Both commands are plan-only by default. The exact plan SHA is
`bd936a94bcdbff83a607acbfca66616c396db526ada50cc89c9e3e312d37007a`; it covers
188 2′×2′ HSC-I warp requests and 1,471 expected visit products, estimates
8.45 GiB, and enforces a 160 MiB per-request cap and a 12 GiB total cap. **Do not
add `--submit` until the user explicitly approves this exact scope and hash.**
The approval-gated command after that approval is:

```bash
python3 download_icl_supported94_warp_noise.py --submit \
  --confirm-plan-sha bd936a94bcdbff83a607acbfca66616c396db526ada50cc89c9e3e312d37007a
```

The runner is sequential and checkpoints each request for resume. It stores no
credentials and keeps `science_measurement=false`.


### Approved 94-pair warp acquisition is complete

The user-approved plan `bd936a94bcdbff83a607acbfca66616c396db526ada50cc89c9e3e312d37007a` has completed: 184/188 archives arrived, totaling 9.26 GiB; four repeat HTTP 404s are documented in the assessment. Do not rerun the acquisition unless a new, exact request scope is separately approved. The offline variance comparison is `analyze_icl_supported94_warp_noise.py`; results and scientific limits are in `icl_supported94_warp_noise_assessment_2026-11-06.md`. The final output is under `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_supported94_warp_noise_analysis_affine_03_2026-11-06/`. The legacy plan commands above describe the completed historical process. `science_measurement=false` remains.


### Same-visit split-stack covariance pilot

`analyze_icl_visit_split_covariance_pilot.py` makes inverse-variance half-stacks from the already downloaded HSC-I visit warps for one supported `i2` target/control pair. Its normalized split-difference covariance is still a pilot: one visit per center has no valid pixels, half-stack depths differ, and nearest-neighbor mapping does not model HSC resampling covariance. The 200–300 kpc normalized RMS before plane fitting is 1.093/1.166 (target/control); the ±8-pixel cumulative lag sums after the local plane estimator are 0.948/1.187. See `icl_visit_split_covariance_pilot_assessment_2026-11-06.md`. It makes no noise correction and leaves `science_measurement=false`.


### Consolidated report includes the supported-94 warp audit

`build_icl_wl_current_evidence_report.py` now validates the approved plan/download hashes and incorporates the 184-center coadd/visit-variance comparison. Its latest 23-input output is under `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_17/`. The report records four repeat-404 gaps, 9.26 GiB received, the 1.00669 center-median variance ratio, and 81.9% median strict coverage at 200–300 kpc. The nearest-neighbor mapping remains approximate; exact HSC resampling covariance and source-Poisson/PSF/sky terms remain open.


### Strict warp-mask association sensitivity

Run `analyze_icl_supported94_warp_mask_association.py` to reproduce the saved 200–300 kpc ICL flux under strict warp-valid pixel support. It checks baseline flux/coverage reproduction, retains the original 40% coverage floor, and compares correlation before/after masking on identical retained pairs. Result: 81 pairs (23 i, 58 i2) from 90 complete warp pairs, with paired correlation-change intervals including zero. This is a selected-subset sensitivity, not a science correction. See `icl_supported94_warp_mask_association_assessment_2026-11-06.md`.


### Supported-sample visit-split covariance audit

`analyze_icl_visit_split_covariance_sample.py` applies the same estimator to all 90 complete supported pairs in the approved warp set (24 `i`, 66 `i2`) and saves per-center diagnostics plus 17×17 kernels. Empty lag cells are treated as unsupported, not zero. A separate joint bootstrap runner adjusts four metrics within each filter using connected 5° blocks. The simultaneous intervals do not resolve target/control differences in split RMS or ±8 covariance. The small positive i2 zero-lag and coverage differences, differing visit sets, and nearest-neighbor resampling remain audit findings; no correction is applied. See `icl_visit_split_covariance_sample_assessment_2026-11-06.md`.


### Heteroscedastic white-noise null

`simulate_icl_visit_split_covariance_null.py` draws noise on the native warp pixels from the supplied visit variance planes and propagates it through the saved mask, split, nearest-neighbor map, inverse-variance stack, local plane fit, clipping, and ±8 covariance estimator (100 realizations for the cluster-101 i2 pair). The null reproduces an estimator-induced negative lag tail but not the observed RMS or covariance excess. This leaves correlated warp noise, residual sky/source structure, and variance calibration unresolved; it does not provide an ICL correction. See `icl_visit_split_covariance_null_assessment_2026-11-06.md`.


### Target/control covariance by annulus

`analyze_icl_annular_coadd_covariance.py` compares masked HSC coadd lag kernels in the 200–300 proper-kpc ICL annulus and a 350–500 kpc outer annulus for the fixed supported 94-pair sample. Paired target-control and radial difference-in-differences 95% intervals include zero in both `i` and `i2`. Strict patch requirements leave only 8 `i` and 16 `i2` complete pairs in the ICL annulus, so the test is underpowered and cannot justify transferring the outer kernel inward. The outer annulus may contain cluster-envelope light and is not treated as blank sky. No correction; `science_measurement=false`. See `icl_annular_coadd_covariance_assessment_2026-11-06.md` and output `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_annular_coadd_covariance_2026-11-06_01/`. Runner SHA-256 `16ad75b7f9ac3f2e0cd2e1759a84037f4c9519307e982a03e32174526b80f7d3`; summary SHA-256 `c03b1c60fcc3415aee615c8685700b42633fa769df0ac4a8173013bafe2f93b0`.


### Familywise ICL–satellite–WL association audit

`adjust_icl_satellite_wl_triple_overlap_max_t.py` tests all 29 conditional relations in the selected 27-lens i2 ICL/satellite/WL overlap with a restricted Webb six-point wild-cluster bootstrap and shared 5-degree spatial-block multipliers. The two originally pointwise patterns become max-T p=0.672 (WL cross vs satellite, 0.4–0.7 Mpc) and p=0.923 (tangential, 2–3 Mpc); none of 29 survives familywise adjustment. This selected 27-lens sample and approximate 20-block partition remain exploratory; no detection or no-effect conclusion. `science_measurement=false`. Assessment `icl_satellite_wl_triple_overlap_max_t_assessment_2026-11-06.md`; runner SHA-256 `162c90585875d0811f6921407c77ec72837f6ecfb5fd9e571bf2c5bbde5e156d`; summary SHA-256 `d5e2d3aa05be32b07be5fa62b8de467aea5314411daba49cb733444850126b0e`; output `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_wl_triple_overlap_maxt_2026-11-06_01/`.


### Latest consolidated ICL/WL evidence report

`build_icl_wl_current_evidence_report.py` now hash-validates and includes the annular coadd covariance comparison and 29-way ICL–satellite–WL wild-cluster max-T audit, alongside the previously consolidated evidence. Report: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_19/icl_wl_current_evidence_report.md`; report SHA-256 `5a4d2f7517e3fa78625c19ec96e51dbe8920c0b3c1bcf6ef3293edb6a7bf5426`; summary SHA-256 `ff2d720b09d6f0112670c956dad9053d440b02508f6ce2f4484d6a28b12d941f`; builder SHA-256 `af028a910fe9e5308d4b382029544f9968818e9dd77c97ca08ba1cbd3298418a`. It reads summaries and pinned audit tables only, with no HSC request or image-pixel access. The approved random-point variance query remains pending: no HSC credential environment variables or live submitter process were present in this Lageunha session.


### ICL–satellite–WL common-support selection audit

`audit_icl_wl_overlap_selection.py` compares the 27 i2 lenses with verified random-subtracted WL support against the other 41 of 68 i2 ICL–satellite-supported lenses. Only 9 of 33 connected 5-degree blocks contain both groups; those blocks hold 10 selected and 16 excluded lenses, while 17 selected and 25 excluded lenses are in single-status blocks. This is a spatial positivity failure, so do not use IPW to transport the 27-lens relation to all 68. Within the nine mixed blocks, four-variable max-T balance finds no adjusted redshift/richness/ICL/satellite difference, but has low power. Shape-pair counts differ by construction (median 4,342 vs 0). Assessment `icl_wl_overlap_selection_i2_assessment_2026-11-06.md`; output `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_overlap_selection_i2_2026-11-06_03/`; runner SHA-256 `b2d474d5ecc234e3b3fcb9bea21739f72af720d895d730017c1dfba111895a94`; summary SHA-256 `fef9cdbfdd595f770d1a1a933d29618443524896a318298266a6a1134cb92a98`. `science_measurement=false`.

The current consolidated report, including annular covariance, max-T association, and common-support selection audits, is `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_20/icl_wl_current_evidence_report.md` (SHA-256 `dae6062a10ca70f0b2eeb6f2fdfbcae76e52665a1bf2ffac31eea9cfa9a7826a`).

### Lanczos-5 warp-to-coadd variance sensitivity

`analyze_icl_supported94_warp_noise_lanczos5.py` applies a separable Lanczos-5 kernel to map the already downloaded visit-warp variance planes onto coadd pixels. It analyzed 184/188 approved centers (four repeat HTTP 404 gaps). Median 200–300 proper-kpc coadd/predicted variance ratios are 1.00707/1.00760 for i target/control and 1.00639/1.00611 for i2; paired target-minus-control block-bootstrap intervals include zero in both filters. Against nearest-neighbor mapping, the median center-level ratio changes are −0.000245 (i) and −0.000161 (i2), with median strict shared-mask coverage 69.0–71.9%. Warp-pixel errors are treated as diagonal, so this does not model covariance already induced by HSC CCD-to-warp resampling or reproduce exact coadd weights, clipping, and mask interpolation. No noise correction; `science_measurement=false`. Assessment `icl_supported94_warp_noise_lanczos5_assessment_2026-11-06.md`; summary SHA-256 `16c8afbcc8287444860bacdcfd01038105bd51d5e82d61e44d8b82d8853e9f74`; per-center table SHA-256 `a0322cc88e63768d172ccc8da52ce5148bef1a4407232b53e93e00bc9ae0f06f`; runner SHA-256 `0df14f6c6be1d858a38437c3ae7300adce81fba859047c8d29d8fe4d07232db9`. No new HSC request or download.


### Selection into the i2 ICL–satellite–WL common sample

`audit_icl_wl_overlap_selection.py` audited WL attrition from the 68 i2 systems with ICL and satellite support: 60 are in the matched 620-lens parent, 35 have target shape rows, 31 have positive random-shape counts, and 27 have verified random-subtracted WL support. These 68 systems occupy 33 connected 5-degree blocks, but only 9 blocks contain both selected and excluded lenses (10 selected and 16 excluded); 17 selected and 25 excluded lenses are in single-status blocks. This fails spatial positivity for IPW transport from 27 to 68. Four-variable within-block max-T finds no adjusted redshift, richness, ICL, or satellite imbalance, but the nine mixed blocks provide little power. Target shape-pair counts (median 4,342 vs 0) are a direct eligibility mechanism. Do not generalize the 27-lens relation using IPW. Assessment `icl_wl_overlap_selection_i2_assessment_2026-11-06.md`; final output `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_overlap_selection_i2_2026-11-06_03/`; runner SHA-256 `b2d474d5ecc234e3b3fcb9bea21739f72af720d895d730017c1dfba111895a94`; summary SHA-256 `fef9cdbfdd595f770d1a1a933d29618443524896a318298266a6a1134cb92a98`.

### Full-783-pair annular covariance scale sensitivity

`analyze_icl_full783_annular_covariance.py` extends the delivered-coadd lag-kernel audit to all 783 same-filter target/control pairs. The 64×64 configuration has complete annular support for 17/302 i and 48/481 i2 pairs; all three within-filter familywise max-T intervals include zero, and support is strongly low-redshift selected. A 32×32 sensitivity (max 64 patches, minimum 8, 80% valid pixels, ±8-pixel lag kernel) raises support to 215/302 i and 402/481 i2, but remains redshift selected (SMD −0.65/−0.39). In 32×32 results, the i cumulative lag-sum interval and i2 zero-lag interval exclude zero after familywise adjustment; annulus coverage differs in both bands. This scale-sensitive, mask-selected delivered-coadd diagnostic may include astrophysical residual structure and is not a pure noise estimate or correction. `science_measurement=false`. Assessment `icl_full783_annular_covariance_32px_assessment_2026-11-06.md`; 32px result `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_annular_coadd_covariance_full783_32px_2026-11-06_01/`; support audit `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_annular_covariance_full783_32px_support_audit_2026-11-06_02/`.

### Approved random-point variance export complete

The exact approved plan `eeef98ef2210337f714e36d8a4afcd352fe312a536c90e8e236ff9f843f53e52` completed all 16 HSC CAS batches: 40,048 rows (11,138,863 bytes) at 1,587/1,598 centers. The independent audit verifies 787/798 matched pairs have valid random-point `i_pix_variance` at both ends; per-center point counts have q16/median/q84 11/22/39. This is sparse random-point metadata, not a coadd variance map and not a complete noise correction. No image, shape, photo-z, or object-catalog data were returned; `science_measurement=false`. Assessment `icl_random_pixel_variance_export_assessment_2026-11-06.md`; audit output `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_random_pixel_variance_audit_2026-11-06_01/`.

### Current consolidated report after random-point variance export

The latest evidence synthesis is `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_23/icl_wl_current_evidence_report.md` (SHA-256 `14ebf8e0965e12f9dfc5aabad4bc359d8e670428de2919c643ca3830f2d1936d`), with summary SHA-256 `5fe134a256ea35e4a7780a796f9b6c38d6cc460432dcddad956c314b318a5cfc`. Builder `build_icl_wl_current_evidence_report.py` SHA-256 `adc7c9d3e86823720643644ea23fe332de32e2602062a3dd94aaa4d3ec08b87d`. It pins the 783-pair 64/32-pixel covariance audits and the completed random-point variance export/audit. The synthesis remains exploratory; `science_measurement=false`.

### Random-point variance sensitivity of the ICL–satellite relation

`analyze_icl_random_point_variance_satellite_relation.py` joins audited local HSC random-point medians to the fixed 94-pair ICL/satellite sample (26 `i`, 68 `i2`). It adds target/control log-ratio or pair-mean terms for pixel variance, sky RMS, input count, and local PSF trace one at a time to the redshift/richness partial-correlation model. A shared 5,000-draw connected-5° max-T family spans all 16 filter/metric adjustments. In `i2`, adjusted r ranges 0.260–0.295 around the baseline 0.293; all adjusted-change simultaneous intervals include zero and all max-T p-values are ≥0.951. This does not establish absence of noise or selection effects: random-point measurement uncertainty is not propagated and the profile sample is selected. No correction; `science_measurement=false`. Assessment `icl_random_point_variance_satellite_relation_assessment_2026-11-06.md`; result `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_random_point_variance_satellite_relation_2026-11-06_01/`.

### S20A/S23B candidate-coverage SQL field correction

An offline review found both earlier expansion queries referenced `mosaic.band`, while the HSC PDR3 mosaic example uses `mosaic.filter01`. The generators now select `m.filter01 AS band` and filter on `m.filter01 = 'HSC-I'`. This is the likely, not server-confirmed, cause of earlier HTTP 406 errors because the response body was not retained. The corrected 2,596-candidate S20A plan has 23,364 probes in 52 batches (plan SHA-256 `96ad3b9078beff082d3dce4de1f9bb7608791915ffc7a46a656d9e6cf30e1149`). The corrected S23B plan has 1,333 candidates, 11,997 probes, 27 batches, and at most 47,988 metadata rows (plan SHA-256 `ee6d2b11893dab0e9ceb0c0d2a3c80c2dee27b600683839d2be294441578ef05`). The offline submitter scope checks passed for both; neither plan was submitted. The S23B plan requires approval of its exact SHA. Both plans request only 14 patch-footprint/depth/seeing metadata columns and exclude images, source rows, shapes, and photo-z. See `icl_expansion_coverage_filter_column_assessment_2026-11-06.md`; `science_measurement=false`.

## 2026-11-06 20/20 forced-profile comparison and support selection

`analyze_icl_satellite_forced_coverage_both_sides.py` compares detected-only
satellite profiles with the existing forced-position profiles at 20% coverage
on both targets and rotated controls. The selection table retains all 224
ICL-supported pairs, including 64 without a supported satellite profile. Of
the supported pairs, 50/71 `i` and 110/153 `i2` pass the 20/20 floor. Selected
pairs have higher median 200–300 kpc ICL residual than excluded pairs in both
filters (3.060 vs 1.116 nJy arcsec⁻² in `i`; 4.446 vs 1.853 in `i2`), despite
similar median redshift and richness. This demonstrates outcome-related
selection and prevents interpreting the supported subset as representative.

Paired connected-5° block-bootstrap decompositions of changes in partial `r`
(controlling redshift and log richness) have intervals containing zero in both
filters. In `i2`, newly supported pairs account for most of the point-estimate
shift (−0.164; 95% interval −0.339 to +0.012), while forced measurements on the
common detected subset change `r` by −0.003 (−0.161 to +0.166). In `i`, the
forced-measurement component is −0.175 (−0.451 to +0.113), also unresolved.
This is a post-hoc sensitivity, not a full-sample correction or a detection;
`science_measurement=false` remains. Assessment:
[`icl_satellite_forced_coverage_both_sides_assessment_2026-11-06.md`](icl_satellite_forced_coverage_both_sides_assessment_2026-11-06.md).
Results: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_forced_coverage_both_sides_2026-11-06_03/`.
No HSC query or download was made. The separately approved 67,071-row
frame-metadata export was already complete at plan SHA-256
`49aa4b7601dd42243209d7b3b394a4dd3b4cbf00d87d28cbb9467553c29f4ce3`; no
duplicate submission was sent.

### Latest consolidated report (2026-11-06, version 27)

The current hash-verified ICL/WL report includes the paired 20/20 forced-position
coverage sensitivity and its outcome-related selection audit. It pins 77 inputs,
including the 224-pair selection ledger and six underlying source tables. Report:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_27/icl_wl_current_evidence_report.md`
(SHA-256 `af88c214eadad945c5177e701f4d7496f95222b297ad84e0be2ba91f595f8b0f`);
summary SHA-256 `42780af71c7392f90f80d53057a8efefe5f9062a84e8f08c595dbbc045eb89b5`;
builder SHA-256 `6a253f424e4e3e3688bbc4698f73c2744a41a6c76413f092944b20d66e8f182a`.
No HSC request/download or image-pixel read was made. `science_measurement=false`.

### Forced unmatched-position radial sensitivity (2026-11-06)

`analyze_icl_forced_position_radial_relations.py` relates the saved 200–300
proper-kpc ICL residual to mean light at forced unmatched CAMIRA positions in
radial bins. It uses the post-hoc 20/20% target/control coverage profiles,
equal position weighting within cluster/bin, redshift and log-richness
adjustment, connected 5-degree blocks, and shared restricted Webb wild-cluster
bootstrap-t weights with Holm and max-T correction. None of nine supported
filter/radius tests has max-T p<0.05. The i 50–75 kpc nominal p=0.030 becomes
max-T p=0.059 and Holm p=0.274. These positions lack unique SEP object
identities, so the response is unmatched-position light rather than confirmed
satellite light. Post-hoc outcome-related selection, limited inner-bin support,
and unpropagated measurement/noise/PSF/sky uncertainty keep
`science_measurement=false`.

Assessment: `icl_forced_position_radial_relations_assessment_2026-11-06.md`.
The integrated report builder hash-validates the result summary, three result
tables, and five source artifacts before including this analysis.

Integrated report version 30 contains 91 hash-verified inputs:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_30/icl_wl_current_evidence_report.md`
(report SHA-256 `62a41b863e1ee2d2150c8647829626972fb3a99d5b123cd993289cc130a65515`;
summary SHA-256 `8d166119df6d21a1d51d3f26a9514e38c6bda8aee727b33d7b4edb1d6e1d2b63`;
builder SHA-256 `0028c1420146454be8b8ce765e957ea133bf41b6e47695c882d7a23330cec713`).
All input/report/builder hashes were independently verified. The report includes
the S23B HTTP 401 submission audit; no CAS job or metadata rows were returned.
`science_measurement=false`.

### S23B footprint metadata screen completed (2026-11-06)

The approved 1,333-candidate, 11,997-probe HSC-I metadata query completed
after fixing the CAS endpoint and replacing unsupported `public.skymap` patch
geometry with `mosaic.areacube && radec_to_coord(...)`. All 27 batches returned
13,575 hash-verified rows. Of 1,333 candidates, 1,043 have metadata matches at
the center and all eight 1-Mpc boundary points; all 1,043 have positive
image/background-product availability flags at each point. This is only a
footprint and QA-metadata screen: actual WCS/mask coverage, image retrieval,
and matched controls remain unverified. Assessment:
`icl_s23b_coverage_metadata_assessment_2026-11-06.md`.

Integrated report version 31 contains 124 hash-verified inputs:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_31/icl_wl_current_evidence_report.md`
(report SHA-256 `3d6aebd01cb2727ca21d2617ef2de1fd17a85c526199315510710f167a1acf15`;
summary SHA-256 `2a0b4edfa9a6536e1f19491365cc5c397a49106b6e1a81d547895b9d9a2f5f9e`;
builder SHA-256 `c7fd66ed582e948b002528bd0b72661608dc5cc2ecf56b209219cf53656ba47f`).
All report/input/builder hashes were independently verified. No images,
object/source rows, shapes, or photo-z were returned; `science_measurement=false`.

### Spatial support overlap for the 20/20 forced-profile selection

The block audit finds selected and excluded pairs in 9/18 connected paired-sky
blocks in `i` and 21/40 in `i2`. Within mixed blocks, equal-block selected-minus-
excluded mean ICL residuals are +1.80 nJy arcsec⁻² (`i`, 95% block interval
−3.26 to +6.66) and +3.73 (`i2`, +0.55 to +6.93). The i2 contrast remains
positive within the overlapping spatial support; i is unresolved. This is a
descriptive block bootstrap, unadjusted for redshift/richness and unable to
transport to single-status blocks. The 20% cut is post-hoc; no correction or
science promotion follows. Assessment:
[`icl_forced_coverage_spatial_selection_assessment_2026-11-06.md`](icl_forced_coverage_spatial_selection_assessment_2026-11-06.md).
Runner `audit_icl_forced_coverage_selection_blocks.py`; output
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_forced_coverage_selection_blocks_2026-11-06_02/`
(summary SHA-256 `3eb63f02507fa83cf1eea377ebc3574bd2b66e86f204bb76face995dcf5f5742`).

### Covariate-adjusted ICL difference by 20/20 forced-profile support

Within mixed connected 5-degree blocks, a pair-weighted fixed-effect model
controls redshift and log richness. The selected-minus-excluded ICL residual
is +1.47 nJy arcsec⁻² in `i` (9 blocks; 95% block-bootstrap interval −2.10 to
+4.79) and +2.13 in `i2` (21 blocks; −0.64 to +5.91). Both intervals include
zero; adjustment attenuates the unadjusted i2 contrast. This is not a correction
or a population estimate because support is limited to mixed blocks, the
threshold is post-hoc, and ICL measurement uncertainty is omitted. Assessment:
[`icl_forced_coverage_adjusted_selection_assessment_2026-11-06.md`](icl_forced_coverage_adjusted_selection_assessment_2026-11-06.md).
Runner `estimate_icl_forced_coverage_selection_adjusted.py`; output
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_forced_coverage_selection_adjusted_2026-11-06_01/`.

The latest 82-input integrated report is version 29 at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_29/icl_wl_current_evidence_report.md`
(report SHA-256 `8350d58462b4eabd228920d970de9852d5b66272e58b629a52577a280dff8d10`;
summary SHA-256 `a3bc38d6191c84c99d58d9391285903670e7b5723674ddfc9627fa74287031a1`).

## S23B 대조 후보 메타데이터 계획 (2026-11-06)

9점 footprint/QA를 통과한 S23B 표적 1,043곳에 대해 HSC PDR3 random table에서
대조 후보 위치를 최대 200개씩 찾는 오프라인 계획을 만들었다. 30분각 탐색,
표적당 최대 200행, 총 상한 208,600행, 위치/ID와 I-band scalar QA/depth 29열,
105 SQL batch, 추정 80 MB다. 영상, 검출원, shape, photo-z 자료는 요청하지 않는다.
생성기 `prepare_s23b_icl_control_candidate_plan.py`; 계획과 SQL 해시를 오프라인
검증했다. 사용자가 승인한 exact-plan 제출과 행 검증은
`submit_s23b_icl_control_candidates.py`가 담당한다 (`--check`는 네트워크 없이 통과).
계획 SHA-256
`143a1bd0a3a3c341576ce420f77ba0089fbfe80fdf58fbfd7afe30ea8e3f6714`; 파일은
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_control_candidate_plan_2026-11-06_01/`
아래에 있다. 이 계획은 신규 대조 후보의 깊이·시잉·오염 매칭을 위한 것이며,
실제 대조장 검증은 아니다. CAS 메타데이터 조회는 진행 중이고 영상 수신은
별도 승인 대상이다. `science_measurement=false`.

### 2026-11-07 S23B random-position export follow-up

세션 3123이 승인된 105개 조회를 순차 제출한다. 체크포인트 기준 batch
000–079의 160,000행이 반환·검증됐고 batch 080이 제출된 상태다. 자격증명은
저장하지 않으며, 영상·검출원·shape·photo-z 자료는 조회 범위에 없다.
완료 후 로컬 감사와 잠정 대조 위치 선별은
`audit_s23b_icl_random_control_candidates.py`가 수행한다. 이 스크립트는
행과 해시를 다시 검사하고 HSC 위치 QA, depth/seeing/SFD 소광 매칭, S20A와
S23B CAMIRA 중심 주변의 0.5/1/2 proper-Mpc 제외 민감도를 평가한다. 산출되는
대조 위치는 영상 WCS·mask·sky·variance·PSF 검증 전까지 science-ready로
간주하지 않는다. 이 postprocessor는 메타데이터 조회 완료 후 실행한다.
`science_measurement=false`.

### 2026-11-07 S23B target image receipt

사용자가 승인한 표적 이미지 범위는 S23B CAMIRA 1,043곳의 HSC PDR3
`pdr3_wide` HSC-I-family `coadd/bg` IMAGE/MASK/VARIANCE plane, 1,000 proper-kpc
half-width다. 예상 58,311,064,944 bytes (54.31 GiB), 30% disk margin
75,804,384,428 bytes. 준비 스크립트
`prepare_s23b_icl_target_image_manifest.py`가 만든 manifest SHA-256은
`e487582ba40cd3e58da67704242db4d7219e381fb08f4b33fed47ce4ff268b44`.
기존 `download_pilot.py`를 이용하며, 확장 sample을 위해 frame 상한을 5,000으로
올리고 파일별 FITS/HDU/WCS 검사와 checkpoint를 유지했다. 세션 5633이 수신 중이며
최신 기록은 121/1,043 FITS (6,826,320,000 bytes) 검증 완료다. 표적 중심은 CAMIRA
centroid이며 확인된 BCG 좌표가 아니다. 이 manifest에는 대조장이 없다.
`science_measurement=false`.
대조장 후보가 metadata audit에서 확정되면 `prepare_s23b_icl_control_image_manifest.py`
가 동일한 HSC-I coadd/bg image/mask/variance 범위의 정확한 수량·용량·해시를
로컬에서 준비한다. 해당 control cutout 세트는 별도 범위로 기록되며, 현재
진행 중인 target download와 섞지 않는다.
전체 표적 수신 완료 후 `audit_s23b_icl_target_images.py`가 모든 FITS hash,
IMAGE/MASK/VARIANCE plane, filter, WCS 정렬과 cutout footprint를 독립적으로
재검증한 다음에 스태킹 단계로 넘긴다.

## 2026-11-06 radial ICL contrast against the outer profile

The saved six-annulus 50–500 proper-kpc stack was reconstructed for the
complete-case 2.5-Kron/40%-coverage sample. The 200–300 kpc target-minus-control
surface brightness exceeds the area-weighted 300–500 kpc profile by 2.780
nJy arcsec⁻² in raw `i` (20 pairs/10 connected blocks; 95% block interval
[0.825, 5.115]) and 2.494 after outer-plane correction ([0.640, 4.774]). In
`i2`, the contrasts are 2.359 (48 pairs/30 blocks; [−0.110, 4.219]) and
2.159 ([−0.296, 4.047]). Both `i2` intervals include zero. The outer profile
can contain real cluster light, so the contrast measures a radial difference
and does not estimate absolute ICL. The selected sample and 40% coverage floor
are exploratory. Source-Poisson, exact coadd covariance, upstream sky transfer,
and extended PSF wings remain unresolved. Assessment:
[`icl_radial_outer_baseline_contrast_assessment_2026-11-06.md`](icl_radial_outer_baseline_contrast_assessment_2026-11-06.md).
Runner [`analyze_icl_radial_outer_baseline_contrast.py`](analyze_icl_radial_outer_baseline_contrast.py).
The result is included in the 933-input integrated report at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_32/`.
Keep `science_measurement=false`.


### S23B paired ICL radial contrasts (2026-11-07)

Run the connected-sky-block paired inner-minus-outer contrast on the already
saved S23B profiles and spatial audit:

```bash
python3 analysis/icl_wgl_membership/hsc_debug/analyze_s23b_icl_relative_radial_contrasts.py \
  --replicas 5000 \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_relative_radial_contrasts_2026-11-07_02
```

The runner compares 50–300, 100–300, and 200–300 proper-kpc target/control
ΔSB against the 300–800 kpc contrast, using 2.5/5/10-degree block widths and
two grid phases. It consumes saved tables/manifests only, performs no HSC
requests or image-pixel reads, and refuses to overwrite an existing output.
The outer annulus may contain real cluster light; results are relative radial
contrasts, not absolute ICL or blank-sky measurements. Keep
`science_measurement=false`; the member catalog and satellite association
remain pending separate download approval.


### S23B matched-control image transfer completed (2026-11-07)

The approved v2 manifest is pinned to SHA-256
`ae9eb406cfc179d9980031bb67f85f8ee71c4580ba59e97713278823db5c06df`. Transfer
and offline audit are complete: 1,041/1,041 FITS files, 58,271,158,080 bytes,
with filter counts `i=248`, `i2=793`. Audit summary:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_control_image_audit_2026-11-07_01/summary.json`.
The audit verifies file checksums, IMAGE/MASK/VARIANCE layout, HSC-I filter,
WCS alignment, and unique target mapping; it does not validate usable pixels
or residual sky. The target-image audit is complete for 1,043 frames. Automated
mask support succeeded for 1,039 targets; four failed with too few valid
pixels. The mask-only common-support audit retained 1,022 same-filter pairs
and found no pair with zero common support inside 1,000 kpc. An exploratory
2D stack exists for these pairs, while `images_ready_for_science=false` and
`science_measurement=false` remain in force. Sky, extinction, PSF wings,
coadd covariance, and CAMIRA member completeness remain open gates.


### S23B richness relation for paired ICL contrasts (2026-11-07)

Run the pinned no-download, no-image-pixel analysis:

```bash
python3 analysis/icl_wgl_membership/hsc_debug/analyze_s23b_icl_richness_relation.py \
  --replicas 5000 \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_richness_relation_2026-11-07_01
```

The model regresses each pair's inner-minus-300–800 kpc ICL contrast on
log10 CAMIRA richness, redshift, and target-control SFD E(B−V) mismatch, then
resamples connected sky blocks under six RA-Dec layouts. The E(B−V) term is a
covariate only; it does not correct flux. The primary 5° intervals include
zero; two i2 slopes are positive only for one of the six layouts. Treat this
as a layout-sensitive richness characterization, not a satellite relation or
a non-detection. CAMIRA richness does not replace member positions. Keep
`science_measurement=false`.


### S23B richness common-support sensitivity (2026-11-07)

Run the pinned sensitivity comparing the primary richness model with a model
that additionally includes the log count of common target/control mask pixels:

```bash
python3 analysis/icl_wgl_membership/hsc_debug/analyze_s23b_icl_richness_coverage_sensitivity.py \
  --replicas 5000 \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_richness_coverage_sensitivity_2026-11-07_01
```

The two models use the same connected-block bootstrap draws; report the paired
slope change as a sensitivity, not a causal correction. The primary 50–300 kpc
slope changes by −0.16 [−5.15, 6.15] in `i` and −0.12 [−1.97, 1.03] in `i2`;
all change intervals include zero. Keep `science_measurement=false`.


### S23B Galactic-extinction sensitivity (2026-11-07)

Run the pinned sensitivity on the completed S23B target/control radial table:

```bash
python3 analysis/icl_wgl_membership/hsc_debug/analyze_s23b_icl_extinction_sensitivity.py \
  --replicas 5000 \
  --seed 20261108 \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_extinction_sensitivity_2026-11-07_04
```

The runner checks the target/control/stack/profile/spatial-block hashes and
reproduces the existing raw contrasts before applying separate `i`/`i2` SFD
factors to target and matched-control means. Its output directory already
exists, so the runner refuses to overwrite it; use a fresh output suffix for
another run. For the 50–300 minus 300–800 kpc contrast, the nominal correction
is +0.369 [0.190, 0.516] nJy arcsec⁻² in `i` and +0.444 [0.337, 0.569] in
`i2` at 5° phase zero. All six spatial-layout intervals for this correction
remain positive. The 300–800 kpc reference may contain cluster light, and
extinction multiplication also rescales residual sky. This remains an
exploratory sensitivity; keep `science_measurement=false`.
## S23B CAMIRA member intake and direct ICL–satellite candidate relation (2026-11-07)

The user-approved one-file member catalog has now been received and audited;
this supersedes earlier README statements that its download was pending.
Only camira_s23b_wide_sm_member.dat.gz was downloaded (10,626,590 bytes;
SHA-256 66529e02c6ae11ebdccc3fb1d4756b639364f372f9d0f09a721786a19906945a).
The importer verifies the gzip stream, joins 12,162 cluster groups and all
1,043 targets, and WCS-checks 45,634 target member positions. There are
39,166 in-footprint positions (85.83%); footprint inclusion does not imply
that the pixel is usable or unmasked. Columns 7–8 remain uninterpreted.

Offline runners:

    python3 analysis/icl_wgl_membership/hsc_debug/import_s23b_camira_member_catalog.py
    python3 analysis/icl_wgl_membership/hsc_debug/analyze_s23b_icl_satellite_relation.py

The association runner reads saved profiles and catalog positions only; it
does not contact HSC or read science-image pixels. It models the common-support
50–300 minus 300–800 kpc relative ICL contrast against unweighted candidate
member surface-density contrast, with richness, redshift, target/control
SFD mismatch, and common-pixel support as covariates. It uses 5,000 connected
RA-Dec block bootstrap resamples in six spatial layouts.

At 5° phase zero, the adjusted slopes are 5.72 [−0.06, 13.63] for HSC-i
(232 pairs/26 blocks) and 4.24 [0.99, 7.80] for HSC-i2 (756 pairs/42 blocks),
in nJy arcsec⁻² per dex. Intervals exclude zero in only 3/6 and 4/6 layouts.
With a 10% support floor in both annuli, the HSC-i2 interval also crosses
zero. This remains exploratory, not a detection; science_measurement=false.

Assessment:
s23b_icl_satellite_relation_assessment_2026-11-07.md

Data outputs (kept outside Git):
/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_catalog_2026-11-07_01/import_2026-11-07_02/
/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_satellite_relation_2026-11-07_03/

Keep science_measurement=false. Remaining gates include valid-pixel/mask
review, calibrated member purity/completeness, residual sky, extinction,
extended PSF wings, coadd covariance, and center/mask sensitivity.

Rerun note: the default output directories above already contain these completed
runs, and both runners refuse to overwrite existing outputs. For a repeat run,
choose a fresh output path; also provide a fresh plan path to the importer.

### Source-aware member-centered outer-light stack and ICL relation (2026-11-07)

The approved member catalog and existing HSC target images have now been used
for an exploratory source-aware member profile stack. No new download or
network access was used. Run the full profile stack with
`measure_s23b_camira_member_outer_profiles.py`, then fit the spatial-block
adjusted ICL relation with `analyze_s23b_member_outer_icl_relation.py`.

The stack includes 979 eligible target/filter frames (226 `i`, 753 `i2`),
16,503 candidate positions at 50–500 proper kpc, and 13,269 retained
candidate-by-bin measurements. Each uniquely SEP-matched member's segment is
preserved. Other detections and bad/bright pixels are masked, and unmatched
member positions receive fallback masks. Target/control annuli require 60%
coverage and at least three same-radius rotated controls.

At 50–75 proper kpc, the median candidate-minus-control surface brightness
is 0.68 nJy arcsec⁻² in `i` (q16–q84 −3.77 to 5.64; 62 clusters) and 1.46 in
`i2` (−4.23 to 8.94; 365 clusters). The adjusted outer-light versus relative
ICL slope is 0.086 [−0.098, 0.334] for `i` and −0.041 [−0.145, 0.090] for
`i2` at the 5° phase-zero layout. All six spatial-layout intervals cross
zero in both filters. The 5%, 10%, and 20% common-support sensitivity fits
also have no available interval excluding zero. This is an unresolved
exploratory relation; keep `science_measurement=false` and
`images_ready_for_science=false`.

Full assessment:
`s23b_camira_member_outer_profile_relation_assessment_2026-11-07.md`.

Outputs (outside Git):
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_outer_profiles_full_2026-11-07_01/`
and
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_outer_icl_relation_2026-11-07_01/`.

### Raw member-field-8 weighting sensitivity for outer satellite profiles

`analyze_s23b_member_outer_raw_field8_sensitivity.py` compares the current
unweighted candidate profiles with a profile weighted by the member file's
eighth raw field. The headerless field's file-specific meaning remains
unverified, so do not label it `w_mem` or membership probability. The analysis
reproduced the unweighted relation inputs, used paired 5,000-draw connected
sky-block resamples, and read no image pixels or network resources.

At 50–75 proper kpc, weighted candidate-minus-control medians are 1.03
nJy arcsec⁻² in `i` (62 clusters; q16–q84 −3.87 to 6.36) and 1.99 in `i2`
(365 clusters; −4.43 to 10.89). The adjusted ICL relation slopes in the 5°
phase-zero layout are 0.017 [−0.247, 0.334] for `i` and −0.013
[−0.118, 0.129] for `i2`; all weighted-slope intervals cross zero in all six
layouts. The paired `i2` weighting shift excludes zero in only 3/6 layouts.
This is a sensitivity check, not a detection. Keep
`science_measurement=false` and `images_ready_for_science=false`.

Assessment:
`s23b_camira_member_outer_raw_field8_sensitivity_assessment_2026-11-07.md`.

Output (outside Git):
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_outer_field8_sensitivity_2026-11-07_01/`
(summary SHA-256
`7156397a478932a9d1fb133b88a758a814bfc0543708719cd36f2a21c1b14875`).

### S23B raw member-field-8 weighting sensitivity (2026-11-07)

Run `analyze_s23b_icl_satellite_field8_weight_sensitivity.py` for the
headerless eighth member-file field sensitivity. The field is numerically
probability-like but its schema is unverified; do not label it w_mem. At
5-degree phase zero the adjusted slope is 5.13 [−0.80,13.05] for HSC-i and
5.01 [2.44,7.62] nJy arcsec⁻² dex⁻¹ for HSC-i2. HSC-i2 intervals exclude
zero across all six spatial layouts, but both filters cross zero at the 10%
common-pixel support floor. Treat this as exploratory and keep
science_measurement=false.

Assessment:
`s23b_icl_satellite_field8_weight_sensitivity_assessment_2026-11-07.md`

Outputs:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_satellite_field8_weight_sensitivity_2026-11-07_01/`
(summary SHA-256
`ec2f6a325ac0a6e349b6f0f52cd106c544d3d811594c35568f5712e3273b2500`).

### S23B CAMIRA member-position mask audit (2026-11-07)

Run `audit_s23b_camira_member_mask_overlap.py` to check imported member
coordinates against existing target mask planes. It verifies 1,039 mask
checksums and mask WCS coordinates, and reads no science-image pixels. The
combined Kron masks flag 98.4–99.1% of in-footprint member centers because
they include detected sources; no segmentation labels were saved, so the
candidate's own source cannot be restored separately. These products are not
ready for member-centered satellite profiles. Four mask products are missing
because the preparation step found too few valid pixels.

Assessment:
`s23b_icl_satellite_member_mask_audit_assessment_2026-11-07.md`

Output:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_satellite_member_mask_audit_2026-11-07_02/`
(summary SHA-256
`98f4d4e5abe043457a19cba229c2b7c0a7f31e3e07ceeb25d7b5e84d1d713710`).

### Existing WL to S23B CAMIRA positional-overlap audit (2026-11-07)

Run `audit_wl_s23b_camira_crossmatch.py` to test whether the current 800-lens
S19A weak-lensing sample overlaps the 1,043-target S23B CAMIRA member-profile
sample. The runner reconstructs the 435-lens common WL support and verifies
its saved ID hash before cross-matching sky coordinates. It performs no CAS
query, source-row transfer, download, or image-pixel read.

There are no matches within 60 arcsec, including among the 427 S23B targets
with member outer-profile rows. The nearest S23B center to any WL lens is
364.846 arcsec away. Numeric ID equality occurs for 66 pairs (29 among the
profile targets), but none is positional; do not join on ID alone. The current
WL and S23B profile products cannot support a same-cluster three-way
weak-lensing–ICL–satellite relation.

Assessment:
`wl_s23b_camira_crossmatch_assessment_2026-11-07.md`.

Output (outside Git):
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_camira_overlap_audit_2026-11-07_01/`.

Next gate: a local count-only HSC S19A plan is ready for the 427 S23B
member-profile targets: 18 batches, 3 proper-Mpc cones, 427 integer counts,
122,198 SQL bytes. It has not been submitted. Plan SHA-256:
`b713eebded489c280fd9ce1b5840bc8d7d3e42d029035f8dbd414b72fb31e984`.
Submission requires separate explicit approval. After counts, validate support
and matched random centers before planning any shape-row export. Keep
`science_measurement=false`.
The plan is verified by `submit_s23b_wl_shape_counts.py --plan`; its submit
mode requires an interactive `lageunha` terminal and the exact approved plan
hash. No query was submitted.

### S23B member-field schema consistency audit (2026-11-07)

`audit_s23b_camira_member_raw_field_schema.py` joined all 637,146 approved
member rows to the 12,162-cluster catalog. Raw field 8 is wholly in
[0.001000069, 1.0] and every row passes the published `w_mem >= 10^-3`
threshold. This strongly suggests it is `w_mem`, but the public README does
not specify the headerless member-file layout. Keep `field8_raw` as the stored
name until confirmed from file-specific documentation. The previous field-8
profile weighting remains sensitivity-only; `science_measurement=false`.

Assessment: `s23b_camira_member_raw_field_schema_assessment_2026-11-07.md`.
Output (outside Git):
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_schema_audit_2026-11-07_04/`.

### S23B field-8 weighted satellite-density sensitivity

`analyze_s23b_icl_satellite_field8_weighted_density.py` repeats the direct
ICL–candidate-density model using the approved member file's raw eighth field
as a sensitivity weight. Literature supports, but does not prove from the
file-specific column order, that this is `w_mem`; keep it named `field8_raw`.
The weighted and unweighted fits use paired connected sky-block resamples.
The i2 primary slope changes from 4.24 [0.99, 7.80] to 5.01 [2.44, 7.62],
but their paired difference is +0.77 [−1.05, 2.46]. Under the 10% common-pixel
floor the primary i2 interval includes zero, and the i result remains
unresolved. This is a sensitivity analysis, not a detection;
`science_measurement=false`.

Assessment: `s23b_icl_satellite_field8_weighted_density_assessment_2026-11-07.md`.
Output (outside Git):
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_satellite_field8_weighted_density_2026-11-07_02/`.

## 2026-11-07 S23B relative-contrast covariance propagation

`propagate_s23b_icl_relative_contrast_covariance.py` replays the existing
1,022-pair stack's bilinear warp, saved K=2.5 masks, pair-common support, and
three signed inner-minus-outer annular contrasts. It contracts the transpose
warp weights with delivered VARIANCE and the measured per-frame 9×9 residual
lag kernel. All six propagated pair counts match the saved 5-degree block
support; the finite-lag conditional pixel-noise σ is 1.089–1.098 times the
delivered-VARIANCE diagonal σ. This is not a total uncertainty or an error-bar
correction: the outer-field kernels are assumed stationary in cluster annuli,
and long-range sky, source Poisson, PSF, masking, and other systematics are
not included. Keep `science_measurement=false` and
`images_ready_for_science=false`.

- Assessment: `s23b_icl_relative_covariance_propagation_assessment_2026-11-07.md`
- Results (outside Git): `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_relative_covariance_propagation_2026-11-07_02/`
- Summary SHA-256: `9c1ef1740490f7cbc4811a9bd2573d69405f5482d7ce754a09642108507fc94e`
- No data request, commit, or push was made for this propagation.

## 2026-11-07 CAMIRA member profile injection recovery

`run_s23b_member_outer_profile_injection_recovery.py` replays the source-aware
50–75 kpc CAMIRA member profile for one seeded-random, already-profiled member
in each of 427 target/filter frames, then reruns detection, matching, masking,
and rotated-control measurement after injecting an exponential halo at three
amplitudes. All baseline replays match exactly. Of 1,281 injected frame rows,
1,221 retain adequate support. Median recovered/input ratios are 0.838–0.896
for `i` and 0.807–0.819 for `i2`; paired amplitude differences have
cluster-bootstrap intervals including zero. Detection/segmentation changes
are frequent, so the measurement is conditional on prior detection and mask
support. The injected model has no PSF convolution or source Poisson noise,
and only target images are injected. This is exploratory transfer calibration,
not a completeness estimate or evidence for an ICL–satellite association.
Keep `science_measurement=false` and `images_ready_for_science=false`.

- Assessment: `s23b_camira_member_outer_injection_recovery_assessment_2026-11-07.md`
- Results: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_outer_injection_recovery_2026-11-07_01/`
- Summary SHA-256: `e9d4cced59f0b0f8ea3f189989886f28af74a339cc5a9b3421f8bca187defc77`

## S23B member-profile PSF request plan (2026-11-07)

The local planner prepared an unsubmitted HSC PDR3 coadd model-PSF request for
the same 427 CAMIRA member coordinates used in the exploratory injection
recovery (62 `i`, 365 `i2`). Existing audited PSF banks do not cover these
positions exactly; the nearest old PSF is 69.6869 arcsec away. The review cap
is 8,608,320 bytes, estimated from the largest PSF FITS file already audited;
stop if the received data would exceed that cap. No images, shape rows, or
photo-z are requested. No network request was made and
`science_measurement=false` remains in force.

- Planner: `prepare_s23b_member_profile_psf_plan.py`
- Assessment: `s23b_member_profile_psf_request_assessment_2026-11-07.md`
- Prepared plan (outside Git): `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_profile_psf_plan_2026-11-07_03/`
- Plan SHA-256: `03f974ad24ead32eabd103ef7a6e0d6b52317aa2647023a5f7402003a2e3b25d`
- Coordinate list SHA-256: `7327d411bacab15954538c64cb44deae9141955377325bc459a1e6a162db7b31`

The existing batch auditor and member-profile injection runner now accept this
request schema. Their PSF path pins the exact approved plan hash, checks the
audited file mapping, and retains signed model-PSF pixels. No PSF retrieval or
PSF-convolved run has been made.

## Existing S19A shape-row reuse for S23B targets (2026-11-07)

The 800-lens S19A source-row export was verified locally: 3,069,544 pair rows,
2,929,303 unique source IDs, consistent repeated source fields, and all rows
inside their original 3 proper-Mpc query cones. Its old query cones intersect
43 of the 427 S23B member-profile apertures, while source rows are present in
only 21. No target aperture is fully covered; estimated union coverage has
median 0 and maximum 60.4%. Partial reuse cannot replace S23B-specific support
counts or source rows. The separate 427-center count-only plan is prepared but
not submitted and requires its own plan-specific approval.

- Assessment: wl_s23b_existing_export_reuse_assessment_2026-11-07.md
- Runner: audit_wl_s23b_reuse_existing_s19a_shape_rows.py
- Summary: /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_existing_export_overlap_2026-11-07_01/summary.json
- Summary SHA-256: 61026a5d82d47303724f576922415c451589f2516f39337adee3abe021c2c54a
- Network use: none; science_measurement=false

## 2026-11-07 Integrated S23B evidence report

`build_s23b_icl_wl_evidence_report.py` checks the saved ICL, member-profile,
Gaia PSF, and S19A/S23B WL summaries and writes a Markdown assessment plus a
hash-pinned JSON provenance record. It reads summaries and plans only; it does
not use the network or read image pixels. The current report is at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_wl_evidence_2026-11-07_01/`.

The S23B shape-count plan remains unsubmitted pending its own approval:
`b713eebded489c280fd9ce1b5840bc8d7d3e42d029035f8dbd414b72fb31e984`.

The follow-up `audit_wl_s23b_reuse_existing_s19a_sources.py` also verified
the already approved 3,100-center random source-row export (9,010,324 rows,
124 batches) and combined it with the 800-lens source export. Source fields
agree for duplicate IDs. Rows from the union occur in 32/427 S23B apertures;
395 have none. Combined geometric query-cone coverage has median 0, p90
0.310, maximum 0.994; 2 apertures reach 95%, and no single old cone contains
a full target aperture. This still does not provide common target support.
The area estimate is an upper bound, and counts are diagnostics only; no n(z),
lensing efficiency, or B-mode polarity calibration was applied.

- Combined audit results: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_existing_source_union_2026-11-07_01/`
- Summary SHA-256: `876bed400f93a55c391c89128aa9e9214ed489e664ea0e618ea3df6920f91ed7`
- Per-target table SHA-256: `c6cc5467d031e76d6f867e222f6a1a40b0b76e86b36129f5fd7901853261cede`
- `science_measurement=false`; no network query or transfer was made.

## 2026-11-07 S23B satellite wild-cluster check

The new local runner `analyze_s23b_icl_satellite_relation_wild.py` tested the
adjusted satellite-density coefficient with 10,000 Webb six-point restricted
wild-cluster bootstrap-t draws. A common 5-degree connected target/control
partition supplied 49 blocks across both filters and all support floors. All
eight filter/support groups were eligible. Three raw p-values were below 0.05;
one `i2` >=20%-support result remained below 0.05 after Holm (slope 4.21
nJy arcsec^-2 dex^-1, raw p=0.0057, Holm p=0.0456), while no result passed the
complete-null max-T diagnostic (minimum p=0.2666). This is not robust detection:
support floors reuse the same images, max-T is calibrated under the complete
null only, and sky/PSF/mask/covariance uncertainties remain. Keep
`science_measurement=false`.

- Assessment: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_satellite_relation_wild_2026-11-07_01/assessment.md`
- Summary SHA-256: `d95f894c6cf54e92d906063fdc920c24e80bdf9ad4194d2ff63dacd8bdd42494`
- Input pair table SHA-256: `98ed9e0022bfa68d5921ead62dce55be7df1cfa169e6a777d66d1e9760f43994`

## 2026-11-09 S23B count-only query and audit

The approved plan SHA-256
`b713eebded489c280fd9ce1b5840bc8d7d3e42d029035f8dbd414b72fb31e984` completed
all 18 batches. The resumed runner retrieved job `1070884`, completed jobs
`1071485` and `1071508`, and verified all 427 integer counts. The local audit
checked the receipts and hashes, SQL center IDs, and mapping to the saved
target manifest. No source rows, coordinate rows, shape values, photo-z, or
images were returned, and credentials were not saved.

The 3 proper-Mpc cones contain 356,191 summed lens-source candidate pairs.
Counts are zero for 329/427 centers (77%); 82 centers have >=1,000 entries,
62 have >=3,000, and 25 have >=5,000. Median count is zero and maximum is
8,959. Descriptive Spearman rho is 0.060 versus redshift and -0.145 versus
CAMIRA richness; no p-values were computed. A local sky-region cross-tab shows
substantial footprint variation in zero counts. Before source-row planning,
define and justify a common-support lens/control sample.

The released S19A weak-lensing table already applies its release weak-lensing
flag and star mask. Project-specific source-redshift, tomographic, B-mode, and
radial cuts are not applied. The totals are availability counts and candidate
pairs, not unique sources, final effective source counts, a lensing signal-to-
noise estimate, or a weak-lensing measurement. Keep `science_measurement=false`.

- Audit output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_shape_count_audit_2026-11-07_01/`
- Assessment SHA-256: `f2ce0c5f22c2f6f894386f0c682d08aaf21236d5dd56188cf896b65064e6d672`
- Per-center table SHA-256: `d650564b4870d4af93ce04b9ebd82475a99366aa68c2dfe236e7f070826c0f7f`
- Auditor: `audit_s23b_wl_shape_counts.py`

The local common-support and row-plan stage was completed on 2026-11-10. The
400-control count-only query and any later source-row submission remain
separate approval gates.

## 2026-11-10 S23B weak-lensing common-support plan

The approved 427-center, 3-proper-Mpc count audit has 329 zero-count centers;
82 centers have at least 1,000 entries in the released S19A shape subset before
project-specific shape, source-redshift, B-mode, or radial cuts. These are raw
availability counts, not final effective source counts or a measurement.

The provisional sample selects 80 of those 82 lenses and assigns five globally
unique HSC random-table controls per lens (400 controls). The other two lenses
do not have five candidates in the saved pool that pass all screens. Controls
pass the saved center flags and depth/seeing/SFD-EBV match; the previous CAMIRA
2-Mpc overlap count is zero; and each is more than 6 proper Mpc from its lens,
so paired 3-Mpc count cones do not overlap. The CAMIRA screen is incomplete and
does not establish a clean large-scale environment.

The 80 selected lenses contain exactly 340,256 verified count-only candidate
pairs, with per-lens min/median/max 1,040/4,323/8,959. The old ICL control
positions are too close to their paired lens for this use: 114/427 are within
3 Mpc and 280/427 within 6 Mpc. Other S23B CAMIRA centers lie within projected
3 Mpc of 110 old controls at any redshift and 23 at |dz|<=0.05. This is a
proximity screen, not a shear-contamination model. The previous S19A
lens/random source-row union has rows in only 29 selected apertures and none
reaches 95% estimated geometric coverage.

The next count-only plan covers 400 control centers in 16 batches of 25 and
returns three columns only: `random_center_id`, `matched_lens_id`, and integer
`n_shape_sources`. It passed the existing random-count plan validator but was
not submitted. Its exact plan-file SHA-256 is
`fe19ac426272d4005e3ee74a902dfdb52d1fe83ecb2cc5ec901c10b880d28e8f`; submit
only after separate explicit approval for that hash.

The prepared source-row plan has 80 lens centers in 4 batches and 400 controls
in 16. The expected lens rows from verified counts are 340,256; exact control
rows require the count-only result. A different completed sample gives only a
planning proxy of 1.05/1.83/2.23 million combined pairs (q16/median/q84) and
0.25/0.55/1.00 GiB uncompressed CSV (low/point/high). The lens and control
schemas carry 22 and 23 columns respectively, including shapes, calibration,
PSF moments, raw B-mode mask, and HSC Y3 source-redshift bin. Full object
photo-z PDFs and images are excluded.

Before a measurement, apply identical source selections to lenses and controls,
verify the raw B-mode polarity, provide calibrated source n(z) and lensing
efficiency, and audit footprint, PSF, control overlap, and spatial covariance.
Keep `science_measurement=false`. Neither count-only nor source-row queries
were submitted in this stage; no network or data transfer was used.

- Planner: `prepare_s23b_wl_common_support_plan.py`
- Assessment and plans (outside Git):
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_wl_common_support_plan_2026-11-10_01/`
- Row-plan file SHA-256: `3ec3c8e77cf9cf38fdd02320718c6e38ee3d23e8948584255cfdc27bfe90cfee`

### 2026-11-10 approved control-count execution and post-count row plan

The approved count-only plan
`fe19ac426272d4005e3ee74a902dfdb52d1fe83ecb2cc5ec901c10b880d28e8f` completed
16/16 batches and verified all 400 control-center rows. It returned only
`random_center_id`, `matched_lens_id`, and integer `n_shape_sources`; no shape
values, source coordinates, photo-z, or images were returned, and credentials
were not saved.

The controls sum to 1,572,453 raw source counts (18 zero and 382 positive).
The 80 selected lenses sum to 340,256, so the future export is expected to
contain 1,912,709 lens/control pairs. The per-lens mean-of-five control count
divided by its lens count has q16/median/q84 0.690/0.936/1.295. Seven lenses
have at least one zero-count control, and none has all five zero. Project source
cuts are not included in these raw counts.

The finalized row plan has 4 lens batches and 16 control batches. Estimated
uncompressed CSV volume is 0.46/0.58/0.86 GiB (low/point/high), from observed
count totals and historical bytes per row. The row query and transfer remain
unsubmitted and need separate approval for exact plan-file SHA-256
`7bfb998149f30a50b5ee7a537bb703ed302b7ea269a3296c35204a836e75aabc`.
`science_measurement=false` remains in force.

- Finalizer: `finalize_s23b_wl_common_support_plan.py`
- Count audit and post-count plans:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_wl_common_support_plan_2026-11-10_01/`

### 2026-11-10 approved S23B source-row export

The user approved the exact post-count row-plan file SHA-256
`7bfb998149f30a50b5ee7a537bb703ed302b7ea269a3296c35204a836e75aabc`.
The plan covers 80 lenses and 400 matched controls in 4+16 serial query
batches, with 1,912,709 expected lens/control source pairs. It returns the
approved 22-column lens and 23-column control schemas; images and full photo-z
PDFs are excluded. The estimated uncompressed CSV volume is 0.46/0.58/0.86
GiB (low/point/high), not an observed transfer size.

`submit_s23b_wl_common_shape_rows.py --check` passed on `lageunha` with no
network access. It validates the immutable plan hash, both count-only sources,
all per-center counts, and every SQL batch. Submission requires the explicit
approval flag and the same exact hash; it uses one hidden password prompt and
records job IDs, row counts, byte counts, and file checksums without saving
credentials. At this checkpoint, no source-row query or transfer has started;
`science_measurement=false`.

### 2026-11-10 approved S23B source-row export complete

The approved plan was submitted from `lageunha`; all 20/20 batches completed.
The four lens batches returned 340,256 lens-source pairs (102,103,795 bytes),
and the 16 control batches returned 1,572,453 control-source pairs
(523,062,067 bytes). The combined export contains 1,912,709 center/source
pairs in 625,165,862 bytes (0.582 GiB). Pair rows are not a count of globally
unique source objects because a source may appear in multiple apertures.

Every batch passed its exact count-only row total, unique center/object-ID
pair, schema, numeric-value, coordinate, and shape-field checks. An independent
post-run audit rechecked all 20 saved file sizes and SHA-256 receipts. No images
or full photo-z PDFs were included, credentials were not saved, and
`science_measurement=false` remains in force. Results and `row_status.json` are
outside Git at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_wl_common_support_plan_2026-11-10_01/source_rows_post_count_01/`.

The approved local source-row quality and common-support audit is complete;
results are recorded below. Calibrated source n(z) uncertainty, PSF, and
spatial covariance remain interpretation gates before a measurement.

### 2026-11-10 approved local S23B source-row QA complete

`audit_s23b_wl_common_shape_rows.py` read only the 20 already-received CSVs and
the previously downloaded HSC Y3 SACC n(z) file. It made no network request,
CAS query, or additional download, and did not read image pixels or compute a
shear profile. It rechecked all batch sizes, SHA-256 receipts, source-pair
uniqueness, approved center IDs, per-center count-only totals, and the 80-lens /
400-control mapping.

The same cuts were applied to both roles: `|e|<2`, `0<=sigma_e<=0.4`,
`mag_a10<=25.5`, `resolution>=0.3`, HSC Y3 source bins 1–4, verified B-mode
mask pass value 1, positive fiducial `n(z)` lensing efficiency, and
0.1–3 proper Mpc. The additional `m>-1` and positive PSF-moment-trace checks
are input-validity gates for effective-support summaries, not HSC shape cuts.

| Stage | Lenses | Five-control centers |
|---|---:|---:|
| Input center/source pairs | 340,256 | 1,572,453 |
| HSC shape + magnitude + resolution + z-bin + B-mode | 240,068 | 1,107,161 |
| Also positive fiducial `n(z)` efficiency | 240,068 | 1,107,161 |
| Also 0.1–3 Mpc and calibration-input valid | 239,565 | 1,105,779 |

The final pair-row retention is 70.4% for lenses and 70.3% for controls.
Those rows are not globally unique sources. Two lens apertures have zero
retained rows; 24/400 control centers have zero, affecting nine paired lens
sets. None of the 80 sets has all five controls at zero. Among lens sets with
nonzero lens support, the q16/median/q84 mean-five-control/lens ratios are
0.690/0.925/1.311 by
retained count, 0.686/0.929/1.311 by HSC weight sum, and 0.690/0.956/1.387 by
the geometry-weight proxy sum.

Median retained `|e|`, `sigma_e`, HSC weight, resolution, and PSF moment trace
are similar between lenses and controls: 0.529/0.532, 0.1056/0.1053,
5.379/5.384, 0.622/0.620, and 0.1127/0.1118, respectively. Median
center-level catalog-weight `N_eff` is 2,964.6 for lenses and 2,840.7 for
controls; geometry-proxy `N_eff` is 1,910.6 and 1,970.4. `N_eff` uses
`(sum w)^2/sum(w^2)`; the geometry proxy uses `weight*(D_l beta)^2` and is not
the calibrated estimator weight.

The B-mode polarity cross-check against the published GAMA09H rectangle
matched exactly: all 48,765 raw-zero pairs were inside RA 132.5–140 degrees,
Dec 1.6–5 degrees, and all 1,863,944 raw-one pairs were outside. There were
zero opposite-value pairs. This confirms raw 0 as excluded for this export,
consistent with the earlier same-release coordinate audit; the official HSC
page itself describes the mask purpose but not the stored-value polarity.

This passes a catalog-row support and distribution QA gate and supports
preparing a calibrated estimator plan. It is not a weak-lensing measurement:
the SACC n(z) is fiducial and its uncertainty is not propagated; the proxy
omits per-object `(1+m)`, responsivity, and selection-bias response; and
control contamination, PSF spatial matching, source-reuse/field covariance,
and injection/recovery remain unresolved. Keep `science_measurement=false`.

- Runner: `audit_s23b_wl_common_shape_rows.py`
- Results: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_wl_common_support_plan_2026-11-10_01/quality_audit_2026-11-10_01/`
- Summary SHA-256: `9ef472c1ae0888e187fce05f3449b1614d1da7c99e7a400403b92f8aa7e2dc54`
- Assessment SHA-256: `4fd4d752f1138810702ec3c1034fa3181f2758505bc9e95deaaed39687d4d753`
- Per-center TSV SHA-256: `e630442cbaeb054763897a46d68d3b4baea7c0884e80ac22f46ece906583f48e`
- Paired-lens TSV SHA-256: `68e2bf3cf4e7dd515400881f8f5a508c66d3cb019fecff32e2b81100de2afe14`
- Radial TSV SHA-256: `6d25f945d2dc567b26d741cba9329b7ddcacd55d60a24341d2f03ca0ce6fb087`
- Next: prepare a pre-specified calibrated, random-subtracted estimator and
  covariance plan, including how to handle zero-support lenses/controls;
  await approval before running a shear profile. Keep `science_measurement=false`.

## 2026-11-11 radial/tomographic support and zero-control sensitivity

The follow-up pass extended the verified 80-lens/400-control source-row QA to
all seven proper-radius bins and HSC Y3 source bins 1–4. It created 13,440
center × radius × tomographic-bin rows and 56 role/radius/bin aggregate cells.
The per-center tomographic counts sum exactly to the earlier radial totals,
and the full radial/tomographic totals sum to the 239,565 lens and 1,105,779
control rows after the common cuts.

The inner 0.1–0.2 Mpc bin is the limiting region. Lens aperture coverage by
z-bin 1/2/3/4 is 65/63/59/52 of 80; control coverage is 285/303/293/262 of
400. In the outer 2–3 Mpc bin, coverage is 78/78/78/78 lenses and
376/376/376/376 controls. All cell counts are center/source pairs, so repeated
sources remain correlated across apertures.

Zero-support sensitivity, among lens sets with nonzero lens rows:

| Minimum positive controls | Eligible lens sets | Positive control centers used | Retained lens pairs | Median control/lens count ratio |
|---:|---:|---:|---:|---:|
| 3/5 | 74/80 (92.5%) | 367 | 233,064 | 0.940 |
| 4/5 | 74/80 (92.5%) | 367 | 233,064 | 0.940 |
| 5/5 | 71/80 (88.8%) | 355 | 224,375 | 0.950 |

The 3/5 and 4/5 rules select the same set. The 5/5 rule removes three more
lens sets; their support-ratio distributions remain similar. These are
selection sensitivities, not a chosen final policy. No shear profile or mass
fit was calculated; `science_measurement=false` remains.

- Runner: `audit_s23b_wl_common_shape_rows.py`
- Results: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_wl_common_support_plan_2026-11-10_01/quality_audit_2026-11-11_01/`
- Summary SHA-256: `eb73aa0365813d880d3edf60b12e4039f2575f8a97eab8868fc918121a34d445`
- Assessment SHA-256: `35369bbc74fa576570860af09d95d29c33c32cf8b1dfb96c30bb5949d5ddcf95`
- Per-center radial/tomographic TSV SHA-256: `f8d604e91c43705cb24152c7e98553efd8838933034a8657bcc603c2f8c11bed`
- Aggregate radial/tomographic TSV SHA-256: `ae84048ee707ffe531ecff80ede270544a7475846079db7520ee59bd62706911`
- Zero-support policy TSV SHA-256: `ba89948681cb0560408030a4df66a232a2989cd5415c34ca1cb99c59d8bc26da`

## 2026-11-11 S23B random-subtracted weak-lensing profile

The approved 80-lens/400-control source rows were passed through the existing
calibrated HSC S19A estimator with fiducial HSC Y3 tomographic `n(z)`. The
primary 4/5-control, bins-2–4 profile retains 74 lens sets and 367 controls.
Tangential `DeltaSigma` declines from `23.845 ± 3.432` in 0.1–0.2 proper Mpc
to `0.762 ± 0.400` in 2–3 proper Mpc, in units of `10^13 Msun Mpc^-2`.
Cross-null behavior changes with the control threshold and jackknife width;
the result is exploratory and `science_measurement=false`.

The calculation used no new query, download, or image pixels. Its pooled
random-center weighting, fiducial `n(z)` uncertainty, survey-region covariance,
random-aperture PSF/depth/mask/sky matching, reduced shear, magnification,
member dilution, and selection-response uncertainty remain open. See the
[full assessment](wl_s23b_random_subtracted_profile_assessment_2026-11-11.md)
and reproducible runner
[`analyze_s23b_wl_random_subtracted_profile.py`](analyze_s23b_wl_random_subtracted_profile.py).
Outputs outside Git:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_random_subtracted_profile_2026-11-11_01/`
(summary SHA-256 `2c899b0090dd94a34f1fd6b0cb111e3d1ada696da2178f98e9fd5b659e699edf`).
# 2026-11-11 가중 방식 및 DESI `n(z)` 민감도

기존 승인된 80-lens/400-control HSC shape rows, 기존 SACC 기준과 DESI/HSC
분포만 사용해 pooled 대조와 세트 동일 가중, 12개 DESI `n(z)` 곡선을 비교했다.
새 조회·다운로드·영상 읽기는 없었다. 렌즈 세트 동일 가중은 중간 반경에서
pooled 대비 최대 약 23% 변했고, 대조 세트 동일 가중은 특정 구간에서 24%까지
달랐다. DESI 중앙 `n(z)` 곡선 변화는 약 −2–4%다. 가중 민감도 공분산과 joint
`n(z)` posterior가 없고 영상/대조장 일치 및 cross-shear 문제가 남아 있어
`science_measurement=false`를 유지한다.

- 실행기: [`analyze_s23b_wl_weight_and_desi_nz_sensitivity.py`](analyze_s23b_wl_weight_and_desi_nz_sensitivity.py)
- 평가: [`wl_s23b_weight_and_desi_nz_sensitivity_assessment_2026-11-11.md`](wl_s23b_weight_and_desi_nz_sensitivity_assessment_2026-11-11.md)
- 산출물 (Git 외): `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_weight_and_desi_nz_sensitivity_2026-11-11_03/`
- summary SHA-256: `9c66b281e0920b79a5c7854e40081b873621867cfaca2a9991d5e93e031f1aa5`

## 2026-11-11 반경별 matched-control 지지와 고정 표본 프로파일

기존 QA 산출표만 재집계해 source bins 2–4의 반경별 matched-control 지지를
확인했다. 전역 74개 렌즈 세트 중 양쪽 모두 반경 지지가 있는 세트는
0.1–0.2부터 2–3 proper Mpc까지 60/67/69/73/74/74/74개다. 그중 반경마다
4개 이상 대조 중심도 지지되는 수는 52/59/64/67/70/72/74개다. 일곱 반경 모두에
렌즈와 대조 ≥4/5를 요구하는 고정 교집합은 52개다.

이 고정 52개 세트의 equal-lens `DeltaSigma_t`는 `10^13 Msun Mpc^-2` 단위로
23.554, 6.048, 4.003, 3.521, 2.521, 0.675, 1.255다. 반경별 지원 표본(60–74개)
profile과 비교하면 일부 구간에서 −55% 또는 +62% 차이가 난다. 52개 표본의
적색편이 중앙값은 전역 74개와 비슷하지만, aperture당 유지 source 행 중앙값은
3,415 대 3,119로 더 높다. 공간 공분산은 계산하지 않았고 per-set 산포는 오차가
아니다. 따라서 이는 지지·표본선택 민감도이며 `science_measurement=false`다.
동일 가중은 렌즈 세트 바깥 평균에만 적용하며, 각 세트 안의 positive random
aperture들은 source-pair lensing weight로 계속 pooled된다.

- 지지 감사기: [`audit_s23b_wl_matched_radial_support.py`](audit_s23b_wl_matched_radial_support.py)
- 고정 표본 재집계기: [`analyze_s23b_wl_fixed_common_support_profile.py`](analyze_s23b_wl_fixed_common_support_profile.py)
- 평가 산출물 (Git 외): `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_matched_radial_support_2026-11-11_02/` 및 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_fixed_common_support_profile_2026-11-11_02/`
- summary SHA-256: `c51a4d9268935a0919d0abd3da2f10b6cd16ee7e35f1e6e23a8f3fe8c9e51234` / `4d270b81737e5bcb6188941f9f26cd6255614fd7ae8d1e6f2977cf51a4686930`
- 다음: 고정 표본 estimator에 matched endpoint 공간 공분산을 붙이고 source 재사용·블록 폭·HSC footprint 민감도를 평가한다. 현 자료로만 수행하며 과학 측정 승격은 별도 관문을 통과한 뒤 판단한다.

## 2026-11-11 fixed-52 paired endpoint spatial covariance

기존 승인된 source rows만 사용해 고정 52개 lens-set의 lens-minus-matched-random
프로파일 공분산을 계산했다. 각 block deletion에서 source 또는 lens/random-center endpoint가
있는 pair를 함께 제거하고 S19A calibration을 다시 계산했다. 두 폭 모두에서 각 replicate에
동일한 7-bin lens-set 집합을 사용하고, annulus마다 최소 4/5 random center를 요구했다.

2.5°/5° grid 점유 셀은 64/28개, covariance rank는 모두 14/14이며 PSD다. deletion 후
최소 공통 set 수는 49/52 및 47/52다. `DeltaSigma_t` (`10^13 Msun Mpc^-2`)는
0.1–0.2 Mpc에서 `23.554 ± 5.286 / 4.410`, 2–3 Mpc에서 `1.255 ± 0.556 / 0.535`다.
전체 52-set point profile은 선행 산출물과 최대 0.03125 Msun Mpc⁻² 안에서 일치한다.
공분산은 근사 equal-area grid에 조건부이며 실측 HSC 영역 공분산은 아니다. `n(z)`,
photo-z dilution, target/control 영상 일치, member contamination, selection response가 남아
`science_measurement=false`를 유지한다.

- 실행기: [`analyze_s23b_wl_fixed_common_support_spatial_covariance.py`](analyze_s23b_wl_fixed_common_support_spatial_covariance.py)
- 상세 평가: [`wl_s23b_fixed_common_support_spatial_covariance_assessment_2026-11-11.md`](wl_s23b_fixed_common_support_spatial_covariance_assessment_2026-11-11.md)
- 결과 (Git 외): `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_fixed_common_support_spatial_covariance_2026-11-11_02/`
- summary SHA-256: `acdd4da69e4d7fe3dbb52d8584b3996bd836a41196da59aecfe3fcdf3e1fc6e9`
- `_01`은 5° 표만 deletion-support ledger에 남겨 superseded했다. `_02`에는 2.5° 64개와
  5° 28개 block deletion이 모두 기록되어 있다.

## Fixed-sample random-center weighting covariance

[`analyze_s23b_wl_fixed_equal_control_spatial_covariance.py`](analyze_s23b_wl_fixed_equal_control_spatial_covariance.py)
compares source-pair-pooled random subtraction with an equal mean of individually
calibrated positive random centers on the same fixed 52-lens sample. It uses paired
2.5°/5° endpoint-cell covariance and reports cross-component conditional diagnostics.
It reads only the approved shape-row export and pinned products, and keeps
`science_measurement=false`.

- Assessment: [`wl_s23b_fixed_equal_control_spatial_covariance_assessment_2026-11-11.md`](wl_s23b_fixed_equal_control_spatial_covariance_assessment_2026-11-11.md)
- Output (outside Git): `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_equal_random_center_spatial_covariance_2026-11-11_02/`
- Summary SHA-256: `00c8427ef6810b78f0900da1d515e3fca2ed77233bd012ebbeddbf5bc5e5c033`

## Official HSC PDR3 Wide tract/patch geometry

[`hsc_pdr3_tract_patch_geometry.py`](hsc_pdr3_tract_patch_geometry.py) parses
the four official HSC PDR3 Wide tract/patch tables already present at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/control_patches_2026-10-14_02/public_geometry`.
It pins byte counts and SHA-256 digests, normalizes unwrapped RA modulo 360°,
and uses patch polygons because the Wide tables list many patch parent tracts
without a full tract polygon in the field's tract section. Same-field patch
overlaps are resolved deterministically to the nearest containing patch center;
cross-field ambiguity is flagged. The parser makes no network request or
download. Official source links are in its module docstring.

## 2026-11-11 HSC S19A estimator calibration-source audit

정적 검토에서 Astrafek의 responsivity, per-object m/c 보정, edge-selection 보정식이
HSC PDR3 문서와 해당 문서에서 연결한 HSC 보정 유틸리티의 정의와 일치했다.
현재 export는 선택 경계값에 놓인 원천이 없어 strict/inclusive 경계 차이도 없다.
이는 구현식 대조이지 측정 검증은 아니다. photo-z 편향·희석 보정과 n(z) 불확실성,
selection-response 오차, 검증된 HSC 공간 공분산 및 대조 영상의 depth/PSF/sky
일치가 남아 science_measurement=false를 유지한다.

- 상세 기록: [wl_s23b_hsc_calibration_conformance_assessment_2026-11-11.md](wl_s23b_hsc_calibration_conformance_assessment_2026-11-11.md)
- 근거: [HSC S19A PDR3 문서](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s19a-shape-catalog-pdr3/)와 문서에서 연결된 [HSC calibration utility](https://github.com/mr-superonion/hsc-y3-calib/blob/main/utils.py)

## 2026-11-11 HSC 공식 기하 기반 WL tract 영향도 audit

고정 52개 렌즈·260개 대조장 표본을 공식 HSC PDR3 Wide patch polygon의
parent tract별로 하나씩 제외해 다시 계산했다. 102개 점유 tract를 제외할 때
공통 렌즈 표본은 48–52개로 남았다. 최대 절대 변화는 안쪽 0.1–0.2 proper-Mpc
구간에서 발생했다. equal-control-center profile의 ΔΣ 변화는
1.727 × 10¹³ M☉ Mpc⁻² (SPRING tract 9609), cross 성분 변화는
2.705 × 10¹³ M☉ Mpc⁻² (SPRING tract 9612)였다.

이 값들은 leave-one-tract 영향도이며 독립 오차나 공분산이 아니다. 인접 tract는
독립 표본이 아니므로 검출 유의도나 halo mass를 해석하지 않는다. 20개 승인
shape-row CSV와 고정 HSC Y3 n(z)만 재사용했고, 원영상·CAS·신규 다운로드는
없었다. 계속 `science_measurement=false`를 유지한다.

- 상세 결과와 재현 해시: [`wl_s23b_fixed52_exact_tract_influence_assessment_2026-11-11.md`](wl_s23b_fixed52_exact_tract_influence_assessment_2026-11-11.md)
- 출력 경로: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_fixed52_exact_tract_influence_2026-11-11_02/`
- summary SHA-256: `7924dfdbe27d9435099bc0d4990c192e28506cf9a50c9ebf0e8ca584d7859e86`

LagEunha에서 수동 CPU 재현:

```bash
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_MAX_THREADS=1 \
  /home/kjhan/miniconda3/bin/python3 analysis/icl_wgl_membership/hsc_debug/analyze_s23b_fixed52_tract_region_influence.py \
  --geometry-dir /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/control_patches_2026-10-14_02/public_geometry \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_fixed52_exact_tract_influence_2026-11-11_02
```

## 2026-11-11 HSC 선택응답 계수 오차 민감도

현재 HSC S19A 적층기는 공식 HSC edge-selection 보정계수의 점추정치를 쓰고
계수 오차는 결과 프로필에 전파하지 않는다. 이를 위해 기존 승인 표본의 고정
52개 렌즈·260개 대조장 무작위 차감 프로필을 기준으로 네 공개 계수 각각을
±1σ 변화시킨 전체 보정 적층을 재계산하는
[`analyze_s23b_wl_selection_coefficient_sensitivity.py`](analyze_s23b_wl_selection_coefficient_sensitivity.py)를
추가했다. 계수는 한 번에 하나씩만 바꾸고, 기준 프로필 재현을 통과하지 않으면
산출을 중단한다.

이 결과는 계수별 반응이지 결합 systematic covariance가 아니다. HSC 유틸리티는
각 multiplicative/additive 보정 내 두 선택경계 기여를 독립으로 합산하지만,
두 보정 사이의 공분산은 제공하지 않는다. 따라서 이 단계는 전체 오차나 검출
유의도를 산출하지 않고 `science_measurement=false`를 유지한다. 새 조회나 데이터
수신 없이 기존 승인 CSV와 n(z)를 사용한다.

```bash
sbatch analysis/icl_wgl_membership/hsc_debug/run_s23b_wl_selection_coefficient_sensitivity.sbatch
```

기본 결과는
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_fixed52_selection_coefficient_sensitivity_2026-11-11_01`에
저장된다.

`lageunha` CPU 실행이 완료됐다. 기준 profile은 고정 52-lens 결과를 최대
`0.03125 Msun Mpc^-2` 차이로 재현했다. 각 계수를 단독으로 ±1σ 변화시킨 profile
반응은 2.5°·5° equal-area spatial jackknife diagonal error의 최대 1.21%
(tangential)와 0.83% (cross)였다. 이 결과는 단독 계수 반응이며 결합 systematic
covariance가 아니다. 상세 수치, 입력 hash와 해석 제한은
[`wl_s23b_selection_coefficient_sensitivity_assessment_2026-11-11.md`](wl_s23b_selection_coefficient_sensitivity_assessment_2026-11-11.md)에
기록했다. `science_measurement=false`를 유지한다.

## 2026-11-11 WL 고정 지지 표본 선택 감사

[`audit_s23b_wl_fixed_support_sample_selection.py`](audit_s23b_wl_fixed_support_sample_selection.py)는
승인된 80개 렌즈 manifest와 기존 radial-support 결과만 읽어 전 구간 공통 지지가
있는 52개 렌즈와 같은 74개 전역 지원 표본의 나머지 22개 렌즈를 비교한다.
80개 중 6개는 74개 전역 지원 집합에 들어오지 않는다. 원본 shape-row CSV,
이미지, CAS 또는 네트워크는 사용하지 않았다.

비교 세부표는
[`wl_s23b_wl_fixed_support_sample_selection_assessment_2026-11-11.md`](wl_s23b_wl_fixed_support_sample_selection_assessment_2026-11-11.md)에
기록했다. 고정 52개의 원시 shape-source 중앙값은 나머지 22개보다 높다
(4,673 대 2,902.5, SMD 1.121). 이 수는 보정된 깊이 지도가 아닌 catalog 가용성
지표다. 52개 WL 표본의 일반화에는 선택효과를 고려해야 한다. 이 감사는 기술
비교이며 가설검정이나 shear 측정이 아니다. `science_measurement=false`를 유지한다.

평가서와 요약은
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_wl_fixed_support_selection_2026-11-11_01/`에
저장했다. `summary.json` SHA-256은
`226d92d196ecafa097ed669ca38d65c00f9a7808bd0c2356e654e5388dcbe76b`이다.

## 2026-11-11 fixed-52 ICL–satellite–weak-lensing relation and max-T audit

The fixed-52 ICL/WL overlap is 25 clusters in `i` and 27 in `i2`; the filter
sets are disjoint. Requiring at least 10% common pixel support in both ICL
annuli leaves 18 and 22 clusters, respectively. The four direct ICL–satellite
partial-correlation intervals all include zero. Across the 112 pointwise WL
association intervals, 8 excluded zero before familywise correction. A
shared-block wild-cluster max-T analysis of all 116 tests found 10 raw
`p<0.05`, but none remained below 0.05 after max-T adjustment; the smallest
adjusted p-value is 0.5154. This is an exploratory result for a selected
52-cluster sample. The 5-degree blocks are approximate, and unresolved
photo-z, shear-response, PSF, masking, magnification, member-dilution, and
source-reuse limitations remain. No physical detection is claimed;
`science_measurement=false`.

- Overlap audit: [`s23b_icl_wl_fixed_support_overlap_assessment_2026-11-11.md`](s23b_icl_wl_fixed_support_overlap_assessment_2026-11-11.md)
- Relation report: [`s23b_icl_satellite_wl_fixed52_relation_assessment_2026-11-11.md`](s23b_icl_satellite_wl_fixed52_relation_assessment_2026-11-11.md)
- Max-T runner/report: [`adjust_s23b_icl_satellite_wl_fixed52_max_t.py`](adjust_s23b_icl_satellite_wl_fixed52_max_t.py), [`s23b_icl_satellite_wl_fixed52_max_t_assessment_2026-11-11.md`](s23b_icl_satellite_wl_fixed52_max_t_assessment_2026-11-11.md)
- Relation summary SHA-256: `a7a9f30107bb800986d5d1aa4bf4e2012497092170e3bc27bf246587f7840e2f`
- Max-T summary SHA-256: `e1bbe2282adca533791dd3b4a3a82e048e99d74c4c2e1a82b27a85bb69a807c2`
- Max-T metrics SHA-256: `c29878527b2bbfdb379e4e4bfdc2fbe4fdee6d38bef67e4b3a204b4a38545010`

## 2026-11-11 S16A shape-catalog coverage count (in progress)

[`prepare_s23b_s16a_shape_count_plan.py`](prepare_s23b_s16a_shape_count_plan.py)
prepares 42 count-only SQL batches for all 1,043 pinned S23B cluster centers,
within 3 proper Mpc. The response would contain only center ID and the integer
number of distinct S16A weak-lensing shape objects (1,043 rows; estimated
maximum CSV size about 31.4 kB). It would not return source IDs, coordinates,
shape values, photo-z, images, or PSF data. The exact plan hash was approved
and submitted through the guarded CAS runner. S16A coverage is a ceiling only
and must not be combined with S19A without a separate release/calibration
audit. At this snapshot, the first five batches returned 125 verified counts;
50 centers have zero sources and the count sum so far is 83,561. The sixth
batch is running.

`submit_s23b_s16a_shape_counts.py` validates the reviewed plan in `--plan`
mode without network access. `--submit` requires the exact approved plan
SHA-256, an interactive Lageunha terminal, and one hidden HSC password prompt.
The adapter submits the 42 count-only batches sequentially and stores the two
aggregate columns for each center.

- Assessment: [`s23b_s16a_shape_count_plan_assessment_2026-11-11.md`](s23b_s16a_shape_count_plan_assessment_2026-11-11.md)
- Submitter: [`submit_s23b_s16a_shape_counts.py`](submit_s23b_s16a_shape_counts.py)
- Prepared plan (outside Git): `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_s16a_shape_count_plan_2026-11-11_01/`
- Running checkpoint (outside Git): `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_s16a_shape_count_plan_2026-11-11_01/results/count_status.json`
- Plan SHA-256: `c426129b150d440c22fca97d1cf2c722414a7f708132a6267563d77ccda997e0`
