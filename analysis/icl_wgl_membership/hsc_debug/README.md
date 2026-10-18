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

## 2026-10-18 약한렌즈 shape-catalog 겹침 사전계산

`prepare_wl_shape_preflight.py`는 800개 은하단의 3 Mpc 원뿔에 대해 HSC S19A
shape source 개수만 확인할 CAS SQL 배치를 만든다. 실행 시 네트워크 조회나 자료
다운로드를 하지 않으며, 준비된 쿼리도 아직 제출하지 않았다.

```bash
python3 analysis/icl_wgl_membership/hsc_debug/prepare_wl_shape_preflight.py \
  --manifest /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/sample_manifest.json \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_shape_preflight_2026-10-18_01
```

32개 SQL 배치와 `preflight_manifest.json`이 생성된다. 렌즈별 source 행 수와 예상
전송량을 확인한 뒤 자료 전송 여부를 정하고, 사전계산 단계에서는 shape 행을
가져오지 않는다.
