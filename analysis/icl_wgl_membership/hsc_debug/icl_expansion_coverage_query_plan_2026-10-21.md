# S20A 확장 표본 HSC PDR3 패치 커버리지 쿼리 계획 — 2026-10-21

## 목적과 표본

현재 800개 이외의 로컬 S20A CAMIRA Wide 목록에서 `0.20 <= z_cl < 0.70`인 후보를 고르고, 기존 표본에서 6 arcmin 이상 떨어지도록 richness 내림차순으로 greedy 선택했다. 재현된 후보 수는 2,596개이며 기존 오프라인 확장 감사와 일치한다. 현 조건은 BCG 분광 적색편이를 요구하지 않고 cluster photo-z를 사용하므로 현재 분석 표본과 같은 모집단으로 가정하지 않는다. 후보 2,596개의 중심을 130,148개 로컬 CAMIRA member position과 richness/redshift/좌표로 연결해 검증했다.

Astrafek은 각 중심과 z_cl에서 계산한 proper 1,000 kpc 원주 8방향을 합해 **23,364 probe positions**를 만들었다. PDR3 Wide coadd 메타데이터에서 HSC-I mosaic patch, tract/patch ID, seeing, zeropoint, patch QA의 i-band depth/seeing 및 image-product 존재 플래그를 반환하도록 50개 cluster씩 52개 SQL 파일로 나눴다. 공식 PDR3 Schema Browser에서 `pdr3_wide.mosaic.band`가 HSC-I를 사용하고 i/i2를 병합한다는 설명과 요청 열의 존재를 확인했고, HSC stored-function 문서의 `patch_contains`/`public.skymap` 조인 예를 사용했다.

이 준비물은 **제출되지 않았다**. HSC 조회, 이미지, source catalog, shape, 개별 photo-z는 수신하지 않았다. 패치 원주 9점 검사는 metadata 기반 coverage screen이며 1-Mpc disk의 모든 픽셀이 덮이는지, 마스크가 usable한지, matched-control 영역이 있는지를 증명하지 않는다. PDR3 mosaic 메타데이터는 coadd 존재와 geometry/patch QA를 점검하는 데만 쓰며 실제 ICL stacking에는 global-sky image와 별도 control/quality 확인이 필요하다.

## 승인된 조회 범위와 제출 상태

- PDR3 Wide `mosaic` + `public.skymap` + `patch_qa` metadata query 52개 배치.
- 후보 2,596개 × 중심/1-Mpc 원주 8점, 총 23,364 probes.
- 반환 열: candidate ID, probe ID, `skymap_id`, tract, patch, band, image dimensions, seeing, zeropoint, i-band PSF depth/seeing, 두 image-product 플래그.
- 예상 결과: probe당 경계에서 1–4개 patch가 매치될 수 있어 약 23k–93k metadata rows, 대략 수 MB에서 15 MB 미만. SQL은 patch coverage가 있는 행만 내보낸다.
- 영상, source/object catalog rows, shape, photo-z, mask pixels, PSF FITS는 요청하지 않는다.

2026-10-21 사용자가 위 52개 metadata 배치 제출을 승인했다. 다만 기존에 승인된 random-center CAS 제출 프로세스가 lageunha에서 여전히 실행 중이다. 2026-10-21 01:35 KST 재확인에서도 PID 1049265는 살아 있었고, 세 번째 coordinate 배치가 `submitted`인 채 체크포인트는 2026-10-19 02:07 UTC 이후 갱신되지 않았다. 프로세스 상태는 `hrtimer_nanosleep`이었고 해당 관찰 순간 열린 네트워크 소켓은 없었다. 60초 뒤에도 같은 프로세스·체크포인트 상태가 유지됐다. 상태 파일에는 `error_type=HTTPError`가 남아 있다. CAS 동시 제출을 피하기 위해 metadata 조회는 아직 제출하지 않았다. 기존 프로세스를 종료·재시작하지 않고, 해당 작업이 완료된 뒤 승인된 범위를 실행한다.

제출기 `submit_icl_expansion_coverage_metadata.py --check`는 네트워크 없이 52개 SQL, 계획 해시, 열 범위 및 데이터 제외 플래그를 검증한다. 실제 제출은 lageunha의 기존 HSC 인증 세션에서 `--submit`으로 한 번에 한 배치씩 수행하며, 결과는 `cas_coverage_metadata/`에 저장한다. 체크포인트 완료 후 `validate_icl_expansion_coverage_metadata.py`가 각 결과의 해시와 스키마를 검증하고 9점 coverage screen을 만든다. 이 후처리 결과도 pixel-level coverage나 과학 측정으로 해석하지 않는다.

## 재현과 해시

```bash
python3 analysis/icl_wgl_membership/hsc_debug/prepare_icl_expansion_coverage_plan.py \
  --output-dir /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_expansion_coverage_plan_2026-10-21_01
```

산출물 경로: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_expansion_coverage_plan_2026-10-21_01/`.
`coverage_query_plan.json` SHA-256 `9da220c05f1f24c409a602fcc4ddffc8cbd5bc33f459d335cc8587f1cd235006`; 52개 SQL·후보 중심·CAMIRA member position·probe 파일 55개의 해시는 manifest의 `output_sha256`에 기록하고 모두 재검증했다. 생성 스크립트 SHA-256 `745161f2568fefb6654d539298c43482d453ae39d05465a37ca33e14f288595b`. Focused tests: 3 passed.

공식 참고: [PDR3 Database](https://hsc-release.mtk.nao.ac.jp/doc/index.php/database__pdr3/), [PDR3 stored functions](https://hsc-release.mtk.nao.ac.jp/rsrc/pdr3/stored_functions.html), [PDR3 Schema Browser](https://hsc-release.mtk.nao.ac.jp/schema/).
