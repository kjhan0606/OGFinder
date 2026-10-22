# S23B 후보 HSC-I coverage metadata 조회안 — 2026-10-22

## 준비한 정확한 범위

S23B CAMIRA의 이미 수신된 후보표 중 1,333개 추가 후보에 대해 PDR3 Wide `pdr3_wide.mosaic`, `public.skymap`, `pdr3_wide.patch_qa`만 검사하도록 준비했다. 각 후보 중심과 proper 1,000 kpc 원주 8개 방향, 총 **11,997 probe 위치**에서 HSC-I mosaic patch coverage 및 patch QA를 묻는다. 최대 50 후보씩 **27 SQL 배치**로 나눈다. 후보 ID는 `S23B:` 접두어를 붙여 S20A ID와 혼동하지 않는다.

CAS가 반환하도록 준비한 14개 열은 `candidate_id`, `probe_id`, `skymap_id`, `tract`, `patch`, `band`, `naxis1`, `naxis2`, `seeing`, `zeropt`, `imag_psf_depth`, `iseeing`, `imeas_exists`, `idet_bkgd_exists`다. 각 probe가 최대 네 patch에 대응한다고 둔 상한은 47,988개 metadata 행이다. 실제 행 수는 CAS 결과로만 알 수 있다.

질의는 이미지, pixel mask, object/source rows, shape, photo-z, PSF FITS를 요청하지 않는다. 9점 screen은 패치 위치 및 기본 QA를 확인할 뿐 1-Mpc 원 전체의 유효 coverage나 matched-control의 존재를 증명하지 않는다. 기존 control 영상에 중심이 들어가는 후보 7개도 query inventory에서 임의로 빼지 않고 표시 대상으로 유지한다. 결과 후 별도 WCS/mask/깊이 감사와 오염 없는 control 설계가 필요하다.

## 기존 조회안과의 구분

2,596개 S20A 후보, 23,364 probe를 담은 기존 계획은 2026-10-21 사용자 승인 범위지만 `prepared_not_submitted`이며 CAS 제출 사실이 없다. 이 계획은 cluster-only `0.20 ≤ z < 0.70` 표본으로, 현재 분석 selection과 같지 않다. S23B query는 다른 catalog release와 다른 좌표 목록이어서 기존 계획의 승인이나 잠재적인 coverage 결과로 대체할 수 없다. 새 S23B 범위는 **사용자 승인 대기**이며 이 문서 작성 중 CAS 연결·질의·자료 반출은 수행하지 않았다.

## 재현

- 생성기: [`prepare_s23b_icl_coverage_plan.py`](prepare_s23b_icl_coverage_plan.py)
- focused tests: [`test_prepare_s23b_icl_coverage_plan.py`](test_prepare_s23b_icl_coverage_plan.py), 2 passed
- 산출 디렉터리: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_coverage_plan_2026-10-22_01/`
- `coverage_query_plan.json` SHA-256: `e615ee7c2414109af6d3c156104d71766355d3b4e2a5dba3636a0142995557d7`
- 29개 artifact hash와 27개 SQL의 query allowlist를 재검증했다. 준비기 결과는 `prepared_pending_user_approval`, `cas_query_submitted=false`, `data_downloaded=false`다.

공통 제출기의 오프라인 `--check`는 기존 S20A 계획과 이 S23B 계획 모두에서 SQL·probe·열 범위·파일 해시를 확인한다. S23B는 추가로 계획 JSON SHA-256을 고정하고, 실제 제출 시 동일한 해시를 `--approved-plan-sha256`으로 지정해야 한다. 결과 후처리기는 S23B release-prefixed ID를 보존해 9점 coverage screen을 낸다.

계획 생성 명령:

```bash
python3 analysis/icl_wgl_membership/hsc_debug/prepare_s23b_icl_coverage_plan.py
```

검증은 `submit_icl_expansion_coverage_metadata.py --check --plan <위 산출 디렉터리>`다. 위 hash를 승인한 뒤 실제 조회할 때만 `--submit --approved-plan-sha256 e615ee7c2414109af6d3c156104d71766355d3b4e2a5dba3636a0142995557d7`을 추가한다. 제출기는 lageunha interactive 세션을 요구하며 비밀번호를 화면에 출력하지 않는다.

CAS 제출은 사용자의 이 정확한 1,333-candidate / 11,997-probe / 14-column metadata 범위 승인 후에만 진행한다. 이미지 다운로드는 이 승인의 범위에 포함되지 않는다.
