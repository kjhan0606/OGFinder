# HSC 800개 표본의 소광·대조장 준비

CPU 실행 위치는 lageunha, 자료 위치는 `/gpfs/kjhan/Astrafek/data`이다.
`prepare_ancillary.py`는 기존 영상과 download manifest를 수정하지 않는다.
HSC 인증을 요구하거나 저장하지 않으며 HSC 영상 다운로드 코드를 대체하지 않는다.

```bash
python3 analysis/icl_wgl_membership/hsc_debug/prepare_ancillary.py \
  --sample /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13 \
  --work /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/ancillary_2026-10-13_01 \
  --cache /gpfs/kjhan/Astrafek/data/catalog/sfd98
```

출력 폴더는 새 이름이어야 한다. 위 이름은 이미 실행했으므로 재실행 시 새 이름을 쓴다.
원본 SFD 두 반구 지도는 NERSC 공개 배포본을 사용하고 파일 크기·SHA-256을 기록한다.
좌표는 ICRS에서 은하좌표로 변환한 뒤 ZEA WCS와 bilinear 보간으로 조회한다.
원본 SFD E(B−V)에 추가 0.86 배율을 적용하지 않는다.

`extinction_candidates.tsv`는 중심 E(B−V), 영상 영역 9지점의 최솟값·최댓값,
헤더의 i/i2 필터, 명목 A_i와 광량 보정 배율을 담는다. 계수 1.641/1.617의 출처는
[HSC 공식 hscPipe7 문서](https://hsc.mtk.nao.ac.jp/pipedoc/pipedoc_7_e/colorterms.html)다.
이 값은 PDR3의 실제 혼합 passband와 동일하다고 확인된 값이 아니다.
**PDR3 카탈로그 a_i 및 유효 passband와 대조하기 전에는 활성 측광 교정에 편입하지 않는다.**
9지점 표본은 변화량 진단이며 픽셀별 소광 지도나 cirrus 제거가 아니다.

`hsc_catalogue_requests.json`에는 표본별 두 SQL을 저장한다.

- 중심 3 arcsec 이내의 primary forced 천체와 a_i 조회. 복수 대응은 그대로 반환하며
  BCG 대응·거리·유효값을 검토한다. 대응이 없으면 임의로 먼 천체로 치환하지 않는다.
- 중심 30 arcmin 이내 공식 random 목록에서 i-band 중심 밝은별·edge 플래그가 명시적으로
  false인 후보를 최대 200개 조회. 고정 문자열을 사용한 해시 순서로 재현성을 유지한다.
  이 후보는 아직 matched control이 아니며 면적 무작위 표본이라고 주장하지 않는다.

`patch_qa.sql`은 공식 patch 관측 메타데이터를 조회한다. 표·함수는
[PDR3 random 안내](https://hsc-release.mtk.nao.ac.jp/doc/index.php/random-points__pdr3/),
[PDR3 데이터베이스 안내](https://hsc-release.mtk.nao.ac.jp/doc/index.php/database__pdr3/),
[좌표 검색 함수](https://hsc-release.mtk.nao.ac.jp/rsrc/pdr3/stored_functions.html)를 따른다.
SQL은 인증된 CAS에서 아직 실행하지 않았으므로 실제 스키마·출력 검증은 미완료다.
현재 환경에 인증이 남아 있지 않으며 사용자의 암호를 재요청하지 않았다.

## 사용자 인증 후 공식 목록 받기

사용자가 직접 인증하겠다고 지시한 후 `fetch_ancillary_catalogues.py`를 추가했다.
기존 cutout 다운로드는 그대로 두고 공식
[hscReleaseQuery API](https://hsc-gitlab.mtk.nao.ac.jp/ssp-software/data-access-tools/tree/master/pdr3/hscReleaseQuery)
규약으로 CAS 목록만 조회한다. syntax에서 실행:

```bash
ssh -t lageunha 'cd /home/kjhan/BACKUP/Astrafek/OGFinder && python3 analysis/icl_wgl_membership/hsc_debug/fetch_ancillary_catalogues.py --user kjhan0606'
```

실행당 숨김 암호 입력 1회다. `/`와 셸 특수문자는 그대로 처리하며 인수·환경변수·파일에
암호를 넣지 않는다. 기본 10개 은하단을 묶어 소광 80회, random 80회, patch QA 1회,
총 161개 쿼리를 순차 수행한다. 메일 알림은 끈다. 터미널은 완료까지 유지한다.
결과는 `sample_800_2026-10-13/catalogues_2026-10-13_01`에 CSV와 체크포인트로 저장된다.
동일 명령을 재실행하면 저장 파일 해시를 확인하고 완료된 결과를 건너뛰며, 기록된 CAS
작업 ID를 재사용한다. 제출 직후 연결이 끊겨 작업 ID를 받지 못했다면 중복 제출하지 않고
중단한다. 그 경우 CAS 작업 목록을 확인해야 한다. 실패한 SQL은 실제 CAS 스키마 검토가
필요하다. `--plan`은 인증 없이 계획만 출력한다.
모의 CAS의 다운로드/재개/변조/인증 실패/불확실한 제출/CSV 검증 5건을 통과했으며,
실제 인증 조회는 사용자 실행 전이므로 통과했다고 주장하지 않는다.

후속 실서버 확인: 사용자 실행 후 patch_qa 56,994행(29,154,612 bytes)은 완료·해시 확인됐다.
첫 소광 배치에서 CAS가 `# target_cluster_id,...` 헤더를 반환하여 기존 열 검사에서
중단된 것을 확인했다. 첫 열의 `# ` 접두사와 BOM을 처리하도록 고쳤으며, 중복 열과
불완전 행 거부는 유지한다. 실제 소광 파일은 10개 대상의 11행으로 파싱되지만 과거
버전의 전송 완료 표시가 없어 `.part`를 완료 파일로 승격하지 않았다. 같은 명령 재실행 시
기존 CAS ID 1059397에서 다시 받고 patch_qa(ID 1059396)는 재사용한다. 새 버전은
중단/오류를 최상위 상태에도 기록한다. 실제 헤더와 재개를 포함한 회귀 10건 통과.
이 확인 시점의 전체 작업은 161개 중 완료 1개이며, 정상 전체 종료가 아니다.

## 예비 검증 표본

`prepare_validation_sample.py`로 165개 면적 조건 충족 집합 중 40개와, 조건 미충족이면서
중심 밝은별 마스크 비율이 높은 321개 집합 중 10개를 별도로 선정했다.
z_cl 0.05 구간과 richness 15–20/20–30/30–50/50+로 층화하며 원본 800개는 보존한다.
할당 수와 조건부 포함 확률·seed·입력 해시는 selection_audit.json에 기록한다.
이 표본은 기술 검증용이며 800개 전체를 대표하는 과학 표본이 아니다. 면적 조건 충족도
PSF·하늘·오염 검증을 대신하지 않는다. 두 집합의 측정을 혼합해 평균하지 않는다.
산출물은 `sample_800_2026-10-13/validation_sample_2026-10-13_01`이다.
`validation_manifest.json`의 원본 영상/맵은 절대 경로로 연결하고 미검증 소광과 대조장은
비워 둔다. 측정 실행은 아직 하지 않았다.

대조장 연결 전 필수 확인은 `control_selection.json`에 기록한다. patch/tract 중복 제거,
전체 cutout footprint, 깊이·PSF·소광·유효면적 대응, 모든 CAMIRA 은하단과의 중첩,
대조장끼리의 중복을 검사한다. 은하단 배제 반경을 바꾸었을 때 생기는 환경 선택 편향도
비교한다. 대상 은하단의 적색편이와 동일한 각크기·마스크·교정을 쓰고, 각 대조장 자체의
소광값을 적용한다. 초기 800개 대조장 목표는 관측 배경 공분산에 따라 조정한다.

공식 목록 확보 후 검증된 좌표로 기존 `download_pilot.py`의 입력 manifest를 만들어
동일 coadd/bg 영상 경로를 재사용한다. 현재는 확정 대조 좌표나 대조 영상이 없으며,
임의의 좌표 이동을 matched control로 대신하지 않는다. PSF 날개·sky 전달·주입 회수와
관측 스태킹은 아직 미완료이고 `science_measurement=false`를 유지한다.


## 검출·마스크·PSF 입력 후속

50장 준비와 독립 산출물 검사 완료. BCG 미일치 39장, PSF 후보 부족 3장을
검토 대상으로 유지한다. 방법·결과·경로는 [detection_psf_README.md](detection_psf_README.md).
과학 스태킹은 아직 실행하지 않았다. 별도 CAS 수집은 최신 확인 시 running/random_008이다.

## 2026-10-14 후속 상태

위 수집 중 상태는 당시 이력이다. CAS 161개 조회는 이후 모두 완료됐고,
공식 random/patch QA의 검증된 CSV를 사용해 787개 잠정 대조 좌표와 은하단
배제 반경 민감도를 계산했다. 실제 검증 영상 50장의 외곽 잔차 진단도 완료했다.
대조 영상은 아직 없으며 대상 patch 대용값·전체 영상 품질 검증이 남아 있다.
방법·한계·최종 결과는 [control_sky_README.md](control_sky_README.md)를 따른다.
