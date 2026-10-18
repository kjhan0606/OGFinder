# 반지름별 마스크와 PSF 외부 대조 검토

2026-10-13 후속. 50장 기술 검증 표본의 면적·방향별 가림과 deblend 민감도를
검사한다. BCG 위치 후보 35장을 과학 표본으로 승인하는 작업이 아니다.

## 마스크 비교

`review_mask_coverage.py`는 저장된 보완 마스크와 원본 영상에서 재계산한
deblend contrast 0.001 결과를 비교한다. 기준 contrast는 0.005다. 두 설정에
동일한 품질 마스크·BCG 판정·CAMIRA 멤버 처리·Kron 2.5/3.5/5.0을 적용한다.
기존 BCG 차단 15장은 대조 설정에서도 차단을 유지하며, 통과 35장의 대조 설정
분할 번호는 이전 검토 결과와 일치하는지 확인한다. 기본 검출 경로를 바꾸지 않는다.

0–500 kpc를 25 kpc 간격으로 나누고 각 고리를 8개 방향 구역으로 나눈다.
`grid_coverage`는 영상 안의 고리 픽셀 수를 분모로 쓰며, `geometric_coverage`는
연속 원환 면적을 분모로 쓴다. 후자는 픽셀화로 생기는 미세한 1 초과값을 1로 제한한다.
최소 방향별 면적과 50% 미만 방향 수를 함께 기록한다. 50%는 기존 기술 검증
기준이며 과학 완전성을 보장하는 컷이 아니다. 이 계산에는 플럭스 측정이 없다.

중심 가림 원인은 원본 품질 마스크, 다른 소스 분할 영역, 이웃 Kron 타원,
검출에 대응하지 않는 CAMIRA 멤버의 원형 마스크로 나눈다. 원인은 서로 겹칠 수 있다.
타원/원 대조는 연속 카탈로그 위치, 마스크 픽셀 조회는 반올림된 픽셀 위치를 써서
경계에서는 미세한 차이가 가능하다. 가려진 중심을 강제로 비우지 않는다.

상위 경로는 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13`이며,
`mask_coverage_2026-10-13_02/radial_coverage.tsv`와 `summary.json`에 저장한다.
01 실행은 WCS 반환 스칼라의 자료형 오류로 첫 영상 처리 중 중단됐다. Python float로
변환해 수정하고 02에서 새로 실행했으며 기존 자료를 덮어쓰지 않았다.

## Gaia 외부 근거

공식 [ESA Gaia TAP 사용법](https://www.cosmos.esa.int/web/gaia-users/archive/programmatic-access)의
익명 테이블 업로드 방식을 사용한다. `match_psf_gaia.py`는 PSF 후보 569개 좌표를
Gaia DR3 `gaiadr3.gaia_source`에 반경 2 arcsec로 대조한다. 원본 ADQL, 업로드
VOTable, 응답 CSV, 입력/응답 해시, 비동기 작업 URL을 보관한다. HSC 암호가 필요 없다.
원본 좌표는 Gaia 기준시기이며 HSC coadd 유효 관측 시기로 전파하지 않았다.
필드 정의는 공식 [Gaia 조회 예시/스키마](https://www.cosmos.esa.int/web/gaia-users/archive/use-cases)를 확인했다.

동기 조회 `psf_gaia_2026-10-13_01`은 HTTP 408로 중단됐고,
`psf_gaia_2026-10-13_02`는 `--async-query`로 조회한다. 과학 승인 없이 조회 근거만 남긴다.
`summarize_psf_gaia.py`는 다음 조건을 모두 충족하는 후보에만 별이라는 외부 근거가
있다는 검토용 표식을 붙인다: 반경 2 arcsec에서 유일한 대응, 거리 ≤1 arcsec,
0<RUWE≤1.4, Gaia star probability≥0.99, duplicated/QSO/galaxy 후보 flag가 false,
양의 시차 또는 고유운동 한 성분의 |S/N|≥5. 이는 연구자가 정한 보수적인 예비
우선순위 조건이며 보정된 별 순도·완전성 모형이 아니다. 시차에는 양의 S/N만 사용한다.
고유운동 두 성분의 공분산을 무시한 결합 S/N을 만들지 않는다.

위치 대응 없음은 은하 판정이 아니다. Gaia 깊이·고유운동·시기 차이로 누락될 수 있다.
근거를 얻어도 HSC 관측시기, 포화·변광·이웃 오염, 넓은 PSF 날개 확인이 남는다.
`psf_approved=false`, `science_measurement=false`를 유지하며 대체 PSF를 생성하지 않는다.


## 면적 비교와 HSC 카탈로그 결과

50/50장·7개 설정(품질 마스크와 분리 2종×Kron 3종)·고리 20개, 총 7,000행을
저장했다. 범위/단조성 검사와 저장 마스크 200개(50장×기준 4종)의 독립 픽셀 재계수가
모두 통과했다. 새 합성 테스트 5건도 통과했다. `artifact_audit.json`,
`gate_review.json`, `coverage_comparison.png`에 검사와 요약을 저장했다.

면적 집합 40장에서 BCG 후보와 모든 0–500 kpc 고리의 면적 ≥50%를 함께 만족한 것은
기준 contrast 0.005/Kron 2.5의 5724·619·5814 세 장이다. contrast 0.001에서는
619·5814 두 장이다. Kron 3.5/5.0에서는 두 contrast 모두 0장이며 스트레스 10장도
전 설정에서 0장이다. 이 기준 통과도 과학 승인을 뜻하지 않는다.

기준 Kron 2.5 중심 가림 원인(중복 허용)은 품질 마스크 18장, 다른 분할 소스 10장,
이웃 타원 20장, 멤버 추가 원형 마스크 20장이다. 연속 위치 대조와 픽셀 조회의 경계
차이를 고려해야 한다. contrast 변화에 따른 최대 고리 면적 변화가 10%p를 넘는
대상은 5724 하나이며 최대 약 72.7%p다. `mask_case_5724_2026-10-13_01`의
이미지 비교와 재계산에서 BCG의 정점·광도 중심·모멘트·Kron은 동일하지만 이웃 검출
목록/마스크가 달라지는 것을 확인했다. 이를 BCG 중심 오차로만 설명할 수 없다.
조건을 완화하거나 넓은 마스크를 임의로 축소해 표본을 늘리지 않았다.

별도 HSC 조회는 161개 모두 complete다: patch QA 56,994행, 소광 800개 대상
1,212행, random 후보 800개 대상 160,000행. 직전 감사 158개 파일과 새 완료 3개를
연결해 크기·SHA-256·행 수 검사를 완료했다. `catalogue_integrity_2026-10-13_01.json`
및 `catalogue_integrity_2026-10-13_final.json`에 기록했다. random은 좌표 후보이며
검증된 대조장 800개 또는 대조 영상 다운로드 완료를 뜻하지 않는다.

다음 과학 검증은 5724와 중심 이웃 겹침 사례의 분할/광량 분배 검토, 보수적 마스크와
이웃 광도 모형의 주입·회수 비교다. 일부 픽셀 면적이 늘었다는 이유만으로 그 설정을
채택하지 않는다. 관측 PSF·하늘 전달·대조장 확인 전에는 실제 스태킹을 승인하지 않는다.


## 완료된 Gaia 대조와 검증

ESA 공간 결합은 제한 시간 내 끝나지 않아 CDS XMatch의 `vizier:I/355/gaiadr3`에
동일 569개 좌표를 2 arcsec/all-match로 대조했다. 대응 73행을 받았고, 그 source_id
73개를 ESA Gaia DR3에서 인덱스 조회해 누락된 QSO/galaxy 후보 flags와 별 확률,
천체측량 정보를 확보했다. 지연된 ESA 공간 결합은 ABORTED로 정리했다.
추가 계정·암호 입력은 없었다. 앞으로 이 사례의 재조회는 거대한 공간 결합보다
CDS 위치 대조→ESA source_id 확인 순서를 우선한다.

CDS 원본·요청 매개변수·자료 정의 ReadMe·ESA ADQL/CSV는
`psf_cds_2026-10-13_01`에 있다. 정규화된 대조와 분류는 그 아래 `evidence/`의
`gaia_matches.csv`, `evidence_summary.json`, `normalization_audit.json`이다.
CDS 정의는 [I/355 ReadMe](https://cdsarc.u-strasbg.fr/viz-bin/ReadMe/I/355?format=html&tex=true),
조회 규약은 [CDS XMatch 클라이언트 문서](https://astroquery.readthedocs.io/en/latest/api/astroquery.xmatch.XMatchClass.html)를 참고했다.
CDS 열의 PSS/PQSO/PGal 확률을 Boolean 후보 flags로 대체하지 않았다.

원본 후보 ID/좌표와 CDS·ESA source_id를 대조하고, 천구 거리를 독립 계산했다.
CDS가 반환한 거리와 ESA 2016.0 좌표에서 재계산한 거리의 최대 차이는
4.92e-7 arcsec다. J2000 전파 좌표를 잘못 연결하지 않았다. 원본 응답 해시,
73개 고유 source_id, 반경·입력 ID를 검사했다.

569개 중 73개가 각각 유일하게 대응했다. 위의 보수적 별 근거 조건은 32개가
통과했고 23개 은하단 영상에 분포한다. 대응했지만 근거 불충분은 41개,
Gaia 위치 대응 없음은 496개다. 여러 후보에 공유되는 Gaia source_id는 0개다.
32개는 **PSF 검토 우선 후보**이며 검증된 PSF나 깨끗한 별 표본으로 승인하지 않는다.
Gaia 대응이 없는 희미한 소스는 HSC 다중대역 형태·색 등 추가 자료로 확인해야 한다.


후속 조건부 주입/회수 432장면을 완료했다. [neighbor_injection_README.md](neighbor_injection_README.md)에 모형 가정, 실패 건수, 전달·오차 공분산과 하늘 민감도를 기록했다. 관측 검증 승인과 구분한다.
