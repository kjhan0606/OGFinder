# HSC 검출·마스크·PSF 입력 예비 준비

2026-10-13, `prepare_detection_psf.py`를 lageunha CPU 단일 스레드로 실행했다.
기존 SEP·Astropy·NumPy·SciPy 환경을 사용하며 Web-local 배포에 의존성을 추가하지 않는다.
원본 영상과 활성 분석 manifest는 수정하지 않는다. 인증이나 GPU는 필요하지 않다.

## 실행과 산출물

```bash
ssh lageunha 'cd /home/kjhan/BACKUP/Astrafek/OGFinder && OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 MPLBACKEND=Agg python3 analysis/icl_wgl_membership/hsc_debug/prepare_detection_psf.py --manifest /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/validation_sample_2026-10-13_01/validation_manifest.json --work /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/detection_psf_2026-10-13_01'
```

실행 완료. 재실행에는 새 `--work` 경로가 필요하며 기존 결과를 덮어쓰지 않는다.
결과 기준 경로는
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/detection_psf_2026-10-13_01`이다.
각 `cluster_<ID>/`에 검출표 `detections.tsv`, 멤버 좌표 대조표 `members.tsv`,
분할지도와 품질/확대 마스크 `masks.fits.gz`, 유효면적 `mask_radial_coverage.tsv`,
PSF 후보표·원본 우표·원본 방사 평균, 전체 보기 `review.png`, BCG 확대
`bcg_review.png`, 상태 `preparation.json`을 저장한다.
입력/코드 해시·버전은 `summary.json`, 독립 검사는 `artifact_audit.json`,
후속 검토 대상은 `review_queue.tsv`에 있다.

## 방법과 한계

- SEP 1.5 sigma, 최소 5픽셀, deblend 32/.005, 배경 mesh 64를 사용한다.
  배경 제거는 검출·Kron 계산에만 적용하며 과학 영상에는 적용하지 않는다.
  S/N은 분산합에 의한 명목값이며 coadd 상관잡음 보정 전이다.
- BCG 위치 1 arcsec 이내 검출이 정확히 하나일 때만 위치 후보로 둔다.
  고유 후보도 검토가 필요하다. 나머지 분할 영역·Kron 2.5/3.5/5.0 배수를 가리고,
  검출과 대응하지 않는 시야 내 CAMIRA 멤버에는 반경 2 arcsec 원을 추가한다.
  이 설정은 마스크 완전성을 보장하지 않는다.
- 검출 시 BRIGHT_OBJECT bit만 별도로 제외하여 가려진 소스를 살펴보되, 출력에는
  원래 품질 bit와 BRIGHT_OBJECT를 보존한다. BCG 위 불량 픽셀이나 이웃 마스크도
  지우지 않는다. BCG 미확정 영상은 모든 소스를 가린 진단용 결과다.
- PSF 후보는 명목 S/N≥50, flags=0, b/a≥0.8, 모멘트 FWHM 0.35–1.5 arcsec,
  크기 분포 최빈구간·이웃 거리·중심부 마스크 조건으로 선정한다. CAMIRA 멤버와
  위치 대응 소스는 제외한다. 프레임당 최대 25개이며 외부 자료로 확인한 별은 아니다.
- 101×101 우표는 원본 선형 플럭스 단위를 유지하고 다른 소스/불량 픽셀은 NaN이다.
  중심 재정렬·정규화·하늘 제거를 하지 않았다. 방사 평균도 원본 진단값이다.
  약 8.4 arcsec 반경만 다루므로 넓은 PSF 날개 검증을 대신하지 않는다.
  후보 부족 시 Gaussian 등 대체 PSF를 생성하지 않는다.
- 모든 결과는 `science_measurement=false`, `mask_ready_for_science=false`다.
  소광 적용·대조장 차감·표면밝기 측정·ICL 스태킹은 수행하지 않았다.

## 실제 실행 결과

50/50개 준비 완료, 오류 0, 약 208초. 두 집합을 합쳐 과학 평균을 내지 않는다.

| 집합 | 영상 | 검출 소스 | PSF 후보 | 후보 ≥4개 영상 | BCG 위치 고유 후보 |
|---|---:|---:|---:|---:|---:|
| 면적 조건 충족 | 40 | 125,786 | 462 | 38 | 9 |
| 중심 마스크 스트레스 | 10 | 31,198 | 107 | 9 | 2 |

BCG 위치 1 arcsec 내 검출 없음 39개, 고유 후보 11개, 복수 후보 0개다.
미일치 39개의 최근접 검출 거리 중앙값은 2.879 arcsec다. 5104 확대 영상에서는
BCG가 보이지만 검출 중심·분할 영역이 어긋난 것을 확인했다. 검출/분리와 품질 마스크의
영향을 조사해야 하며, 이를 BCG 부재나 카탈로그 오류로 판정하지 않는다.
원본 품질 마스크가 BCG 중심 픽셀을 덮는 영상은 18개(면적 8, 스트레스 10)다.
PSF 후보 부족은 3806(1개), 6940(3개), 3700(2개)다.

lageunha 합성 회귀 3건 및 독립 산출물 검사 50/50개 통과. 검사 항목은 마스크
HDU/크기, 원본 마스크 보존, 확대 단조성, BCG 상태, 검출/우표 수, 표본 좌표
WCS 왕복 오차 <1e-5 픽셀, 우표 크기/미처리 상태, 유효면적 범위, 검토 이미지 존재다.
이는 파일·계산 경로 검증이며 과학 승인 검증은 아니다.

후속 순서는 BCG 중심 검출·분리 검토, 별 분류 확인과 경험적 PSF/넓은 날개 준비,
마스크 민감도·주입 회수, 공식 random/patch QA 기반 대조장 선정·영상 준비다.
사용자 실행 카탈로그 수집은 별도 진행 중이며 이 준비에 인증을 다시 요구하지 않는다.


후속 BCG 정점/분할 대조 및 멤버 마스크 보완 결과는 [bcg_review_README.md](bcg_review_README.md)에 기록했다. 기존 결과 폴더는 보존했다.
