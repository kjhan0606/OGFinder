# HSC ICL–위성 연관성 및 배경은하 렌즈 연구 현황 — 2026-10-20 갱신 (초판 2026-10-19)

## 현재 판정

Astrafek에는 HSC PDR3 영상·대조장·PSF 입력을 사용한 ICL 스택, 위성 외곽광 연관성 분석, 배경은하 약한렌즈 스택과 각 단계의 재현 스크립트가 구현되어 있다. 지금 결과는 모두 탐색 단계다. ICL–위성 관계는 다중비교와 sky 처리 민감도를 넘지 못했고, 렌즈 B-mode 마스크 방향은 외부 영역과 좌표 대조로 검증되었다. 렌즈 신호에는 아직 무작위 중심 차감과 대응 matched-field 공분산이 없다. 따라서 두 갈래 모두 `science_measurement=false`를 유지한다.

## 보유 자료와 실행 이력

- HSC 표본은 CAMIRA 은하단 800개다. 표적 영상 800개와 확보된 매칭 대조장 798개에서 동일 필터로 쌓을 수 있는 표적–대조장 쌍은 783개(i 302, i2 481)다. HSC PDR3는 `i`/`i2` 응답곡선이 비슷하지만 동일하지 않고 적색 천체에서는 보정 차이가 더 커질 수 있다고 명시하므로 두 필터는 섞지 않는다. 픽셀별 색이 없는 ICL 영상에 객체별 유효응답 보정을 추정해 적용하지 않고, 현재 Astrafek 스택은 각 target/control 쌍의 동일 필터만 비교해 필터별로 보고한다 ([HSC PDR3 처리 설명](https://hsc-release.mtk.nao.ac.jp/doc/index.php/processing__pdr3/), [필터 응답 차이](https://hsc-release.mtk.nao.ac.jp/doc/index.php/filterfrac-2/)).
- ICL 분석의 공통지원 2-D 스택과 K=2.5, 3.5, 5.0 독립지원 비교는 Astrafek의 고정 스크립트로 생성·검증되었다. 독립지원 추정치는 픽셀 교집합의 영향을 줄이는 민감도 분석이며 기본 측정값으로 채택하지 않았다.
- 2026-10-20에 Astrafek 독립 검증기를 보완해 세 마스크 설정의 스택 원자료 연결, 쌍별 프로파일, 출력 지도 차감식, 마스크 정의 및 실행 코드 해시를 재검증했다. 검증 산출물은 각 스택 폴더의 `independent_verification_2026-10-20.json`이며, 이 재검증은 과학 측정 게이트를 바꾸지 않는다.
- 재검사한 표적–대조 2-D 스택의 50–100 kpc 표적 초과 평균은 i에서 `41.39 ± 2.90`, i2에서 `44.59 ± 2.32` nJy/arcsec²이고, 100–200 kpc에서는 각각 `12.92 ± 1.49`, `10.75 ± 1.01`이다. ±는 은하단 쌍 산포에서 계산한 기술적 SEM이며 전체 공간 공분산·sky·PSF 계통오차가 아니다. 해당 고리의 중앙 공통 픽셀 지원률은 i 15%/14%, i2 19%/18%로 제한적이다. 스택 요약 SHA-256은 `9c804cfb12ff2b86a5bdc8fdbfea94aa1dc473711e6f7ac6f4d233cc6b853173`이며 결과는 계속 `science_measurement=false`다.
- 2026-10-20의 추가 400–500 kpc 바깥 고리 기준 민감도에서는 i 기준광이 내부 50–100, 100–200, 200–300 kpc 고리별로 각각 `1.658`, `1.668`, `1.783` nJy/arcsec²였다. HSC tract-center 근사 연결영역을 하나씩 제외해도 i 기준광은 양수(`1.079–2.196`)로 남았고 i2 전체 평균은 0 부근이며 제외 범위는 0을 가로질렀다. 따라서 i 오프셋이 한 연결영역에만 집중됐다는 설명은 지지되지 않지만, 잔류 sky와 실제 외곽 ICL은 구분되지 않는다. 세 점유 영역의 leave-one-region-out은 공분산이 아닌 영향도다. 상세 수치·방법·입력 해시는 [바깥 고리 지역 영향도 평가](icl_outer_baseline_region_influence_assessment_2026-10-20.md)에 있고 재현기는 [`audit_icl_outer_baseline_regions.py`](audit_icl_outer_baseline_regions.py)다.
- 동일한 지원 선택 표본에서 400–500 kpc 기준광과 내부 target-minus-control 신호의 쌍별 공변동도 redshift, richness, 면적지원을 조정한 뒤 100–300 kpc에서 양수였다(i partial `r=0.338/0.455`, i2 `0.355/0.421`); 50–100 kpc i 구간은 0을 포함했다. 이 여섯 탐색 비교는 외곽 기준광과 안쪽 광량이 독립이 아님을 보일 뿐, 공변동을 잔류 sky와 실제 ICL 중 어느 것으로 설명할지 판별하지 않는다. [바깥 기준광 공변 평가](icl_outer_baseline_coherence_assessment_2026-10-20.md)에 5,000회 연결 블록 재표집·입력 해시를 기록했다.
- 위성 연관성에는 반경 구간별 94쌍(i 26, i2 68)이 있다. 주입·회수 검증은 위치별 HSC-I `pdr3_wide` coadd 모델 PSF 744개와 승인된 시험 PSF 1개를 사용해 57개 중심 및 688개 위성 위치에서 실행되었다. PSF 감사 매니페스트 SHA-256은 `6375bc58df6d5a15d837cf5b13d3db79ebdc6944dbf7510164f465b6cbe72f7e`다. 이 약 7 arcsec 코어는 넓은 ICL PSF 날개를 검증하지 않는다.

양수 성분 모멘트로 설명적으로 잰 745개 모델 코어의 가우스 등가 장축 FWHM 중앙값은 `1.084 arcsec`(q16–q84 `1.015–1.147`)였고, 중심 57곳과 위성 688곳의 중앙값은 비슷했다. 위치별 폭에는 변동이 있어 기존 개별 커널 주입은 타당하지만, 이 지표는 seeing 보정이나 확장 날개 검증이 아니다. [PSF 코어 분포 평가](psf_core_variation_assessment_2026-10-20.md)에 방법과 파일 해시를 기록했다.

동일한 57개 프레임에서 저장 ICL 200–300 kpc target-minus-control 신호와 중심 PSF 코어 폭의 양의 공변동을 관찰했다. redshift·richness·공통 면적지원을 조정한 `r`은 i `0.500`(13개), i2 `0.247`(44개)였지만 연결공간 블록 재표집 범위는 넓었다. 이는 소표본의 코어 PSF 민감도 진단이며 원인이나 보정값을 정하지 않는다. [ICL–PSF 코어 공변동 평가](icl_psf_core_dependence_assessment_2026-10-20.md)에 마스크 coverage 차이와 입력 해시를 기록했다.

근사 tract-region 평균효과를 추가한 민감도에서도 z·richness·공통 coverage 보정 후 `r`은 i `0.484`, i2 `0.231`로 남았지만, i 표본은 13개뿐이고 영역 고정효과 후 자유도 8이다. 네 연결영역 중 i는 2개, i2는 3개에 걸쳐 있었고 57개 matched target/control은 각각 같은 영역으로 배정됐다. 이는 근사한 지역 평균 차이만으로 양의 관계가 사라지지 않았음을 보이는 점검이지 검출 또는 유의성 검정은 아니다. [PSF–ICL 영역 민감도 보고서](icl_psf_region_fixed_effects_assessment_2026-10-20.md)에 해시와 제한을 기록했다.
- 기존 영상의 8,078개 compact 후보와 Gaia 진단은 target-specific extended PSF 날개 표본을 제공하지 못했다. Gaia 대조는 569개 후보 위치에서 73개 행을 반환했고, 32개 별 후보 중 전체 네 배경/마스크 조건 통과는 4개였다. 2026-10-20 사용자 승인으로 나머지 7,509개 위치도 2초각 Gaia DR3 대조를 완료해 944개 후보–source 쌍을 얻었다. 그중 480개는 288개 프레임에 분포한 항성 근거 후보, 464개는 항성 근거 불충분, 6,565개는 무매칭으로 분류됐다. 480개 모두 Gaia 변광 플래그가 `NOT_AVAILABLE`이고, 추가로 승인된 Gaia DR3 `vari_summary`·`vari_classifier_result` 조회도 480개 source 모두 일치 행이 없었다. 이는 변광 분석 자료가 없다는 뜻이며 비변광 확인은 아니다. 480개 후보의 원본 HSC 영상 진단은 288개 프레임의 333개 우표, 1,920개 방법 조합을 처리해 완료했다. 358개 조합이 진단 프로필을 생성했고, 137개 후보가 한 조합 이상, 42개가 네 조합 모두 유효했다. 64–96픽셀 고리의 기하 지원을 통과한 274개 방법 프로필은 양수 118개·0 이하 156개로 부호가 섞였다. 이를 후보별로 묶으면 외곽 지원 후보 113개의 profile/core 중앙값은 `−3.07×10⁻⁷`(q16–q84 `−1.85×10⁻⁶–1.29×10⁻⁶`)이며, 지원 설정 전체 양수/전체 비양수/혼합은 각각 38/58/17개다. 네 설정 외곽 지원을 모두 확보한 후보는 30개뿐이다. Gaia `ref_epoch=2016.0` 고유운동만 적용하고 HSC 후보 좌표를 고정한 ±10년 시나리오에서는 480개 모두 2″ 매칭 반경 내였으며 최대 분리는 1.572″였다. 실제 HSC coadd epoch는 확인되지 않았지만 우표 중심은 HSC 검출 위치다. 이 민감도는 정확한 epoch별 Gaia cross-match 확인을 대체하지 않는다. 반복 방법은 독립 표본이 아니며 배경/마스크 선택에 민감하므로 날개 검출 증거로 사용할 수 없다. 영상 기반 별 검증도 남아 있어 PSF 과학 승인은 없다. 정량 요약·해시·한계는 [Gaia PSF 후보 확장 기록](extended_psf_README.md)에 있다.
- 기존 Gaia 별 우표 중 333개(227 표적 프레임)의 matched-control 위치 지원을 추가로 검사했다. target cluster center 기준 각거리 오프셋을 control center에 WCS로 옮긴 결과 319개(217 프레임)는 동일 `i`/`i2` 하위 필터와 401×401 전체 footprint를 만족했고, 14개는 하위 필터 불일치로 제외됐다. 이는 다음 paired background-profile 진단의 입력 적합성일 뿐 PSF wing 검출이나 ICL 보정은 아니다. 보고서·체크섬은 [matched-control 지원 평가](empirical_psf_matched_control_support_assessment_2026-10-20.md)에 있다.
- 약한렌즈 표적 쪽은 CAS 원자료 3,069,544 lens–source 행(800 렌즈)을 가져왔다. 품질·거리·반경 컷 뒤 2,153,004 쌍, 609 렌즈, 2,051,991개 고유 배경원으로 줄었다. 기존 결과는 표적 중심 스택만 포함한다.

## ICL과 위성 외곽광

독립지원 비교에서 50–200 kpc 구간의 최대 같은 쌍 추정치 이동은 K=5, i2의 175–200 kpc에서 32.9%였다. 따라서 마스크 지원 방식이 외곽 프로파일의 크기를 실질적으로 바꿀 수 있다. 94쌍 분석에서 i2 전체지원 부분상관은 `r≈0.293`이지만 2.5°, 5°, 10° 공간 블록 폭에서 신뢰구간 해석이 달라진다. 별도 가족단위 검정은 48개 설정 중 추론지원 기준을 통과한 30개 어디에서도 Holm 또는 max-T 보정 뒤 유의성을 보이지 않았다.

2026-10-20 원자료 해시를 다시 확인한 표본 보존 감사에서 독립 50–75 kpc 위성 스택은 i 70개/42개 은하단 `2.51 [-0.30, 5.50]`, i2 308개/140개 `1.19 [-0.24, 2.49]` nJy/arcsec²로, 두 필터 모두 은하단 부트스트랩 구간이 0을 포함했다. 결합된 ICL–위성 표본은 94개 cluster/filter 쌍(i 26, i2 68)이다. i2 부분 Pearson `r=0.293`의 5° 공간 블록 구간 `[0.042, 0.497]`은 양수지만, 일반 쌍 부트스트랩 구간 `[-0.007, 0.526]`과 Spearman 공간 블록 구간 `[-0.088, 0.301]`은 0을 포함한다. i 부분 Pearson은 `r=0.079`, 공간 블록 구간 `[-0.186, 0.479]`다. 따라서 선택된 i2 부분상관은 흥미롭지만 순위 통계와 재표집 방식에 민감하며, 검출로 주장할 수 없다. 재현 가능한 해시 연결 감사는 `icl_satellite_sample_retention_2026-10-20_02/summary.json`에 있으며 SHA-256은 `7952f53bd29289dc182b24abdd7b0356b9a09e40c4e15055a0ea1f6aa1ca354f`다.

연결 target/control 블록 하나씩을 제외하는 추가 영향도 검사에서는 i2 `r`이 세 블록 폭 모두에서 양수였으나 `0.166–0.177`까지 낮아졌다. i는 2.5°와 10° 폭에서 한 블록 제외 시 부호가 바뀌었다. 이는 표본 내 단일 블록 영향 진단일 뿐 신뢰구간이나 검출 검정이 아니며, 기존 다중비교 결과를 바꾸지 않는다. 상세 결과와 영향 블록 ID는 [블록 제외 영향도 평가](satellite_icl_leave_one_block_out_assessment_2026-10-20.md)에 기록했다.

94쌍의 부분상관에 patch-QA 연결영역 고정효과를 추가해 영역 간 평균 오프셋 민감도를 점검했다. i2는 `r=0.29264`에서 `0.29105`로 유지됐고, i는 `0.07898`에서 `0.13860`으로 바뀌었다. 각 target/control 쌍의 양 끝점이 이 근사 영역 분할에서 항상 같은 연결요소에 속해 두 효과는 선형 중복이며 한 번만 적합했다. 결과는 지역 평균 오프셋만으로 i2 상관이 생기지는 않았음을 시사하는 탐색 진단이다. 소표본·선택 민감성, sky/PSF·마스크 계통효과와 공간 공분산은 해결되지 않았고 `science_measurement=false`를 유지한다. 재현 및 제한은 [영역 고정효과 평가](satellite_icl_region_fixed_effects_assessment_2026-10-20.md)에 있다.

동일 표본의 은하단별 leave-one-out 영향도에서는 전체 지원 i2 `r`이 항상 양수였지만 `0.151–0.338`까지 변했고, i는 한 은하단 제외에서 부호가 바뀌었다. 중심 적격 부분표본은 i2 15쌍, i 8쌍뿐이며 공간 블록 구간이 모두 0을 포함한다. 따라서 현재 상관은 검출 근거가 아니다. 이는 영향도 검사이지 신뢰구간·검출 검정이 아니다. 세부 결과는 [cluster 영향도 평가](satellite_icl_cluster_influence_assessment_2026-10-20.md)에 기록했다.

redshift·richness 조정 순위 부분상관도 함께 확인했다. 전체 지원 i2는 0.139, i는 0.248로, Pearson 계수와 크기 차이가 있고 i의 단일점 영향도도 달랐다. 이는 이 선택 표본에서 관계 형태가 상관 정의에 민감함을 보여주며, 독립적인 검출 근거는 아니다.

14개 반경별 상관과 paired raw-to-corrected 차이에 대해 같은 연결 5° target/control 블록을 일곱 반경에 걸쳐 공동 재표집한 max-T family-wise 구간을 추가했다. raw, 외곽 차감, 두 값의 차이로 나눈 6개 필터·통계 family의 42개 구간 모두 0을 포함한다. 미보정 i 5–10 kpc 외곽 차감 상관 `rho=0.212`의 동시 95% 구간은 `[-0.054, 0.479]`다. 따라서 관측된 반경별 ICL–위성 관계는 현재 표본에서 다중비교 보정 뒤 검출되지 않는다. [반경별 family-wise 평가](satellite_radial_icl_association_familywise_assessment_2026-10-20.md)에 전체 구간, 해시, 재현기를 기록했다.

94개 ICL–위성 결합 행은 200–300 kpc ICL 지원 parent 224개 중 i 37%, i2 44%만 차지한다. 선택 여부에 따라 면적지원, redshift/richness 및 근사 HSC 영역별 결합률이 달라 현재 관계를 전체 CAMIRA 모집단에 일반화할 수 없다. 표준화 평균차 감사는 [연관표본 선택 평가](icl_satellite_association_selection_assessment_2026-10-20.md)에 있다. 작은 표본에서 선택효과를 보정할 수 없어 탐색 상태를 유지한다.

이후 같은 94행에서 ICL 공통 면적지원률을 redshift·richness와 함께 조정하면 i2 `r`은 `0.293`에서 `0.320`, 표적·대조 지원률을 별도로 넣으면 `0.305`로 이동했다. 하지만 같은 5° 블록의 조정 전후 차이 백분위 범위는 두 모델 모두 0을 포함했고, i 26행 결과도 매우 넓었다. 면적지원 조정은 제외된 행이나 모집단 선택을 복구하지 못한다. 결과는 [면적지원 조정 민감도](satellite_icl_coverage_adjustment_assessment_2026-10-20.md)에 있으며 검출 판단은 바뀌지 않는다.

이 94개는 독립 50–75 kpc 위성 프로파일 182개와 ICL 200–300 kpc 지원 행 224개의 정확한 교집합이다. 교집합 밖에는 ICL 전용 i/i2 45/85행과 위성 전용 16/72행이 있어, 연관분석 규모는 두 지원 관문의 겹침에 의해 결정된다. 추후 표본 확장 시 한쪽 관문만 완화하지 않고 양쪽 영상·마스크 지원을 함께 재평가해야 한다.

같은 patch-QA 트랙 연결요소로 matched target–control ICL 프로파일도 점검했다. 783쌍 중 세 영역에 487, 244, 52쌍이 분포하고 1.6–3.0° 인접 기준에서 분할은 안정적이다. 가장 큰 영역을 제외하면 50–100 kpc 초과광은 i에서 1.8%, i2에서 9.2% 감소했고, 100–200 kpc는 i에서 10.2% 증가, i2에서 9.8% 감소했다. 200 kpc 바깥 저신호 프로파일은 영역 제외에 따라 부호가 바뀌는 고리가 있다. 50–200 kpc의 pair별 중앙 공통 픽셀 지원은 14.5–19.0%에 그친다. 이는 저지원 외곽광의 지역 민감도를 보여줄 뿐, 신뢰구간이나 ICL 검출 검정이 아니다. [ICL tract-region 영향도 평가](icl_hsc_tract_region_influence_assessment_2026-10-20.md)에 재현 스크립트, 데이터 매칭, 해시, 전체 고리 결과를 남겼다.

같은 HSC 트랙 연결요소로 위성–ICL 상관 자체도 재계산했다. 94행의 target과 matched control 중 어느 한 점이라도 제외 영역에 있으면 그 행을 빼고 redshift·richness 회귀를 다시 적합했다. 전체 i2 partial Pearson `r=0.293`은 영역 1을 빼면 `r=0.048`로 감소하고, Spearman `rho`는 `+0.111`에서 `−0.048`로 바뀐다(68행 중 29행 제외). 따라서 앞서 직사각 공간 셀 삭제에서 양수였던 partial Pearson 부호만으로 연관성이 지역에 강건하다고 볼 수 없다. i는 전체 26행의 `r=0.079`이고 영역 1 삭제 후 8행만 남아 `r=−0.378`이지만, 이 잔여 표본은 너무 작아 부호를 해석하지 않는다. 세 영역 영향도는 공분산·신뢰구간이 아니며, 관계의 원인이 천체물리적 환경인지 선택·영상 계통인지 구분하지 못한다. 상세 결과는 [위성–ICL tract-region 영향도](satellite_icl_tract_region_influence_assessment_2026-10-20.md)에 기록했다.

은하수 소광을 픽셀별로 적용했을 때 기존 중심 소광 보정과의 ICL 차이는 사실상 0이었다. 반면 국소 평면 제거는 고정된 74쌍에서 i2 부분상관을 0.239에서 0.124로 낮추고 구간이 0을 포함하게 했다. 100–500 kpc 확산광 주입 시험에서는 100 kpc 신호의 99.1–99.5%, 250 kpc의 86.0–90.8%, 500 kpc의 63.0–71.1%가 Astrafek의 downstream 국소 평면 제거 뒤 남았다. 500 kpc 원형/납작 프로파일 모두 표적·대조장 회수율이 비슷했고 쌍차 부트스트랩 구간은 0을 포함했다. 이 시험은 고정 마스크와 전달된 coadd/bg를 사용해 상류 HSC sky 처리, PSF·광자 잡음·마스크 재검출을 재현하지 않으므로 관측 프로파일 보정값으로 쓰지 않는다. 수치와 실행 해시는 [확장 규모 평면 주입 평가](icl_plane_injection_extended_assessment_2026-10-20.md)에 있다. 공개 확장 PSF 참조를 적용한 i 대조군 13개 시험에서는 단위 diffuse 주입의 중앙 21.5%가 ICL 측정 고리에 들어오고, 위성 외피는 35.3%가 들어왔다. 이 값도 특정 입력·마스크의 전달 민감도이지 관측값에 곱할 보정계수가 아니다.

추가로 같은 최신 코드와 57개 위치별 HSC-I 모델 PSF에서 100/500 kpc diffuse 주입의 검출·마스크 재생성 응답을 짝비교했다. 평면 보정 ICL 응답/원시 annulus 응답비의 중앙값은 i `1.002 → 0.644`, i2 `0.997 → 0.633`이었고, 필터별 paired 500–100 차이 중앙값은 `−0.341`과 `−0.362`였다. 후속 200/300 kpc 주입에서 retention 중앙값은 i `0.913/0.816`, i2 `0.918/0.812`였다. 모두 q16–q84 프레임 간 산포이지 신뢰구간이 아니며 HSC 상류 sky·PSF 날개·coadd 공분산을 포함하지 않는다. i2 한 프레임은 300/500 kpc diffuse-only 시나리오에서 위성 측정 지원을 잃었지만 별도 ICL 응답은 보존됐다. 이 민감도는 국소 평면의 광범위 확산광 감쇠가 마스크 재검출 뒤에도 나타남을 보일 뿐, 관측 ICL 보정이나 ICL–위성 검출 근거가 아니다. 100/500 결과는 [PSF·마스크 재검출 규모 비교](icl_psf_redetection_scale_comparison_assessment_2026-10-20.md), 4개 규모의 결과·해시·제한은 [확산광 규모 격자 평가](icl_psf_diffuse_scale_grid_assessment_2026-10-20.md)에 기록했다.

같은 57개 프레임의 네 규모를 동일한 연결 5° target/control 블록으로 5,000회 공동 재표집했다. 양수 raw-response 공통 표본(i 12개/9블록, i2 44개/24블록)에서 paired median retention 변화는 100 kpc 대비 200/300/500 kpc에서 i `−0.064/−0.178/−0.349`, i2 `−0.075/−0.182/−0.362`였고 95% 블록 부트스트랩 구간은 모두 0 아래였다. i2 100 kpc raw 주입 응답이 `0.00338`인 cluster 4687은 차감 후 음수여서 평균 ratio를 크게 왜곡했으므로 median을 주 통계량으로 삼고, 분모 없는 `corrected−raw` 응답 차이도 함께 기록했다. 이는 선택 주입표본에서의 방법 전달 추론이며 물리적 ICL 검출이나 관측 보정이 아니다. 해시·구간·분모 감사는 [확산광 규모 블록 부트스트랩 평가](icl_psf_diffuse_scale_block_bootstrap_assessment_2026-10-20.md)에 있다.

공분산 감사는 position별 PSF 주입 표본 57개(i 13, i2 44)의 표적·매칭 대조 coadd에서 9×9 단거리 잔차 커널을 측정하고 실제 signed aperture 가중치에 전파했다. 큰 top-hat 근사 분산 비는 표적 1.59, 대조장 1.54–1.55였다. 정확한 ICL 및 위성 가중치와 양쪽 ICL 영상의 분산을 사용하면 대각 분산만 쓴 경우보다 ICL 잡음 표준편차가 i에서 1.234배, i2에서 1.240배였고, 잡음만으로 계산한 ICL–위성 상관 중앙값은 각각 −0.0073, −0.00091이었다. 이는 해당 57개 프레임에 조건부인 단거리·정상 커널 민감도다. source Poisson 잡음, 상류 sky 전달, 위치별 완전 공분산은 포함하지 않으며 공간 블록 오차에 단순 곱하면 같은 변동을 중복 셀 수 있다. 따라서 관측 연관성을 보정하거나 검출 여부를 결정하는 수치로 사용하지 않는다.

57쌍의 잡음 공분산을 같은 표본의 측정 ICL·위성값으로 전파한 10,000회 Monte Carlo에서는 i2 부분표본 44개의 측정 편상관 `r=0.409`가 finite-lag 잡음 범위 `[0.320, 0.465]`를 보였다. 하지만 전체 68쌍 i2 추정값은 `r≈0.293`이고 공간 블록 구간은 블록 폭에 민감하다. 공분산이 있는 표본만의 값이 전체보다 큰 점은 선택 효과를 경고한다. i 부분표본은 finite-lag 잡음 범위가 0을 포함한다. 이 조건부 잡음 실험은 spatial-block 범위와 합치지 않으며 과학 검출을 뜻하지 않는다. 세부 결과는 [측정 잡음 Monte Carlo 평가](icl_satellite_measurement_noise_mc_assessment_2026-10-19.md)에 기록했다.

현재 CAMIRA 후보를 더 넓히는 오프라인 검토에서는 기존 선택 아래 새로 분리되는 후보가 57개뿐이고 새 5도 하늘 셀은 없었다. 적색이동 범위를 넓혀도 분광 확인을 유지하면 새 5도 셀은 2개였다. 따라서 같은 발자국에서 영상 수만 늘리는 것으로 독립 하늘 영역의 부족을 해결하기 어렵다.

2026-10-20에 승인된 S23B CAMIRA Wide 별마스크 표 12,162행을 기존 S20A 목록 및 완결된 800개 PDR3 영상 표본과 좌표 대조했다. 1초각 이내 S20A 중복은 4,492행이며, 현 표적 위치와 겹치는 607행 중 559행은 현재 적색편이/BCG 컷을 통과한다. 기존 800개에 6분각 간격을 적용한 S20A 비중복 후보 풀은 탐욕 선택 기준 1,333개, 신규 5도 좌표 셀은 22개다. 이는 향후 독립 하늘표본을 늘릴 가능성을 보여주지만, 이 카탈로그는 최종연도 S23B/PDR4 쪽 자료인 반면 현재 영상은 PDR3다. 신규 후보의 PDR3 coadd·마스크·PSF·깊이·matched-control 지원은 확인되지 않았으며, 5도 셀 수는 target/control 연결 공분산 블록 수가 아니다. 따라서 새 후보의 영상 취득 가능성이나 기존 33개 독립 블록을 넘어서는 유효 블록 수를 아직 주장하지 않는다. 세부 수치·체크섬과 재현 절차는 [S23B 카탈로그 감사](s23b_camira_catalog_audit_2026-10-20.md)에 기록했으며 `science_measurement=false`를 유지한다.

2026-10-20에 기존 CAMIRA 후보 풀과 power 근사를 기계적으로 연결했다. 탐색적 `r=0.293`, 30개 비교 보정의 80% 검정력 계획은 독립 연결 블록 180개를 요구하며 현재 값은 33개다. 가장 넓은 기존 목록 후보 시나리오는 2,596개를 추가하지만 새 target 5도 셀은 6개, 전체 84개다. 현재 pair-block 알고리즘은 target 셀을 공유하는 행을 연결하므로, 84개는 해당 후보 표본이 만들 수 있는 연결 블록 수의 상한이며 계획치보다 최소 96개 적다. 실제 영상 지원·매칭 대조장·공분산을 적용한 유효 블록 수는 이보다 낮을 수 있다. 이 판단은 사후 탐색효과를 사용한 Gaussian Fisher-z 설계 근사와 현재 카탈로그에 한정된다. HSC 공식 PDR3 CAMIRA 자료는 Wide와 Deep/UD 목록 모두를 PDR3(S20A)에서 구성했으며, 현재 Astrafek 입력은 이미 CAMIRA S20A Wide 별마스크 목록이다 ([PDR3 CAMIRA catalog](https://hsc-release.mtk.nao.ac.jp/doc/index.php/camira-cluster-catalog__pdr3/)). 따라서 Wide 목록을 다시 받아도 새 독립 릴리스가 아니다. Deep/UD는 명목 면적이 각각 27 및 3.5 deg²로 작아 주 분석의 독립 블록 부족을 해결하지 못한다 ([HSC survey](https://hsc-release.mtk.nao.ac.jp/doc/index.php/survey__pdr3/)). 이 면적은 설계 개요이지 현재 표본의 실제 유효 면적·독립 공분산 블록을 뜻하지 않는다. 재현 스크립트와 해시 연결 산출물은 `assess_icl_catalog_expansion_feasibility.py` 및 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_catalog_expansion_feasibility_2026-10-20_04/`에 있다. 신규 HSC 질의·좌표 반출·영상 다운로드는 없으며, `science_measurement=false`다.

같은 `r=0.293` 설계 근사를 가설적 검정 수 1·6·30과 효과 `r=0.2, 0.293, 0.4`로 확장했다. 단일 사전등록 검정이면 80% 검정력에 92개 블록이 필요해 현재 33개와 기존 후보 상한 84개 모두에 미달한다. 30개 family에서는 180개가 필요하다. 단일 검정도 효과가 `r=0.2`라면 196개가 필요하고, 84블록 상한에서 근사 검정력은 `0.437`이다. 반대로 탐색치 `r=0.293`을 실제 효과로 간주해도 84블록에서 단일 검정 검정력은 `0.765`에 그친다. 이 결과는 검정 수를 사후에 줄여 기존 분석을 정당화하지 않으며, 향후 하나의 주평가변수는 새로운 자료를 보기 전에 정해야 한다. 재현기 [`assess_icl_power_family_sensitivity.py`](assess_icl_power_family_sensitivity.py), 테스트 `test_assess_icl_power_family_sensitivity.py`, 입력 해시가 든 결과는 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_power_family_sensitivity_2026-10-20_01/`에 있다. 요약 SHA-256은 `40a472ab0483fd3cccbb5ba0730f7207ba6773c343067eb9213c3adf7caa1883`; 스크립트 SHA-256은 `8279553f57c0c42733e747e760855a7f00f8e022f8a240278ef5424d1b1c8bf3`다. 모든 경우는 Fisher-z 근사이고 `science_measurement=false`다.

후속 독립 표본을 위한 검토용 설계 후보는 i2 `icl_delta_200_300`과 `satellite_delta_50_75`의 단일 양방향 partial-Pearson 주검정, redshift/richness 조정, 그리고 보수적인 `r=0.2` 계획값에서 196 usable 연결블록을 제안한다. 반경·필터 선택과 추정량은 현재 탐색 결과에서 나온 것이므로 이 설계는 사전등록이 아니며 독립 표본에만 적용할 수 있다. 지원·매칭 탈락률, 경험 PSF 날개, diffuse+satellite 공동 주입 검증을 선행조건으로 기록했다. 문서 [`icl_satellite_confirmatory_replication_design_candidate_2026-10-20.md`](icl_satellite_confirmatory_replication_design_candidate_2026-10-20.md)는 연구자 검토용 후보일 뿐 `science_measurement=false` 상태를 바꾸지 않는다.

**해석:** 양의 ICL–위성 관계는 흥미로운 가설이지만 검출로 보고할 수 없다. 다음 ICL 관문은 반경·필터·마스크·통계량을 분석 전에 고정하고, 실제 검출·마스크 재생성을 포함하는 diffuse+위성 공동 주입으로 sky·PSF 날개·마스크 누출 전달함수를 평가하는 것이다. 독립 통계력을 높이려면 단순히 기존 CAMIRA 선택을 완화하기보다 새 독립 하늘 영역과 완전도·렌즈 보정을 함께 설계해야 한다.

## 배경은하 약한렌즈

기존 표적 전용 `DeltaSigma` 프로파일은 `0.1–0.2 Mpc`에서 약 `1.61×10^14 M☉ Mpc⁻²`, `2–3 Mpc`에서 약 `1.53×10^13 M☉ Mpc⁻²`로 감소한다. 52개 공간 셀 삭제 공분산은 독립 쌍 SEM보다 낫지만 셀이 연결된 HSC 관측영역과 일치하지 않을 수 있다. 이 표적 전용 프로파일에는 무작위 중심 차감, 대조장 분산, `n(z)`·shear 보정 불확실성, reduced-shear 보정이 빠져 있어 검출 유의도나 최종 질량으로 해석하면 안 된다.

2026-10-20 네 HSC Y3 source-bin 약한렌즈 profile을 동일 52개 sky-cell jackknife로 비교했다. 0.1–0.7 Mpc의 세 구간에서 `DeltaSigma`가 source bin에 따라 달랐고, 반경 일곱 구간 Holm 보정 뒤에도 GLS 일관성 검정이 남았다. 이 결과는 SACC `n(z)`와 source selection/calibration의 bin별 transfer를 추가 점검해야 함을 보인다. cluster-only 진단이며 random 차감, `n(z)` 불확실성 및 연결 관측영역 공분산을 포함하지 않으므로 질량 또는 렌즈 검출로 해석하지 않는다. 세부 프로파일과 공분산은 [tomographic 일관성 평가](wl_tomographic_consistency_assessment_2026-10-20.md)에 기록했다.

이어 같은 SACC 입력에서 source-bin 3·4의 평균 적색이동을 DESI 논문의 주변 68% 이동 구간 중앙값과 네 모서리로 평행 이동한 다섯 시나리오를 실행했다. 독립 검증기는 다섯 56×56 공분산 모두의 대칭성과 순위 51, 입력 체크섬, 요약 envelope를 통과시켰다. bin 간 일관성 불일치는 모든 시나리오에서 0.1–0.7 Mpc에 남았고, 반경 일곱 구간 Holm 보정 p값은 각각 0.1–0.2 Mpc에서 `7.6–8.6×10⁻⁵`, 0.2–0.4 Mpc에서 `9.2×10⁻⁹–1.7×10⁻⁸`, 0.4–0.7 Mpc에서 `7.4×10⁻⁴–1.04×10⁻³`였다. 이는 이 근사 이동 범위로 불일치가 해소되지 않음을 보이며 실제 보정 posterior 검증의 우선순위를 높인다. 상수 이동은 DESI의 보정 `n(z)` posterior 전파가 아니고, cluster-only 프로파일에 random 차감 및 연결 관측영역 공분산도 없으므로 렌즈 검출이나 질량 제약이 아니다. 상세 수치와 파일 해시는 [n(z) 이동 격자 평가](wl_published_nz_shift_grid_assessment_2026-10-20.md)에 기록했다.

같은 승인 lens–source 행에서 source object ID 재사용의 영향을 5도 52개 블록으로 평가했다. 선택된 2,153,004쌍 중 고유 source 2,051,991개 가운데 99,671개(4.86%)가 여러 lens 쌍에 재사용됐고, 그 source가 차지하는 행은 200,684개(9.32%)였다. lens endpoint, source endpoint, 양 endpoint 공동 삭제 공분산은 모두 56차원·순위 51이었다. joint-endpoint 대 lens-only 프로파일 오차비의 중앙값은 tangential 1.008(범위 0.973–1.046), cross 1.006(0.988–1.046)으로, 이 블록 구성에서 source 재사용이 marginal error를 크게 바꾸지는 않았다. 세 삭제 방식 모두 source-bin 일관성 Holm p값이 0.1–0.7 Mpc에서 `6.0×10⁻⁹–4.7×10⁻⁴`였고 더 큰 반경에서는 유의하지 않았다. 이는 source-bin별 calibration/`n(z)` 차이를 가리키는 진단이지 렌즈 검출이 아니다. Astrafek 독립 검증기 [`validate_wl_source_reuse_covariance.py`](validate_wl_source_reuse_covariance.py)가 대칭성·PSD·순위·프로파일 오차·GLS/Holm·재사용 집계를 검사해 통과했다. 요약 SHA-256은 `ac578d5818aae383586b62855e293bbfb5a1427e8a53e252d2569e6e4302b428`, 검증 JSON SHA-256은 `3f565f508a40ca1d98aef4d7f58b1152303694546811bb2cee0602801c870fd0`이다. 작은 음의 최소 고유값은 최대 행렬 원소 대비 `8×10⁻¹⁷` 이하로 수치 반올림 수준이며 PSD 상대 허용치는 `10⁻¹²`다. 이 민감도는 source-level 전체 공분산, 연결된 HSC 관측영역 공분산, random subtraction, 보정 `n(z)` posterior를 대체하지 않으며 `science_measurement=false`를 유지한다.

후속 문헌 검토에서 2026-08-25 개정된 HSC PDR3 S19A 형상표본 전용 DESI DR1/DR2 클러스터링 적색이동 보정을 확인했다. 2026-10-20 사용자 승인으로 고정 commit `abc70af29910e49fce17a63ef3083160d515fb6c`의 두 scale-cut NPZ와 JSON sidecar만 수신했다(4개, 총 549,818 bytes; 파일별 크기와 SHA-256은 GPFS manifest에 기록). 자료는 네 source bin의 1,000점 `n(z)` mean/median 및 점별 lower/upper와 magnification-corrected 곡선을 포함하지만 joint posterior draw는 포함하지 않는다. 같은 승인 lens–source 행과 52-cell jackknife로 실행한 paired 민감도 분석에서 DESI central curve는 SACC 기준 대비 첫 source bin의 median `DeltaSigma`를 +12.5–14.8%, bins 2–4를 약 −1.1–3.5% 바꿨다. 토모그래피 내부 불일치는 작아졌지만 0.1–0.7 Mpc에서 계속 남았다. 따라서 DESI 보정이 불일치를 해소하거나 질량 신호를 확정하지 않는다. 분석 방법·독립 검증·정량표·입력 및 출력 해시는 [DESI `n(z)` 분포 민감도 평가](wl_desi_nz_distribution_sensitivity_assessment_2026-10-20.md)에 기록했다. 결과는 cluster-only이며 random 차감과 posterior 불확실성 전파가 없어 `science_measurement=false`다.

약한렌즈 adapter의 shear responsivity, 개별 `m/c` 보정, HSC 경계 선택편향 계수를 논문과 HSC가 연결한 S19A 보정 코드에 대조했다. 소스 수준 수식 검토에서 구현은 일치한다. 이는 추정기 공식 검증일 뿐 새 `n(z)`, 표적–random 차감, 대응 공분산 검증은 아니다. 특히 토모그래픽 내부 bin 불일치가 남아 있어 현재 프로파일은 여전히 탐색 결과다. 세부 대조와 endpoint 주의사항은 [S19A 보정식 검토](wl_s19a_shape_calibration_formula_review_2026-10-20.md)에 기록했다.

HSC 문서는 `b_mode_mask`가 큰 B-mode 구역을 제거하는 데 쓰인다고 설명하지만 0/1 방향은 적지 않는다. Astrafek은 형상 품질 및 tomography 컷 뒤 좌표를 HSC Y3가 제외한 GAMA09H 사각영역(RA 132.5–140°, Dec 1.6–5°)과 직접 대조했다. 안쪽의 77,021개 lens–source 행은 모두 raw 0이고 바깥쪽 2,157,489개는 모두 raw 1이었다. shear 결과를 보지 않고 raw 0을 제외값으로 판별했다. 이 극성을 사용한 HSC Y3 four-bin SACC `n(z)`, 개체별 shear 보정, selection-bias 보정을 적용한 cluster-only 스택은 2,153,004쌍과 609개 렌즈를 유지했다. 52개 연결 sky 셀 jackknife는 별도로 계산됐지만, matched-random 차감과 source-overlap 보정은 아직 없다. 이 해석은 HSC가 공개한 마스크 설명과 Y3 시스템atics 연구의 제외영역 정의에 근거한다 ([HSC shape-catalog 문서](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s19a-shape-catalog-pdr3/), [HSC Y3 PSF systematics 연구](https://academic.oup.com/mnras/article/525/2/2441/7202336)).

52개 균등 셀 잭나이프가 연결된 HSC 영역 변동을 대표하는지 확인하기 위해 patch-QA 트랙 중심 연결성으로 leave-one-region 영향도도 계산했다. 트랙 중심 인접 임계값 1.6–3.0°에서 전체 840개 트랙은 네 연결요소, 표본은 세 요소(렌즈 499, 246, 55개)로 안정적으로 나뉘었다. 가장 큰 요소를 뺄 때 표적 `DeltaSigma`는 0.1–0.7 Mpc에서 최대 16.4% 이동했고, 가장 안쪽 cross 성분은 `+14.3 × 10^12 M_sun Mpc^-2` 이동했다. 이는 표본의 지역 영향이 무시할 수 없다는 점을 보여주는 영향도 진단이지 공분산이나 신뢰구간이 아니다. 세 영역만으로는 안정적인 전체 공분산을 추정할 수 없고 최근접 트랙 중심 및 근사 인접 그래프의 경계 불확실성도 남는다. 입력·요약 해시 및 재현 방법은 [tract-region 영향도 평가](wl_hsc_tract_region_influence_assessment_2026-10-20.md)에 있다. 따라서 기존 52셀 결과와 함께 보되 무작위 중심 차감 및 대응 field 공분산 관문을 대체하지 않는다.

무작위 중심 후보 좌표의 사용자 승인 CAS 내보내기는 800 렌즈, 예상 780,191 행, 32개 배치 계획이다. 2026-10-20 21:42 KST SSH 재확인에서 lageunha 제출 PID `1049265`는 약 34시간 35분째 살아 있었고 `hrtimer_nanosleep` 상태였다. 실제 좌표 체크포인트는 계획의 `cas_coordinates/coordinate_status.json`이며 2/32 배치(33,514행) 완료, 세 번째 CAS 작업 `1062312` 제출 상태였다. 체크포인트에는 최상위 `error_type=HTTPError`가 있지만 시각·대상 요청이 없어 CAS 작업의 terminal 실패로 해석할 수 없다. 제출기는 구버전 코드를 메모리에 올린 채 살아 있으므로 그대로 두고 중복 제출하지 않는다. `1062312`는 Slurm 작업 ID가 아니다.

완료된 앞의 두 배치만 분리해 부분 기하 파일럿을 실행했다. 승인 계획·배치 해시와 렌즈별 예상 행 수를 검증한 33,514좌표 중 2,129개(6.35%)만 CAMIRA 3 Mpc 조리개 중첩과 정확 좌표 중복을 함께 고려한 뒤 남았다. 첫 50개 렌즈 중 31개에만 후보가 있고 19개는 0개였다. 계획 순서의 첫 50개는 대표 표본이 아니므로 이 비율을 800개에 외삽하지 않는다. 이 결과는 초기 후보 풀이 렌즈별로 불균일할 수 있음을 보여준다. 재현 및 상세 제한은 [부분 좌표 기하 파일럿](wl_random_coordinate_geometry_pilot_2026-10-20.md)에 기록했다. 전체 기하 스크리닝은 좌표 내보내기가 32/32 완료되고 각 배치의 행 수·해시가 검증된 뒤에 시작한다.

2026-10-20에 Astrafek의 random-stack 실행기에 paired endpoint jackknife를 추가했다. 검증된 lens/random shape 행에서 같은 5도 등면적 셀을 만들고, lens 또는 source 끝점 중 하나라도 삭제 셀에 놓인 pair를 두 표본에서 함께 제외한 뒤 각 복제에서 전체 토모그래픽 보정·스택을 다시 계산한다. 실행기는 lens 행에서 재계산한 프로파일이 저장된 cluster-only 산출물과 일치하는지도 검사한다. 합성자료와 기존 표적 산출물 provenance 검사는 통과했지만 random shape 행이 아직 없어 paired covariance는 실제 자료로 실행하지 않았다. 셀 폭 민감도와 연결된 HSC 영역의 영향도, photo-z posterior 및 reduced-shear 불확실성은 남아 있으므로 과학 게이트는 계속 닫혀 있다.

후보 풀이 부족할 때를 위해 기존 2R–30′ query와 겹치지 않는 30′–60′ annulus의 800-lens count-only SQL 계획을 로컬에서 준비하고 validator로 확인했다. HSC CAS에는 제출하지 않았다. 실행 전 사용자 승인이 필요하며 기존 coordinate export가 끝날 때까지 병렬 제출하지 않는다. 범위·계획 SHA·후속 승인 관문은 [무작위 중심 보충 고리 계획](wl_random_coordinate_expansion_plan_2026-10-20.md)에 기록했다.

별도 오프라인 감사는 ICL용 HSC 무작위 대조장 798곳을 약한렌즈 random 중심으로 재사용할 수 있는지 평가했다. 768곳이 CAMIRA 3 Mpc 측정 조리개와 겹쳤고, 남은 30곳은 모두 서로 3 Mpc 조리개가 분리되어 있었다. 따라서 이 기존 대조장 표본만으로는 800개 렌즈의 대표 random stack을 만들 수 없다. 감사는 이미지 및 target manifest와 hash-verified 7,319개 CAMIRA 목록만 사용했으며 새 HSC query나 shape 행 반출은 하지 않았다. 재현기는 [`audit_wl_matched_control_centers.py`](audit_wl_matched_control_centers.py)이고, 결과는 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_icl_matched_control_geometry_2026-10-20_04/summary.json` (SHA-256 `12eec246699c24a6f194e271758bfdc09d85660b1eeedcc73a54c7aec1109cc9`)에 있다. 이 30곳도 shape-source footprint와 선택함수를 검증하기 전에는 과학용 random 중심이 아니다.

현재 인증된 좌표 내보내기 프로세스가 CAS를 계속 확인 중이다. 그 프로세스를 유지하고, 32개 배치가 모두 행 수·해시 검증을 통과한 뒤 Astrafek 전체 기하 스크리너를 실행해 CAMIRA 3 Mpc 중첩·중복 좌표·렌즈별 지원량을 검사한다. 그 뒤 전체 조리개 영상 범위, 마스크, 깊이, PSF, sky, 소광과 공간 블록을 맞춰 무작위 중심을 선택한다. 모양 카탈로그 행 수와 전송량을 count-only로 산정한 뒤 무작위 원천의 shape 행을 반출하려면 별도 사용자 승인을 받아야 한다. 그 데이터 없이는 표적–무작위 차감과 대응 공분산을 완료할 수 없다.

교차 성분의 1.5–2.0 Mpc 구간은 `3.1 ± 1.6 M☉ pc⁻²`로 개별 bin에서 약 2σ 벗어나므로, 무작위 차감과 대응 공분산을 적용한 뒤 null-test 대상으로 확인한다. 다른 bin은 이 값에 맞춰 재선택하지 않는다.

## Astrafek 재현 경로와 근거 문서

실행기, Slurm 스크립트, 검증기, 요약 생성기는 `analysis/icl_wgl_membership/hsc_debug/`에 보관한다. 대형 HSC 파일과 실행 결과는 Git 밖의 `/gpfs/kjhan/Astrafek/data/`에 두고, 코드·SQL 계획·체크섬·요약을 저장소에 기록한다.

- [표적–대조장 2-D ICL 스택](icl_2d_target_control_stack_README.md), [독립지원 민감도](icl_2d_independent_support_assessment_2026-10-19.md)
- [ICL–위성 다중비교 검정](icl_satellite_familywise_2026-10-18.md), [반경별 블록 폭 민감도](icl_satellite_perbin_block_scale_assessment_2026-10-19.md), [ICL tract-region 영향도](icl_hsc_tract_region_influence_assessment_2026-10-20.md), [위성–ICL 상관 tract-region 영향도](satellite_icl_tract_region_influence_assessment_2026-10-20.md)
- [소광 민감도](icl_pixel_extinction_assessment_2026-10-19.md), [평면 주입 전달](icl_plane_injection_assessment_2026-10-19.md), [위치별 PSF 및 공동 주입](psf_convolved_joint_injection_assessment_2026-10-19.md), [coadd 잔차 공분산](icl_coadd_residual_covariance_assessment_2026-10-19.md), [ICL 바깥 고리 지역 영향도](icl_outer_baseline_region_influence_assessment_2026-10-20.md), [ICL 바깥 기준광 공변동](icl_outer_baseline_coherence_assessment_2026-10-20.md), [ICL–위성 lag 공분산](satellite_icl_lag_covariance_assessment_2026-10-19.md), [표적·대조 공분산](satellite_icl_target_control_covariance_assessment_2026-10-19.md), [측정 잡음 전파](icl_satellite_measurement_noise_mc_assessment_2026-10-19.md), [PSF 파일 감사](psf_picker_batch_2026-10-19.md)
- [CAMIRA 표본 확장 검토](camira_expansion_pool_assessment_2026-10-19.md), [B-mode 극성 좌표 검증과 calibrated stack](wl_bmode_coordinate_validation_2026-10-19.md), [표적 전용 약한렌즈 공분산](wl_cluster_spatial_jackknife_assessment_2026-10-19.md), [tract-region 영향도](wl_hsc_tract_region_influence_assessment_2026-10-20.md), [tomographic 일관성 평가](wl_tomographic_consistency_assessment_2026-10-20.md), [무작위 좌표 기하 파일럿](wl_random_coordinate_geometry_pilot_2026-10-20.md), [무작위 중심 보충 고리 계획](wl_random_coordinate_expansion_plan_2026-10-20.md), [ICL 대조장 random 적합성 검사기](audit_wl_matched_control_centers.py)

이 문서는 결과를 통합해 읽기 위한 과학 현황 기록이다. 위의 탐색 분석을 논문용 측정으로 승격하지 않으며 새 영상·shape·photo-z 반출을 승인하지 않는다.

### Matched-control PSF profile diagnostic — 2026-10-20

기존 target/control coadd의 같은 각거리 위치에서 source-masked radial profile을 비교했지만, 외곽 bin 지원은 최대 12 프레임뿐이고 sky annulus에 민감했다. 대조 위치는 별이 없는 placebo이며 target의 ICL curvature도 남으므로 이 차이는 PSF wing 추정치가 아니다. 해석·제외 사유·출력 해시는 [matched-control profile 평가](empirical_psf_matched_control_profiles_assessment_2026-10-20.md)에 있으며, `psf_approved=false`, `science_measurement=false`다.

Gaia 경험 PSF의 epoch 범위를 보강하기 위해 333개 별 우표가 걸친 227개 target coadd의 모든 FITS HDU header를 검사했다. 표준 관측 DATE/MJD 키는 0/227이지만, 별 우표 중심은 Gaia 좌표가 아니라 HSC 검출 `X_IMAGE/Y_IMAGE`이므로 이 누락이 영상 중심정렬을 막지는 않는다. 이전 ±2/5/10년 Gaia 이동은 cross-match 민감도이지 실제 epoch 확인은 아니다. 정확한 Gaia 후보 식별, 변광·고립도, 포화 및 wing/sky 전달 검증은 남아 있다. 재현기·해시는 [coadd epoch 감사](empirical_psf_hsc_coadd_epoch_audit_2026-10-20.md)에 기록했다.

333개 후보의 포화·혼잡도·마스크 지지율을 실제 HSC detection 중심 및 저장된 stamp mask로 감사했다. 12픽셀 중심 source-mask 개방률 99% 문턱을 통과한 표본은 Kron 2.5 290개/201프레임, Kron 3.5 242개/179프레임이다. 하지만 이 중심 적격 별 중 각각 다수가 48–96픽셀 고리의 20% 이상을 mask로 가려 외곽 PSF wing 추정에는 부족하다. 12픽셀 내 다른 검출은 없지만 96픽셀 내 검출 이웃은 전원에게 있다. 따라서 기존 프로필 스택의 중심 정규화는 가능하나 확장 wing 보정은 승인할 수 없다. 전체 정량 결과와 해시는 [별 우표 품질 감사](empirical_psf_star_stamp_quality_assessment_2026-10-20.md)에 기록했으며 `psf_approved=false`, `science_measurement=false`를 유지한다.

약한렌즈 random-center 후보 count-only 자료도 전량 감사했다. 800개 렌즈 모두 각 렌즈당 계획 목표 200보다 많은 random point 후보를 갖고(최소 364, 중앙값 1,000), 원시 후보 고갈은 현재 주요 병목이 아니다. 후보 수의 렌즈 redshift 의존성은 per-lens 유효 annulus area로 나누면 Spearman `rho=0.595`에서 `0.080`으로 줄어든다. 다만 첫 50개 geometry pilot에서 중첩·중복 제거 후 6.35%만 남아 실제 제약은 모든 cluster aperture와의 기하 중첩 및 이후 영상/형상 지원일 수 있다. 승인 좌표 export는 2/32 batch에서 아직 진행 중이라 전체 표본의 잔존 지원은 미확인이다. 재현 코드와 입력·출력 해시는 [random-coordinate candidate count 평가](wl_random_coordinate_candidate_count_assessment_2026-10-20.md)에 있으며 추가 CAS 요청은 없었다.

경험 PSF 별 스택의 반경 지원을 HSC cluster 적색편이와 같은 `H0=70`, `Omega_m=0.3` 우주론으로 환산했다. 333개 우표에서 실제 쌓인 96픽셀 끝점은 중앙값 약 88 kpc이고 0/227 cluster frame이 ICL annulus 시작점 200 kpc까지 닿는다. 전체 우표 inscribed footprint도 200 kpc에 닿는 프레임이 52/227뿐이며 300 kpc에는 하나도 없다. 따라서 이 PSF 스택은 50–75 kpc satellite outskirts의 inner validation에는 일부 관련되지만 ICL 200–300 kpc 고리의 wing 보정값을 제공하지 않는다. 이 coverage gap을 닫기 전까지 PSF 보정된 ICL–satellite 상관은 과학 측정으로 승격하지 않는다. [물리 반경 감사](empirical_psf_physical_radius_coverage_assessment_2026-10-20.md)에 입력 provenance와 정량 결과를 보존했다.

### 기존 full-frame target/control 영상의 300 kpc PSF-wing 파일럿 — 2026-10-20

96픽셀 우표의 반경이 ICL 200–300 kpc보다 짧아, 이미 보유한 800개 표적 영상과
798개 매칭 대조 영상에서 Gaia 별 중심 프로필을 full-frame으로 측정했다. 외곽장,
양쪽 영상의 300 kpc 기하 지원, 동일 하위필터, HSC 검출 무결성, Kron 2.5 중심 개방률
조건으로 18개 별/프레임 쌍을 선별했다. Gaia variability는 모두 `NOT_AVAILABLE`이며
epoch 검증도 남아 있다.

25 kpc 고리와 기존 두 국소 sky annulus 설정으로 864개 중복 없는 프로필 행을
얻었다. 그러나 반경별 target-control 쌍은 최대 4개뿐이고 250–300 kpc에는 없으며,
200–225 kpc에는 1쌍이다. 프레임 단위 재표집과 독립 산술 검증은 통과했지만 이 지원은
ICL 외곽 PSF 보정이나 위성–ICL 분석에 부족하다. `_01`은 control 행 중복 버그가 있어
사용하지 않으며 `_02`가 수정 결과다. 판정은 `science_measurement=false`,
`psf_approved=false`다. 상세 방법·표본·해시는 [외곽장 PSF-wing 평가](empirical_psf_outer_field_wings_assessment_2026-10-20.md)에 있다.

2026-10-20 후속으로 여섯 ICL 고리와 일곱 위성 고리의 교집합 84개를 전량 계산하고,
420개의 target/control cluster-filter 쌍을 공통 5도 공간블록에 묶어 반경 가족검정을
했다. 같은 표본의 equal-satellite 및 CAMIRA column-8 weighted 프로필을 함께 포함한
168개 가족에서 5,000회 Webb wild-cluster bootstrap, Holm 및 max-T 보정 후 유의한 셀은
없었다. 주 외곽 조합(200–300 × 50–75 kpc)에서 i2 equal-weight `r=0.293`, raw
`p=0.030`; CAMIRA-weighted `r=0.241`, raw `p=0.055`였지만 168개 Holm 값은 모두
`≥1.000`이고 최소 max-T `p=0.255`였다. 이는 반경/필터 탐색과 estimator 선택을
보정한 뒤 현재 표본에서 연관성을 검출하지 못했다는 뜻이지, 연관 부재의 증명은 아니다.
CAMIRA `weight`도 불편한 독립 멤버십 확률이 아니며 선택편향을 교정하지 않는다. 미검증
PSF 날개·sky 전달·위성 완전도·coadd 공분산은 여전히 핵심 관문이라 `science_measurement`
상태를 바꾸지 않는다. 재현과 독립 검증은 [CAMIRA 가중 반경행렬 평가](icl_satellite_weighted_radial_matrix_familywise_assessment_2026-10-20.md)에
기록했다.

weak-lensing의 28성분 cross-shear null을 source-bin과 반경 shell로 분해했다. 전체
결합 불일치는 유의(χ²=74.75/28, p=3.93×10⁻⁶)하지만, 각 shell의 4-bin 검정은
7개 반경 보정 후 모두 비유의(최소 Holm p=0.102)이고, 각 source-bin의 7-shell
검정도 4-bin 보정 후 비유의(최소 Holm p=0.310)다. 28개 marginal 성분도 최대
|residual/error|=2.11, 최소 Holm p=0.966으로 특정 요소를 고립하지 못했다. 이는
전체 null을 해소하지 않으며 random-lens 차감 전 체계오차 경고로 남긴다. 상세와
재현 해시는 [cross-shear null 분해](wl_tomographic_consistency_assessment_2026-10-20.md#후속-cross-shear-null-반경source-bin-분해)에 기록했다.

같은 endpoint(200–300 × 50–75 kpc)와 같은 공간블록에서 CAMIRA 가중 estimator 변화량을
10,000회 paired bootstrap으로 직접 비교했다. equal 대비 weighted 변화는 i에서 −0.003
(95% 구간 −0.225…0.092), i2에서 −0.052 (−0.139…0.061)로 둘 다 0을 포함한다.
따라서 가중으로 점 추정이 낮아지는 방향은 관측되지만 estimator 차이가 분리됐다고
말할 수 없다. 이 endpoint는 앞선 탐색에서 선택됐으므로 확증 검정으로 해석하지 않으며,
`science_measurement=false`다. 상세 및 재현 해시는 [paired CAMIRA 가중 민감도](icl_satellite_weighted_radial_matrix_familywise_assessment_2026-10-20.md#동일-표본의-estimator-차이-paired-민감도--2026-10-20)에 있다.

2026-10-21 PSF 후속에서는 600 kpc 기준 91개 별 쌍/78개 프레임을 측정했다. 첫 집계가
`i`와 `i2`를 혼합한 사실을 발견해 해당 `_01` 출력을 폐기하고, 4,368 프로필 행을
필터별 48개 cell로 다시 집계한 `_02`만 남겼다. 독립 validator가 저장된 프로필과
10,000회 frame bootstrap을 재현했다. `_02`에서 `i`는 29쌍/25프레임, `i2`는
62쌍/53프레임이다. 엄격 기준으로 `i`의 50–75 kpc 지원은 0프레임, `i2`는 4–6이며,
ICL 200–300 kpc에서는 `i2` 지원이 0이고 `i`도 250 kpc 위는 0이다. 20% 면적/6개
sector 완화 후에도 필터별 outer-bin 최대 지원은 1–6프레임뿐이다. 따라서 넓은 표본도
PSF 보정을 허용하지 못한다. 상세·해시는 [필터분리 PSF 지원 평가](empirical_psf_outer_frame_bootstrap_assessment_2026-10-21.md)에
기록했고 `psf_approved=false`, `science_measurement=false`를 유지한다.

사용자 승인으로 공개 HSC PDR3 확장 PSF의 HSC-I 대칭·비대칭 모델만
Zenodo v4.0.0에서 받았다(2개, 19,618,560 bytes; 각 크기와 MD5 및 SHA-256 검증).
4001² FITS 모델은 5.6′까지 뻗어 기존 별 우표가 도달하지 않는 ICL 반경의
점광원 누출 기준을 제공한다. `H0=70`, `Omega_m=0.3`에서 점광원 총광량의
0.222–0.242%가 200–300 kpc 고리에 들고, 대칭/비대칭의 원형 고리 광량은
거의 같다. 독립 validator가 입력 체크섬과 적색편이별 물리 고리 적분을
재현했다. HSC PDR3는 i/i2 노출을 공합성하므로 public HSC-I는 coadd-level
참조로 유망하지만, local frame PSF·mask·sky transfer 검증은 미완료다.
그 43×43 코어는 기존 위치별 PSF Picker 코어 중앙값보다 21.5% 넓고 `R90`도
19.1% 크므로, 전체 public PSF를 위치별 모델 대신 쓸 수 없다. 외곽 wing과
위치별 core의 1.5–3.3″ 중첩 비교에서는 744개 지점의 scale 중앙값 0.739
(5–95백분위 0.467–1.045), 모양 잔차 RMS 중앙값 12.6%를 얻었다. i와 i2는
각각 158/586 지점으로 따로 요약했다. 이 범위는 접합 검토의 출발점이지
3.3″ 밖 날개의 국소 정규화 증거는 아니며, 위치별 core를 접합하는 감도
분석이 다음 검증 단계다. 따라서 현 결과는
PSF 보정이나 ICL 측정이 아니며 `psf_approved=false`, `science_measurement=false`를
유지한다. 데이터와 감사 해시는 [공개 확장 PSF 평가](extended_psf_public_reference_assessment_2026-10-21.md)에,
도구는 [`audit_public_extended_psf_reference.py`](audit_public_extended_psf_reference.py),
[`validate_public_extended_psf_reference_audit.py`](validate_public_extended_psf_reference_audit.py),
[`compare_public_extended_psf_to_picker_cores.py`](compare_public_extended_psf_to_picker_cores.py),
[`audit_public_psf_picker_core_overlap.py`](audit_public_psf_picker_core_overlap.py)에 있다.

공개 HSC-I 확장모델을 i2 프레임에 프록시로 적용한 주입 민감도도 추가했다.
승인된 core PSF 위치가 있는 44개 i2 프레임에 대칭·비대칭 모델을 각각 적용했고,
44개 전부 오류 없이 완료됐다. 독립 검증기가 저장된 transfer 응답을 시나리오 값에서
재계산해 두 변형 모두 최대 산술 오차 0을 확인했다. 매끈한 확산 주입의 ICL 고리
응답 중앙값은 0.224, 위성 외곽광 주입은 0.485였지만 프레임별 범위가 각각
0.004–0.318, −0.047–1.424로 넓다. 두 공개 모델의 paired 차이는 이 반경 측정에서
최대 절댓값 9.54×10⁻⁷로 사실상 무시 가능하지만, 이는 local i2 PSF 일치나 wing 보정을
검증하지 않는다. 따라서 프록시 감도 결과는 추정 보정계수로 쓰지 않고
`science_measurement=false`를 유지한다. 입력·결과·검증 해시와 제한은 [공개 PSF 평가의
i2 감도 절](extended_psf_public_reference_assessment_2026-10-21.md#exploratory-hsc-i-proxy-sensitivity-on-i2-frames)에,
재현기·검증기는 [`extended_psf_injection.py`](extended_psf_injection.py),
[`validate_extended_psf_injection_batch.py`](validate_extended_psf_injection_batch.py),
[`compare_extended_psf_injection_variants.py`](compare_extended_psf_injection_variants.py)에 기록했다.

### 전체 표적·대조 PSF 코어와 공개 확장 PSF 중첩 평가 — 2026-10-21

새로 확보한 코어를 포함한 800 target/798 control, 총 1,598개 위치에서
공개 HSC-I 확장 PSF와 위치별 코어의 1.5–3.3 arcsec 프로필 진폭·모양을
비교했다. 표적은 기존 stack의 `i`/`i2`, 대조장은 각 대조 manifest의 실제
`coadd_header_filter`를 사용했다. 정확히 같은 필터인 783쌍(i 302, i2 481)만
paired 비교했고 다른 필터 15쌍은 제외했다. target-control full-window scale
차이 중앙값은 i −0.0239, i2 −0.0033이며, 진폭 scale은 두 필터 모두 1.5–2.4에서
2.4–3.3 arcsec로 감소한다. 따라서 단일 amplitude로 공개 외곽 PSF를 붙일 근거가
없고, 넓은 물리 반경의 wing과 sky 전달도 여전히 검증되지 않았다. 초기 `_01`
출력은 대조장 필터를 표적에서 상속해 paired 비교에서 폐기하고, 실제 필터를
사용한 `_02`만 권위 결과로 삼았다. Astrafek 분석기와 독립 FITS 재계산기가
1,598개 scale을 재현했고 7개 focused test가 통과했다. 상세 방법·결과·입력 및
산출물 해시는 [전체 표적/대조 공개 PSF 중첩 평가](extended_psf_target_control_overlap_assessment_2026-10-21.md)에
기록했다. 이 진단은 보정이나 측정이 아니며 `psf_approved=false`,
`science_measurement=false`를 유지한다. 다음 검증은 기존 target/control mask에서
splicing 반경·진폭을 바꾼 점광원 및 diffuse-profile 주입회수 민감도다.

### HSC `i` 필터 공개 PSF 주입회수 민감도 — 2026-10-21

동일한 13개 support-selected `i` frame에서 공개 HSC-I 대칭·비대칭 모델을
각각 적용하고 실제 source/ICL/satellite mask 재검출을 포함해 주입회수했다.
독립 validator가 각 run 13/13 응답 산술을 오차 0으로 재현했고, paired 모델
차이는 최대 4.77×10⁻⁷로 사실상 무시 가능했다. 그러나 대칭 모델의 smooth diffuse
→ ICL 200–300 kpc 회수 중앙값은 0.215 (프레임 범위 −2.734–0.271), 위성 외곽광
→ 위성 추정량은 0.294 (−3.024–1.895)로 프레임 간 응답이 넓고 음수도 포함한다.
따라서 단일 transfer correction은 지지되지 않는다. 이는 공통 공개 PSF의
모델민감도이며 위치별 PSF, sky transfer, source photon noise 또는 coadd
공분산을 검증하지 않는다. 결과는 `science_measurement=false`, `psf_approved=false`.
상세 수치·해시·제한은 [i-filter 주입회수 평가](extended_psf_i_injection_recovery_assessment_2026-10-21.md)에
있다. `compare_extended_psf_injection_variants.py` 비교기를 i/i2 선택형으로
일반화하고 8개 focused tests를 통과시켰다.
