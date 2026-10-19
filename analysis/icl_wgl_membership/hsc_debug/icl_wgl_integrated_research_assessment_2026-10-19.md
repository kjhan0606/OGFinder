# HSC ICL–위성 연관성 및 배경은하 렌즈 연구 현황 — 2026-10-20 갱신 (초판 2026-10-19)

## 현재 판정

Astrafek에는 HSC PDR3 영상·대조장·PSF 입력을 사용한 ICL 스택, 위성 외곽광 연관성 분석, 배경은하 약한렌즈 스택과 각 단계의 재현 스크립트가 구현되어 있다. 지금 결과는 모두 탐색 단계다. ICL–위성 관계는 다중비교와 sky 처리 민감도를 넘지 못했고, 렌즈 B-mode 마스크 방향은 외부 영역과 좌표 대조로 검증되었다. 렌즈 신호에는 아직 무작위 중심 차감과 대응 matched-field 공분산이 없다. 따라서 두 갈래 모두 `science_measurement=false`를 유지한다.

## 보유 자료와 실행 이력

- HSC 표본은 CAMIRA 은하단 800개다. 표적 영상 800개와 확보된 매칭 대조장 798개에서 동일 필터로 쌓을 수 있는 표적–대조장 쌍은 783개(i 302, i2 481)다. HSC PDR3는 `i`/`i2` 응답곡선이 비슷하지만 동일하지 않고 적색 천체에서는 보정 차이가 더 커질 수 있다고 명시하므로 두 필터는 섞지 않는다. 픽셀별 색이 없는 ICL 영상에 객체별 유효응답 보정을 추정해 적용하지 않고, 현재 Astrafek 스택은 각 target/control 쌍의 동일 필터만 비교해 필터별로 보고한다 ([HSC PDR3 처리 설명](https://hsc-release.mtk.nao.ac.jp/doc/index.php/processing__pdr3/), [필터 응답 차이](https://hsc-release.mtk.nao.ac.jp/doc/index.php/filterfrac-2/)).
- ICL 분석의 공통지원 2-D 스택과 K=2.5, 3.5, 5.0 독립지원 비교는 Astrafek의 고정 스크립트로 생성·검증되었다. 독립지원 추정치는 픽셀 교집합의 영향을 줄이는 민감도 분석이며 기본 측정값으로 채택하지 않았다.
- 2026-10-20에 Astrafek 독립 검증기를 보완해 세 마스크 설정의 스택 원자료 연결, 쌍별 프로파일, 출력 지도 차감식, 마스크 정의 및 실행 코드 해시를 재검증했다. 검증 산출물은 각 스택 폴더의 `independent_verification_2026-10-20.json`이며, 이 재검증은 과학 측정 게이트를 바꾸지 않는다.
- 재검사한 표적–대조 2-D 스택의 50–100 kpc 표적 초과 평균은 i에서 `41.39 ± 2.90`, i2에서 `44.59 ± 2.32` nJy/arcsec²이고, 100–200 kpc에서는 각각 `12.92 ± 1.49`, `10.75 ± 1.01`이다. ±는 은하단 쌍 산포에서 계산한 기술적 SEM이며 전체 공간 공분산·sky·PSF 계통오차가 아니다. 해당 고리의 중앙 공통 픽셀 지원률은 i 15%/14%, i2 19%/18%로 제한적이다. 스택 요약 SHA-256은 `9c804cfb12ff2b86a5bdc8fdbfea94aa1dc473711e6f7ac6f4d233cc6b853173`이며 결과는 계속 `science_measurement=false`다.
- 위성 연관성에는 반경 구간별 94쌍(i 26, i2 68)이 있다. 주입·회수 검증은 위치별 HSC-I `pdr3_wide` coadd 모델 PSF 744개와 승인된 시험 PSF 1개를 사용해 57개 중심 및 688개 위성 위치에서 실행되었다. PSF 감사 매니페스트 SHA-256은 `6375bc58df6d5a15d837cf5b13d3db79ebdc6944dbf7510164f465b6cbe72f7e`다. 이 약 7 arcsec 코어는 넓은 ICL PSF 날개를 검증하지 않는다.
- 기존 영상의 8,078개 compact 후보와 Gaia 진단은 target-specific extended PSF 날개 표본을 제공하지 못했다. Gaia 대조는 569개 후보 위치에서 73개 행을 반환했고, 32개 별 후보 중 전체 네 배경/마스크 조건 통과는 4개였다. 전체 미조회 후보 7,509개의 Gaia DR3 위치 대조 계획을 오프라인 준비했지만 외부 좌표 업로드와 조회는 승인 전이므로 제출하지 않았다. 계획·전송 범위·체크섬은 [Gaia PSF 후보 확장 계획](extended_psf_README.md)에 기록했다.
- 약한렌즈 표적 쪽은 CAS 원자료 3,069,544 lens–source 행(800 렌즈)을 가져왔다. 품질·거리·반경 컷 뒤 2,153,004 쌍, 609 렌즈, 2,051,991개 고유 배경원으로 줄었다. 기존 결과는 표적 중심 스택만 포함한다.

## ICL과 위성 외곽광

독립지원 비교에서 50–200 kpc 구간의 최대 같은 쌍 추정치 이동은 K=5, i2의 175–200 kpc에서 32.9%였다. 따라서 마스크 지원 방식이 외곽 프로파일의 크기를 실질적으로 바꿀 수 있다. 94쌍 분석에서 i2 전체지원 부분상관은 `r≈0.293`이지만 2.5°, 5°, 10° 공간 블록 폭에서 신뢰구간 해석이 달라진다. 별도 가족단위 검정은 48개 설정 중 추론지원 기준을 통과한 30개 어디에서도 Holm 또는 max-T 보정 뒤 유의성을 보이지 않았다.

2026-10-20 원자료 해시를 다시 확인한 표본 보존 감사에서 독립 50–75 kpc 위성 스택은 i 70개/42개 은하단 `2.51 [-0.30, 5.50]`, i2 308개/140개 `1.19 [-0.24, 2.49]` nJy/arcsec²로, 두 필터 모두 은하단 부트스트랩 구간이 0을 포함했다. 결합된 ICL–위성 표본은 94개 cluster/filter 쌍(i 26, i2 68)이다. i2 부분 Pearson `r=0.293`의 5° 공간 블록 구간 `[0.042, 0.497]`은 양수지만, 일반 쌍 부트스트랩 구간 `[-0.007, 0.526]`과 Spearman 공간 블록 구간 `[-0.088, 0.301]`은 0을 포함한다. i 부분 Pearson은 `r=0.079`, 공간 블록 구간 `[-0.186, 0.479]`다. 따라서 선택된 i2 부분상관은 흥미롭지만 순위 통계와 재표집 방식에 민감하며, 검출로 주장할 수 없다. 재현 가능한 해시 연결 감사는 `icl_satellite_sample_retention_2026-10-20_02/summary.json`에 있으며 SHA-256은 `7952f53bd29289dc182b24abdd7b0356b9a09e40c4e15055a0ea1f6aa1ca354f`다.

연결 target/control 블록 하나씩을 제외하는 추가 영향도 검사에서는 i2 `r`이 세 블록 폭 모두에서 양수였으나 `0.166–0.177`까지 낮아졌다. i는 2.5°와 10° 폭에서 한 블록 제외 시 부호가 바뀌었다. 이는 표본 내 단일 블록 영향 진단일 뿐 신뢰구간이나 검출 검정이 아니며, 기존 다중비교 결과를 바꾸지 않는다. 상세 결과와 영향 블록 ID는 [블록 제외 영향도 평가](satellite_icl_leave_one_block_out_assessment_2026-10-20.md)에 기록했다.

94쌍의 부분상관에 patch-QA 연결영역 고정효과를 추가해 영역 간 평균 오프셋 민감도를 점검했다. i2는 `r=0.29264`에서 `0.29105`로 유지됐고, i는 `0.07898`에서 `0.13860`으로 바뀌었다. 각 target/control 쌍의 양 끝점이 이 근사 영역 분할에서 항상 같은 연결요소에 속해 두 효과는 선형 중복이며 한 번만 적합했다. 결과는 지역 평균 오프셋만으로 i2 상관이 생기지는 않았음을 시사하는 탐색 진단이다. 소표본·선택 민감성, sky/PSF·마스크 계통효과와 공간 공분산은 해결되지 않았고 `science_measurement=false`를 유지한다. 재현 및 제한은 [영역 고정효과 평가](satellite_icl_region_fixed_effects_assessment_2026-10-20.md)에 있다.

동일 표본의 은하단별 leave-one-out 영향도에서는 전체 지원 i2 `r`이 항상 양수였지만 `0.151–0.338`까지 변했고, i는 한 은하단 제외에서 부호가 바뀌었다. 중심 적격 부분표본은 i2 15쌍, i 8쌍뿐이며 공간 블록 구간이 모두 0을 포함한다. 따라서 현재 상관은 검출 근거가 아니다. 이는 영향도 검사이지 신뢰구간·검출 검정이 아니다. 세부 결과는 [cluster 영향도 평가](satellite_icl_cluster_influence_assessment_2026-10-20.md)에 기록했다.

redshift·richness 조정 순위 부분상관도 함께 확인했다. 전체 지원 i2는 0.139, i는 0.248로, Pearson 계수와 크기 차이가 있고 i의 단일점 영향도도 달랐다. 이는 이 선택 표본에서 관계 형태가 상관 정의에 민감함을 보여주며, 독립적인 검출 근거는 아니다.

같은 patch-QA 트랙 연결요소로 matched target–control ICL 프로파일도 점검했다. 783쌍 중 세 영역에 487, 244, 52쌍이 분포하고 1.6–3.0° 인접 기준에서 분할은 안정적이다. 가장 큰 영역을 제외하면 50–100 kpc 초과광은 i에서 1.8%, i2에서 9.2% 감소했고, 100–200 kpc는 i에서 10.2% 증가, i2에서 9.8% 감소했다. 200 kpc 바깥 저신호 프로파일은 영역 제외에 따라 부호가 바뀌는 고리가 있다. 50–200 kpc의 pair별 중앙 공통 픽셀 지원은 14.5–19.0%에 그친다. 이는 저지원 외곽광의 지역 민감도를 보여줄 뿐, 신뢰구간이나 ICL 검출 검정이 아니다. [ICL tract-region 영향도 평가](icl_hsc_tract_region_influence_assessment_2026-10-20.md)에 재현 스크립트, 데이터 매칭, 해시, 전체 고리 결과를 남겼다.

같은 HSC 트랙 연결요소로 위성–ICL 상관 자체도 재계산했다. 94행의 target과 matched control 중 어느 한 점이라도 제외 영역에 있으면 그 행을 빼고 redshift·richness 회귀를 다시 적합했다. 전체 i2 partial Pearson `r=0.293`은 영역 1을 빼면 `r=0.048`로 감소하고, Spearman `rho`는 `+0.111`에서 `−0.048`로 바뀐다(68행 중 29행 제외). 따라서 앞서 직사각 공간 셀 삭제에서 양수였던 partial Pearson 부호만으로 연관성이 지역에 강건하다고 볼 수 없다. i는 전체 26행의 `r=0.079`이고 영역 1 삭제 후 8행만 남아 `r=−0.378`이지만, 이 잔여 표본은 너무 작아 부호를 해석하지 않는다. 세 영역 영향도는 공분산·신뢰구간이 아니며, 관계의 원인이 천체물리적 환경인지 선택·영상 계통인지 구분하지 못한다. 상세 결과는 [위성–ICL tract-region 영향도](satellite_icl_tract_region_influence_assessment_2026-10-20.md)에 기록했다.

은하수 소광을 픽셀별로 적용했을 때 기존 중심 소광 보정과의 ICL 차이는 사실상 0이었다. 반면 국소 평면 제거는 고정된 74쌍에서 i2 부분상관을 0.239에서 0.124로 낮추고 구간이 0을 포함하게 했다. 평면 주입 시험에서는 250 kpc 규모의 주입광 중 약 9–14%가 측정에서 빠졌지만, 해당 시험은 검출 마스크를 고정했고 상류 HSC sky 제거를 재현하지 않았다. 공개 확장 PSF 참조를 적용한 i 대조군 13개 시험에서는 단위 diffuse 주입의 중앙 21.5%가 ICL 측정 고리에 들어오고, 위성 외피는 35.3%가 들어왔다. 이 값은 특정 입력·마스크의 전달 민감도이지 관측값에 곱할 보정계수가 아니다.

공분산 감사는 position별 PSF 주입 표본 57개(i 13, i2 44)의 표적·매칭 대조 coadd에서 9×9 단거리 잔차 커널을 측정하고 실제 signed aperture 가중치에 전파했다. 큰 top-hat 근사 분산 비는 표적 1.59, 대조장 1.54–1.55였다. 정확한 ICL 및 위성 가중치와 양쪽 ICL 영상의 분산을 사용하면 대각 분산만 쓴 경우보다 ICL 잡음 표준편차가 i에서 1.234배, i2에서 1.240배였고, 잡음만으로 계산한 ICL–위성 상관 중앙값은 각각 −0.0073, −0.00091이었다. 이는 해당 57개 프레임에 조건부인 단거리·정상 커널 민감도다. source Poisson 잡음, 상류 sky 전달, 위치별 완전 공분산은 포함하지 않으며 공간 블록 오차에 단순 곱하면 같은 변동을 중복 셀 수 있다. 따라서 관측 연관성을 보정하거나 검출 여부를 결정하는 수치로 사용하지 않는다.

57쌍의 잡음 공분산을 같은 표본의 측정 ICL·위성값으로 전파한 10,000회 Monte Carlo에서는 i2 부분표본 44개의 측정 편상관 `r=0.409`가 finite-lag 잡음 범위 `[0.320, 0.465]`를 보였다. 하지만 전체 68쌍 i2 추정값은 `r≈0.293`이고 공간 블록 구간은 블록 폭에 민감하다. 공분산이 있는 표본만의 값이 전체보다 큰 점은 선택 효과를 경고한다. i 부분표본은 finite-lag 잡음 범위가 0을 포함한다. 이 조건부 잡음 실험은 spatial-block 범위와 합치지 않으며 과학 검출을 뜻하지 않는다. 세부 결과는 [측정 잡음 Monte Carlo 평가](icl_satellite_measurement_noise_mc_assessment_2026-10-19.md)에 기록했다.

현재 CAMIRA 후보를 더 넓히는 오프라인 검토에서는 기존 선택 아래 새로 분리되는 후보가 57개뿐이고 새 5도 하늘 셀은 없었다. 적색이동 범위를 넓혀도 분광 확인을 유지하면 새 5도 셀은 2개였다. 따라서 같은 발자국에서 영상 수만 늘리는 것으로 독립 하늘 영역의 부족을 해결하기 어렵다.

**해석:** 양의 ICL–위성 관계는 흥미로운 가설이지만 검출로 보고할 수 없다. 다음 ICL 관문은 반경·필터·마스크·통계량을 분석 전에 고정하고, 실제 검출·마스크 재생성을 포함하는 diffuse+위성 공동 주입으로 sky·PSF 날개·마스크 누출 전달함수를 평가하는 것이다. 독립 통계력을 높이려면 단순히 기존 CAMIRA 선택을 완화하기보다 새 독립 하늘 영역과 완전도·렌즈 보정을 함께 설계해야 한다.

## 배경은하 약한렌즈

기존 표적 전용 `DeltaSigma` 프로파일은 `0.1–0.2 Mpc`에서 약 `1.61×10^14 M☉ Mpc⁻²`, `2–3 Mpc`에서 약 `1.53×10^13 M☉ Mpc⁻²`로 감소한다. 52개 공간 셀 삭제 공분산은 독립 쌍 SEM보다 낫지만 셀이 연결된 HSC 관측영역과 일치하지 않을 수 있다. 이 표적 전용 프로파일에는 무작위 중심 차감, 대조장 분산, `n(z)`·shear 보정 불확실성, reduced-shear 보정이 빠져 있어 검출 유의도나 최종 질량으로 해석하면 안 된다.

2026-10-20 네 HSC Y3 source-bin 약한렌즈 profile을 동일 52개 sky-cell jackknife로 비교했다. 0.1–0.7 Mpc의 세 구간에서 `DeltaSigma`가 source bin에 따라 달랐고, 반경 일곱 구간 Holm 보정 뒤에도 GLS 일관성 검정이 남았다. 이 결과는 SACC `n(z)`와 source selection/calibration의 bin별 transfer를 추가 점검해야 함을 보인다. cluster-only 진단이며 random 차감, `n(z)` 불확실성 및 연결 관측영역 공분산을 포함하지 않으므로 질량 또는 렌즈 검출로 해석하지 않는다. 세부 프로파일과 공분산은 [tomographic 일관성 평가](wl_tomographic_consistency_assessment_2026-10-20.md)에 기록했다.

이어 같은 SACC 입력에서 source-bin 3·4의 평균 적색이동을 DESI 논문의 주변 68% 이동 구간 중앙값과 네 모서리로 평행 이동한 다섯 시나리오를 실행했다. 독립 검증기는 다섯 56×56 공분산 모두의 대칭성과 순위 51, 입력 체크섬, 요약 envelope를 통과시켰다. bin 간 일관성 불일치는 모든 시나리오에서 0.1–0.7 Mpc에 남았고, 반경 일곱 구간 Holm 보정 p값은 각각 0.1–0.2 Mpc에서 `7.6–8.6×10⁻⁵`, 0.2–0.4 Mpc에서 `9.2×10⁻⁹–1.7×10⁻⁸`, 0.4–0.7 Mpc에서 `7.4×10⁻⁴–1.04×10⁻³`였다. 이는 이 근사 이동 범위로 불일치가 해소되지 않음을 보이며 실제 보정 posterior 검증의 우선순위를 높인다. 상수 이동은 DESI의 보정 `n(z)` posterior 전파가 아니고, cluster-only 프로파일에 random 차감 및 연결 관측영역 공분산도 없으므로 렌즈 검출이나 질량 제약이 아니다. 상세 수치와 파일 해시는 [n(z) 이동 격자 평가](wl_published_nz_shift_grid_assessment_2026-10-20.md)에 기록했다.

같은 승인 lens–source 행에서 source object ID 재사용의 영향을 5도 52개 블록으로 평가했다. 선택된 2,153,004쌍 중 고유 source 2,051,991개 가운데 99,671개(4.86%)가 여러 lens 쌍에 재사용됐고, 그 source가 차지하는 행은 200,684개(9.32%)였다. lens endpoint, source endpoint, 양 endpoint 공동 삭제 공분산은 모두 56차원·순위 51이었다. joint-endpoint 대 lens-only 프로파일 오차비의 중앙값은 tangential 1.008(범위 0.973–1.046), cross 1.006(0.988–1.046)으로, 이 블록 구성에서 source 재사용이 marginal error를 크게 바꾸지는 않았다. 세 삭제 방식 모두 source-bin 일관성 Holm p값이 0.1–0.7 Mpc에서 `6.0×10⁻⁹–4.7×10⁻⁴`였고 더 큰 반경에서는 유의하지 않았다. 이는 source-bin별 calibration/`n(z)` 차이를 가리키는 진단이지 렌즈 검출이 아니다. Astrafek 독립 검증기 [`validate_wl_source_reuse_covariance.py`](validate_wl_source_reuse_covariance.py)가 대칭성·PSD·순위·프로파일 오차·GLS/Holm·재사용 집계를 검사해 통과했다. 요약 SHA-256은 `ac578d5818aae383586b62855e293bbfb5a1427e8a53e252d2569e6e4302b428`, 검증 JSON SHA-256은 `3f565f508a40ca1d98aef4d7f58b1152303694546811bb2cee0602801c870fd0`이다. 작은 음의 최소 고유값은 최대 행렬 원소 대비 `8×10⁻¹⁷` 이하로 수치 반올림 수준이며 PSD 상대 허용치는 `10⁻¹²`다. 이 민감도는 source-level 전체 공분산, 연결된 HSC 관측영역 공분산, random subtraction, 보정 `n(z)` posterior를 대체하지 않으며 `science_measurement=false`를 유지한다.

후속 문헌 검토에서 2026-08-25 개정된 HSC PDR3 S19A 형상표본 전용 DESI DR1/DR2 클러스터링 적색이동 보정이 확인됐다. 이는 현재 사용 중인 S19A 모양 카탈로그와 직접 같은 표본을 보정하며, 기존 HSC Y3 보정에 없던 `z>1.2` ELG/QSO 추적자를 포함한다. 논문은 1번 bin이 약간 낮은 적색이동 쪽으로, 2번은 광도측정 보정과 일치, 3·4번은 기존 Rau+2022 결과보다 높은 적색이동 쪽으로 이동한다고 보고한다. 3·4번의 상대 이동은 논문 정의로 각각 `Δz₃=−0.043^{+0.022}_{−0.020}`, `Δz₄=−0.050±0.012`다. 이 결과는 현재 source-bin 불일치의 보정 민감도를 다시 계산해야 한다는 근거이지, 불일치를 설명하거나 질량 신호를 확정하는 근거는 아니다. 새 posterior `n(z)` 결과 파일은 저자 공개 `desi-y3-hsc/results/distributions`에 있다고 논문 코드 저장소 README가 명시한다. 본 저장소나 데이터 폴더로 그 파일을 아직 가져오지 않았다. 다음 단계에서는 사용자 승인 뒤 해당 보정 결과만 수신하고, 파일 형식·라이선스·체크섬을 확인한 뒤 SACC 기준 결과와 동일한 paired jackknife에서 비교한다. 출처와 범위는 [적색이동 보정 업데이트](wl_pdr3_desi_nz_calibration_update_2026-10-20.md)에 기록했다.

약한렌즈 adapter의 shear responsivity, 개별 `m/c` 보정, HSC 경계 선택편향 계수를 논문과 HSC가 연결한 S19A 보정 코드에 대조했다. 소스 수준 수식 검토에서 구현은 일치한다. 이는 추정기 공식 검증일 뿐 새 `n(z)`, 표적–random 차감, 대응 공분산 검증은 아니다. 특히 토모그래픽 내부 bin 불일치가 남아 있어 현재 프로파일은 여전히 탐색 결과다. 세부 대조와 endpoint 주의사항은 [S19A 보정식 검토](wl_s19a_shape_calibration_formula_review_2026-10-20.md)에 기록했다.

HSC 문서는 `b_mode_mask`가 큰 B-mode 구역을 제거하는 데 쓰인다고 설명하지만 0/1 방향은 적지 않는다. Astrafek은 형상 품질 및 tomography 컷 뒤 좌표를 HSC Y3가 제외한 GAMA09H 사각영역(RA 132.5–140°, Dec 1.6–5°)과 직접 대조했다. 안쪽의 77,021개 lens–source 행은 모두 raw 0이고 바깥쪽 2,157,489개는 모두 raw 1이었다. shear 결과를 보지 않고 raw 0을 제외값으로 판별했다. 이 극성을 사용한 HSC Y3 four-bin SACC `n(z)`, 개체별 shear 보정, selection-bias 보정을 적용한 cluster-only 스택은 2,153,004쌍과 609개 렌즈를 유지했다. 52개 연결 sky 셀 jackknife는 별도로 계산됐지만, matched-random 차감과 source-overlap 보정은 아직 없다. 이 해석은 HSC가 공개한 마스크 설명과 Y3 시스템atics 연구의 제외영역 정의에 근거한다 ([HSC shape-catalog 문서](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s19a-shape-catalog-pdr3/), [HSC Y3 PSF systematics 연구](https://academic.oup.com/mnras/article/525/2/2441/7202336)).

52개 균등 셀 잭나이프가 연결된 HSC 영역 변동을 대표하는지 확인하기 위해 patch-QA 트랙 중심 연결성으로 leave-one-region 영향도도 계산했다. 트랙 중심 인접 임계값 1.6–3.0°에서 전체 840개 트랙은 네 연결요소, 표본은 세 요소(렌즈 499, 246, 55개)로 안정적으로 나뉘었다. 가장 큰 요소를 뺄 때 표적 `DeltaSigma`는 0.1–0.7 Mpc에서 최대 16.4% 이동했고, 가장 안쪽 cross 성분은 `+14.3 × 10^12 M_sun Mpc^-2` 이동했다. 이는 표본의 지역 영향이 무시할 수 없다는 점을 보여주는 영향도 진단이지 공분산이나 신뢰구간이 아니다. 세 영역만으로는 안정적인 전체 공분산을 추정할 수 없고 최근접 트랙 중심 및 근사 인접 그래프의 경계 불확실성도 남는다. 입력·요약 해시 및 재현 방법은 [tract-region 영향도 평가](wl_hsc_tract_region_influence_assessment_2026-10-20.md)에 있다. 따라서 기존 52셀 결과와 함께 보되 무작위 중심 차감 및 대응 field 공분산 관문을 대체하지 않는다.

무작위 중심 후보 좌표의 사용자 승인 CAS 내보내기는 800 렌즈, 예상 780,191 행, 32개 배치 계획이다. 2026-10-20 16:57 KST SSH 확인에서 lageunha 제출 스크립트 PID `1049265`가 약 29시간 49분째 살아 있었다. GPFS 체크포인트는 2/32 배치(33,514행) 완료, 세 번째 CAS 작업 `1062312` 제출 상태이며 파일 수정 시각은 프로세스 시작 5초 뒤인 2026-10-19 11:07:53 KST에 멈춰 있다. 체크포인트의 `error_type=HTTPError`는 현재 실패로 단정할 수 없다. 재개 코드는 기존 상태에서 전역 오류 필드를 지우지 않고 `running`으로 바꾸며, 매 상태 폴링 결과를 파일에 저장하지 않는다. 따라서 현재 CAS 상태와 시작 후 배치 진행 여부는 확인되지 않았다. 살아 있는 제출기를 유지하고 중복 제출하지 않는다. `1062312`는 Slurm 작업 ID가 아니므로 `squeue`로 확인할 수 없다. 전체 기하 스크리닝은 좌표 내보내기가 32/32 완료되고 각 배치의 행 수·해시가 검증된 뒤에 시작한다.

별도 오프라인 감사는 ICL용 HSC 무작위 대조장 798곳을 약한렌즈 random 중심으로 재사용할 수 있는지 평가했다. 768곳이 CAMIRA 3 Mpc 측정 조리개와 겹쳤고, 남은 30곳은 모두 서로 3 Mpc 조리개가 분리되어 있었다. 따라서 이 기존 대조장 표본만으로는 800개 렌즈의 대표 random stack을 만들 수 없다. 감사는 이미지 및 target manifest와 hash-verified 7,319개 CAMIRA 목록만 사용했으며 새 HSC query나 shape 행 반출은 하지 않았다. 재현기는 [`audit_wl_matched_control_centers.py`](audit_wl_matched_control_centers.py)이고, 결과는 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_icl_matched_control_geometry_2026-10-20_04/summary.json` (SHA-256 `12eec246699c24a6f194e271758bfdc09d85660b1eeedcc73a54c7aec1109cc9`)에 있다. 이 30곳도 shape-source footprint와 선택함수를 검증하기 전에는 과학용 random 중심이 아니다.

현재 인증된 좌표 내보내기 프로세스가 CAS를 계속 확인 중이다. 그 프로세스를 유지하고, 32개 배치가 모두 행 수·해시 검증을 통과한 뒤 Astrafek 전체 기하 스크리너를 실행해 CAMIRA 3 Mpc 중첩·중복 좌표·렌즈별 지원량을 검사한다. 그 뒤 전체 조리개 영상 범위, 마스크, 깊이, PSF, sky, 소광과 공간 블록을 맞춰 무작위 중심을 선택한다. 모양 카탈로그 행 수와 전송량을 count-only로 산정한 뒤 무작위 원천의 shape 행을 반출하려면 별도 사용자 승인을 받아야 한다. 그 데이터 없이는 표적–무작위 차감과 대응 공분산을 완료할 수 없다.

교차 성분의 1.5–2.0 Mpc 구간은 `3.1 ± 1.6 M☉ pc⁻²`로 개별 bin에서 약 2σ 벗어나므로, 무작위 차감과 대응 공분산을 적용한 뒤 null-test 대상으로 확인한다. 다른 bin은 이 값에 맞춰 재선택하지 않는다.

## Astrafek 재현 경로와 근거 문서

실행기, Slurm 스크립트, 검증기, 요약 생성기는 `analysis/icl_wgl_membership/hsc_debug/`에 보관한다. 대형 HSC 파일과 실행 결과는 Git 밖의 `/gpfs/kjhan/Astrafek/data/`에 두고, 코드·SQL 계획·체크섬·요약을 저장소에 기록한다.

- [표적–대조장 2-D ICL 스택](icl_2d_target_control_stack_README.md), [독립지원 민감도](icl_2d_independent_support_assessment_2026-10-19.md)
- [ICL–위성 다중비교 검정](icl_satellite_familywise_2026-10-18.md), [반경별 블록 폭 민감도](icl_satellite_perbin_block_scale_assessment_2026-10-19.md), [ICL tract-region 영향도](icl_hsc_tract_region_influence_assessment_2026-10-20.md), [위성–ICL 상관 tract-region 영향도](satellite_icl_tract_region_influence_assessment_2026-10-20.md)
- [소광 민감도](icl_pixel_extinction_assessment_2026-10-19.md), [평면 주입 전달](icl_plane_injection_assessment_2026-10-19.md), [위치별 PSF 및 공동 주입](psf_convolved_joint_injection_assessment_2026-10-19.md), [coadd 잔차 공분산](icl_coadd_residual_covariance_assessment_2026-10-19.md), [ICL–위성 lag 공분산](satellite_icl_lag_covariance_assessment_2026-10-19.md), [표적·대조 공분산](satellite_icl_target_control_covariance_assessment_2026-10-19.md), [측정 잡음 전파](icl_satellite_measurement_noise_mc_assessment_2026-10-19.md), [PSF 파일 감사](psf_picker_batch_2026-10-19.md)
- [CAMIRA 표본 확장 검토](camira_expansion_pool_assessment_2026-10-19.md), [B-mode 극성 좌표 검증과 calibrated stack](wl_bmode_coordinate_validation_2026-10-19.md), [표적 전용 약한렌즈 공분산](wl_cluster_spatial_jackknife_assessment_2026-10-19.md), [tract-region 영향도](wl_hsc_tract_region_influence_assessment_2026-10-20.md), [tomographic 일관성 평가](wl_tomographic_consistency_assessment_2026-10-20.md), [무작위 좌표 기하 파일럿](wl_random_coordinate_geometry_pilot_2026-10-20.md), [ICL 대조장 random 적합성 검사기](audit_wl_matched_control_centers.py)

이 문서는 결과를 통합해 읽기 위한 과학 현황 기록이다. 위의 탐색 분석을 논문용 측정으로 승격하지 않으며 새 영상·shape·photo-z 반출을 승인하지 않는다.
