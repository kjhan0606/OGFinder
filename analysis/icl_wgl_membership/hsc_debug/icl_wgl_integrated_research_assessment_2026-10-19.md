# HSC ICL–위성 연관성 및 배경은하 렌즈 연구 현황 — 2026-10-19

## 현재 판정

Astrafek에는 HSC PDR3 영상·대조장·PSF 입력을 사용한 ICL 스택, 위성 외곽광 연관성 분석, 배경은하 약한렌즈 스택과 각 단계의 재현 스크립트가 구현되어 있다. 지금 결과는 모두 탐색 단계다. ICL–위성 관계는 다중비교와 sky 처리 민감도를 넘지 못했고, 렌즈 B-mode 마스크 방향은 외부 영역과 좌표 대조로 검증되었다. 렌즈 신호에는 아직 무작위 중심 차감과 대응 matched-field 공분산이 없다. 따라서 두 갈래 모두 `science_measurement=false`를 유지한다.

## 보유 자료와 실행 이력

- HSC 표본은 CAMIRA 은하단 800개다. 표적 영상 800개와 확보된 매칭 대조장 798개에서 동일 필터로 쌓을 수 있는 표적–대조장 쌍은 783개(i 302, i2 481)다.
- ICL 분석의 공통지원 2-D 스택과 K=2.5, 3.5, 5.0 독립지원 비교는 Astrafek의 고정 스크립트로 생성·검증되었다. 독립지원 추정치는 픽셀 교집합의 영향을 줄이는 민감도 분석이며 기본 측정값으로 채택하지 않았다.
- 위성 연관성에는 반경 구간별 94쌍(i 26, i2 68)이 있다. 주입·회수 검증은 위치별 HSC-I `pdr3_wide` coadd 모델 PSF 744개와 승인된 시험 PSF 1개를 사용해 57개 중심 및 688개 위성 위치에서 실행되었다. PSF 감사 매니페스트 SHA-256은 `6375bc58df6d5a15d837cf5b13d3db79ebdc6944dbf7510164f465b6cbe72f7e`다. 이 약 7 arcsec 코어는 넓은 ICL PSF 날개를 검증하지 않는다.
- 약한렌즈 표적 쪽은 CAS 원자료 3,069,544 lens–source 행(800 렌즈)을 가져왔다. 품질·거리·반경 컷 뒤 2,153,004 쌍, 609 렌즈, 2,051,991개 고유 배경원으로 줄었다. 기존 결과는 표적 중심 스택만 포함한다.

## ICL과 위성 외곽광

독립지원 비교에서 50–200 kpc 구간의 최대 같은 쌍 추정치 이동은 K=5, i2의 175–200 kpc에서 32.9%였다. 따라서 마스크 지원 방식이 외곽 프로파일의 크기를 실질적으로 바꿀 수 있다. 94쌍 분석에서 i2 전체지원 부분상관은 `r≈0.293`이지만 2.5°, 5°, 10° 공간 블록 폭에서 신뢰구간 해석이 달라진다. 별도 가족단위 검정은 48개 설정 중 추론지원 기준을 통과한 30개 어디에서도 Holm 또는 max-T 보정 뒤 유의성을 보이지 않았다.

은하수 소광을 픽셀별로 적용했을 때 기존 중심 소광 보정과의 ICL 차이는 사실상 0이었다. 반면 국소 평면 제거는 고정된 74쌍에서 i2 부분상관을 0.239에서 0.124로 낮추고 구간이 0을 포함하게 했다. 평면 주입 시험에서는 250 kpc 규모의 주입광 중 약 9–14%가 측정에서 빠졌지만, 해당 시험은 검출 마스크를 고정했고 상류 HSC sky 제거를 재현하지 않았다. 공개 확장 PSF 참조를 적용한 i 대조군 13개 시험에서는 단위 diffuse 주입의 중앙 21.5%가 ICL 측정 고리에 들어오고, 위성 외피는 35.3%가 들어왔다. 이 값은 특정 입력·마스크의 전달 민감도이지 관측값에 곱할 보정계수가 아니다.

공분산 감사는 position별 PSF 주입 표본 57개(i 13, i2 44)의 표적·매칭 대조 coadd에서 9×9 단거리 잔차 커널을 측정하고 실제 signed aperture 가중치에 전파했다. 큰 top-hat 근사 분산 비는 표적 1.59, 대조장 1.54–1.55였다. 정확한 ICL 및 위성 가중치와 양쪽 ICL 영상의 분산을 사용하면 대각 분산만 쓴 경우보다 ICL 잡음 표준편차가 i에서 1.234배, i2에서 1.240배였고, 잡음만으로 계산한 ICL–위성 상관 중앙값은 각각 −0.0073, −0.00091이었다. 이는 해당 57개 프레임에 조건부인 단거리·정상 커널 민감도다. source Poisson 잡음, 상류 sky 전달, 위치별 완전 공분산은 포함하지 않으며 공간 블록 오차에 단순 곱하면 같은 변동을 중복 셀 수 있다. 따라서 관측 연관성을 보정하거나 검출 여부를 결정하는 수치로 사용하지 않는다.

57쌍의 잡음 공분산을 같은 표본의 측정 ICL·위성값으로 전파한 10,000회 Monte Carlo에서는 i2 부분표본 44개의 측정 편상관 `r=0.409`가 finite-lag 잡음 범위 `[0.320, 0.465]`를 보였다. 하지만 전체 68쌍 i2 추정값은 `r≈0.293`이고 공간 블록 구간은 블록 폭에 민감하다. 공분산이 있는 표본만의 값이 전체보다 큰 점은 선택 효과를 경고한다. i 부분표본은 finite-lag 잡음 범위가 0을 포함한다. 이 조건부 잡음 실험은 spatial-block 범위와 합치지 않으며 과학 검출을 뜻하지 않는다. 세부 결과는 [측정 잡음 Monte Carlo 평가](icl_satellite_measurement_noise_mc_assessment_2026-10-19.md)에 기록했다.

현재 CAMIRA 후보를 더 넓히는 오프라인 검토에서는 기존 선택 아래 새로 분리되는 후보가 57개뿐이고 새 5도 하늘 셀은 없었다. 적색이동 범위를 넓혀도 분광 확인을 유지하면 새 5도 셀은 2개였다. 따라서 같은 발자국에서 영상 수만 늘리는 것으로 독립 하늘 영역의 부족을 해결하기 어렵다.

**해석:** 양의 ICL–위성 관계는 흥미로운 가설이지만 검출로 보고할 수 없다. 다음 ICL 관문은 반경·필터·마스크·통계량을 분석 전에 고정하고, 실제 검출·마스크 재생성을 포함하는 diffuse+위성 공동 주입으로 sky·PSF 날개·마스크 누출 전달함수를 평가하는 것이다. 독립 통계력을 높이려면 단순히 기존 CAMIRA 선택을 완화하기보다 새 독립 하늘 영역과 완전도·렌즈 보정을 함께 설계해야 한다.

## 배경은하 약한렌즈

기존 표적 전용 `DeltaSigma` 프로파일은 `0.1–0.2 Mpc`에서 약 `1.61×10^14 M☉ Mpc⁻²`, `2–3 Mpc`에서 약 `1.53×10^13 M☉ Mpc⁻²`로 감소한다. 52개 공간 셀 삭제 공분산은 독립 쌍 SEM보다 낫지만 셀이 연결된 HSC 관측영역과 일치하지 않을 수 있다. 이 표적 전용 프로파일에는 무작위 중심 차감, 대조장 분산, `n(z)`·shear 보정 불확실성, reduced-shear 보정이 빠져 있어 검출 유의도나 최종 질량으로 해석하면 안 된다.

HSC 문서는 `b_mode_mask`가 큰 B-mode 구역을 제거하는 데 쓰인다고 설명하지만 0/1 방향은 적지 않는다. Astrafek은 형상 품질 및 tomography 컷 뒤 좌표를 HSC Y3가 제외한 GAMA09H 사각영역(RA 132.5–140°, Dec 1.6–5°)과 직접 대조했다. 안쪽의 77,021개 lens–source 행은 모두 raw 0이고 바깥쪽 2,157,489개는 모두 raw 1이었다. shear 결과를 보지 않고 raw 0을 제외값으로 판별했다. 이 극성을 사용한 HSC Y3 four-bin SACC `n(z)`, 개체별 shear 보정, selection-bias 보정을 적용한 cluster-only 스택은 2,153,004쌍과 609개 렌즈를 유지했다. 52개 연결 sky 셀 jackknife는 별도로 계산됐지만, matched-random 차감과 source-overlap 보정은 아직 없다. 이 해석은 HSC가 공개한 마스크 설명과 Y3 시스템atics 연구의 제외영역 정의에 근거한다 ([HSC shape-catalog 문서](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s19a-shape-catalog-pdr3/), [HSC Y3 PSF systematics 연구](https://academic.oup.com/mnras/article/525/2/2441/7202336)).

무작위 중심 후보 좌표의 사용자 승인 CAS 내보내기는 800 렌즈, 예상 780,191 행, 32개 배치 계획이다. GPFS 체크포인트는 `running`이며 2/32 배치가 완료되고 세 번째 배치 `1062312`가 제출된 상태다. 원격 `lageunha`에서 이 계획을 실행하는 제출 스크립트 PID `1049265`가 실제로 살아 있고 CAS 응답을 기다리는 것을 확인했다. `1062312`는 Slurm 작업 ID가 아니므로 `squeue`로 확인할 수 없다. 실행 중인 제출자를 유지하고 이 배치를 다시 제출하지 않는다. 전체 기하 스크리닝은 좌표 내보내기가 32/32 완료되고 각 배치의 행 수·해시가 검증된 뒤에 시작한다.

현재 인증된 좌표 내보내기 프로세스가 CAS를 계속 확인 중이다. 그 프로세스를 유지하고, 32개 배치가 모두 행 수·해시 검증을 통과한 뒤 Astrafek 전체 기하 스크리너를 실행해 CAMIRA 3 Mpc 중첩·중복 좌표·렌즈별 지원량을 검사한다. 그 뒤 전체 조리개 영상 범위, 마스크, 깊이, PSF, sky, 소광과 공간 블록을 맞춰 무작위 중심을 선택한다. 모양 카탈로그 행 수와 전송량을 count-only로 산정한 뒤 무작위 원천의 shape 행을 반출하려면 별도 사용자 승인을 받아야 한다. 그 데이터 없이는 표적–무작위 차감과 대응 공분산을 완료할 수 없다.

교차 성분의 1.5–2.0 Mpc 구간은 `3.1 ± 1.6 M☉ pc⁻²`로 개별 bin에서 약 2σ 벗어나므로, 무작위 차감과 대응 공분산을 적용한 뒤 null-test 대상으로 확인한다. 다른 bin은 이 값에 맞춰 재선택하지 않는다.

## Astrafek 재현 경로와 근거 문서

실행기, Slurm 스크립트, 검증기, 요약 생성기는 `analysis/icl_wgl_membership/hsc_debug/`에 보관한다. 대형 HSC 파일과 실행 결과는 Git 밖의 `/gpfs/kjhan/Astrafek/data/`에 두고, 코드·SQL 계획·체크섬·요약을 저장소에 기록한다.

- [표적–대조장 2-D ICL 스택](icl_2d_target_control_stack_README.md), [독립지원 민감도](icl_2d_independent_support_assessment_2026-10-19.md)
- [ICL–위성 다중비교 검정](icl_satellite_familywise_2026-10-18.md), [반경별 블록 폭 민감도](icl_satellite_perbin_block_scale_assessment_2026-10-19.md)
- [소광 민감도](icl_pixel_extinction_assessment_2026-10-19.md), [평면 주입 전달](icl_plane_injection_assessment_2026-10-19.md), [위치별 PSF 및 공동 주입](psf_convolved_joint_injection_assessment_2026-10-19.md), [coadd 잔차 공분산](icl_coadd_residual_covariance_assessment_2026-10-19.md), [ICL–위성 lag 공분산](satellite_icl_lag_covariance_assessment_2026-10-19.md), [표적·대조 공분산](satellite_icl_target_control_covariance_assessment_2026-10-19.md), [측정 잡음 전파](icl_satellite_measurement_noise_mc_assessment_2026-10-19.md), [PSF 파일 감사](psf_picker_batch_2026-10-19.md)
- [CAMIRA 표본 확장 검토](camira_expansion_pool_assessment_2026-10-19.md), [B-mode 극성 좌표 검증과 calibrated stack](wl_bmode_coordinate_validation_2026-10-19.md), [표적 전용 약한렌즈 공분산](wl_cluster_spatial_jackknife_assessment_2026-10-19.md), [무작위 좌표 기하 파일럿](wl_random_coordinate_geometry_pilot_2026-10-20.md)

이 문서는 결과를 통합해 읽기 위한 과학 현황 기록이다. 위의 탐색 분석을 논문용 측정으로 승격하지 않으며 새 영상·shape·photo-z 반출을 승인하지 않는다.
