# Prospective ICL–위성 연관 확인 설계안

**상태:** 검토용 설계 후보. 사용자 승인이나 사전등록이 아니며 현재 표본의 탐색 결과를 확인 결과로 바꾸지 않는다. `science_measurement=false`.

## 왜 독립 확인 표본이 필요한가

현재 Astrafek 탐색 분석의 후보 주 비교는 HSC `i2`에서 ICL target-minus-control 평균광도 `icl_delta_200_300`과 위성 외곽광 target-minus-control 평균 `satellite_delta_50_75`의 관계다. redshift와 CAMIRA richness를 선형 조정한 partial Pearson 상관은 68개 cluster/filter 쌍, 33개 연결 target/control 공간블록에서 `r=0.293`이다. 일반 pair bootstrap은 0을 포함하고, 비모수 순위 결과도 0을 포함한다. 이 표본은 반경·필터·지원 선택을 탐색한 뒤 얻었으므로 확인 검정에 재사용할 수 없다.

후보풀 확장 검토에서 현재 S20A CAMIRA 발자국의 가장 넓은 오프라인 표본은 target 5도 셀 84개를 점유할 수 있었다. 각 셀을 하나의 연결블록으로 가정해도 상한은 84개이며 영상 지원과 대조장 검사를 거치면 usable block 수는 줄 수 있다. 따라서 동일 발자국에서 cluster 행만 추가하는 방법으로 충분한 확인 표본을 확보한다고 전제하지 않는다.

## 사전 고정할 후보 추정량

- **주 표본:** 별도 HSC 하늘영역의 사전에 정의한 CAMIRA 은하단. 기존 800개 표본과 겹치는 cluster, 5도 target cell 및 사용 영상은 확인 표본에서 제외한다. 선택 기준·적색이동·richness 범위와 HSC `i2` 영상/대조장 적격 기준은 측정값과 상관을 보기 전에 고정한다.
- **주 설명변수:** 200–300 kpc 고리의 표적–매칭 대조장 차등 표면광도, `icl_delta_200_300`.
- **주 반응변수:** 50–75 kpc 위성 외곽 고리의 표적–매칭 대조장 차등 표면광도, `satellite_delta_50_75`.
- **주 연관 추정량:** cluster 단위 두 광도에 대해 redshift 및 CAMIRA richness를 선형 조정한 partial Pearson `r`. 양방향 검정, family-wise 단일 주검정 `α=0.05`. 이 추정량과 i2 선택은 현 탐색결과에서 유래하므로 독립 확인 표본에만 적용한다.
- **공간 불확실도:** 매칭 target/control 쌍 끝점이 속한 공간 단위를 공동으로 재표집한다. 기존 5도 연결블록 구현을 주 분석으로 고정하고, 2.5도·10도 폭은 민감도 분석으로만 보고한다. 최종 블록 정의가 새 footprint에서 연결영역을 대표하는지는 영상 수집 전에 모의·트랙 연결 자료로 평가한다.

## 표본 크기 목표

최소 계획 시나리오는 partial `r=0.2`, 두 nuisance 변수, 양방향 `α=0.05`, 근사 검정력 80%다. Astrafek Fisher-z 설계 계산상 **196개의 usable 독립 연결 공간블록**이 필요하다. `r=0.293`을 그대로 가정하면 92개가 필요하지만, 이 효과는 사후 선택된 탐색값일 수 있어 최소 표본 목표로 사용하지 않는다. 후보풀 등록 수는 지원·매칭·마스크 탈락률을 독립 pilot에서 측정한 뒤 `196 / 예상 usable fraction`으로 정하고, 블록별 실제 수를 기록한다. cluster 수만으로 196블록 목표를 대체하지 않는다.

검정가족을 한 개로 두는 것은 미래 표본에 대한 사전 설계일 뿐이다. 현재 30개 탐색 비교의 보정을 소급해서 없애지 않는다. 보조 반경·필터·통계량은 명시적으로 2차 결과로 두고, 해석 가능한 불확실도와 다중비교 보정을 함께 보고한다.

## 필수 품질 관문

1. 표본 선정·i2 target/control 매칭·각 고리의 픽셀/면적 지원을 결과변수 계산 전에 고정한다.
2. 실제 detection과 mask 재생성을 포함한 diffuse+satellite 공동 주입으로 크기별 sky transfer, 마스크 누출, 표적/대조장 차등응답을 시험한다.
3. 선택된 별 표본만으로 충분하다고 가정하지 않고, ICL 측정 반경까지 도달하는 경험 PSF 날개를 영상 기반으로 검증한다.
4. 모든 1차·2차 스택, 상관, 블록 재표집 입력과 산출물을 hash-linked manifest로 저장한다.
5. 위 게이트 중 하나라도 실패하면 결과를 방법 민감도 또는 탐색 결과로 보고하고 물리적 연관을 주장하지 않는다.

약한렌즈는 별도 분석 갈래다. 이 문서는 무작위 렌즈 차감, DESI 보정 `n(z)`, shear 보정 불확실도나 대응 공분산을 승인·대체하지 않는다.

## 재현 근거

- 표본·공분산 결과: [`satellite_icl_spatial_block_assessment_2026-10-19.md`](satellite_icl_spatial_block_assessment_2026-10-19.md)
- 블록·효과·다중검정 민감도 실행기: [`assess_icl_power_family_sensitivity.py`](assess_icl_power_family_sensitivity.py)
- 현재 출력: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_power_family_sensitivity_2026-10-20_01/`
- 현재 수치의 실행 상태: planning approximation; `science_measurement=false`.
