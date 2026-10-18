# ICL–위성 외곽광 연관성의 공간 블록 다중비교 검증

평가일: 2026-10-18 (Asia/Seoul)<br>
판정: **다중비교 후 유의한 연관성 없음; 탐색 분석만 가능**<br>
상태: `science_measurement=false`

## 목적과 입력

기존 800개 표본 탐색에서 i2의 일부 조정 부분상관 구간이 0을 벗어났지만,
같은 영상에서 만든 24개 Kron/면적 조건을 사후 비교한 결과였다. 컷 하나를 골라
유의성을 주장하지 않도록, 위성 50–75 kpc 외곽 초과광을 ICL 200–300 kpc 초과광,
적색이동, richness에 회귀하고 모든 필터·컷·중심 적격 집합에 대해 공간 블록 wild
bootstrap과 다중비교 보정을 수행했다.

입력은 `icl_satellite_cut_sensitivity_2026-10-18_01/matched_cluster_values.tsv`다.
기존 요약이 저장한 SHA-256과 일치하는 649개 행을 사용했다. 결과는
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_familywise_2026-10-18_02/`
에 있다. 계산기는
[`familywise_block_wild_association.py`](familywise_block_wild_association.py)다.

## 방법

대상과 대조장의 5도 하늘 셀을 행 전체에서 연결해 51개 공간 블록을 만들었다.
회귀의 귀무모형은 ICL 변수를 빼고 절편·적색이동·richness를 포함한다. 귀무모형
잔차에 HC2 leverage 보정을 적용하고, Webb 6점 분포의 공통 wild weight를 각 공간
블록에 부여했다. 같은 bootstrap 복제의 블록 weight를 모든 설정에서 공유해 컷
간 의존성을 보존했다. 검정 통계량은 공간 블록 CR1 학생화 t이며, 10,000회 복제의
최대 절대 t로 single-step max-T 보정값을 계산했다. 원시 bootstrap p 값에 대한
Holm 보정도 함께 저장했다. 48개 필터×Kron×면적×선택 설정을 모두 표에 남겼다.

최소 5개 은하단과 6개 이상의 독립 연결 블록을 가진 30개 설정만 추론 가족에
포함했다. 공간 블록이 부족한 5개와 은하단 수가 부족한 13개는 결과표에 사유와
함께 보존했다. 이 최소 블록 규칙은 관측 효과와 무관하게 실행 전에 고정했다.

## 결과

- 48개 설정을 보고했고, 30개가 추론 지원 기준을 통과했다.
- 지원 통과 30개 중 원시 wild-bootstrap `p < 0.05`는 6개였다.
- Holm 보정 `p < 0.05`는 0개, max-T 보정 `p < 0.05`도 0개였다.
- 가장 작은 원시 p는 i2, Kron 2.5, 50% 면적 조건의 전체 공통 표본이었다
  (`n=35`, 24개 공간 블록, 부분상관 `r=0.267`, `t=2.692`, raw `p=0.00930`,
  Holm `p=0.279`, max-T `p=0.636`).
- 기준 탐색 조건인 i2, Kron 2.5, 40%에서는 `n=53`, 28개 블록,
  `r=0.234`, raw `p=0.0334`, Holm `p=0.969`, max-T `p=0.621`이었다.

따라서 기존에 관찰한 양의 조정 상관은 컷·선택 집합 다중비교를 반영하면 안정적인
검출로 남지 않는다. 이는 관계가 없음을 증명하지 않는다. 현재 표본·측정 오차와
공간 지원으로는 관계를 식별하지 못했다는 뜻이다.

## 해석 한계와 다음 단계

max-T는 완전 귀무가설 아래 single-step 진단이며, 설정별 표본이 달라 subset
pivotality를 확인하지 못했으므로 강한 familywise 보장을 주장하지 않는다. 5도
연결 블록은 하늘의 공간 의존성을 근사하지만 coadd 픽셀 공분산, 장거리 sky 잔차,
PSF 날개, 소광, 멤버 질량·완전도, 마스크 선택 편향을 포함하지 않는다. 24개
설정은 같은 영상의 재분석으로 독립 재현이 아니다. 인과 해석이나 과학 검출 주장은
하지 않는다.

다음 과학 관문은 실제 위성 외곽광과 ICL을 함께 주입·회수해 광량 분배 누출을
측정하고, 잔류 천체 마스크·PSF 날개·상류 sky 전달을 보정한 뒤, 본 분석 전에
표본·마스크·반경·가설을 고정하는 것이다. 현재 결과만으로 특정 컷을 채택하지 않는다.

재현 산출물:

- `familywise_association.tsv` — 48개 설정의 지원량·계수·원시 및 보정 p 값
- `maxT_null.tsv` — 10,000회 최대 절대 t 귀무분포
- `summary.json` — 입력·프로그램·출력 해시, 설정, 상태, 제한 사항

방법론 참고: [Cameron, Gelbach & Miller (2008), wild cluster bootstrap](https://doi.org/10.1162/rest.90.3.414);
[MacKinnon & Webb (2017), wild bootstrap with unequal cluster sizes](https://doi.org/10.1002/jae.2508).
