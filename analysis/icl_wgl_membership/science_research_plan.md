# 장기 연구주제: HSC ICL·위성 외곽구조·은하-은하 약한중력렌즈

등록일: 2026-10-13 (Asia/Seoul)<br>
상태: **사용자 승인 장기 연구주제 / 관측 타당성 검증 단계**<br>
프로젝트: Astrafex Desktop / OGFinder<br>
연결: [실행 계획](implementation_plan.md), [분석 도구 설계](research_software_design.md)

도구의 상위 요구사항은 [목적 기반 파이프라인 설계](../../docs/science_workflow_designer.md)다.
이 연구는 구상성단 측광·LSBG 분석과 함께 범용 도구에서 구성하는 과학 recipe로 유지한다.

## 1. 핵심 질문과 범위

같은 별질량·색·적색편이 조건에서 위성은하의 외곽 별빛은 은하단 중심 거리와
ICL의 양·분포에 따라 달라지는가? 그 관계는 검출 누락, 배경 제거, PSF 및
광량 분배에 따른 측정 공분산을 보정한 뒤에도 남는가? 같은 은하단 표본에서
배경 은하의 전단을 스택해 질량 분포를 제약하고, 그 결과가 광도 측정의
질량·환경 해석을 어떻게 바꾸는가?

사용자는 위성은하 외곽구조 연구에 ICL과의 관계를 포함하고 장기 주제로
등록하도록 지시했다. 은하단 중심의 확산광 스택과 위성은하 중심의 별빛 스택을
같은 표본에서 연결하는 것이 주 분석이다. 사용자의 추가 지시에 따라 배경 은하
약한중력렌즈(galaxy-galaxy lensing)도 독립적인 공동 연구 갈래로 진행한다.
ICL 영상 스택과 약한렌즈 전단 스택은 서로 다른 자료와 추정량을 쓰며, 검증된
은하단 ID·적색편이·선택함수를 통해 마지막에 연결한다. 렌즈 질량은 BCG+확산광을
순수 ICL로 분해하는 대체물이 아니다.

연구 주제명은 **“HSC 은하단 ICL·위성 외곽구조와 배경 은하 약한중력렌즈의
공동 분석”**으로 둔다. 최초 측정량은 BCG+확산광이다. BCG 분리
방법을 명시하고 검증하기 전에는 이 값을 순수 ICL 또는 ICL 비율로 부르지 않는다.

## 2. 검증할 가설과 해석 한계

- H1: 통제된 위성은하 표본에서 외곽 광도 프로파일이 은하단 중심 거리에 따라
  달라진다. 중심부 외곽광 감소는 검증할 방향성 가설이며 전제하지 않는다.
- H2: 은하단 특성과 관측 조건을 통제한 뒤에도 확산광 지표와 위성 외곽구조
  잔차 사이에 관계가 남는다. 양·음·무상관을 모두 허용한다.
- H0: 겉보기 관계가 혼잡도, 검출 누락, 미검출 은하, PSF, 하늘 잔차 또는
  동일 픽셀의 광량을 두 성분에 나누는 과정만으로 설명된다.

조석 박탈은 가능한 물리적 설명이다. 병합, 은하단 유입 전 처리, 별의 종족,
궤도·투영, 생존 위성만 관측하는 선택도 고려한다. 현재 생존한 위성의 외곽광
부족분을 누적 ICL의 공급량과 직접 같다고 놓지 않는다. 상관은 인과 또는
개별 별의 이동 증거가 아니다. 검출되지 않는 경우에도 검증된 감도와 상한을 낸다.

## 3. 표본과 추가 입력

800개는 다운로드 목표이자 입력 후보 수이며 과학 분석 통과 수가 아니다.
기존 213개를 포함한 `0.30 <= z_cl < 0.55` 표본의 선택·층화 기록은
[실행 계획](implementation_plan.md)에 있다. 입력 위치는
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/sample_manifest.json`이다.
이 문서는 다운로드 완료 여부를 판정하지 않는다.

BCG 분광 선택, 기존 213개 고정, richness 우선 근접쌍 제거, 층화 추출을
선택함수에 기록한다. CAMIRA 멤버는 전체 위성 모집단의 완전 목록이 아니다.
새 검출 및 다중대역 카탈로그 대응으로 멤버십·전경·배경 오염을 추정한다.
원본 멤버 카탈로그의 정의되지 않은 열을 등급이나 별질량으로 해석하지 않는다.

필요한 추가 자료는 HSC 다중대역 측광·측광 적색편이와 품질 정보, 별질량 또는
명시적인 대체 변수, 위치별 PSF와 날개, 깊이·소광, matched random 및 필드
비교 표본이다. 약한렌즈에는 HSC S19A re-Gaussianization shape catalog, 해당
shape weight·multiplicative/additive calibration·B-mode mask, object_id에 연결한
source photo-z PDF, 무작위 렌즈/대조점을 추가한다. HSC PDR3 shape catalog는
S19A 내부 릴리스 기반 제품으로 PDR3의 일반 측광 테이블과 구별하며, S19A Wide
shape/측광 테이블의 실제 하늘 겹침을 먼저 센다([HSC S19A shape catalog 안내](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s19a-shape-catalog-pdr3/)).
처음에는 분해 가능한 밝은 위성부터 사용한다. DES 영상은
필수 입력이 아니며, 독립 검증이 필요한 경우에 별도로 검토한다.

## 4. 두 종류의 스태킹과 측정 정의

| 갈래 | 중심·좌표 | 주 측정량 | 처리 원칙 |
|---|---|---|---|
| 은하단 확산광 | BCG 중심, 고유거리 R (kpc) | BCG+확산광 표면밝기, 검증된 구간의 광도·기울기 | 위성·전경·배경·결함을 가리고 survey random으로 하늘 잔차 추정 |
| 위성 외곽구조 | 각 위성 중심, 고유거리 r (kpc) | 위성 평균 프로파일, 외곽 구간 광도와 비교 표본 대비 잔차 | 목표 위성을 보존하고 이웃·BCG·주변 확산광을 모형화/통제 |
| 공동 관계 | 같은 은하단 ID와 통제된 위성 집단 | 확산광 지표와 위성 외곽 잔차의 조건부 관계 | 두 측정의 공분산·선택·측정오차를 함께 전파 |
| 배경 은하 약한렌즈 | 은하단 중심, 고유거리 R (Mpc) | 배경 source shape에서 얻은 tangential shear 또는 ΔΣ(R), 질량모형 | S19A shape 응답·편향 보정, source P(z) 희석 보정, random-lens 및 B-mode null |

약한렌즈는 영상 픽셀을 평균하는 것이 아니라 shape catalog의 은하단–배경원
쌍으로 전단 프로파일을 스택한다. HSC S19A 카탈로그의 공식 shape 정의와 보정,
photo-z 희석 보정 절차를 사용하고, PDR3 일반 `meas`/shape 측정을 대신 끼워
넣지 않는다. 800개 후보와 shape catalog의 중첩률, 렌즈별 유효 배경원 수 및
적색편이별 유효 β를 산출한 뒤 분석 반경·질량 적합도를 고정한다.

R과 r을 혼용하지 않는다. 분석 가능한 반경은 PSF와 주입/회수 성능으로
예비 표본에서 정하고 본 분석 전에 고정한다. 현재 0–500 kpc 은하단 프로파일
설정은 초기 비교 범위이며 위성은하 자체의 반경으로 사용하지 않는다.

관측 HSC-I 표면밝기 스택은 좁은 적색편이 및 관측 조건 구간에서 시작한다.
넓은 적색편이 비교에는 소광, i/i2 응답, K 보정, 표면밝기 감광 및 광도 단위를
일관되게 처리한다. 1D 원본 고리 평균을 주 추정량으로, 면적 보존 2D 스택을
진단으로 사용한다. PSF 전방모형 또는 공통 해상도 처리와 유도 공분산을 기록한다.

위성은하를 은하단 중심 거리별로 나누고 별질량·색·적색편이·해상도를 맞춘다.
은하단별 위성 평균을 먼저 만든 뒤 같은 은하단 가중으로 결합하는 것을 주
추정량으로 설계한다. 위성별 가중 평균은 별도 감도 분석이며 서로 다른 모집단
평균임을 표시한다. 반경별 참여 은하단·위성·유효 면적과 가중을 저장한다.
공동 비교의 양쪽은 같은 유효 은하단 집합을 사용한다.

외곽/전체 광량비는 공통 분모가 만드는 상관을 점검하는 보조량으로 둔다.
우선 독립적으로 정의한 외곽 광도와 질량·색 조건부 잔차를 사용한다. 관측된
현재 별질량을 맞추는 선택 자체도 박탈 이전 질량을 맞추는 것과 다름을 기록한다.

## 5. 가짜 상관을 배제하기 위한 필수 실험

1. **공동 주입/회수:** 위성 외곽과 확산광 진폭을 독립적으로 바꾼 합성 장면,
   무상관 장면, 알려진 상관 장면을 실제 관측 조건에 넣는다. 검출·분리·마스크·
   배경·스택·관계 추정 전체를 다시 수행해 두 측정 사이의 누출을 측정한다.
2. **광량 분배 민감도:** 위성 마스크를 넓히거나 위성 모형을 빼는 방법을
   비교한다. 같은 빛이 위성에서 빠지고 ICL로 들어가 생기는 인위적인 반상관을
   반드시 검증한다. BCG/ICL/위성 외곽의 분해 불확실성을 공동으로 전파한다.
3. **두 종류의 대조 위치:** survey matched random은 관측 하늘 잔차를,
   같은 은하단 내 비슷한 R의 대조 위치는 위성 주변 BCG·확산광의 영향을
   검증한다. 후자는 진짜 확산광과 상관된 이웃을 포함하므로 빈 하늘로 취급하지
   않는다. 국소 배경 고리가 위성 외곽을 지우는 효과도 주입으로 확인한다.
4. **선택과 PSF:** 밝기·크기·표면밝기·R별 검출 및 분리 완전도, 멤버 오염,
   밝은 별 마스크, PSF 날개를 검증한다. 정보가 사라진 영역을 큰 역확률 가중으로
   억지 복원하지 않고 분석 범위에서 제외하거나 식별 불가로 보고한다.
5. **공동 공분산:** 은하단 단위 재표집과 하늘 영역 jackknife로 두 갈래를
   함께 재계산한다. 같은 은하단의 위성을 독립 은하단처럼 세지 않는다.
   두 평균 스택 한 쌍으로 상관계수를 계산하지 않는다. 사전 정의 집단별 차이
   또는 은하단별 잠재 측정량을 둔 오차 포함 모형으로 관계를 추정한다.
6. **잡음에 의한 집단 선택:** 측정된 ICL로 높은/낮은 집단을 나눌 때는 분류
   오차와 선택 편향을 전파한다. 분리 가능한 노출·픽셀 또는 교차 적합은 보조
   검증이며 공통 하늘·PSF 오차까지 독립으로 만든다고 가정하지 않는다.
7. **약한렌즈 null·교정:** shape weight와 S19A shear-response/bias 항을 적용하고
   source P(z)로 foreground/member dilution을 보정한다. B-mode, 회전 shape,
   무작위 렌즈, 적색편이 순서 교환, source density/마스크 null을 확인한다.
   적색편이 겹침이 부족하거나 shape footprint 밖인 렌즈는 탈락 사유와 수를
   남기고, 광도 스택 표본과 다른 유효 표본을 같은 모집단인 것처럼 합치지 않는다.

코애드에 넣은 인공광은 그 이전 처리에서 이미 사라진 신호를 검증하지 못한다.
서베이 배경 처리의 전달 성능은 처리 단계가 다른 제품 또는 노출 단계 실험으로
별도 확인한다. 공동 측정이 식별되지 않으면 성분별 물리 해석을 보류한다.

## 6. 단계와 완료 기준

| 단계 | 산출물 | 완료 기준 | 상태 |
|---|---|---|---|
| A. 자료·선행연구 | 표본 계보, 입력/교정 목록, 선행연구 비교표 | 누락 입력 및 연구 차별점 명시 | 일부 자료·초기 문헌 확인 |
| B1. ICL 예비 측정 | matched-control 차감 프로파일·위성 위치 밀도 | 고정 반경·대역별로 재현 가능하고 BCG/sky/PSF 게이트 통과 | 반경별 지원 224쌍(137 strict쌍과 행별 일치); Kron/면적 24조건 민감도에서 큰 마스크는 지지가 급감하고 flux도 이동. 과학 게이트 미통과 |
| B2. 위성 외곽 예비 측정 | HSC 다중대역 위성 측광·위성 중심 스택 | 멤버 완전도·질량/색 맞춤·광량 누출 실험 통과 | 전체 800장 탐색 스택 5,832행(1,768개 위성); 7-bin 공통 은하단 i 29개/i2 103개. 외곽 평균은 0과 양립; 공통 표본 ICL 조인 21/53개, 상관 미확정·과학 게이트 미통과 |
| B3. 약한렌즈 예비 측정 | S19A shape source 겹침표·전단 profile | 유효 source 수, 보정 n(z), 응답 보정 및 null 시험 통과 | 800 렌즈·3,069,544 shape 행 검증 완료; 2,234,510개 컷 통과 쌍의 raw `e_t/e_x` 진단도 완료. 보정 n(z), B-mode polarity, random/null·공분산 미완료; science gate 미통과 |
| C. 공동 검증 | 주입/회수, null, 공분산·감도 보고서 | 무상관 주입의 가짜 관계와 알려진 신호 회수 오차가 사전 오차 예산 내 | 확산광 코애드 주입·대조장 공분산 실행; joint satellite–ICL 및 lensing null 미완료 |
| D. 분석 정의 고정 | 표본/반경/단위/가중/모형/판정 기준 버전 | 본 분석 전 숫자 기준과 주·보조 가설 고정 | 미실행 |
| E. 전체 표본 | 두 스택, 공동 관계, 상한·민감도 | 효과와 체계오차, 참여 표본 및 선택함수 함께 보고 | 미실행 |
| F. 물리 해석 | 관측을 재현한 시뮬레이션 비교, 논문 | 박탈·병합·선택 설명을 구분 가능한 범위에서 평가 | 미실행 |

30–50개는 예비 설계 규모다. 800개의 통계력이나 검출 유의도는 아직 보장하지
않는다. 단계 B/C에서 편향·산포를 구해 필요한 유효 은하단 수와 세분화 수준을
정한다. 잔류 편향 허용치 등 수치 기준은 효과를 본 뒤 유리하게 바꾸지 않는다.
`science_measurement=false`는 주장별 검증을 통과하기 전까지 유지한다.

2026-10-15 실행 상세는 [주입·회수와 공분산 검증](hsc_debug/prestack_validation_README.md)에 있다.
대조장 798개, 확산광 300개 주입, 기존 은하단 50개를 사용했다. 원본 고리 추출의
50–300 kpc 회수는 지원 집합에서 기술 기준 내였지만, 전체 50–500 kpc의 공통
은하단−대조장 쌍이 부족해 전체 스태킹 공분산은 승인하지 않았다.
이 결과를 위성 외곽과 ICL의 공동 관계 검증 완료로 간주하지 않는다.

2026-10-16에는 같은 50쌍에서 마스크 크기와 양쪽 면적 기준을 교차한 공통 지원
민감도를 추가했다. 원래 Kron 3.5/50% 정의에서는 50–300 kpc 네 구간을 모두
지원하는 쌍이 1개였다. Kron 2.5/40% 탐색 조건은 24쌍·21개 공동 하늘 블록,
rank-4 기술 공분산 진단을 만들었고, 여섯 주입 모형은 모형별 28개 대조장에서
네 구간 중앙 회수 편향 절대값 최대 5.90%였다. 이 표본은 기존 40개 면적 적격
후보에서만 왔으며 전체 800개 대표 표본이 아니다. 40% 컷과 작은 마스크를
과학 기준으로 채택하지 않았고 실제 검출원 외곽빛 잔류·상류 PSF/sky 전달도
검증되지 않았다. 결과와 입력/출력 해시는
`sample_800_2026-10-13/prestack_pair_support_2026-10-16_02/`에 기록했다.
다음은 이 후보의 잔류 천체광 민감도를 기존 matched sample 전체에 확장한 뒤,
본 분석 전 유효 면적·마스크·반경 정의를 고정하는 일이다. 공동 위성–ICL 주입과
공분산은 여전히 미완료다.

2026-10-16 전체 표본의 paired 면적 지원 집계도 완료했다. 800개 은하단 중 798개가
대조장과 연결되며, Kron 2.5/40%에서 공통 면적 기준을 충족한 수는 50–300 kpc
138쌍(동일 필터 137쌍, 공동 5° 블록 44개), 50–500 kpc 69쌍(필터 일치 68쌍,
블록 37개)이다. 이는 사전 정의된 최소 표본·블록 수 스크린을 넘긴 것이며 밝기
측정이나 flux covariance 검증은 아니다. Kron 3.5/50% 기준에서는 같은 범위가
각각 4쌍과 1쌍뿐이다. 50–300 kpc 후보를 전체 표본 예비 측정으로 시험할 수는
있지만 40% 면적 컷/작은 마스크를 채택하거나 과학적 검출로 해석할 근거는 부족하다.
실제 BCG/잔류원 마스크, 필터 응답·소광, PSF/상류 sky 전달과 800개 분포에 맞춘
주입 회수 및 공간 공분산을 확인하기 전까지 영상 스택/위성–ICL 관계는 보류한다.
해당 후보는 i 44개/i2 94개이며 richness 50+가 없어, 필터 응답 보정과 표본 선택
효과를 별도로 검토해야 한다.
세부 수치와 산출물은 `hsc_debug/prestack_validation_README.md` 및
`sample_800_2026-10-13/prestack_full_pair_area_support_2026-10-16_01/`에 있다.

2026-10-18에는 저장된 지원 인덱스와 해시가 확인되는 기존 50개 검증 표본에서
CAMIRA 멤버의 투영 위치 밀도와 matched-control 차감 광도 프로파일을 비교했다.
Kron 2.5/공통 면적 40%/50–300 kpc 조건의 24쌍은 모두 대응 필터가 같았지만
i 6쌍과 i2 18쌍으로 나뉘며, 자동 유일 중심 후보는 8개뿐이다. i2 18쌍의 반경별
Spearman ρ는 0.51, 0.44, 0.44, 0.46이며 네 구간의 명목 p 값은 각각 0.030,
0.066, 0.065, 0.055다. richness와 적색편이를 선형 제거한 Pearson 상관은
0.45–0.50 범위다. 이 작은, 면적 선택된 검증 집합의 탐색 통계로서 p 값은
공간 공분산·네 구간 다중비교를 반영하지 않는다. 이 상관 계산에 투입한 멤버 변수는
CAMIRA 위치 count뿐이다. 별도로 보존된 `camira_s20a_wide_sm_member.dat.gz`를
현재 800개 표본에 대조하니 manifest에 기록된 SHA-256과 일치했고 800/800 은하단의
좌표·적색이동·richness가 맞았으며 멤버 행 36,513개가 연결됐다. 문헌의 CAMIRA 멤버
카탈로그는 별질량 `M_*`와 멤버 가중치 `w`를 제공한다. 로컬 파일은 머리글이 없고
기존 `06_fetch_members.py`와 그 해시 일치 실행 요약 `logs/member_summary.json`은
열 8을 `weight`, 열 7을 미정 `field7`로 기록한다. 분석은 이 원본 기록과 카탈로그
SHA를 확인하고 열 8을 CAMIRA 가중치로 사용한다. 열 7 값 범위(9.986–12.320)는 log
별질량과 부합하지만, 파일 제작 기록에 정의가 없어 질량 후보로만 남긴다. 특히 CAMIRA
가중치는 질량·공간 필터를 포함하므로 단순한 비편향 위성 밀도로
쓸 수 없다. 이 파일에는 위성별 보정 측광이 없고 현 분석도 위성 광도를 사용하지
않는다. 따라서 이 결과는 위성 외곽광–ICL 관계의 측정이 아니다. 산출물은
`sample_800_2026-10-13/satellite_association_pilot_2026-10-18_02/`에 있다.

전체 800개 대조에서도 CAMIRA 멤버 행이 800/800 클러스터에 연결됐고, 위치·z·richness
검사가 통과했다. 다만 위성별 cmodel/다중대역 측광이나 실제 영상 유래 외곽광은 없다.
후속 단계는 해당 무헤더 파일의 열 7 정의를 확인하고 CAMIRA 가중 규약을 고려해, 위성 질량
구간별 표본 및 HSC 영상에서 측정한 위성 중심 잔광 프로필을 구성하는 것이다.
CAMIRA 멤버의 `M_*`, `w` 제공과 필터 의미는
[Oguri et al.의 CAMIRA member catalog 설명](https://doi.org/10.3847/1538-4357/ad76ad),
[CAMIRA richness·질량 필터 설명](https://academic.oup.com/pasj/article-abstract/70/SP1/S20/3885861)을
참고하되, 이 문헌 설명만으로 로컬 파일 열 번호를 확정하지 않는다.

2026-10-18에는 800개 대상 영상의 전체 반경 프로필 계산도 완료했다(Slurm
427234). 800개 영상·마스크 해시를 계산 중 검증했고 798개 대조장 프로필과
대응시켰다. 입력 필터는 i 309개/i2 491개, 대조장 필터는 i 312개/i2 486개이며
동일 필터 대상-대조장 쌍은 783개다. Kron 2.5/공통면적 40%/50–300 kpc 조건의
필터별 지원쌍은 i 44개, i2 93개로 총 137개뿐이고, 각각 15개와 37개의 공동
5도 하늘 블록을 갖는다. 주 지원 컷에서 4×4 공간 부트스트랩 공분산은 두 필터
모두 유한·전계수였지만, 이는 측정된 공분산의 수치 안정성만 보이며 누락된
PSF·sky·마스크 체계오차를 포함하지 않는다.

이 137쌍의 차감 표면밝기는 두 대역 모두 50–100 kpc에서 200–300 kpc로 갈수록
감소했다. 다만 50–300 kpc 스택에서 K2.5 중심 픽셀이 마스크된 대상은 90/137이고,
800개 전체에서 자동 유일 중심 후보는 258/800뿐이다. K2.5/50%/50–300 kpc로
강화하면 지원쌍은 i 15개/i2 26개로 줄며, Kron 3.5/50%에서는 각 대역 2쌍뿐이라
공분산 추정이 성립하지 않는다. 50–500 kpc에서 i 표본도 20쌍/10블록으로
6개 빈 공분산 표본량 기준을 통과하지 못한다. 따라서 결과는 전체 표본 계산의
재현성·지원·오류 예산을 평가하는 탐색 산출물이지, 과학용 ICL 검출이나 최종
반경 선택으로 채택하지 않는다. 상세 수치와 판정은
[`stack_quality_assessment_2026-10-18.md`](stack_quality_assessment_2026-10-18.md),
표·공분산·해시 기록은 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_flux_profiles_2026-10-18_02/`
에 있다.

같은 137쌍에 대해 5,000회 공동 5도 하늘 블록 부트스트랩으로 멤버 위치 밀도–차감광
상관의 공간 민감도를 확인했다. i 대역의 네 Spearman 95% 구간은 모두 0을 포함했다.
i2의 Spearman 구간은 50–100, 100–150, 150–200, 200–300 kpc에서 각각
[0.080, 0.588], [0.105, 0.456], [0.054, 0.392], [0.154, 0.496]으로 양수였지만,
richness·z를 함께 선형 제거한 Pearson 구간은 100–150 kpc의 [0.058, 0.365]만
0을 제외했다. i 대역도 100–150 kpc에서 부분상관 구간 [0.125, 0.573]이 양수였으나
44쌍의 네 구간 탐색 중 하나다. 이 구간은 다중비교·마스크/sky/PSF 체계오차를
교정하지 않은 공간 안정성 진단이며 위성 외곽 별빛–ICL 검출이 아니다. 재현 산출물은
`sample_800_2026-10-13/satellite_association_full_2026-10-18_01/`이다.

2026-10-18에는 저장된 SEP 분할·검출 산출물을 전체 800장에서 독립 감사했다. 고유 ID,
표 행 수, FITS 차원·이진 마스크·중첩 관계와 이미지-마스크 WCS는 800/800 통과했고,
WCS 왕복 오차 최댓값은 1.1e-7 픽셀이었다. 원본 HSC 품질 비트 재계산은 필터·적색이동
층화 16장 표본에서만 대조해 16/16 일치했다. 이 감사는 마스크 파일 무결성을 보일 뿐,
BCG 중심 258/800 유일 후보, PSF 날개 및 천체원 마스크 완전성을 인증하지 않는다.

기존 800장 영상에서 CAMIRA 위치를 검출표에 1 arcsec 이내 일대일 매칭해 위성 중심
5–75 kpc 표면밝기 프로필을 계산했다. 엄격한 7개 반경 동시 커버리지 조건을 적용한
권위본은 9,468개 후보 중 90개(0.95%) 프로필(i 12개/11개 은하단, i2 78개/49개 은하단)을
보존했다. 반경별 유효 표본을 따로 두는 보완 분석을 같은 800장에 적용하자, 60% 대상·대조
커버리지와 반경당 3개 대조 기준을 만족하는 위성-반경 측정은 5,832개가 됐다. 고유 위성은
최소 한 구간에서 1,768개였다. 50–75 kpc 외곽 구간은 i 70개 위성/42개 은하단,
i2 308개/140개 은하단이다. 각 구간의 측정 은하단 수가 서로 다르므로 구간 사이 변화는
동일 표본의 프로필 기울기로 해석할 수 없다. 해당 외곽 평균은 i 2.51 [-0.30, 5.50],
i2 1.19 [-0.24, 2.49] nJy/arcsec² (은하단 부트스트랩 95%)로 두 대역 모두 0을 포함한다.
공분산, PSF 날개, sky transfer, CAMIRA 멤버 가중·질량·색, 위성 완전도, 소광·대역 응답은
반영하지 않았다.

같은 필터에서 7개 반경 모두 지원되는 은하단 교집합은 i 29개, i2 103개다. 은하단 목록을
고정하고 5도 하늘 블록 공분산을 계산한 공통 표본의 50–75 kpc 평균은 i 2.86
[-1.15, 7.95], i2 1.00 [-0.21, 2.07] nJy/arcsec²이며 두 구간 모두 0을 포함한다.
프로필 모양은 같은 은하단 표본으로 비교할 수 있지만, 반경에 따라 사용되는 위성 자체와
위성 수는 달라져 같은 위성을 추적한 프로필은 아니다. 7×7 공간 블록 공분산은 두 필터에서
모두 full-rank·양의 정부호였으나, PSF/sky 계통오차와 픽셀 공분산은 해결되지 않았다.

기존 ICL 스택은 50–300 kpc 네 구간이 모두 면적 기준을 통과해야 해 137쌍만 남았다.
반경별로 독립 선택한 보완 분석은 전체 800개 대상·798개 대조장에서 1,355개 쌍-반경
행을 보존했고, 200–300 kpc에는 224쌍(i 71, i2 153)이 남았다. 이 반경의 평균
차감 표면밝기는 i 3.23 [2.15, 4.23], i2 3.79 [2.68, 4.80] nJy/arcsec²이며,
괄호는 5도 연결 하늘 블록 부트스트랩 95% 구간이다. 기존 strict 137쌍의 200–300 kpc
값은 행별로 전부 정확히 재현했다. 구간별 표본이 다르므로 이 평균들을 공통 표본의
방사 프로필로 읽어서는 안 된다. sky/PSF/마스크 교정 전 탐색 결과라 과학 게이트는
통과하지 않았다.

전체 sample에서 200–300 kpc annulus를 Kron 인자 2.5/3.5/5.0, target·control 면적 문턱
25/35/40/50%로 교차한 24조건 민감도도 수행했다. K2.5/Kron 3.5의 지지가 상대적으로
충분한 25–40% 컷에서는 평균 차감값이 두 대역에서 양수였지만, Kron 5.0에서는 선택 수가
크게 줄고 평균이 작아지거나 오차 구간이 0을 포함했다. 5.0/40%의 i는 1쌍, 5.0/50%는
i 0쌍/i2 1쌍이라 공간 구간을 안정적으로 추정할 수 없다. 면적 문턱 변화에도 평균이
단조롭지 않아 표본 선택과 마스크 변화의 영향을 분리하지 못했다. 이 민감도는 science cut을
선정하지 않으며 `science_measurement=false`를 유지한다.

위성의 7-bin 공통 은하단 표본과 200–300 kpc ICL annulus를 조인하면 i 21개, i2 53개다.
5도 블록 bootstrap Spearman ρ는 i 0.312 (p=0.169; 95% [-0.122, 0.620]), i2 0.166
(p=0.234; [-0.061, 0.369])이며 두 구간 모두 0을 포함한다. z·richness 조정 i2 partial
Pearson r은 0.234 [0.024, 0.421]이지만 비조정 상관·중심 민감도·다른 지원 정의에서는
안정적으로 재현되지 않는다. 중심 후보 고유·비가림 i2 표본은 10개뿐이다. 따라서 현재
자료는 위성 외곽광–ICL 물리적 상관의 증거로 해석할 수 없다.

같은 공통 위성 표본에서 ICL의 Kron 마스크 2.5/3.5/5.0과 target/control coverage
25/35/40/50%를 교차한 12개 설정을 i와 i2에서 각각 평가해 24개 상관 민감도 조건을 검사했다.
전체표본 20/24 조건에서 추정
가능했고, 비조정 Spearman 공간 블록 구간이 0을 제외한 것은 i2의 50% coverage 두 조건뿐이다.
z·richness 조정 partial 상관은 i2 네 조건에서 구간이 0을 제외했으나, 서로 연관된 사후
비교이며 선택·마스크 민감도와 독립 재현이 해결되지 않았다. 이는 i2 사전 고정 후속 분석을
검토할 근거이지 물리적 상관의 검출은 아니다. 산출물은
`sample_800_2026-10-13/icl_satellite_cut_sensitivity_2026-10-18_01/`이며
`science_measurement=false`를 유지한다.

같은 48개 컷·필터·중심선택 설정에 대해 5도 연결 하늘 블록 wild-cluster bootstrap을
10,000회 수행하고, 공통 block weight로 단일 단계 max-T와 Holm 보정을 계산했다.
48개 설정을 모두 보존했으며, 은하단 5개·공간 블록 6개 이상인 30개만 추론 지원
기준을 통과했다. 30개 중 원시 p<0.05는 6개였지만 Holm 또는 max-T 보정 후 p<0.05는
없었다. 기준 i2/Kron 2.5/40% 설정도 raw p=0.0334, Holm p=0.969, max-T p=0.621이었다.
이는 다중비교 아래 독립적으로 재현된 연관성을 보이지 않는다. max-T는 완전 귀무
가설의 탐색 진단이며 strong familywise 보장은 주장하지 않는다. 블록은 픽셀/공분산,
PSF 날개, sky 전달, 소광 및 멤버 완전도를 포함하지 않는다. 상세 방법·수치는
[공간 블록 다중비교 검증](hsc_debug/icl_satellite_familywise_2026-10-18.md)과
`icl_satellite_familywise_2026-10-18_02/`에 있고, 과학 측정 플래그는 계속 거짓이다.

정량·정성 평가와 산출물 경로는 [위성 외곽광 전체표본 평가](hsc_debug/satellite_outer_profile_assessment_2026-10-18.md)에
기록했다. 권위 데이터는 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/` 아래
`satellite_detection_psf_2026-10-18_01/`,
`satellite_outer_profile_full_2026-10-18_03/`,
`satellite_outer_profile_perbin_full_2026-10-18_01/`,
`satellite_common_support_2026-10-18_01/`,
`icl_annulus_support_2026-10-18_01/`,
`icl_annulus_sensitivity_2026-10-18_02/`,
`satellite_outer_icl_association_perbin_2026-10-18_03/`,
`icl_satellite_cut_sensitivity_2026-10-18_01/`에 있다. 이전 `_01`·`_02` 완전 프로필과
수정 전 `_01` 상관 결과는 사용하지 않는다.

2026-10-18 현재 약한렌즈 shape source 자료는 없다. HSC 공식 문서는 S19A
re-Gaussianization shape catalog를 PDR3 incremental product로 설명하며,
`s19a_wide.weaklensing_hsm_regauss`를 `meas2`·`meas`와 `object_id`로 연결하고
photo-z PDF를 별도 이용하도록 안내한다. 이에 따라 다음은 800개 렌즈 좌표의
CAS 중첩 row count와 추정 전송 크기를 먼저 산출하는 일이다. source 자료 취득은
기존 HSC 다운로드 승인 원칙에 따라 크기·열·행을 제시하고 별도 사용자 승인을
받은 뒤 진행한다. shape stack·질량 측정은 아직 하지 않았다.

S19A 분석에서는 shape catalog의 sky-coordinate `e1/e2`, derived weight, per-object
`m/c1/c2`, ellipticity RMS·measurement uncertainty, B-mode mask와 `hsc_y3_zbin`을
보존한다. 기본 자료 경로는 HSC Y3 네 tomographic source-bin에 해당하는 SACC의 fiducial
`n(z)` 곡선을 사용하며, shape 행마다 full photo-z PDF를 받을 필요가 없다. 이 fiducial
곡선은 posterior 표본 자체가 아니므로 posterior/systematic 불확실성은 별도로 전파한다. HSC 문서에서 zbin 0은
다중봉우리 P(z) 때문에 배제된 source이고 1–4가 네 tomographic bin이다. raw
ellipticity를 이미 보정된 shear로 간주하지 않는다. HSC estimator의 weighted responsivity와
multiplicative/additive 및 selection-bias 보정, photo-z dilution correction을 적용한
뒤 lens-source pair profile을 만든다. 자세한 catalog columns와 PDF 처리 요구는
[HSC S19A shape catalog 안내](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s19a-shape-catalog-pdr3/)
와 [HSC 3-year shear catalog 방법론](https://arxiv.org/abs/2107.00136)을 따른다.
이를 위한 `hsc_debug/s19a_shape_stack.py` adapter를 추가했다. 현재 full P(z) 입력과
HSC Y3 tomographic `n(z)` 입력을 모두 지원한다. 후자는 계획된 source-row의 `source_zbin`을
사용해 pair별 `<D_ls/D_s>`를 계산하므로 134.4 GiB full-PDF tar가 필요하지 않다.
빈별 HSC responsivity와 per-object `m/c`를 반영해
접선·교차 전단 및 weak-shear 근사 `DeltaSigma` 프로필을 계산한다. 선택편향은 HSC
Y3/S19A edge-density 식으로 반경 빈마다 자동 계산하며, pair 가중치는
`shape_weight * SigmaCritInv**2`, additive 항은 가중 평균 PSF ellipticity를 쓴다.
식과 계수는 HSC가 연결한
[hsc-y3-calib 보정 코드](https://github.com/mr-superonion/hsc-y3-calib)의
`get_sel_bias`를 따른다. 경계 밀도 식의 수치 검사, 합성 주입 회수, PDF·tomographic
`n(z)`·SACC 입력·ID·mask·선택컷 검사를 포함한 15개 시험이 통과했다.
`load_hsc_y3_nz_sacc`는 공식 SACC의 `wl_0`–`wl_3` tracer 격자와 분포를
검사한다. 승인된 파일을 내려받아 실제 파일 로딩과 중복 테이블 검사를 통과했고,
현재 `source_zbin` 1–4 분포와 shape 선택 후 빈별 행 수를 대조했다
([검증 보고서](hsc_debug/hsc_y3_nz_sacc_validation_2026-10-19.md)). 이 함수와
fiducial 분포만으로 n(z) 불확실성이 전파되는 것은 아니다.
HSC 안내에서 B-mode 영역 마스크는 분석별 선택 사항으로 설명되므로 adapter는 기본
적용과 명시적 미적용 민감도 실행을 모두 지원한다. 실자료 비교 때 두 결과의 source 수,
shear profile, null 진단을 함께 기록하며 마스크 정책을 결과를 본 뒤 선택하지 않는다.
고정 `n(z)` 적용은 HSC Y3 사진·clustering-redshift 기반 보정 경로지만,
fiducial `n(z)` 보정 및 posterior/systematic 불확실성을 현재 stack에서 전파하지 않는다. 실제 S19A
source 행과 보정 `n(z)` 입력도 아직 없다.
`subtract_random_lens_profile`은 cluster−random 점 추정량을 계산하고 독립 shape-noise
오차를 임의로 합치지 않는다. `paired_jackknife_covariance`는 같은 하늘 블록을
cluster/random 양쪽에서 제거한 replicate의 차이로 공분산을 계산한다. 아직 실제 random
pair 행이나 공통 하늘-block replicate는 없다. reduced-shear 적합, 실제 random-lens
차감, source 중복 및 하늘 블록 공분산 추정, 선택 보정 오차 전파는 미완료다.
그러므로 이 adapter의 출력도 항상
`science_measurement=false`이며 일반 `stack_shear_profiles`나 합성 metacal 시험이
서베이 교정을 대신하지 않는다.

같은 날짜에 800개 렌즈의 각거리와 고유반경 3 Mpc를 사용한 기하학적 사전 계산에서
렌즈 원뿔 반경은 7.71–11.43 arcmin (10–90 백분위 7.97–10.56 arcmin), 합산 면적
상한은 58.09 deg²였다.
은하단 쌍 246개에서 두 렌즈 원뿔이 겹치며 143개 렌즈가 하나 이상의 다른 렌즈와
겹친다. 이는 shape 자료의 실측 중첩/겹침 수가 아니라 좌표 기하만의 계산이다.
CAS 조회는 렌즈 ID별 원뿔 결과를 보존하고, 중복 source object_id가 여러 렌즈의
profile에 기여하는 현상을 허용하되 최종 공분산은 하늘 블록과 source 중복을
반영해야 한다. 전부 합쳐 source를 한 번만 남기는 union query는 lens-source
pair 정보를 잃으므로 분석 입력으로 쓰지 않는다.
렌즈별 원뿔의 shape 행 수만 세는 CAS SQL 32개 배치는
`sample_800_2026-10-13/wl_shape_preflight_2026-10-18_01/`에 준비했다. 최초 문서화
시점에는 인증된 CAS 제출이 아직 이뤄지지 않았으나, 이후 승인된 count-only 조회를
`lageunha`에서 완료했다. `hsc_debug/submit_wl_shape_counts.py --plan`은 SQL 해시와
렌즈 ID 800개를 로컬 검증하고, `--submit`은 렌즈별 정수 집계만 요청한다.

2026-10-18 후속: 사용자가 승인한 범위는 **렌즈별 정수 count만**이다. `lageunha`에서
32개 CAS 배치가 모두 완료됐고 각 결과의 SHA-256·25행·요청 렌즈 ID를 검증했다.
합계 800행·고유 렌즈 800개·음수 count 0개이며, 렌즈별 count 합은 3,069,544,
0개인 렌즈는 168개, 중앙값 4,438, 최댓값 9,767이다. 이 합계는 렌즈 원뿔마다
계산한 count라 중복 source와 source 품질/적색편이 컷 적용 전의 lens-source 후보량이며,
고유 source 수가 아니다. 체크포인트와 요약 모두 `source_rows_downloaded=false`,
`science_measurement=false`, `credentials_saved=false`다. 실제 source 행 취득은
예상 전송량과 필요한 열을 제시한 별도 승인 전까지 시작하지 않는다.

count-only 집계의 표본 의존성도 사후 기록했다. 렌즈별 원시 후보 count와 기하 적색편이의
Spearman `rho=-0.427`, richness와는 `rho=+0.048`이었다(기술 요약이며 p 값 검정 없음).
고적색편이 구간 `0.50–0.57`에서 170개 중 55개(32.4%)가 count 0이며 중앙값은
3,601개였다. 이 값은 3 Mpc 원뿔의 기하·survey footprint/depth·catalog 선택이 섞인
후보량이므로 최종 렌즈 효율이나 물리적 적색편이 추세가 아니다. 품질·zbin·radial 및
random-lens 컷 뒤의 유효 source 수를 다시 산출해야 한다. 상세 구간표와 재현 코드는
[count-only 표본 의존성 감사](hsc_debug/wl_shape_count_selection_audit_2026-10-18.md),
산출 요약은
`sample_800_2026-10-13/wl_shape_count_selection_audit_2026-10-18_01/`에 있다.

2026-10-18 source-row 준비: count 결과를 기준으로 동일한 3 Mpc 렌즈 원뿔에서
22개 필드만 반환하는 32개 SQL batch를 로컬 생성했다. 예상치는 3,069,544
lens-source pair 행, 약 0.73–1.37 GiB 미압축 CSV이며 실제 전송 전 추정 범위다.
HSC Y3 `hsc_y3_zbin`과 별도 보정 `n(z)`를 이용한다. adapter의
`stack_s19a_tomographic_pairs`가 해당 입력을 처리하며 bin 0을 제외한다. HSC shape
안내에 따라 bin 0은 다중봉우리 P(z)로 제외되고 1–4는 Y3 cosmic-shear tomographic
bin이다. fiducial `n(z)` calibration/posterior systematic 오차 전파는 후속 분석으로 남는다. 전체 DNNz-Wide PDF
tar는 144,312,176,640 byte(134.4 GiB)이므로 포함하지 않는다. 준비 쿼리는
네트워크 실행하지 않았고 원자료 행은 아직 없다. 열·용량·검증·잔여 과학 게이트는
[shape-row export 검토서](hsc_debug/wl_shape_row_export_review_2026-10-18.md)에 있다.

### 2026-10-18 CAMIRA 후보 질량·가중 민감도

기존 800개 영상의 위성 반경별 프로필 5,832행을 무헤더 CAMIRA S20A 멤버 자료와
위치로 연결해 질량 구간 및 CAMIRA `w` 가중/비가중 평균을 비교했다. 800/800 클러스터의
좌표·richness·CAMIRA `z_cl` 및 멤버 좌표 집합이 맞았고, 프로필이 있는 543개
클러스터의 멤버 인덱스는 저장 좌표를 거쳐 연결했다. 파일 수집 코드와 SHA-일치
`member_summary.json`에서 열 8은 `weight`로 확인돼 CAMIRA `w`로 사용한다. 열 7은
`field7`로만 기록되며 `log10(M*/M_sun)` 해석은 값 범위에 근거한 후보 질량 민감도다.
CAMIRA 가중치는 질량·공간 필터를 포함하므로 비편향 위성 평균으로 취급하지 않는다.

비가중 결과 14개 필터×반경 값은 기존 스택 요약을 최대 `3.9e-8 nJy/arcsec²` 이내에서
재현했다. 50–75 kpc에서 후보 질량 11.0–11.5 구간의 평균은 i `4.63` (21 위성/19개
은하단), i2 `4.88` (62/49) nJy/arcsec²였고, i2 최고 질량 구간은 7/7로 지지가 너무
작다. 이 값들은 마스크·선택·PSF/sky·공분산 미교정 탐색 결과이며 물리적 질량 추세나
ICL 연관성 검출이 아니다. 세부 감사·제한사항은
[CAMIRA 멤버 가중 프로필 감사](hsc_debug/member_weighted_profile_audit_2026-10-18.md),
수치 산출물은 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_member_weight_sensitivity_2026-10-18_04/`에
있다. 계산기와 회귀시험은 `hsc_debug/analyze_member_weighted_satellite_profiles.py`와
`hsc_debug/test_analyze_member_weighted_satellite_profiles.py`다.

확인된 CAMIRA `weight`를 사용해 50–75 kpc 위성 외곽광과 200–300 kpc ICL의
부분상관도 비가중·가중으로 비교했다. K2.5/40% 탐색 분기에서 조정 변수를 z와
richness로 고정하고 연결된 5도 하늘 블록을 10,000회 함께 재표집했다. i의 상관은
0.120에서 가중 후 0.090으로, i2는 0.234에서 0.174로 낮아졌다. 가중−비가중 차이의
95% 구간은 i `[-0.285, 0.093]`, i2 `[-0.187, 0.100]`으로 모두 0을 포함한다.
이는 사후 민감도 분석이며 다중비교 보정이나 확정 검출이 아니다. 상세 표와 한계는
[CAMIRA 가중 위성광–ICL 민감도](hsc_debug/camira_weighted_icl_association_2026-10-18.md),
입력·코드·결과 해시는
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/weighted_satellite_icl_association_2026-10-18_02/summary.json`에 있다.

다음에는 이 정확한 S20A 무헤더 파일의 7열 배포 스키마를 확인하고, CAMIRA 후보
질량을 실제 다중대역 위성 측광과 완전도·오염 보정으로 검증한다. 그 뒤에만 질량별
공동 위성–ICL 주입/회수 및 사전 고정 표본 상관을 진행한다.

2026-10-18에는 필터별 부모 은하단과 ICL 반경 지원 표본을 대조했다. K2.5/40%의
네 구간 모두를 지원하는 strict 표본은 137개(i 44, i2 93), 200–300 kpc 단일
annulus 표본은 224개(i 71, i2 153)다. 각 대역의 부모와 비교한 z SMD는 각각
strict i/i2 `+0.216/+0.132`, annulus i/i2 `+0.089/+0.103`; richness SMD는
strict `−0.194/−0.218`, annulus `−0.108/−0.153`이다. 선택 표본 모두 부모보다
평균 z가 높고 richness가 낮다. strict에서도 중심 마스크 가림이 i 63.6%, i2
66.7%이고, annulus 표본은 i 69.0%, i2 77.8%다. 따라서 반경별 지원 표본을 같은
모집단으로 간주하거나 두 대역을 무보정 결합하지 않는다. 이 기술 감사는 selection
correction이 아니며 `science_measurement=false`를 유지한다. 해시·전체 결과는
[ICL 반경 표본 선택 감사](hsc_debug/icl_sample_selection_audit_2026-10-18.md)와
`sample_800_2026-10-13/icl_sample_selection_audit_2026-10-18_02/`에 있다.

같은 외곽 annulus 224쌍을 자동 유일 중심 후보이면서 K2.5 중심 픽셀이 가려지지 않은
은하단과 나머지로 나눠, 대상·대조장 연결 5도 공간 블록 10,000회 부트스트랩으로
민감도를 확인했다. i 중심 적격 표본은 22개에서 4.53 [2.31, 6.60], 나머지 49개는
2.64 [1.17, 4.00] nJy/arcsec²였다. i2는 적격 34개에서 1.99 [−0.85, 5.03], 나머지
119개에서 4.30 [2.89, 5.58]였다. 연결 블록을 하나씩 제외한 평균 범위는 모든 그룹에서
양수였지만, i2 중심 적격 구간은 bootstrap 구간이 0을 포함한다. 필터마다 중심 적격
분할 결과가 달라 자동 중심 품질에 대한 단일 보정 방향을 정할 수 없다. 사후
부분집합 분석이며 `science_measurement=false`다.
[중심 품질 민감도 보고서](hsc_debug/icl_center_quality_sensitivity_2026-10-18.md)와
해시 포함 산출물 `sample_800_2026-10-13/icl_center_quality_sensitivity_2026-10-18_02/`을
참조한다.

ICL 200–300 kpc와 위성 50–75 kpc 고리의 기하학적 중첩을 줄이기 위해 위성의 은하단
중심 반경 125–375 kpc 범위 여부로 기존 74개 cluster-filter 연관 표본을 층화했다.
5도 연결 하늘 블록 5,000회 bootstrap에서 i2 부분상관은 중첩 가능 표본 45쌍/28블록에서
`r=0.280 [0.042, 0.538]`, 기하 비중첩 표본 33쌍/22블록에서 `r=0.043
[-0.318, 0.345]`였다. 두 집단 모두 있는 같은 25쌍/21블록에서 상관 차이는
`r_overlap-r_no-overlap=0.143 [-0.419, 0.759]`로 0을 포함한다. 따라서 중첩 가능
집단의 더 큰 점추정치는 같은 픽셀 누출 가능성과 일치하지만, 기하층 간 차이를 검출하지
못했다. 이 탐색 분석은 mask-valid 픽셀, PSF/sky 전달, 픽셀 공분산을 측정하지 않으며
전체 연관성을 보정하지 않는다. 상세 결과와 hash는
[기하 중첩 층화 감사](hsc_debug/icl_satellite_geometric_strata_2026-10-18.md)와
`sample_800_2026-10-13/satellite_icl_geometric_strata_2026-10-18_01/`에 저장했다.

후속 실제 픽셀 감사를 800개 기존 영상에서 공통 지원을 갖는 63개 은하단·162개 위성에
수행했다. ICL `MASK_K2P5`와 위성 SEP/Kron·분할 마스크를 각각 적용하자 115개 기하
중첩 중 114개에서 같은 유효 픽셀이 남았다(i 24/30, i2 90/132). 위성 유효 외곽 고리
픽셀 가운데 공통 픽셀 비율의 중앙값은 i 36.0%, i2 33.4%였다. 이는 두 영상 측정이
독립 픽셀 표본이 아님을 직접 확인하지만 광량 공분산을 측정하거나 상관을 교정하지는
않는다. `science_measurement=false`를 유지한다. [실제 픽셀 중첩 보고서](hsc_debug/satellite_icl_pixel_overlap_2026-10-18.md),
해시 포함 산출물 `sample_800_2026-10-13/satellite_icl_pixel_overlap_2026-10-18_02/`.

방위각 대조 고리의 음의 가중치를 포함한 명목 잡음 공분산도 같은 63개 쌍/162개 위성에서
계산했다. HSC 분산맵을 독립 픽셀로 가정할 때 i 명목 픽셀 잡음 상관 중앙값은 −0.0112
([16,84백분위] −0.0209,+0.0080), i2는 −0.0023 (−0.0243,+0.0058)였다. 63쌍 중
36쌍에서 공분산 부호가 음수였다. 이는 위성 control 고리가 단순 target 중첩의 양의
기여를 상쇄할 수 있음을 보이나, HSC coadd 상관잡음과 ICL random-control 항이 빠진
명목값이다. 실제 total covariance가 작거나 음수라고 단정하지 않으며 `science_measurement=false`를
유지한다. [공통 픽셀 공분산 보고서](hsc_debug/satellite_icl_shared_pixel_covariance_2026-10-18.md),
해시 산출물 `sample_800_2026-10-13/satellite_icl_shared_pixel_covariance_2026-10-18_01/`.

## 7. 실행과 장기 유지

연구 설계·코드는 OGFinder, 카탈로그·FITS·실행 산출물은
`/gpfs/kjhan/Astrafek/data`에 둔다. CPU 분석은 사용자 지시대로 `lageunha`에서
수행하며 syntax 로그인 노드에서 장시간 계산하지 않는다. 기존 다운로드
스크립트와 재개 기록을 재사용한다. 인증정보는 연구 설정·로그·문서에 저장하지 않는다.

다음 세션은 이 문서와 [실행 계획](implementation_plan.md)의 최신 항목을 읽고,
자료 상태를 확인한 뒤 미완료 단계를 이어간다. 장기 주제 등록은 완료된 관측
측정이나 백그라운드 자동 실행을 뜻하지 않는다. 2026-10-18 familywise null 검정은
같은 영상을 재사용했으므로 독립 재현을 대신하지 않는다. 다음 우선 작업은 사전 고정한
i2 독립 표본/null 검정, 위성 외곽광–ICL 공동 주입/회수, 잔류 천체 마스크와 PSF/sky
전달 감사다.
마스크·반경 정의와 공동 측정 기준을 고정하기 전에는 민감도 표에서 science cut을 선택하지 않는다.
문헌·데이터 변경 시 날짜와 근거를 기록한다.

## 8. 초기 문헌과 차별화 검토

- [Brown et al., 2024, Horizon-AGN ICL assembly](https://arxiv.org/abs/2409.10607):
  시뮬레이션에서 박탈·병합·유입 전 처리가 ICL에 기여한다. 단일 기원 가정을 피한다.
- [Martin et al., 2026, satellite demographics and ICL](https://arxiv.org/abs/2601.18693):
  유입 위성 집단과 ICL 분포를 연결한 이론 연구다. ICL–위성 연결 자체를 최초라고
  주장하지 않으며 현재 생존 집단과 유입 집단을 구별한다.
- [George et al., 2025, environment and quiescent galaxy sizes](https://arxiv.org/abs/2503.22788):
  HSC 기반 크기·환경 비교의 선행연구다. 외곽 프로파일과 공동 관측 편향 검증의
  추가 기여를 구체화해야 한다.
- [Watkins et al., 2024, low-surface-brightness sky subtraction](https://arxiv.org/abs/2401.12297):
  희미한 광원과 산란광이 하늘 추정에 주는 영향을 검토한다. 배경 처리 검증의 근거다.

위 목록은 초기 검토이며 완전한 문헌 조사나 독창성 확정이 아니다. 후보 차별점은
같은 관측 표본에서 ICL과 위성 외곽을 공동 측정하고 광량 분배·완전도·공분산을
검증하는 데 있다.

### 2026-10-19 weak-lensing source-row QA 보완

승인된 800 렌즈·3,069,544 lens-source 행을 `lageunha`에서 스트리밍 감사했다.
원 논문의 shape-selection 조건 `|e|<2`, `0<=sigma_e<=0.4`, `R2>=0.3`,
`mag_A<=25.5`에 맞춰 어댑터를 수정했다. 기존 구현은 `|e|>=1`을 오류로 막아
88,460개의 정상 범위 후보를 포함한 전체 스택을 중단시킬 수 있었다. 수정 뒤
어댑터 단위·합성 회수 시험 17개가 통과했다. 실제 행 기준 네 shape cut과
tomographic bin 1–4를 통과한 pair는 2,234,510개이며, B-mode와 radial 선택 전
기술 지원량이다. 원자료 해시, counts, 선택 세부값은
[shape-row quality 보고서](hsc_debug/wl_shape_row_quality_assessment_2026-10-19.md)에
있다. `b_mode_mask` 0/1 방향은 여전히 미확정이고 보정 SACC `n(z)`, random-lens
차감, source overlap을 반영한 하늘 블록 공분산도 필요하다. 따라서 보정 전단·질량
결과는 아직 계산하지 않았으며 `science_measurement=false`를 유지한다.

### 2026-10-19 raw tangential ellipticity 진단

`audit_wl_raw_shape_profiles.py`가 이미 반출된 3,069,544 lens-source 행에 HSC
Y3 shape cuts와 tomographic bin 1–4 선택을 적용하고, proper radius 0.1–3 Mpc의
raw shape-weighted `e_t`와 cross `e_x` 프로필을 계산했다. 2,234,510개 쌍이
shape·bin 선택을 통과했다. 무마스크 기술 진단에서 `e_t`는 0.1–0.2 Mpc의
0.03434에서 2–3 Mpc의 0.005649로 감소했고, cross 성분은 두 구간에서 −0.000757,
+0.000075였다. 두 가능한 `b_mode_mask` 극성의 제외 결과도 나란히 만들었으며,
반대 극성의 pair count가 무마스크 count를 정확히 분할함을 검사했다. raw profile은
양의 접선 정렬 패턴과 양립하지만 shape response·선택편향·source `n(z)`와 전경
희석, random-lens 차감, source 중복/공간 공분산을 보정하지 않았으므로 렌즈 신호나
질량 검출로 해석하지 않는다. 결과는
[raw shape profile 평가](hsc_debug/wl_raw_shape_profile_assessment_2026-10-19.md)와
`sample_800_2026-10-13/wl_raw_shape_profiles_2026-10-19_02/`에 저장했다.
`science_measurement=false`를 유지한다.

### 2026-10-19 ICL 하늘 추정기 짝비교

Astrafek 도구 `hsc_debug/audit_icl_sky_estimator_delta.py`로 저장된 74개
matched pair의 pixel-SFD와 plane-first ICL 상관을 같은 공간 블록 재표집에서
비교했다. 사전 지정한 i2 표본 53쌍에서 z·richness 조정 partial r은 0.239에서
0.124로 바뀌었고, 짝지은 차이(plane-first−pixel-SFD)의 95% 구간은 −0.329–0.063으로
0을 포함한다. 따라서 estimator 민감도는 관찰되지만, estimator 간 차이나 물리적
ICL–위성 연관을 확정하지 않는다. 입력 행과 필터 배정은 바꾸지 않았고 i 필터 21쌍은
보조 기술 결과로만 계산했다. 재현 스크립트, 검사, 결과 기록은
[sky-estimator 보고서](hsc_debug/icl_sky_estimator_delta_assessment_2026-10-19.md)와
`sample_800_2026-10-13/icl_sky_estimator_delta_2026-10-19_01/summary.json`에 있다.
다음 과학 게이트는 실제 확산광 주입·회수, 대응 random field, PSF wing 및 mask 누출을
포함한 하늘 추정 비교와 완전 공분산 검증이다. `science_measurement=false`를 유지한다.

### 2026-10-19 matched random-field plane 전달 주입

Astrafek 도구 `hsc_debug/audit_icl_plane_injection.py`와
`hsc_debug/run_icl_plane_injection.sh`가 74개 target/control matched pair의
coadd/bg 영상에 네 가지 확산광 모형을 주입하고, 200–300 kpc 고리의 raw 및
600–900 kpc tile-plane 보정 회수율을 계산하도록 구성됐다. K2.5 마스크와 고정된
대상·대조장 선택을 사용하며, 5,000회 공간 블록 재표집으로 target-control 회수 차이를
평가한다. 이는 전달받은 coadd/bg 위에서 Astrafek의 후단 plane 추정기만 검사하므로
HSC upstream sky 처리, PSF 날개, mask growth와 상관잡음 검증을 대신하지 않는다.
첫 Slurm job 427839는 74쌍 측정 후 요약 그룹의 행 수 검사 오류로 종료됐다.
Astrafek에 shape별 고정 pair 수 regression test를 추가하고 job 427845로 재실행해
완료했다. i2 53쌍에서 plane 보정 평균 회수율은 100 kpc 원형/타원 모형에 99.1%/99.5%,
250 kpc 모형에 86.0%/90.8%였다. target−random-field 회수 차이의 95% 공간 블록 구간은
네 모형 모두 0을 포함했다. 따라서 downstream plane 제거는 확산광 크기에 따라
약 0.5–14%를 흡수하지만, 이 matched sample과 모형에서 target-control 차등 전송은
검출되지 않았다. 이는 관측 보정값이 아니며 HSC upstream sky subtraction, PSF 날개,
mask growth, photon noise 및 전체 공분산을 검증하지 않는다. 원자료는 변경하지 않았고
`science_measurement=false`를 유지한다. 결과·해시는
[plane 주입 보고서](hsc_debug/icl_plane_injection_assessment_2026-10-19.md)와
`sample_800_2026-10-13/icl_plane_injection_2026-10-19_02/`에 저장했다.

### Weak-lensing 교정 입력

약한 중력렌즈용 fiducial n(z) 입력 정합성과 두 B-mode 극성 가정의 예비 stack은 끝났지만,
과학용 calibrated 분석은 아직 완료되지 않았다. HSC S19A shape 문서는
4개 source tomographic bin의 `n(z)`를 포함하는 공식 자료
`hsc_y3_fourier_space_data_vector.sacc`를 안내하고, point/photo-z 기반 추정에는
dilution 및 selection-bias 교정이 필요하다고 명시한다
([HSC S19A shape-catalog 문서](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s19a-shape-catalog-pdr3/)).
사용자의 명시 승인에 따라 이 SACC 파일만 다운로드했다. 파일은 승인된 shape-row
export와 동일한 sample data 경로에 보존했고, 실제 tracer와 `source_zbin`의 빈 수를
대조했다. 추가 영상·원자료 shape row·개별 photo-z PDF 반출은 이 단계에 포함하지
않았다. 2026-10-19에 승인된 800 렌즈/3,069,544 shape 행과 SACC를 사용해 HSC의
shape·selection 보정을 적용한 두 개의 B-mode mask bit 극성 민감도 stack을 계산했다.
raw 0을 masked로 보면 2,153,004 pair/609 lens, raw 1을 masked로 보면 76,864 pair/26 lens가
0.1–3 Mpc에 남는다. 첫 가정의 weak-shear `DeltaSigma`는 0.1–0.2 Mpc에서
161.5±11.5, 2–3 Mpc에서 15.3±0.9 `Msun/pc^2`이고, 두 번째 가정은 각각
58.1±61.1, 12.8±4.9 `Msun/pc^2`이다. 이 오차는 독립 pair shape noise만 나타내며
서베이 공분산이 아니다. 외관상 profile이 더 자연스러워 보인다는 이유로 bit 극성을
선택하지 않았다. 전체 출력과 한계는
[B-mode 극성 민감도 보고서](hsc_debug/wl_bmode_polarity_sensitivity_2026-10-19.md)에
있다. random-lens 차감, source overlap/공간 공분산, n(z) 오차 전파 및 reduced-shear
보정이 남아 `science_measurement=false`를 유지한다.

HSC PDR3 안내에서 random points가 공개 데이터베이스에 있다고 확인했다
([HSC PDR3 available data](https://hsc-release.mtk.nao.ac.jp/doc/index.php/available-data__pdr3/)).
다음은 그 table schema와 800 렌즈 표본에 맞는 footprint/depth/PSF 조건을 확인하는
count-only 질의 초안을 준비하는 단계다. random position row 조회·반출은 기존 다운로드
승인 원칙에 따라 별도 승인을 얻기 전에는 실행하지 않는다.

### 2026-10-19 반경별 ICL–위성 공간 블록 불확실도

기존 94개 반경별 matched cluster/filter 행에 ICL 대상 좌표와 검증된 대조장 추출 좌표를
복구해 connected target/control 5도 공간 블록 재표집을 적용했다. Astrafek 도구
`hsc_debug/audit_satellite_icl_spatial_blocks.py`는 5,000회 일반 pair 재표집과 공간 블록
재표집을 같은 고정 표본에서 비교한다. 주 i2 68쌍/33블록의 redshift·richness 조정 partial
Pearson r은 0.293, 공간 블록 95% 구간은 [0.042, 0.497]이지만, 일반 pair 구간과 보조
Spearman 공간 구간은 0을 포함한다. 작은 중심 적격·i 표본도 구간이 넓다. 선택 표본의
사후 탐색 결과이므로 검출로 해석하지 않고 `science_measurement=false`를 유지한다.
좌표 결합 94/94, 입력 및 결과 해시, 테스트 3개를 확인했다. 세부 결과와 한계는
[공간 블록 평가](hsc_debug/satellite_icl_spatial_block_assessment_2026-10-19.md) 및
`sample_800_2026-10-13/satellite_icl_spatial_blocks_2026-10-20_02/`에 기록했다.
