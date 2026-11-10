# 장기 연구주제: HSC ICL·위성 외곽구조·은하-은하 약한중력렌즈

등록일: 2026-10-13 (Asia/Seoul)<br>
상태: **사용자 승인 장기 연구주제 / 관측 타당성 검증 단계**<br>
프로젝트: Astrafex Desktop / OGFinder<br>
연결: [실행 계획](implementation_plan.md), [분석 도구 설계](research_software_design.md)
최신 ICL·위성·약중력 근거 상태: [2026-11-06 연구 상태 평가](hsc_debug/icl_wgl_research_status_2026-11-06.md)

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
| B3. 약한렌즈 예비 측정 | S19A shape source 겹침표·전단 profile | response·n(z)·random subtraction·B-mode·공분산 및 matched-field systematics 검증 | HSC 응답·선택보정과 4-bin HSC Y3 `n(z)` 적용. 좌표를 이용해 `b_mode_mask=0`의 제외 방향을 독립 검증. 보정된 609-lens cluster-only stack, 공통지원 435 lenses/2,175 matched random 중심의 차감 stack, source-density/estimator-weight boost 및 2.5°/5° covariance/null 분석 완료. Tangential은 0과 불일치하고 cross는 0과 양립하나 공분산·영역 선택 민감성이 남음. Full-aperture random depth/mask/PSF/sky 일치, `n(z)` 불확도·dilution, reduced shear, magnification 및 contiguous-field covariance 미완료; science gate 미통과 |
| C. 공동 검증 | 주입/회수, null, 공분산·감도 보고서 | 무상관 주입의 가짜 관계와 알려진 신호 회수 오차가 사전 오차 예산 내 | 137쌍 확산광 진폭 주입과 target/control coadd·위성–ICL finite-lag covariance 분석 완료. Joint 94-pair IPW+pixel-noise 분석도 수행. 절대 ICL 회수율은 약 23%이고 상류 sky 전달·source-Poisson·확장 PSF 불확실성은 미검증. WL random/cross null과 source reuse/boost covariance 분석 완료하였으나 matched-field systematics와 covariance estimation uncertainty는 미해결; 게이트 미통과 |
| D. 분석 정의 고정 | 표본/반경/단위/가중/모형/판정 기준 버전 | 본 분석 전 숫자 기준과 주·보조 가설 고정 | 미실행 |
| E. 전체 표본 | 두 스택, 공동 관계, 상한·민감도 | 효과와 체계오차, 참여 표본 및 선택함수 함께 보고 | 탐색 스택은 수행: ICL 동일필터 target/control 783쌍, 반경별 ICL 지원 224쌍, 위성 50–75 kpc 지원 94쌍, random-subtracted WL common-support 435 lenses, ICL–위성–WL i2 교집합 27 lenses. IPW·noise·boost·반경·block-width 민감도 기록됨. selection/transfer/calibration을 통과한 confirmatory full-sample 분석은 미실행 |
| F. 물리 해석 | 관측을 재현한 시뮬레이션 비교, 논문 | 박탈·병합·선택 설명을 구분 가능한 범위에서 평가 | 미실행 |

30–50개는 예비 설계 규모다. 800개의 통계력이나 검출 유의도는 아직 보장하지
않는다. 단계 B/C에서 편향·산포를 구해 필요한 유효 은하단 수와 세분화 수준을
정한다. 잔류 편향 허용치 등 수치 기준은 효과를 본 뒤 유리하게 바꾸지 않는다.
`science_measurement=false`는 주장별 검증을 통과하기 전까지 유지한다.

단계표의 약중력 상태는 2026-11-06 기준으로 갱신했다. B-mode 좌표 판정과
random subtraction이 완료된 점은 독립 기록
[`wl_bmode_coordinate_validation_2026-10-19.md`](hsc_debug/wl_bmode_coordinate_validation_2026-10-19.md),
[`wl_source_density_boost_joint_covariance_assessment_2026-11-05.md`](hsc_debug/wl_source_density_boost_joint_covariance_assessment_2026-11-05.md),
[`wl_random_shape_stack_detection_assessment_2026-11-05.md`](hsc_debug/wl_random_shape_stack_detection_assessment_2026-11-05.md)에
남아 있다. Cross null이 통과했다는 사실은 모든 systematics가 닫혔다는 뜻이 아니다.

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
5–75 kpc 표면밝기 프로필을 계산했다. 각 위성과 같은 은하단 중심 투영 반경에
놓인 방위각 회전 위치를 대조 프로필로 삼았다. 32개 후보 방위각 중 검출원 3 arcsec
이내 또는 CAMIRA 멤버 15 kpc 이내인 위치를 제외하고 최대 8개를 사용했으며,
위성 프로필에서 이 대조 프로필 평균을 뺐다. 따라서 평균적인 국소 확산광은 줄지만
비대칭 ICL/은하단 하위구조, 남은 천체원 빛과 PSF 날개까지 분리·제거한 측정은 아니다.
엄격한 7개 반경 동시 커버리지 조건을 적용한
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
공식 random-points 설명에서 `pdr3_wide.random`, 점밀도 100/arcmin²,
`adjust_density`, patch/tract interior 및 HSC-I mask/depth/PSF 열을 확인했다
([HSC random-points 설명](https://hsc-release.mtk.nao.ac.jp/doc/index.php/random-points__pdr3/)).
Astrafek 스크립트 `hsc_debug/prepare_wl_random_control_counts.py`는 승인된 800 렌즈
preflight에서만 32개 count-only SQL batch를 생성한다. 렌즈별 3 Mpc 렌즈 지름보다
떨어진 local annulus (2×각 aperture–30′)의 `adjust_density<0.01`, primary, patch/tract
interior, i-band non-edge/non-bright-object-center random 수만 집계한다. 출력 열은
`lens_id,n_random_candidates`뿐이며 좌표·shape/photo-z 행을 선택하지 않는다.
계획은 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/`
`wl_random_point_count_plan_2026-10-19_01/`에 있으며 `prepared_not_submitted`다.
이는 후보 수 가용성 검사일 뿐 random center 선정이 아니다. 아직 CAS 질의를 제출하지
않았다. 질의 결과가 승인되더라도 좌표 반출은 별도 승인을 받아야 하고, 이후 전체
CAMIRA 중심 제외·완전 aperture footprint·HSC-I depth/PSF·중복/공간 공분산을 검증한다.

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

### 2026-11-03 약중력 random-center 지원 관문

이미 승인된 random 좌표와 기존 렌즈별 shape-count 계획을 재사용해, CAMIRA 중심과
3 proper-Mpc aperture가 겹치지 않는 중심 중 lens-specific 재사용을 허용하는 표본을
평가했다. 620개 렌즈에 중심 5개씩, 총 3,100개 lens–center 연결과 3,020개 고유 좌표를
선정했다. 77개 좌표는 여러 렌즈에 공유되므로 최종 공분산에서 공유 위치를 보존해야 한다.
각 중심의 HSC S19A shape 객체 수만 `COUNT(DISTINCT object_id)`로 집계하는 124개 SQL
배치(배치당 최대 25개)를 만들고 해시·열·집계 표현을 오프라인에서 대조했다. 결과 계획은
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/`
`wl_random_shape_count_plan_lens_reuse_2026-11-03_01/`에 있으며 총 SQL 크기는 1,069,022
bytes다. shape 값·source 행·photo-z·영상은 반환하지 않는다. 이 집계는 기존 800개 렌즈
count 승인과 다른 데이터 반출 범위이므로 사용자 명시 승인이 아직 필요하며, HSC CAS에
연결하거나 질의를 제출하지 않았다. count 결과가 확보되더라도 random shape-row export와
cluster-minus-random stack은 별도 지원·전송량·공분산 검증 뒤 진행한다. 약중력 결과는 계속
cluster-only 진단이고 `science_measurement=false`다. 상세한 재사용 민감도는
[`wl_random_coordinate_reuse_assessment_2026-11-03.md`](hsc_debug/wl_random_coordinate_reuse_assessment_2026-11-03.md),
격자 크기별 잭나이프 민감도는
[`wl_spatial_jackknife_cell_scale_assessment_2026-11-03.md`](hsc_debug/wl_spatial_jackknife_cell_scale_assessment_2026-11-03.md)에
정리했다.

### 2026-11-04 matched random-center weak-lensing stack

승인된 HSC S19A shape 행 반출 124개 배치가 완료되어 3,100 random center,
620개 렌즈 연결, 9,010,324개 center-source 행(2,962,522,918 bytes)을 받았다.
기존 800-lens cluster stack을 우선 재현한 뒤, target과 random 양쪽에 source support가
있는 435개 렌즈와 그에 연결된 2,175 random center(양의 source 수 2,046,
0-count 129)만 사용했다. 동일 표본의 random-subtracted DeltaSigma와 paired
5-degree equal-area spatial jackknife 공분산을 구성했고 51개 하늘 블록이 점유됐다.
접선 프로필은 양수이며 반경과 함께 감소하고, 7-bin cross 성분 null 검정은
chi-square=10.594/7 dof, p=0.157이다. 결과는 탐색용이고 `science_measurement=false`다.
random aperture 전체의 depth/mask/PSF/sky 일치, 연결된 HSC footprint 공분산,
블록 크기 민감도, n(z) 보정 오차와 reduced-shear 보정이 남아 질량 검출로 해석하지
않는다. 원본 summary SHA-256은
`83c5484454033b2652adcba5bd0699f28082fa6983eefa960613da75491bafc8`이며 상세 수치,
제한 및 그림 재현 명령은
[hsc_debug 약중력 평가](hsc_debug/wl_random_shape_stack_assessment_2026-11-04.md)에 있다.

같은 날짜의 ICL coadd 헤더 감사는 기존 target/control 1,598장 모두에서 IMAGE,
MASK, VARIANCE 평면을 확인했으나, synthetic source Poisson variance를 유도할 gain 또는
노출량 헤더는 찾지 못했다. 기존 관측 영상에 임의의 shot-noise 보정을 적용하지 않는다.
ICL–위성 주입/회수 결과는 이 노이즈 한계와 coadd 공분산, matched-control sky transfer,
PSF wing 및 mask 누출 검증이 끝날 때까지 `science_measurement=false`로 유지한다.
재현 스크립트와 해시 기록은
[hsc_debug ICL Poisson 입력 감사](hsc_debug/icl_coadd_poisson_noise_assessment_2026-11-04.md)에 있다.

### 2026-11-05 약중력 하늘 블록 크기 민감도

기존 5도 equal-area paired endpoint jackknife(51 blocks)와 동일한 435개 공통지원
렌즈·2,175개 random 중심·8,905,892 random source pair를 사용해 2.5도 분할(154 blocks)을
추가 계산했다. 접선 ΔΣ 점추정은 저장 정밀도 내에서 같고, 대각 jackknife 오차 비율은
2.5도/5도에서 반경별 0.820–1.067이었다. 1.5–2 Mpc 오차가 18% 감소하는 등 공분산은
격자에 일부 의존한다. 7-bin cross null은 5도 χ²=10.594 (p=0.157), 2.5도 χ²=7.294
(p=0.399)로 바뀐다. 두 14×14 covariance는 full-rank 양의 정부호지만 이는 올바른
coverage를 보장하지 않는다. 분석기 full-recompute 방식으로 2.5도 실행에 약 90분이 걸려
더 넓은 격자 탐색 전 sufficient-statistic 기반 최적화가 필요하다. 무작위 중심의 전체
aperture depth/mask/PSF/sky 일치와 연결된 HSC footprint 공분산이 검증되지 않아
`science_measurement=false`를 유지한다. 수치·해시·재현 정보는
[hsc_debug 블록 민감도 평가](hsc_debug/wl_random_shape_jackknife_sensitivity_assessment_2026-11-05.md)에
정리했다.

### 2026-11-05 ICL–위성 근거 통합과 frame metadata 제출기

12개 hash-bound 진단 JSON을 통합하는 최신 보고서는 74쌍 robustness 분석과 전체 94쌍
결합지원 분석, 지원 선택가중·PSF 주입표본 대표성·frame-count 민감도를 함께 제시한다.
i2 부분상관은 74쌍에서 `r=0.234`, 94쌍에서 `r=0.293`이며,
paired 94-minus-74 구간 `[-0.188, 0.343]`은 0을 포함한다. 따라서 표본 확장이 상관을
유의하게 바꿨다는 근거는 없다. 보고서 생성기
[`build_icl_association_evidence_report.py`](hsc_debug/build_icl_association_evidence_report.py)는
저장 요약만 읽고 영상이나 HSC를 조회하지 않는다. 최신 Astrafek 산출물은
[통합 근거 보고서 v3](hsc_debug/icl_association_evidence_status_2026-11-05_03.md)
(SHA-256 `4443a16955a1e60461cfadf65cfa38a305be154db7fda2622971e17f4c433274`)이다.
800개 표본은 보수적인 다중비교 검정력 목표에 부족하고 관측/완전도·PSF·sky·공분산도
미해결이므로 과학 게이트는 계속 닫혀 있다. 8개 입력을 사용한 이전 통합본은
[초기 ICL association evidence status](hsc_debug/icl_association_evidence_status_2026-11-05.md)다.

ICL source-Poisson 잡음 관문을 위해 기존 target/control 1,598개 coadd 중심의
16-batch count-only HSC PDR3 mosaicframe v3 조회(13개 열, plan SHA-256
`3dac4db631d0dcbee5c5f4509e639f2bd90ca2af99fa7cf5f66c2a9f60bd9014`)는 사용자 승인 뒤
16/16 batch 완료됐다. 표적과 대조장의 frame 수 중앙값은 모두 42이며, 224개 ICL 지원
cluster/filter 쌍의 frame·visit·CCD count offset은 ICL 대비와 뚜렷한 관계를 보이지
않았다. 이는 픽셀별 depth/coadd weight/분산 일치 검증이 아니다. 다음 단계의 67,071개
center-frame gain/exposure 행 반출 계획(SHA-256
`49aa4b7601dd42243209d7b3b394a4dd3b4cbf00d87d28cbb9467553c29f4ce3`)은
16/16 batch 완료 상태로 확인됐으며, 결과와 한계는 2026-11-06 기록에 갱신했다.

### 2026-11-05 ICL–위성 94쌍 지원표본 상관 민감도

기존 74쌍 pixel-extinction 진단이 저장된 전체 94쌍 결합지원 표본을 대표하는지
검토했다. Astrafek의 `analyze_icl_satellite_supported_94.py`가 영상 픽셀이나 HSC를
조회하지 않고 pair 요약 행 및 좌표 manifest만 읽어, raw ICL 200–300 kpc와 위성
50–75 kpc 대비의 redshift·richness 조정 부분상관 및 필터별 5도 연결블록 bootstrap
5,000회를 산출했다. 전체 지원표본은 i 26쌍 `r=0.079` (95% 구간 `[-0.183, 0.477]`),
i2 68쌍 `r=0.293` (`[0.039, 0.500]`)이었다. 동일 블록을 사용한 94-minus-74 paired
차이 구간은 i `[-0.245, 0.092]`, i2 `[-0.188, 0.343]`로 모두 0을 포함했다. 이는
표본범위 민감도일 뿐 선택 편향을 해소하지 않으며, i2 점추정치 상승만으로 신호가
강해졌다고 판단할 수 없다. 보정·PSF·sky·완전도·coadd 공분산 한계가 남아
`science_measurement=false`를 유지한다. 수치·입력 및 산출물 해시는
[94쌍 민감도 평가](hsc_debug/icl_satellite_supported_94_assessment_2026-11-05.md)에
기록했다.

### 2026-11-05 약중력 source-density boost 관문

기존 matched-random subtraction 이후 아직 적용하지 않은 cluster-member source
dilution의 크기를 평가하기 위해, 동일한 HSC S19A 선택을 거친 source pair 수를
435개 공통지원 렌즈와 2,175개 matched random 중심으로 정규화했다. paired delete-one
equal-area sky-cell 재표본화는 5도 52개 및 2.5도 159개 endpoint block으로 수행했다.
count ratio는 0.1–0.2 Mpc에서 `1.765 ± 0.064`, 2–3 Mpc에서 `1.148 ± 0.030`이었다.
초기 `beta^2` 계산에 렌즈 거리 인자가 빠진 것을 발견해 보정했다. 정확한
`w_shear SigmaCritInv^2` response-adjusted ratio는 같은 구간에서 각각
`1.139 ± 0.048`과 `1.126 ± 0.031`으로 count ratio보다 작다. 수밀도 초과만으로
렌즈 추정량 보정치를 정할 수 없다는 점을 확인했다. 오류는 선택한 공간 셀의 jackknife
표준오차다. Astrafek은 네 boost 지표와 일곱 반경 빈, 무작위 차감 tangential/cross
`DeltaSigma` 14개 값의 전체 42×42 공분산을 저장했다. boost–`DeltaSigma` 상관계수는
2.5도 block에서 -0.346–0.325, 5도 block에서 -0.393–0.465였다. 저장된 stack의 전체
프로필과 원래 pair-endpoint 공분산을 재현했다. 다만 coverage·mask·depth·PSF·sky 응답
일치는 검증하지 않는다. `DeltaSigma`에는 비율을 적용하지 않았고 photo-z 희석 보정도
미해결이므로 `science_measurement=false`다. 올바른 `SigmaCritInv` 정의·재현 스크립트·해시는
[joint boost/DeltaSigma 평가](hsc_debug/wl_source_density_boost_joint_covariance_assessment_2026-11-05.md)에
기록했다. dimensionless `beta^2`를 사용한 이전 평가는
[보정 전 기록](hsc_debug/wl_source_density_boost_spatial_assessment_2026-11-05.md)으로 보존하고
인용 금지 표시를 했다.

### 2026-11-05 기존 800개 표본의 ICL–위성 반경별 지원률

저장된 ICL 200–300 kpc 지원 224개 cluster/filter 쌍과 per-annulus 위성 프로필을
교차해 Astrafek의 [`audit_icl_satellite_support_yield.py`](hsc_debug/audit_icl_satellite_support_yield.py)로
반경별 지원 수를 다시 계산했다. 위성 반경이 5–10 kpc일 때 교집합은 198쌍이지만,
50–75 kpc에서는 94쌍( i 26, i2 68)으로 감소한다. 50–75 kpc 쌍 중 기존 중심 마스크
적격은 23쌍( i 8, i2 15)뿐이다. 따라서 더 안쪽 빈의 큰 표본을 외곽 구조 신호로
대체하지 않는다.

위성 후보 9,468개 가운데 6,383개는 어떤 target annulus도 60% coverage 바닥값을 넘지
못했고, 1,306개는 target annulus가 있으나 최소 세 개의 방위각 control이 부족했다.
40% floor로 같은 800개 이미지를 다시 측정하자 50–75 kpc 교집합은 94쌍에서 204쌍
(i 64, i2 140)으로 늘고 중심 마스크 적격은 23쌍에서 52쌍으로 증가했다. 하지만
i2 부분상관은 `r=0.293`에서 `0.249`로 바뀌며 40% 구간 `[-0.002,0.454]`과 paired
floor 차이 `[-0.252,0.169]`가 모두 0을 포함한다. 중심 마스크 적격 부분집합도 40% 기준
i2 30쌍 `r=0.242`, 구간 `[-0.024,0.684]`; 60% 기준 15쌍 구간 `[-0.081,0.785]`로
모두 0을 포함한다. 새 control 프로필이 추가되어 기존
쌍의 표면밝기 값도 바뀐다. 따라서 표본이 커져도 estimator가 안정됐다고 볼 수 없고
40% floor는 채택하지 않았다. 낮은 coverage가 cutout 경계,
source/neighbor mask, detection geometry 중 무엇에서 주로 비롯되는지는 아직 분리하지
못했다. 우선 저장된 mask와 WCS 및 member/detection 위치를 대조해 원인을 진단하고,
cutout 크기가 원인일 때만 별도의 넓은 영상 계획을 세운다. 상세 표·제약·해시는
[support-yield 평가](hsc_debug/icl_satellite_support_yield_assessment_2026-11-05.md)에
기록했다. 이 집계는 이미지 픽셀을 읽거나 HSC에 접속하지 않았고
`science_measurement=false`다.

40% 재측정은 기존 CPU Astrafek 파이프라인을 800개 이미지에 실행했다. 동일한
support/연관 분석 코드와 2,000회 5도 연결블록 paired bootstrap 결과는
[coverage-threshold 민감도 보고서](hsc_debug/icl_satellite_coverage_threshold_assessment_2026-11-05.md)에,
재현 스크립트는 [`analyze_icl_satellite_coverage_threshold.py`](hsc_debug/analyze_icl_satellite_coverage_threshold.py)에
기록했다. 지원률과 부분상관 변화를 함께 보는 그림은
[`plot_icl_satellite_coverage_threshold.py`](hsc_debug/plot_icl_satellite_coverage_threshold.py)가
생성한다. 낮은 coverage 표본은 systematic sensitivity로만 보관하며 baseline stack에
합치지 않는다.

2026-11-05에는 이 42개 공동변수의 시각화도 추가했다. Astrafek의
[`plot_wl_boost_delta_sigma_joint.py`](hsc_debug/plot_wl_boost_delta_sigma_joint.py)는 저장된
v2 JSON에서 네 밀도비, 무작위 차감 `DeltaSigma`/cross 성분, 전체 공간 jackknife 상관행렬을
함께 그린다. [시각화 평가](hsc_debug/wl_source_density_boost_joint_visualization_2026-11-05.md)는
점추정치를 과학 검출로 해석하지 않으며 boost 보정도 적용하지 않는다. 그림 SHA-256은
`2ecfcae2a312386a09504463251c6d328598c287857624a2063c7bfb2d100d2e`다.

## ICL–위성 outer-annulus mask selection (2026-11-05)

The 40% 50–75 kpc diagnostic support is non-random in projected cluster radius.
Among 2,691 candidates in 224 ICL-supported cluster/filter frames, the fraction
retained after both target coverage and control requirements is 20.6–22.7
percentage points higher in the largest than the smallest cluster-centric
radius quartile in both filters, using cluster-ID bootstrap intervals. The
earlier mask decomposition attributes 69.4% of target-coverage failures to
source-related masks, chiefly neighboring-source masks. Neither result is a
science correction. Before interpreting the retained-sample association,
estimate multivariable selection and calibrate it with injection/recovery;
compare only predeclared radial-support/completeness analyses and carry the
selection uncertainty. Do not remove protective masks or adopt a relaxed floor
as the baseline without that calibration. Details and hashes are in
`hsc_debug/icl_satellite_mask_selection_assessment_2026-11-05.md`.

The adjusted cluster-robust logistic diagnostic retains the pattern after
controlling for detection S/N, redshift, and richness: per-SD radius odds
ratios for profile retention are 1.55 (95% CI 1.25–1.92) in `i` and 1.45
(1.28–1.64) in `i2`; per-SD log-richness odds ratios are 0.64 (0.51–0.81) and
0.81 (0.69–0.94), respectively. This descriptive model is not
injection-calibrated and cannot yet be used to weight the ICL–satellite signal.

### 2026-11-05 224→94 joint-support selection weighting sensitivity

Astrafek의 [`analyze_icl_satellite_joint_support_ipw_sensitivity.py`](hsc_debug/analyze_icl_satellite_joint_support_ipw_sensitivity.py)는 224개 ICL-supported cluster/filter pair 중 94개에서 이용 가능한 50–75 kpc 위성 프로파일의 선택확률을 추정했다. 모든 pair에 있는 filter, redshift, log richness, 200–300 kpc ICL coverage, central-mask 상태를 사용한 pooled logistic model에서 coverage는 profile support odds와 강한 양의 관계를 보였다(OR 1.86/SD). 지원표본과 비지원표본의 coverage 표준화 평균 차이는 i `+0.44`, i2 `+0.70`이었다. stabilized IPW 후 redshift·선형 richness 보정 부분상관은 i `0.079→−0.011`, i2 `0.293→0.227`로 낮아졌다. 수정된 2,000회 연결 5° 블록 bootstrap IPW 구간은 각각 `[-0.566, 0.581]`, `[-0.034, 0.416]`이며, 유효 N은 22/26, 60/68이다. 기존 BFGS의 216 precision-loss 응답은 모두 gradient norm 허용기준을 충족해 실패가 아니었고, 최대 gradient norm은 `1.40×10⁻⁶`, 선택모형 stationarity 실패는 0이었다. 중심 mask 제외, coverage 제곱항, filter별 slope를 허용한 세 대안 모델에서도 가중 부분상관은 i `[-0.014, 0.021]`, i2 `[0.226, 0.231]` 범위였다. 위성 신호를 관측하지 못한 130쌍의 조건부 무작위 결측 가정은 검증되지 않았다. 이 결과는 sensitivity only이며 completeness correction이나 과학 검출로 인용하지 않는다 (`science_measurement=false`). [평가와 한계](hsc_debug/icl_satellite_joint_support_ipw_sensitivity_2026-11-05.md), [요약 JSON](/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_joint_support_ipw_2026-11-05_06/summary.json; SHA-256 `8dbfae2169361c9fb6ce5ccfcce684b486162dbf4820e531112cdfc4f39fdbc7`)에 재현 정보를 기록했다.


### 2026-11-05 ICL–약중력 공통 표본 감사

Astrafek [`audit_icl_wl_sample_overlap.py`](hsc_debug/audit_icl_wl_sample_overlap.py)가 승인된 source-row/count exports를 로컬에서 읽어 435개 lens의 기존 random-subtracted 공통지원 ID와 저장된 SHA-256을 정확히 재현했다. 224개 ICL profile-supported pair 가운데 공통 약중력 표본은 `i` 38개, `i2` 70개다. 50–75 kpc 위성 프로필까지 갖춘 94쌍 중에는 `i` 12개, `i2` 27개만 약중력 target/random 양쪽 지원이 있다. 이 교집합은 filter별 확인 분석에 작으므로 질량–위성 상관의 검출로 해석하지 않는다. 우선 38/70 ICL 지원렌즈의 별도 random-subtracted WL 프로필을 exploratory stack으로 만들 수 있지만, 이는 위성 구조 관계 자체를 측정하지 않는다. 94쌍 공동 분석은 저검정력 진단으로 남기고, shape 행이 없는 렌즈의 원인을 footprint/선택과 대조한 뒤 추가 자료 범위를 정한다. 결과와 한계는 [공통 표본 평가](hsc_debug/icl_wl_sample_overlap_assessment_2026-11-05.md), 재현 산출물은 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_sample_overlap_audit_2026-11-05_03/`에 있다. 추가 HSC 조회는 실행하지 않았고 `science_measurement=false`다.


### 2026-11-05 ICL-supported weak-lensing subset stack

Astrafek의 [`run_icl_supported_wl_subset_stack.py`](hsc_debug/run_icl_supported_wl_subset_stack.py)가 기존 검증 shape 행, calibrated tomographic estimator, matched random centers를 사용해 ICL 200–300 kpc profile 지원과 435-lens common WL 지원의 교집합을 filter별 random-subtracted stack했다. `i`는 38렌즈, `i2`는 70렌즈다. 5°와 2.5° endpoint jackknife 사이 tangential 점추정치는 같고 median 오차비(2.5°/5°)는 1.05/1.12다. 하지만 `i` cross 성분은 두 분할에서 모두 0과 불일치(p=`8.0e-4`, `3.4e-11`)하므로 tangential 프로필을 검출로 해석하지 않는다. `i2` cross는 p=`0.42`, `0.31`로 0과 양립하고 tangential 프로필은 유망하지만 탐색적이다. 위성 50–75 kpc 프로필까지 지원되는 더 좁은 표본(`i` 12, `i2` 27)에서도 별도 WL stack을 계산했다. 이 중 `i` cross는 5°/2.5° 모두 0과 불일치(p=`2.5e-274`, `3.5e-22`), `i2` cross는 양쪽 모두 0과 양립(p=`0.282`, `0.575`)한다. joint 표본은 너무 작고 `i` null 이상이 지속돼 위성–질량 관계의 검출로 해석하지 않는다. random aperture의 depth/mask/PSF/sky 일치, contiguous HSC covariance region, photo-z 불확도, boost 및 reduced-shear 보정은 미해결이라 모두 `science_measurement=false`다. 프로필 그림, report 및 covariance null 검사는 [subset 평가](hsc_debug/icl_supported_wl_subset_stack_assessment_2026-11-05.md), 재현 runner와 `compare_icl_supported_wl_blockwidth.py`에 저장했다.

### 2026-11-05 HSC coadd frame-count v3

사용자 승인 v3 CAS 조회(계획 SHA-256 `3dac4db631d0dcbee5c5f4509e639f2bd90ca2af99fa7cf5f66c2a9f60bd9014`)에서 800 target와 798 control 중심의 정수 frame/visit/CCD count-only 결과를 받았다. 1,598곳 모두 HSC-I mosaic patch에 대응하며 중앙값은 target/control 각각 42/42 frame, 9/9 visit, 38/37 CCD다. 798 matched target-control pair에서 세 count의 중앙 차는 모두 0으로, 중심 coadd coverage의 총체적 불균형은 보이지 않는다. 이는 픽셀 깊이/분산 매칭을 뜻하지 않으며 gain, exposure, sky variance, coadd weight는 조회하지 않았다. 후속 frame-level variance/calibration 원자료는 별도 제한 계획과 승인 후에만 가능하다. 재현 auditor와 한계는 [frame-count 평가](hsc_debug/icl_frame_metadata_count_assessment_2026-11-05.md), 산출물은 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_metadata_count_audit_2026-11-05_02/`에 기록했고 `science_measurement=false`로 유지한다.

이후 Astrafek runner `hsc_debug/submit_icl_frame_metadata_rows.py`의 체크포인트에서 위 1,598개 중심·16 batch·67,071행 조회가 모두 완료된 것을 확인하고, `hsc_debug/audit_icl_frame_metadata_row_results.py`로 전체 범위를 다시 감사했다. 21개 요청 열, center/patch identity, 중복 없는 frame membership, 1,598개 중심과 798개 타깃-대조 쌍을 검증했고 exposure·gain·zeropoint·sky RMS·MJD 열 결측은 0이었다. 현재 결과 summary의 plan SHA-256은 `49aa4b7601dd42243209d7b3b394a4dd3b4cbf00d87d28cbb9467553c29f4ce3`, 16/16 batch complete다. 기존 결과 체크포인트를 사용했으며 이번 재검증 과정에서 새 CAS job은 제출하지 않았다. 이미지, shape, photo-z, object catalog 및 per-pixel weight는 받지 않았다. 이 frame 메타데이터만으로 pixel-level coadd variance를 복원할 수 없으므로 주입 잡음 보정이나 과학 측정 승격은 하지 않는다. 검증본은 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_metadata_row_audit_2026-11-06_02/summary.json`에 있다.

HSC frame-count offset이 저장된 ICL target-control 신호와 연관되는지 별도 점검했다. 200–300 kpc supported 224 pair에서 frame/visit/CCD 수 차이 각각과 ICL 표면밝기 차이의 partial Pearson 상관을 filter별로 계산해 redshift, log richness, common aperture coverage를 조정했고, 5° 연결 target/control sky block을 5,000회 재표집했다. i (71 pair/18 block)는 r=`−0.085, −0.061, +0.078`, i2 (153/40)는 `+0.017, −0.041, −0.008`; 여섯 구간이 모두 0을 포함한다. 이는 정수 count offset이 이 지원 표본에서 ICL 대비를 설명한다는 증거를 보이지 않지만, 픽셀별 깊이·분산·coadd weight를 통제한 것은 아니다. 표본 선택 및 작은 i 표본 제한을 두고 `science_measurement=false`를 유지한다. runner `hsc_debug/analyze_icl_frame_count_icl_sensitivity.py`, 평가 `hsc_debug/icl_frame_count_icl_sensitivity_assessment_2026-11-05.md`, summary `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_count_sensitivity_2026-11-05_01/summary.json` (SHA-256 `dd00410649e230359df6a6ea6582be2df6881e247bf81d29347bdca33d569cca`). 분석은 저장된 표만 사용했고 HSC query/image pixel read는 없었다.

### 2026-11-06 ICL annulus random-point variance 지원

기존 target 중심당 200개(30′ 원뿔)씩 저장한 PDR3 random point 160,000개에 대해 `hsc_debug/audit_icl_existing_random_noise_support.py`로 proper-kpc 반경 지원을 확인했다. 800개 target 중 200–300 kpc ICL 고리에는 총 109점만 들어왔고, 98개 target에 적어도 하나, 어느 target도 10개 이상 점을 갖지 않았다. 50–75 kpc 위성 고리에는 3점뿐이었다. 따라서 기존 random 표본의 `i_pix_variance`는 두 목표 반경을 위한 target/control 짝지은 보정에 쓸 수 없다. 평가 문서 `hsc_debug/icl_random_variance_support_assessment_2026-11-06.md`, 스크립트와 요약 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_existing_random_noise_support_2026-11-06_02/`에 남겼다. 첫 `_01` 결과는 각도-물리 반경 변환 오류를 발견해 무효화했고 `_02` 결과만 사용한다.

다음 제한 HSC 데이터 단계로, target 800개와 matched control 798개의 200–300 proper-kpc 고리에서 HSC PDR3 `pdr3_wide.random`의 `i_pix_variance`, `i_sky_std`, `i_inputcount_value`, local PSF moments 및 mask flags를 추출하는 SQL 계획을 Astrafek에 준비했다. `adjust_density < 0.2`의 명목 밀도 20점/arcmin²로 예상 행 수는 46,424개이며 16개 batch다. 계획 SHA-256은 `eeef98ef2210337f714e36d8a4afcd352fe312a536c90e8e236ff9f843f53e52`, SQL aggregate SHA-256은 `5d82df638435ee39029dc51acdcf580ab6e242d21562bc5523ef9c807368f585`이다. 계획 파일은 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_random_pixel_variance_plan_2026-11-06_03/plan.json`이며 준비만 되었고 제출되지 않았다. 제출기 `hsc_debug/submit_icl_random_pixel_variance.py`는 이 계획 해시의 명시 승인을 요구한다. HSC 공식 문서는 random table의 `pix_variance`와 `inputcount_value`를 픽셀 분산 및 해당 위치 기여 노출 수로 설명한다([PDR3 random-point 자료](https://hsc-release.mtk.nao.ac.jp/doc/index.php/random-points__pdr3/)). 결과를 받아도 완전한 coadd 가중치 지도, correlated-noise, sky subtraction 또는 source-Poisson correction을 대체하지 않으며 50–75 kpc 위성 고리의 noise correction에는 충분한 표본밀도가 아니다. `science_measurement=false`를 유지한다.

### 2026-11-05 PSF 주입 표본 대표성 감사

기존 PSF-convolved joint injection 57쌍을 전체 224개 ICL-supported
cluster/filter pair와 비교했다. `audit_joint_injection_psf_sample_representativeness.py`
가 저장된 profile/summary table만 읽어 필터별 redshift, richness, ICL
200–300 kpc coverage와 50–75 kpc satellite support를 비교했다. 57개는
224개의 정확한 부분집합이지만 ICL coverage가 더 높다(SMD +0.800; i
+1.061, i2 +0.736). 위성 50–75 kpc profile이 60% floor를 통과한 비율도
54/57 (94.7%)로 parent의 98/224 (43.8%)보다 높다. 따라서 기존 57개
injection transfer를 224개로 일반화하지 않는다. 재현 스크립트와 요약은
[PSF subset 대표성 감사](hsc_debug/icl_joint_injection_psf_sample_representativeness_assessment_2026-11-05.md),
머신 출력 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_joint_injection_psf_sample_representativeness_2026-11-05_01/`
(summary SHA-256 `e9bd84a6ccb309fc57b405edaca18c0c0c71f91702fda88364987a6e4fb6bbf3`)에
기록했다. `science_measurement=false`다. 224개 지원쌍에 위치별 PSF 코어를
확장하는 것은 주입·회수 검증 범위만 넓히며 관측 검출이나 완전도 보정이 아니다.
이를 준비해 `position_psf.py`에 세 개의 감사된 PSF batch를 결합하는
`SupportedPairsPositionPSFBank`를 추가했고, `run_joint_icl_satellite_batch.py`에
위성 프로필 지원 실패를 보존하는 `--icl-supported-only` 선택을 넣었다.
`verify_supported_pairs_position_psf_bank.py`는 224개 프레임·2,693개 위성 역할과
파일 해시를 수신 뒤 확인한다. 아직 추가 PSF는 받지 않아 이 확장 실행 경로는
실행되지 않았고, 결과는 계속 `science_measurement=false`다.

2026-11-05에 로컬 S23B CAMIRA 표본 풀의 여섯 가지 적색편이/BCG 분광
선택 시나리오를 비교했다. 기존 `0.30–0.55` 및 BCG 분광 일치 기준은
이전 감사의 후보 ID 1,333개를 정확히 재현했다. BCG 분광 조건을 유지한
채 `0.25–0.65` 또는 `0.20–0.70`로 넓히면 catalog-only 추가 후보는
1,983/2,339개지만 새 5도 셀은 25/26개뿐이다. 분광 일치를 요구하지 않는
최광 선택은 4,137개/31개 셀이지만 cluster redshift 교정이 더 필요해
현재 기본 표본으로 삼지 않는다. 5도 셀은 독립 공분산 블록이 아니며,
어떤 후보도 PDR3 영상·마스크·깊이·PSF·전체 aperture 및 matched-control
지원을 입증하지 않는다. 따라서 현재 선택은 재현 기준으로 유지하고,
확장 시료는 사전 지정 redshift 층화와 영상/대조장 coverage 감사를 통과한
뒤 별도 분석한다. Astrafek 재현기
[`audit_s23b_camira_selection_scenarios.py`](hsc_debug/audit_s23b_camira_selection_scenarios.py),
결과 평가 [`icl_s23b_camira_selection_scenarios_2026-11-05.md`](hsc_debug/icl_s23b_camira_selection_scenarios_2026-11-05.md),
machine summary `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_camira_selection_scenarios_2026-11-05_02/summary.json`
(SHA-256 `e0fddda7ecc973327842796278dd01b8e2ae3de33f2d07244590cfb085507ff6`)에 기록했다.
새 HSC 조회나 자료 반출은 하지 않았고 `science_measurement=false`다.

2026-11-05에는 승인된 표적/random shape 행에서 HSC Y3 source z-bin별 exact
S19A estimator-weight boost와 source-bin 1 제외 스택 민감도를 평가했다.
공통 435 렌즈/2,175 random 중심을 유지했고, 네 bin의 pair count 합계와
재구성한 전체 `DeltaSigma`가 기존 저장 스택과 정확히 일치했다. z-bin 1
weight excess는 0.1–0.2 Mpc에서 `4.605±0.274`, 2–3 Mpc에서
`1.260±0.049`였고, bins 2–4는 전 반경에서 약 0.97–1.16이다. z-bin 1과
각 higher bin의 21 radial contrasts는 2.5°와 5° paired jackknife 모두에서
Holm 보정 뒤 유의했으나 normal/jackknife 근사에 조건부인 진단이다.
z-bin 1을 제외한 profile은 양의 tangential 신호와 0과 양립하는 cross
성분을 유지하지만, 0.1–0.4 Mpc tangential amplitude가 각각 약 12.6%,
9.8% 커진다. 이는 낮은 source-bin의 적색이동/선택 문제가 stack에 영향을
준다는 근거다. member contamination, magnification, mask/depth/PSF/sky,
또는 `n(z)` 편향을 구별하지 못하므로 source-density boost를 보정으로
곱하지 않는다. 분석기
[`analyze_wl_tomographic_weight_boost.py`](hsc_debug/analyze_wl_tomographic_weight_boost.py),
평가 [`wl_tomographic_weight_boost_assessment_2026-11-05.md`](hsc_debug/wl_tomographic_weight_boost_assessment_2026-11-05.md),
산출물 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_tomographic_weight_boost_2026-11-05_04/summary.json`
(SHA-256 `47419ac229ffec9aab2a001f629252627126a456812b01cf11908ad44122fe14`)에 기록했다.
무작위 aperture depth/mask/PSF/sky 비교, contiguous footprint 공분산,
photo-z posterior 및 reduced-shear 오차는 남아 `science_measurement=false`다.

ICL–위성 지원 교집합 94쌍을 중심 마스크 적격/부적격으로 추가 분리했다.
적격 표본은 i 8쌍(7개 5° block), i2 15쌍(11개 block)뿐이다. z와
richness를 조정한 부분상관은 i에서 적격 `−0.421` (`95% [-1.000,1.000]`),
부적격 `+0.234` (`[0.014,0.802]`); i2에서 적격 `+0.307`
(`[-0.069,0.792]`), 부적격 `+0.312` (`[-0.070,0.551]`)다. 각 그룹의
5° 연결 block bootstrap 5,000회와 같은 block을 양쪽에 적용한 paired
차이 구간은 i `−0.655` (`95% [-1.665,0.632]`, paired valid 86%), i2
`−0.005` (`[-0.516,0.693]`)였다. 두 paired 구간은 모두 0을 포함하고,
i paired 표본의 invalid draw가 많아 정밀한 마스크 효과 추론은 불가능하다.
작은 적격 표본에서 ICL–위성 상관은 여전히 불확정이다. 마스크를 완화하거나
두 그룹을 동일 모집단으로 합치지 않는다.
분석기
[`analyze_icl_satellite_center_mask_sensitivity.py`](hsc_debug/analyze_icl_satellite_center_mask_sensitivity.py),
평가 [`icl_satellite_center_mask_sensitivity_assessment_2026-11-05.md`](hsc_debug/icl_satellite_center_mask_sensitivity_assessment_2026-11-05.md),
산출물 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_center_mask_sensitivity_2026-11-05_02/summary.json`
(SHA-256 `8b5ad1230b78816a68b8d8ea810f64df32740f1ec199caac86f6be65c7e0e9ee`)에 기록했다.
기존 요약/좌표만 사용했으며 새 HSC 조회나 영상 판독은 없고
`science_measurement=false`다.

2026-11-05 sample-expansion projection은 S23B catalog-only 후보 1,333개를
현재 800개 parent에 더할 때 50–75 kpc `i2` pair가 60% coverage yield의
8.5% 유지 가정에서 약 181개, 40% yield의 17.5% 유지 가정에서 약 373개가
될 것으로 계산했다. 40% 표본의 관측 partial-r은 0.249이며 95% spatial-block
interval `[-0.002,0.454]`로 0을 포함한다. 이 후보 전체는 기존 discovery와
겹치지 않는 5° target 셀이 22개뿐이며 통합 셀 상한도 100개다. 이는 단일
사전 지정 endpoint 80% 계획의 92 independent blocks보다 작고, 30 비교
familywise 계획의 180개보다도 작다. 따라서 S23B 자료는 탐색적 pair 수를
늘릴 수 있지만 확인 표본의 독립 공간력을 채우지 않는다. image/mask/depth/
PSF/member/control 지원은 아직 미검증이다. 재현 분석기
[`assess_s23b_support_yield_and_blocks.py`](hsc_debug/assess_s23b_support_yield_and_blocks.py),
평가 [`icl_s23b_support_yield_and_blocks_assessment_2026-11-05.md`](hsc_debug/icl_s23b_support_yield_and_blocks_assessment_2026-11-05.md),
기계 요약 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_support_yield_and_blocks_2026-11-05_01/summary.json`에 기록했다.

2026-11-06 후속으로 같은 1,333개 후보를 840개 HSC tract-center 연결
성분에 직접 배정했다. 800개 discovery lens는 세 성분에 499, 246, 55개가
있으며 S23B 후보도 같은 성분에 694, 564, 72개가 있다. 네 번째 성분에는
S23B 후보 3개만 있고 discovery lens는 없다. 따라서 22개의 새 5도 셀은
독립 HSC 연결 영역 22개를 뜻하지 않는다. S23B 후보 전체는 탐색적 표본
확장에는 도움이 될 수 있지만 독립된 영역 재현 표본은 거의 늘리지 않는다.
분석기 `hsc_debug/audit_s23b_candidate_tract_region_support.py`, 요약
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_tract_region_support_2026-11-06_01/summary.json`
(SHA-256 `e9a9b9424ff1e8ec9f7014a2ab323645d0bc492c292e6714339167b79733c703`).
분류는 근사 nearest-tract-center 결과이며 candidate 영상은 다운로드하지 않았다.

### 2026-11-05 random-subtracted weak-lensing NFW mass screen

Astrafek runner [`fit_wl_random_stack_nfw_mass.py`](hsc_debug/fit_wl_random_stack_nfw_mass.py)
forward-modelled the saved 435-lens random-subtracted tangential profile using
source redshift bins 2–4, exact lens/source geometry and estimator weights,
and proper 0.2–1.5 Mpc. With fixed concentration `c=3`, the effective
`M200c` is `1.78e14 Msun` (5° endpoint jackknife) and `1.80e14 Msun` (2.5°);
`c=5` gives `1.54e14` and `1.44e14 Msun`, with poor conditional profile
chi-square at both block widths. This dependence is a warning that these
conditional single-mass fits do not calibrate the richness-selected mass
distribution. Reduced shear, concentration distribution, magnification,
miscentering, two-halo structure, photo-z systematics, random-aperture
depth/mask/PSF/sky matching, and covariance uncertainty remain unresolved.
Treat the mass as an exploratory effective scale only. This does not establish
a mass–ICL–satellite relation or change `science_measurement=false`. The
[assessment](hsc_debug/wl_random_stack_nfw_mass_screen_assessment_2026-11-05.md)
and machine summary `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_stack_nfw_mass_2026-11-05_01/summary.json`
(SHA-256 `d3ca78a24fc030f64694247fd45d33f607d89c18fb102b6b5a6879d53c52cd4b`)
record the inputs and fit details. No new HSC query or image read was used.

### 2026-11-05 random-subtracted weak-lensing reduced-shear sensitivity

Astrafek [`fit_wl_random_stack_reduced_shear_sensitivity.py`](hsc_debug/fit_wl_random_stack_reduced_shear_sensitivity.py)
integrated the HSC Y3 tomographic `n(z)` moments through ninth power and fit
the reduced-shear geometric series through order 8 to the existing 435-lens
random-subtracted profile. It retained the same exact pair geometry, source
bins 2–4, radial range 0.2–1.5 proper Mpc, estimator calibration, and saved
5°/2.5° jackknife covariance as the weak-shear NFW screen. The conservative
series-tail bound is 0.012–0.022%. Effective `M200c` shifts down 1.7–5.5%
across fixed `c=3,4,5`; fixed-concentration dependence remains about 20%, and
`c=5` still has poor conditional profile fit. Reduced shear therefore matters
at a few-percent level but does not resolve the larger single-mass/concentration
model uncertainty. This is not a mass calibration or ICL–satellite detection;
`science_measurement=false`. Full equations, fit table and limitations are in
[the assessment](hsc_debug/wl_random_stack_reduced_shear_sensitivity_assessment_2026-11-05.md).
The provenance-complete summary is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_stack_reduced_shear_sensitivity_2026-11-05_03/summary.json`
(SHA-256 `310a03e29069b7799239f21a9ae3c0d0ef0abda5e2d0da93c9e8bdc2b50bc286`).
No HSC request or image read was made.

### 2026-11-05 random-subtracted richness dependence

The cluster-only richness contrast was re-evaluated after matched-random
subtraction using Astrafek
[`analyze_wl_random_subtracted_richness_split.py`](hsc_debug/analyze_wl_random_subtracted_richness_split.py).
From 435 saved common-support lenses, the exact estimator cuts leave 415 with
usable target and random support (211 lower, 204 upper richness). A predeclared
parent split followed by maximum-cardinality redshift matching yields 204
lower/upper pairs. The high-minus-low random-subtracted `DeltaSigma` contrast
is positive from 0.2 to 3 proper Mpc (for example, `+60.16 +/- 18.09` in
`10^12 M_sun Mpc^-2` at 0.2–0.4 Mpc); the inner 0.1–0.2 Mpc bin is consistent
with zero. The conditional seven-bin Wald diagnostic (`chi2=67.84`, rank 7,
nominal `p=4.04e-12`) is not a detection significance because the 50 equal-area
jackknife cells are not validated contiguous HSC regions and systematics are
not propagated. The signal supports an exploratory richness-dependent lensing
profile in this selected sample, not a calibrated richness–mass relation or
ICL–satellite detection. Matched-aperture depth/mask/PSF, source-support
selection, `n(z)`, shear calibration, mass/concentration distribution, and
covariance remain open; `science_measurement=false`. Full results and limits:
[assessment](hsc_debug/wl_random_subtracted_richness_split_assessment_2026-11-05.md).
The summary is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_subtracted_richness_split_2026-11-05_02/summary.json`
(SHA-256 `b38a78f6005dc17d8a12ddd6117debd63a005d0b21ff8dbda05961bd31eb65b7`).
The initial 214-pair contrast is retained only as a superseded diagnostic; it
failed to require estimator-usable target support.

### 2026-11-05 joint two-mass fit to the richness contrast

Astrafek [`fit_wl_richness_split_nfw_mass_contrast.py`](hsc_debug/fit_wl_richness_split_nfw_mass_contrast.py)
jointly fits lower- and higher-richness effective `M200c` to the 204-pair
random-subtracted `DeltaSigma` contrast, using the observed shared 50-block
covariance and four proper-radius bins from 0.2–1.5 Mpc. Across fixed
concentrations `c=3,4,5`, the higher/lower mass ratio is 1.21–1.89 with source
bins 1–4 and 1.39–2.85 with source bin 1 excluded. The conditional profile
intervals are broad (ratio upper limit as high as 28), and the absolute masses
shift strongly with concentration and source-bin selection. The direction
suggests a higher effective lensing mass in the high-richness group in this
selected sample, but does not calibrate the CAMIRA richness–mass relation.
The approximate sky-block covariance and unresolved depth/mask/PSF, photo-z,
shear, reduced-shear and mass-distribution systematics preclude a detection
claim; `science_measurement=false`. Full limits and fit statistics are in
[the NFW contrast assessment](hsc_debug/wl_richness_split_nfw_mass_contrast_assessment_2026-11-05.md).
All-bin result
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_richness_split_nfw_mass_contrast_2026-11-05_04/summary.json`
(SHA-256 `c0e4c92e5ef00fc77ffef7e26e918d3da312049e6d5ee2877edbe4385c32217b`);
source-bin-1-excluded result
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_richness_split_nfw_mass_contrast_2026-11-05_03/summary.json`
(SHA-256 `6d293f41c12a2d44fc4dd7303be238c9bfbcca27570db180f0bf93e3224eaa94`).
No HSC query or image read was made.

### 2026-11-05 ICL–위성 지원 선택·주입 표본 통합 재평가

Astrafek의 `build_icl_association_evidence_report.py`가 기존 94쌍 상관 및
공간·PSF·노이즈 진단에 선택확률 가중 민감도, PSF 주입표본 대표성,
HSC frame-count offset 검사를 해시와 함께 통합했다. 선택확률 가중 후 부분상관은
i에서 `−0.011` (원시 `0.079`), i2에서 `0.227` (원시 `0.293`)이다. 조건부
독립 결측 가정은 검증되지 않았고 블록 불확실도는 0을 포함한다. 57개 PSF 주입
표본의 ICL 200–300 kpc coverage는 224개 지원표본보다 전체 `0.800 SD`
( i2 `0.736 SD`) 높으며, satellite 지원률도 더 높아 주입 전달 결과를 전체에
적용할 수 없다. Frame/visit/CCD 수 차이는 ICL 대비와 뚜렷한 상관을 보이지
않았지만 픽셀 깊이·coadd 가중치·분산을 검증하지 않는다. 따라서 선택편향과
주입 외삽 문제는 여전히 미해결이며 `science_measurement=false`다. 해시를 포함한
[통합 평가 v2](hsc_debug/icl_association_evidence_status_2026-11-05_03.md)를 참조한다.

### 2026-11-05 richness 대비 약전단–reduced-shear 모형 민감도

동일한 204개 적색이동 matched lens pair, source-row estimator weight,
random-subtracted `DeltaSigma`와 저장 50-block covariance를 유지하고, source
HSC Y3 `n(z)`를 적분한 reduced-shear NFW 두 질량 contrast를 추가했다. ICL–위성
연관의 직접 측정은 아니며 weak-lensing branch의 질량 해석 민감도다. 급수는
source-bin `n(z)`의 비영 지지에서 order 16까지 계산했다. 고/저 richness
유효 질량비는 source bin 1–4에서 고정 `c=3,4,5`별 `1.168, 1.200, 1.214`,
source bin 1 제외에서 `1.150, 1.172, 1.198`이다. 기존 약전단 적합의 대응값은
각각 `1.894, 1.396, 1.209`와 `2.847, 1.727, 1.391`이었다. reduced shear는
특히 낮은 concentration 적합의 ratio를 낮추며 concentration 간 차이를 줄인다.
order 8→16에서 mass ratio는 최대 0.3% 변했다. 가장 극단적인 c=3/source-bin
1 제외의 보수적 per-pair 생략항 상한은 6.94%로 잔여 모형 한계로 기록한다.
구간은 같은 저장 covariance의 conditional Δχ²=1 profile interval이며 검출
유의도가 아니다. 한 halo per richness group, selected sample, approximate sky
covariance와 photo-z/calibration/magnification/miscentering/depth-mask-PSF-sky
systematics가 남아 `science_measurement=false`다. 실행기와 SHA, 입력·결과는
[reduced-shear 평가](hsc_debug/wl_richness_split_reduced_shear_assessment_2026-11-05.md)에 있다.

### 2026-11-06 확산광 주입 진폭 선형성과 HSC frame metadata

PSF를 위치별로 컨볼루션하고 source mask를 매번 재생성하는 137개 동일 프레임에서 확산광 주입 진폭 0.5, 1.0, 2.0 nJy arcsec⁻² (100 kpc 기준)을 비교했다. 위성 주입은 1.0으로 고정했다. 진폭당 ICL 고리 응답 중앙값은 `i=0.240, 0.226, 0.224`, `i2=0.225, 0.226, 0.226`이며, 1.0 대비 paired 5° 공간블록 bootstrap 구간은 모든 진폭·필터에서 0을 포함한다. 검증 범위에서 상대 응답은 대체로 선형이지만 절대 회수율은 약 23%이고 마스크도 여러 프레임에서 변한다. 주입원 Poisson 잡음, coadd 픽셀 가중치, HSC 상류 sky 전달, 확장 PSF 날개가 빠져 있으므로 이 응답을 관측 ICL 스택의 flux 보정으로 쓰지 않는다. 결과와 재현기는 [진폭 응답 평가](hsc_debug/icl_diffuse_injection_amplitude_assessment_2026-11-06.md), `analyze_icl_diffuse_injection_amplitude_grid.py`, `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_diffuse_injection_amplitude_grid_supported137_2026-11-06_01/`에 있다. `science_measurement=false` 유지.

사용자가 plan SHA-256 `49aa4b7601dd42243209d7b3b394a4dd3b4cbf00d87d28cbb9467553c29f4ce3`을 승인한 뒤, 1,598개 target/control 중심의 frame metadata 67,071행(21열)을 16개 CAS 배치로 수신하고 독립 오프라인 감사했다. exposure·4 gain·FLUXMAG0·sky RMS·MJD는 모두 결측 없이 반환됐고, frame 원 필터는 `i` 20,341행 및 `i2` 46,730행이다. 798 matched target/control 쌍에서 frame/visit/CCD 수, exposure 합, sky RMS 중앙값, FLUXMAG0 중앙값, amplifier별 gain 중앙 target−control 차는 모두 0이다. 이는 patch-level 표본의 거친 불균형을 배제할 뿐 pixel depth나 coadd variance 측정이 아니다. 원본 frame map을 coadd 픽셀에 연결하는 가중치/증폭기 영역이 없으므로 Poisson 분산을 구성하지 않고, `science_measurement=false`를 유지한다. 독립 auditor `hsc_debug/audit_icl_frame_metadata_row_results.py`와 summary `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_metadata_row_audit_2026-11-06_01/summary.json` (SHA-256 `0f361fcb496fbd7be7c1059f638414f6900535f73a9ac48696030480f400a237`)에 기록했다.

### 2026-11-06 ICL/WL split의 tract 연결성분 분포

parent 800 cluster에서 이미 만든 근사 nearest-tract 성분표로 저장된 split 표본 ID를 대조했다. i 공통표본 38개는 성분 0/1에 18/20개이고, 매칭쌍은 고 ICL 그룹 3/7, 저 ICL 4/6으로 대체로 균형이며, 10쌍 중 5쌍이 성분을 가로지른다. i2 공통표본 70개는 0/1/2에 48/18/4개이고, 29쌍의 고/저 그룹 수는 21/7/1 및 22/4/3, 같은 성분 17쌍·교차 성분 12쌍이다. 따라서 i 강한 cross-null은 단순한 high/low lens 수 불균형만으로 설명되지 않지만, source·random 위치를 성분별로 분해하지 않아 원인 판정은 아니다. 성분은 공식 HSC field가 아니라 근사 footprint 연결이다. Astrafek 분석기 `hsc_debug/audit_wl_icl_split_tract_region_balance.py`; 결과 hash는 i `8cba95a740c5d332b8c313c5687e58289652042497a5cc977061032203657130`, i2 `adb313d390c9def350255dd8719d46ea9c536b1864f142b73283455421c0f6d1`. 추가 HSC 자료 미수신, `science_measurement=false`.

### 2026-11-06 i 대역 ICL split weak-lensing 블록 폭 민감도

처음 i2에서만 실행한 ICL split 약중력 분석을 저장된 i 대역 공통표본으로 확장했다. 유효 공통렌즈 37개 중 200–300 kpc ICL 대비 중앙값 split 후 적색이동·log-richness caliper matching (`|Δz|≤0.03`, `|Δln λ|≤0.25`)은 10쌍만 남겼다. source tomography 2–4, 같은 5개 random aperture/lens와 저장된 HSC source-row를 유지하고 공간블록 폭 2.5°, 5°, 10°를 바꿨다. 매칭된 쌍 목록과 7개 radial bin의 점 프로파일은 세 계산에서 정확히 같다. 그러나 tangential Wald 진단은 각각 χ²=15.79/161.92/101.28로 크게 달라지고 cross-component null도 χ²=24.70/90.97/105.50으로 모두 크다. 이 표본은 10개 lens pair에 블록이 9–24개뿐이고 공분산이 블록 정의에 민감해 p값은 검출 유의도가 아니다. radially mixed 점 프로파일과 강한 cross-null은 i 결과를 신뢰 가능한 ICL–halo-mass 연관으로 해석하지 못하게 한다. 기존 i2 split도 블록 폭에 따라 tangential 명목 p가 `2.42×10⁻⁴`에서 `0.099`로 달라져 두 필터가 견고한 공통 증거를 주지 않는다. 상세 표와 입력·결과 SHA는 [i-band 평가](hsc_debug/wl_icl_i_band_blockwidth_assessment_2026-11-06.md)에 저장했다. image pixel/HSC 추가 조회 없음, `science_measurement=false` 유지.

### 2026-11-06 i split cross-null의 source/랜덤 영역 분해

저장된 승인 shape-row와 random-row를 재검증한 뒤 i-band 10 matched pair의 target와 다섯 random center/lens를 parent tract component별로 분리해 HSC Y3 stack을 다시 계산했다. 기존 radial, source-bin, shape, magnitude, resolution, B-mode, `n(z)` 효율 cuts를 재적용해 target 사용 행 53,399개 및 random 사용 행 241,258개를 재현했다. 사용된 모든 source는 대응 lens/center와 같은 tract component에 매핑되어 영역 경계 넘어선 source 연결은 cross-null 설명이 아니다. high-minus-low random-subtracted cross 프로파일은 두 성분 모두 진동하며 첫 두 radial bin은 같은 부호와 큰 진폭이다. 한 component만 단독으로 결과를 만드는 패턴은 아니지만, 성분당 3–7개 lens/group뿐이고 region covariance가 없으므로 물리·기기 원인을 판별할 수 없다. Astrafek tool `hsc_debug/audit_wl_icl_split_source_region_cross_null.py`, summary SHA-256 `fb90c982da784dff60b270efdfe986d0fbe5b8f3a656e15dd5bc14629810862c`. 추가 HSC query/download 없음; `science_measurement=false`.

### 2026-11-06 target/random cross 분해와 외곽 tangential 후보

성분별 endpoint arithmetic에서 첫 두 radial bin의 cross high−low 값은 주로 target 쪽에 있으며 random-center 차분이 이를 일부 줄이거나 component 1에서는 더 크게 만든다. 저장 covariance를 `R≥1 Mpc`로 잘라 post-hoc 검토했을 때 i cross null nominal p는 2.5°/5°/10° 블록에서 0.307/0.303/0.107로 감소하고, tangential high−low contrast는 각각 p=0.0093/6.6×10⁻⁵/5.3×10⁻⁵였다. 그러나 더 큰 i2 29-pair 표본에서는 같은 cut의 tangential p가 0.172/0.114 (2.5°/5°)이고 cross p는 0.639/0.550으로, outer i tangential 후보가 재현되지 않는다. i cut은 프로파일 확인 뒤 선택했고 10 matched lens pairs에 9–24개 근사 block만 있어 어느 결과도 검출이 아니다. 재현 도구 `hsc_debug/audit_wl_icl_split_target_random_cross_decomposition.py`와 `hsc_debug/audit_wl_icl_split_radial_cut_nulls.py`; 결과 hash `ffc453ac4b78e7f25a354d19469b31411856cc6d4e0ab41274f6940f331fe5af` 및 `583426fef8b60c011ea392dd80c3a7783d1572fa06006d72b729bebea1eade45`. `science_measurement=false` 유지.

### 2026-11-06 연속 ICL 대비와 i2 약중력 스택

저장된 승인 HSC source/random 행만 사용해 ICL/WL 공통 지지 i2 렌즈 70개를 개별 분석했다. 각 렌즈의 target `DeltaSigma`에서 대응 랜덤 중심 5개의 합산 프로파일을 빼고, 연속형 target−control ICL 200–300 kpc 대비와 redshift 및 log-richness를 각 radial-bin 회귀에 넣었다. source tomography 2–4를 사용했고, radial-bin별 유효 렌즈 수는 58, 63, 65, 66, 67, 68, 70이다. 가장 바깥 2–3 proper Mpc bin에서 ICL 1 SD 증가당 기울기는 `−12.3×10^12 M☉ Mpc⁻²` (HC3 명목 p=0.021, approximate 27-cell delete-one 오차 `4.8×10^12`)이나, 일곱 tangential-bin 보정 최소 p/BH q는 0.149다. 나머지 radial 기울기는 부호가 바뀌고 모든 cross 기울기는 명목 p≥0.20이다. 바깥 음의 후보는 기존 i2 median split과 재현되지 않으므로 ICL–halo-mass 검출로 보지 않는다. 표본은 parent 800 중 70개이고, 블록은 검증된 HSC field가 아니며 ICL Poisson/coadd-weight, photo-z, shear, reduced-shear, magnification, centering, depth/mask/PSF 계통오차가 미해결이다. 실행기와 정량 평가는 [연속 ICL–WL 평가](hsc_debug/wl_continuous_icl_relation_assessment_2026-11-06.md), 결과 hash `12c9456799e2671e3a0452c09a52465bdc6434c5b7546659c4bf21ac340f5fea`에 기록했다. 새 HSC 조회는 하지 않았고 `science_measurement=false` 유지.

같은 결과의 방사거리별 tangential/cross 회귀 기울기와 HC3·공간 delete-one 오차를 두 패널로 그렸다. 재현기 `hsc_debug/plot_wl_continuous_icl_relation.py`; PNG `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_continuous_icl_relation_i2_2026-11-06_02/wl_continuous_icl_relation_i2_v2.png` (SHA-256 `f8225b95111ab29a040edf97826eb3a6469fa0f87ed66c81b27078273c82ec08`), PDF SHA-256 `c3e8b5a58838ec02282736fa1f50dd51f9bbf50bffe53743cd20140d21a53cc3`. 그림은 외곽-bin 명목 p가 7-bin 보정 뒤 유의하지 않음을 표시한다.

### 2026-11-06 HSC warp 분산층 pilot 계획 준비

이미 승인·수신한 CCD-frame metadata에서 표적 `target:6966` (`RA=346.791503°`, `Dec=+1.039977°`, tract/patch `9703/306`)의 43 frame행, 9 visit, 40 CCD를 검증했다. DAS Cutout manual의 semi-width/height 규약으로 2′×2′ tangent-plane `pdr3_wide/HSC-I/type=warp` 범위를 만들고 image, mask, variance 세 층을 지정했다. 0.168″/pixel, 세 32-bit plane, 9 full-square visit을 가정한 비압축 평면 용량은 55,212,300 bytes (52.65 MiB)이고 downloader는 70 MiB에서 중지한다. Astrafek planner `hsc_debug/prepare_icl_warp_pilot_plan.py`가 offline plan을 생성하고 승인-gated `hsc_debug/download_icl_warp_pilot.py`는 기본 plan-only이며 승인 SHA가 일치할 때만 요청한다. plan SHA-256 `a436df245a33bda09cd5599e564b2b9302e5d0d4c006b59fcf6fab1c2c33a65e`; 저장 경로 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_warp_pilot_plan_2026-11-06_02/plan.json`. 현재 미제출이며 warp image/mask/variance 층을 새로 받지 않았다. Pilot에서 archive·visit·FITS HDU·variance 정렬을 확인하고 나서만 전체 데이터 요청 범위를 별도 계획한다. 그 전까지 `science_measurement=false`.

### 2026-11-06 i2 ICL–위성–weak-lensing 공통 표본

기존 ICL–위성 지원표본 68개와 연속 weak-lensing 지원표본 70개 중 감사에서 확인된 i2 렌즈 27개를 결합했다. 표적/control 5도 연결 블록을 다시 구성한 5,000회 bootstrap에서 ICL–위성 부분상관은 redshift 및 `ln(richness)` 통제 후 `r=0.067`, 95% 구간 `[-0.144, 0.293]`이다. ICL과 WL의 조건부 연관 14개는 모두 구간이 0을 포함한다. 위성 프로파일과 WL의 14개 조합 중 두 점별 구간만 0을 벗어나지만, 29개 관계를 다중비교 보정하지 않은 선택표본 탐색 결과여서 검출로 해석하지 않는다. 자세한 구간·해시·한계는 [삼중 공통 표본 평가](hsc_debug/icl_satellite_wl_triple_overlap_assessment_2026-11-06.md)에 기록했고, Astrafek 재현기는 `hsc_debug/analyze_icl_satellite_wl_triple_overlap.py`다. 결과 경로 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_wl_triple_overlap_2026-11-06_01/`. 추가 HSC 조회는 하지 않았으며 `science_measurement=false` 유지.

### 2026-11-06 지원 ICL–위성 쌍의 반경별 coadd variance

기존 승인 이미지 중 94개 ICL·위성 지원 target/control 쌍 (i 26, i2 68)의 VARIANCE plane을 0–900 kpc 열 개 고리에서 비교했다. FITS 이미지 해시·HDU·중심 WCS를 검증하고 bad/bright-star mask bits, 비유한·비양수 variance 픽셀을 제외했다. 필터·고리별 5도 연결블록 5,000회 bootstrap에서 20개 점별 비교 중 세 구간만 target/control 차이가 0을 벗어났다 (i 0–50 kpc 분산비 1.058; i2 50–75와 75–100 kpc 0.974/0.952). 다중비교 보정 전의 약한 탐색 신호이고, ICL 대비와 variance 비의 부분상관은 모든 구간에서 bootstrap 구간이 0을 포함한다. 외곽 600–900 kpc의 median 유효 픽셀 비율은 약 0.53이다. 이 coadd variance는 합성원 source-Poisson의 유효 gain을 주지 않으며 관측 프로파일이나 과학 오류를 보정하지 않는다. Astrafek 도구 `hsc_debug/analyze_icl_variance_radial_balance.py`, 평가 `hsc_debug/icl_variance_radial_balance_supported94_assessment_2026-11-06.md`, 산출물 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_variance_radial_balance_supported94_2026-11-06_01/`. 추가 HSC 조회·다운로드 없음, `science_measurement=false` 유지.

### 2026-11-06 전체 783 target/control 표본의 반경별 분산 및 가족단위 검정

같은 필터의 이미지가 모두 있는 전체 783쌍 (i 302, i2 481)에서 delivered coadd VARIANCE plane을 같은 열 개 proper-kpc 고리로 측정했다. 20개 필터·반경 대비를 함께 보정하는 10,000회 max-absolute-t 동시구간은 두 필터 사이 공간블록 연결을 유지했고, 38개 근사 block 중 가장 큰 block은 107쌍이다. 가족단위 95% 동시구간에서 `i2`의 75–100, 100–150, 200–300 kpc 분산비만 같음(1)보다 낮았으며 점 추정은 0.961, 0.955, 0.961; 동시구간은 각각 [0.924,0.999], [0.915,0.996], [0.934,0.988]이다. i 대역에는 차이가 남지 않았다. 이는 선택된 영상의 delivered variance plane 차이일 뿐 전체 잡음·source-Poisson이나 ICL 추정량 공분산을 보정하지 않으며, 38개 근사 block/큰 단일 block의 공간 한계가 있다. 94개 ICL–위성 지원표본으로 제한한 ICL 대비 상관 구간은 모두 0을 포함한다. 상세 결과와 해시는 [전체표본 반경 분산 평가](hsc_debug/icl_variance_radial_balance_full783_assessment_2026-11-06.md). 재현 도구 `hsc_debug/analyze_icl_variance_radial_balance.py --scope all783` 및 `hsc_debug/analyze_icl_variance_radial_familywise.py`. 첫 familywise `_01` 표는 필터별 지표 분리가 누락되어 invalid 처리했고, 고친 `_02` 표를 유효 결과로 사용한다. 추가 HSC 조회·다운로드 없음, `science_measurement=false`.

### 2026-11-06 warp pilot FITS 사후검증기 보강

승인-gated `hsc_debug/download_icl_warp_pilot.py`의 기존 방문·HDU 검사에 이미지/마스크/분산층 shape, 세 층의 WCS 일치, 요청 중심이 영상 내부에 있는지, 0.168″/pixel ±5%, 유한 이미지 비율, 마스크 점유율, 분산층의 유한·양수 비율을 기록하는 검사를 추가했다. 방문 집합이 맞아도 공간·분산 검사에 실패하면 `complete_with_qa_warnings`로 표시한다. 현재 default plan-only 실행에서 plan SHA `a436df245a33bda09cd5599e564b2b9302e5d0d4c006b59fcf6fab1c2c33a65e`와 9개 예상 방문을 확인했으며 네트워크 요청은 제출하지 않았다. 다운로더 SHA-256 `766dbedb2b76951307b6c10a9d4be309f96ba1f509d6e09073c0716a6c71e291`. 실제 warp 자료가 없으므로 사후 FITS 검사는 아직 실행되지 않았고 `science_measurement=false`다.

### 2026-11-06 프레임 메타데이터와 coadd 분산의 짝별 연관

승인된 프레임 메타데이터 집계와 기존 coadd VARIANCE 측정을 결합해 200–300 proper-kpc 고리의 같은 필터 target/control 774쌍(i 299, i2 475)을 오프라인 분석했다. outcome은 `log(target/control median delivered VARIANCE)`이며, 프레임 수·방문 수·FLUXMAG0로 정규화한 sky RMS proxy·FLUXMAG0·4개 median gain의 기하평균 차이를 각각 redshift와 log-richness를 통제해 비교했다. 두 필터의 10개 연관을 한꺼번에 재표집하는 공통 5도 연결블록 5,000회 single-step max-T에서 38개 근사 블록 중 최대 블록은 106쌍이었다. `i` 프레임 수의 점별 r은 −0.147이나 동시 95% 구간은 [−0.312, 0.017], 조정 p=0.088이다. `i2`에서 가장 큰 점 추정인 FLUXMAG0 r=0.182도 동시 구간 [−0.247, 0.611], 조정 p=0.750이다. 열 개 동시구간 모두 0을 포함하므로 patch-level frame 요약이 관측된 분산비를 설명한다는 가족단위 근거는 없다. 프레임 행은 pixel-to-frame/amp 매핑을 제공하지 않고 sky RMS는 coadd 잡음이 아니며, 이 분석은 source-Poisson이나 상관잡음을 보정하지 않는다. `science_measurement=false`를 유지한다. 상세 제한·해시·재현법은 [frame metadata–variance 평가](hsc_debug/icl_frame_metadata_variance_association_assessment_2026-11-06.md), 재현기는 `hsc_debug/analyze_icl_frame_metadata_variance_association.py`, 산출물은 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_metadata_variance_association_2026-11-06_02/`이다. 추가 HSC 조회·다운로드는 하지 않았다.
# 2026-11-06: delivered-variance sensitivity of the ICL–satellite relation

Using only the existing 94-pair supported sample and full-sample variance
table, an offline paired block-bootstrap sensitivity analysis compared the
ICL–satellite partial correlation under controls for redshift/richness,
variance ratio, and variance ratio plus a coarse pair-depth proxy. For `i2`
(N=68), the baseline partial correlation is 0.293 (95% connected-block
interval 0.043–0.498); adjusting for target/control variance ratio gives
0.285, with a paired change of −0.008 (−0.078–0.034). Adding the depth proxy
gives 0.260, change −0.032 (−0.148–0.027). Both changes include zero. The
`i` sample (N=26) is imprecise, with all model intervals spanning zero. This
does not show that the delivered `i2` variance difference explains the
supported-sample relation. It is a selected-sample sensitivity check, not a
causal or calibrated science result; `science_measurement=false` remains.
Detailed method, limitations and hashes are in
`hsc_debug/icl_variance_adjusted_satellite_assoc_assessment_2026-11-06.md`.
Runner: `hsc_debug/analyze_icl_variance_adjusted_satellite_assoc.py`.
Machine output:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_variance_adjusted_satellite_assoc_2026-11-06_01/`.

## 2026-11-06: variance as a predictor of satellite-profile support

Because only 94 of 224 ICL-supported cluster/filter pairs have saved
50–75 kpc satellite profiles, the IPW selection model was extended with the
existing 200–300 kpc target/control delivered variance ratio and mean log
variance. In the `i2` sample, the IPW partial correlation is 0.227 at baseline,
0.266 with variance ratio, and 0.193 with both variance terms. Paired
5-degree block-bootstrap changes are +0.039 (95% interval −0.069–0.162) and
−0.034 (−0.163–0.092), respectively; both include zero. The full model also
reduces the `i` effective sample size to 13.5 of 26 supported pairs. These
results show model sensitivity but do not identify a correction: the
missing-at-random assumption is unverified, and the variance plane is not a
total-noise measure. Keep `science_measurement=false`. Full methods and hashes:
`hsc_debug/icl_satellite_ipw_variance_sensitivity_assessment_2026-11-06.md`.
Runner: `hsc_debug/analyze_icl_satellite_ipw_variance_sensitivity.py`.
Output:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_ipw_variance_sensitivity_2026-11-06_01/`.

## 2026-11-06: combined support-selection and pixel-noise sensitivity

For the 94 of 224 ICL-supported pairs with a satellite profile, Astrafek now
combines the satellite-support IPW model with the saved finite-lag
ICL/satellite pixel covariance. The baseline `i2` noise-corrected IPW partial
correlation is 0.238 (95% connected-block interval −0.023–0.446). Adding the
target/control delivered-variance ratio gives 0.278 (−0.018–0.495), and adding
the coarse mean-depth proxy gives 0.201 (−0.054–0.410). Paired changes from
baseline include zero. This remains a sensitivity analysis: support missing
at random is unverified, and the finite-lag correction omits source Poisson,
PSF wings and residual-sky covariance. `science_measurement=false` remains.
Methods and hashes are in
`hsc_debug/icl_satellite_ipw_pixelnoise_joint_assessment_2026-11-06.md`.
Runner: `hsc_debug/analyze_icl_satellite_ipw_pixelnoise_joint.py`.
Output:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_ipw_pixelnoise_joint_2026-11-06_01/`.

## 2026-11-06: forced profiles at undetected CAMIRA member positions

To test whether the satellite stack could extend beyond SEP-detected members,
3,197 CAMIRA positions at 50–500 proper kpc with no unique SEP match within
1 arcsec and no SEP detection within 3 arcsec were measured in the existing
800 HSC-I cutouts. A profile bin required target coverage of at least 60% and
at least three rotated controls meeting the same threshold. Only 14/3,197
(0.44%) positions passed in any annulus and 9/3,197 (0.28%) passed at
50–75 kpc. Target coverage alone reached 60% at that radius for 8/1,292 `i`
and 10/1,905 `i2` positions; only 4 and 5 passed after control support.
Most positions had eight valid rotated controls, so target-annulus coverage
was the main bottleneck. This forced-position step therefore does not
materially recover the detection-selected sample.

Forced-only ICL association support is insufficient (2 `i` and 4 `i2` pairs).
Combining retained forced profiles with existing detected profiles leaves 26
and 68 pairs, respectively. The exploratory redshift/log-richness adjusted,
connected-5-degree-block partial correlations are `r=0.131` (95% interval
−0.211 to 0.556) for `i` and `r=0.287` (0.037 to 0.490) for `i2`. The `i2`
interval is not a science detection or completeness correction; CAMIRA
positions are not verified individual galaxy centroids, the retained forced
sample is strongly coverage-selected, and PSF/sky/source-Poisson/total
covariance validation remains incomplete. `science_measurement=false` remains.
See [the assessment](hsc_debug/icl_satellite_forced_position_assessment_2026-11-06.md),
runner `hsc_debug/measure_forced_satellite_profiles.py`, and support audit
`hsc_debug/audit_satellite_forced_position_support.py`. Machine output:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_forced_outer_profile_2026-11-06_01/`.
No additional HSC query or download was made.

## 2026-11-06: random-subtracted weak-lensing mass sensitivity to source n(z)

기존 승인 shape 및 random-center 행을 이용해 동일한 435개 공통 렌즈와
2,175개 대응 대조 중심에서 source-bin 1을 제외한 cluster/random `DeltaSigma`
프로필을 다시 계산하고 reduced-shear NFW 적합을 반복했다. 무이동 SACC 기준은
저장된 5° block-tensor 프로필과 최대 상대차 `8.0e-16`으로 일치했다. source-bin
3/4의 공개된 평균 고적색편이 이동을 적용하면 고정 `c=3` 유효 질량이
`1.750e14`에서 `1.706e14 M☉`로 2.60% 낮아졌다. 해당 marginal 68% 경계의
네 조합은 중심 이동값 대비 −0.96%~+1.07%였다. 반면 농도 가정 `c=3,4,5`의
기준 질량은 `1.750,1.604,1.450e14 M☉`로 변했고, `c=4,5` 적합은 χ²/dof
`12.58/3`, `32.54/3`로 좋지 않았다. 이 간이 n(z) 보정보다 NFW 농도/모형 가정이
더 큰 영향을 준다. 사용한 5° jackknife 공분산은 모든 n(z) 시나리오에서
고정했고 spline posterior, miscentering, magnification, two-halo와 연속
HSC field 공분산은 미전파다. 따라서 질량은 여전히 탐색적 유효 stack scale이며
`science_measurement=false`다. [상세 평가](hsc_debug/wl_random_subtracted_nz_mass_sensitivity_assessment_2026-11-06.md),
재현기 `hsc_debug/analyze_wl_random_subtracted_nz_mass_sensitivity.py`, 결과
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_subtracted_nz_mass_sensitivity_2026-11-06_01/summary.json`.
신규 HSC 조회·다운로드 없음.

## 2026-11-06: supported-sample long-lag ICL covariance

To check whether the earlier ±4-pixel residual covariance kernel had
converged, it was remeasured to ±8 pixels for the same 224 target/control
pairs (71 `i`, 153 `i2`) that support the 200–300 kpc ICL profile. Filter-wise
5-degree spatial-block bootstrap (5,000 paired resamples) gives square-kernel
target-minus-control median-factor differences of `i=0.025` (95% −0.002–0.047)
and `i2=0.059` (0.026–0.090) at 4 pixels; at 8 pixels they are `i=0.081`
(−0.001–0.132) and `i2=0.077` (0.030–0.129). Negative covariance shells make
the cumulative sums non-monotonic, so convergence at ±4 is not established.
The `i2` 8-pixel filter-wise interval excludes zero, while `i` includes zero;
these are exploratory, unadjusted intervals. This kernel does not model
spatial variation, lags beyond 8 pixels, source Poisson, pixel coadd weights,
sky transfer, or PSF wings. It does not correct the ICL noise or qualify a
science measurement (`science_measurement=false`). Methods and hashes:
`hsc_debug/icl_supported_long_lag_noise_assessment_2026-11-06.md`; runner
`hsc_debug/audit_icl_supported_long_lag_noise.py`; summary
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_supported_long_lag_noise_2026-11-06_01/summary.json`.

The ±8 kernels were also contracted with the actual annulus masks and delivered
variance planes for these same 224 pairs; reconstructed profile mean and
coverage matched the saved values for every pair. The independent-frame
finite-lag pixel-noise stack sigma is 0.0793 nJy arcsec⁻² in `i` and 0.0745 in
`i2`, only 0.8% and 1.3% above the exact-aperture 9×9 estimate. These are 7.6%
and 7.0% of the respective saved spatial-block 95% interval half-widths. This
suggests finite-lag pixel noise is subdominant to measured sample scatter under
the stationary-kernel/independent-frame approximation; it does not establish
the total uncertainty because source Poisson, sky transfer, PSF, covariance
beyond 8 pixels, and shared-field covariance are absent. The measured profile
was not corrected, and `science_measurement=false` remains. Detailed methods
and hashes: `hsc_debug/icl_long_lag_aperture_covariance_assessment_2026-11-06.md`;
runner `hsc_debug/analyze_icl_long_lag_aperture_covariance.py`; summary
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_aperture_covariance_200_300_2026-11-06_01/summary.json`.

### 2026-11-06 target-only 주입의 ICL–위성 교차누출 진단

기존 94개 ICL·위성 지원 쌍에 감사된 PSF-convolved 주입 응답을 연결했다.
유효 주입 응답은 58쌍(i 14, i2 44)에만 있으며, 36쌍은 PSF 기록 부재 31개와
위성 고리 지원 부족 5개 때문에 연결되지 않았다. 대상 은하단 영상에 넣은
위성 주입은 ICL 200–300 kpc 고리에 i 0.357, i2 0.490의 중앙 응답을 냈다.
확산광 주입의 같은 고리 응답은 각각 0.211, 0.224였다. 위성 프로파일 응답은
산포가 크고 0을 통과한다. 은하단 신호는 target에만 존재하므로 target-only
주입은 matched random control을 고정한 target-minus-control ICL 추정량의 증분
응답을 측정한다. 대조장에 같은 은하단 신호를 넣는 것은 이 응답을 측정하는 데
필요하지 않다.

중앙 응답비를 빼는 단순 누출 stress test에서 부분상관은 i −0.743,
i2 −0.567로 바뀌었다. 이 계산은 보정량이 아니며 물리적 반상관의 증거가
아니다. 유효 응답은 선택된 58쌍뿐이며, 개별 response matrix 산포와 지원 선택을
공동 관계에 전파하지 않았다. 단일 중앙 누출계수는 검증된 보정값이 아니다.
후속 분석은 target-only 주입과 실제 고정 대조장을 이용해 정확한
target-minus-control 추정량을 재계산하고, pair별 matrix 및 공간 블록 불확도를
상관 추정에 전파해야 한다. 주입원 Poisson, coadd 가중치, 확장 PSF 날개,
HSC 상류 sky 전달도 별도로 검증해야 한다. 대조장 null 주입은 가짜 신호 응답을
측정할 수 있지만 대조장에 은하단 광원을 넣는 방식으로 수행하지 않는다.
상세 결과는 [주입 교차누출 평가](hsc_debug/icl_satellite_injection_cross_talk_assessment_2026-11-06.md),
재현기는 `hsc_debug/analyze_icl_satellite_injection_cross_talk.py`, 해석을 수정해
재실행한 수치 산출물은
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_injection_cross_talk_2026-11-06_02/`에 있다.
추가 HSC 조회·다운로드는 없으며 `science_measurement=false`를 유지한다.

58개 유효 쌍에서 개별 2×2 응답행렬을 역산하는 추가 민감도 분석을 수행했다.
i2 부분상관은 관측값 0.399에서 역산값 −0.381로 바뀌었다. 고정 응답행렬을
조건으로 한 5도 블록 bootstrap 구간은 [−0.971, −0.281]이지만, 44쌍 중 21쌍의
확산광 성분과 18쌍의 위성 성분이 음수였고 절대 역산 진폭은 각각 최대 1,840 및
2,211 nJy arcsec⁻²에 달했다. 응답행렬 불확도와 측정 공분산을 전파하지 않은
zero-intercept 선형 역산이므로, 이 구간은 반상관 검출이 아니라 현 전달행렬의
불안정성을 나타낸다. [행렬 역산 평가](hsc_debug/icl_satellite_transfer_matrix_deconvolution_assessment_2026-11-06.md),
Astrafek 재현기 `hsc_debug/analyze_icl_satellite_transfer_matrix_deconvolution.py`,
산출물 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_transfer_matrix_deconvolution_2026-11-06_01/`.
교정 전까지 과학 게이트는 유지한다.

추가로 같은 PSF-convolved 대상 137개를 위성 주입 진폭 0.5, 1.0, 2.0
nJy arcsec⁻²에서 비교했다. 공통 유효 지원은 54쌍(i 12, i2 42)이며 1.0 진폭의
응답은 기존 교차누출 결과와 정확히 일치했다. 그러나 pairwise i2 matrix 역산의
부분상관은 진폭 0.5에서 −0.922, 진폭 2.0에서 −0.075로 달라졌고, 최대 조건수는
각각 2,685와 962였다. 음수·극단 진폭이 다수 발생했으며 bootstrap은 응답행렬을
고정했다. 따라서 음의 구간을 물리적 반상관으로 읽을 수 없고 현재 zero-intercept
pairwise 역산은 불안정하다. [진폭별 역산 평가](hsc_debug/icl_satellite_amplitude_grid_transfer_matrix_assessment_2026-11-06.md),
재현기 `hsc_debug/analyze_icl_satellite_amplitude_grid_transfer_matrix.py`, 결과
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_amplitude_grid_transfer_matrix_2026-11-06_01/`.
`science_measurement=false`를 유지한다.

### 2026-11-06 paired lensing 스택의 연결 HSC tract 영향도

기존 승인된 표적·무작위 shape 행만 사용해 cluster-minus-random weak-lensing
스택의 연결 HSC 영역 영향도를 계산했다. 435개 공통 지지 렌즈와 2,175개
matched random 중심에서 유효한 표적 pair 2,075,062행과 random pair
8,905,892행을 재적재하고, 7개 반경 구간의 저장된 표적·cross profile을
저장 정밀도에서 정확히 재현했다. HSC patch-QA tract 중심 840개를 2도 이웃
연결 성분으로 묶으면 4개 성분이며, 스택과 겹치는 성분은 3개다. 각 성분을
하나씩 제외하며 lens/random 중심과 source 중 제외 영역에 배정된 모든 행을
제거하고 전체 보정 estimator를 다시 계산했다.

가장 큰 영역을 제외하면 표본 렌즈 261/435개가 제거되고, ΔΣ는 0.1–0.2 Mpc에서
28.6%, 0.4–0.7 Mpc에서 26.8% 감소한다. 각각 5도 격자 jackknife 대각 오차의
3.01배와 3.29배다. 다음 영역은 141개 렌즈를 포함하며 최대 shift가 기존
5도 오차의 1.33배다. 나머지 영역은 33개 렌즈를 포함한다. 따라서 균등 면적
격자 공분산은 지배적 HSC 영역 영향도를 충분히 설명하지 못한다. 단 세 점유
영역과 근사 tract-center 그래프는 공분산 추정에 충분하지 않으므로 본 결과는
leave-one-region-out 영향도 진단으로만 사용한다. 양의 감소형 ΔΣ 프로파일은
탐색적이며 질량 측정 또는 렌즈 검출이 아니다. `science_measurement=false`를
유지한다.

tract-center graph를 1.6°, 2.0°, 2.4°, 3.0° adjacency에서 다시 계산한
별도 Astrafek 감사에서는 네 설정 모두 같은 4개 성분과 같은 공통표본 분배
(렌즈 261, 141, 33개; random 중심 1,305, 705, 165개)를 만들었다.
따라서 측정한 leave-one-region shift는 이 임계값 범위의 성분 연결 방식에
민감하지 않다. 하지만 이 점검은 nearest-center 경계와 실제 tract polygon을
검증하지 않으며 region covariance의 유효성을 높이지 않는다. 재현기
`hsc_debug/audit_wl_random_tract_partition_sensitivity.py`, 결과
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_tract_partition_sensitivity_2026-11-06_01/summary.json`.

Astrafek 재현기 `hsc_debug/analyze_wl_random_tract_region_influence.py`, 결과
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_tract_region_influence_2026-11-06_01/summary.json`,
상세 [영역 영향도 평가](hsc_debug/wl_random_tract_region_influence_assessment_2026-11-06.md).
새 HSC CAS 조회나 다운로드는 수행하지 않았다. 후속 lensing 분석은 실제 survey
영역 geometry와 충분한 독립 영역을 이용한 공분산, random-center 전체 조리개
depth/mask/PSF/sky 검증, photo-z 및 shear 보정 불확도, reduced shear와
magnification 검토를 요구한다.

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

### 2026-11-06: ICL 세기별 위성은하 외곽 프로파일 스택

기존 94개 ICL–위성 공통 지지 쌍(i 26개, i2 68개)의 은하단별 위성 초과
표면밝기 프로파일을 200–300 kpc ICL 신호의 필터별 중앙값으로 나눠 일곱
반경 구간에서 비교했다. 적색편이와 log richness를 보정한 high-minus-low
대비의 공간 블록 CR1 95% 구간은 14개 모두 0을 포함했다. 10,000회 제한형
Webb wild-cluster bootstrap의 Holm 보정 p값은 모두 1.0, max-T 보정 p값의
최솟값은 0.520이었다. i2 추정치는 모두 양수지만 선택된 표본의 공간 상관을
고려하면 검출되지 않았고, i 추정치는 불확실하며 반경에 따라 부호가 바뀐다.
따라서 현재 결과는 ICL–위성 연관성을 입증하지 않는다. 분석은 기존 자료만
사용했으며 `science_measurement=false`를 유지한다. 구현·결과·한계는
[`icl_strength_split_satellite_profiles_assessment_2026-11-06.md`](hsc_debug/icl_strength_split_satellite_profiles_assessment_2026-11-06.md),
재현 스크립트는 `hsc_debug/analyze_icl_strength_split_satellite_profiles.py`,
그림 스크립트는 `hsc_debug/plot_icl_strength_split_satellite_profiles.py`에
기록했다. 산출물은
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_icl_strength_stacked_profiles_2026-11-06_01/`에 있다.

### 2026-11-06: 연속형 ICL–위성 프로파일 관계

같은 94개 쌍의 은하단별 위성 프로파일에서 ICL 200–300 kpc 초과 신호를
연속 예측변수로도 사용하고, 필터별 중앙값 분할과 함께 적색편이 및 log
richness를 보정했다. 양 방식·두 필터·일곱 반경에 대한 28개 검정을 같은
공간 블록 Webb bootstrap으로 묶었다. 최저 원 p값은 i2 50–75 kpc의 0.0369지만
max-T 보정 p값은 0.770이다. 전체 가족의 최소 max-T p는 0.190이며 Holm
보정 p값은 모두 1.0이다. i 내부 세 구간의 pointwise CR1 구간은 음수지만
가족 단위 공간 재표집 검정을 통과하지 못했다. 따라서 이 선택 표본에서
ICL–위성 관계를 확증하지 못하며 `science_measurement=false`를 유지한다.
세부 결과는
[`icl_satellite_radial_relations_assessment_2026-11-06.md`](hsc_debug/icl_satellite_radial_relations_assessment_2026-11-06.md),
재현 코드는 `hsc_debug/analyze_icl_satellite_radial_relations.py`, 그림 코드는
`hsc_debug/plot_icl_satellite_radial_relations.py`이다. 산출물은
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_radial_relations_joint_2026-11-06_01/`에 있다.


### 2026-11-06: weak-lensing estimator-weight boost sensitivity

Using only the saved S19A shape rows, calculated the candidate profile formed
by multiplying each random-subtracted DeltaSigma bin by its exact
response-adjusted estimator-weight boost. The product and covariance were
recomputed within paired 2.5- and 5-degree delete-one spatial-block
replicates. The seven-bin candidate rises by 12.5–16.3%, with matching central
values across the two block widths. This quantifies sensitivity to a proposed
boost only; it does not validate member dilution or separate it from
magnification, depth, mask, PSF, sky, or selection effects. The saved science
profile was not changed, `boost_correction_applied=false`, and
`science_measurement=false`. No new HSC data were requested. See
[`hsc_debug/wl_boost_profile_sensitivity_assessment_2026-11-06.md`](hsc_debug/wl_boost_profile_sensitivity_assessment_2026-11-06.md)
and output `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_source_density_boost_spatial_2026-11-05_06/`.


### 2026-11-06: forced-position satellite support versus annulus coverage

An existing-data sensitivity tested the 3,197 CAMIRA positions without SEP
detections. Turning off masks around other unmatched CAMIRA positions did not
recover support under the original 60% target/control annular-coverage cut.
At an exploratory 20% floor with at least three rotated controls, the
50–75 kpc forced stack increased from 9 to 65 position profiles and the
combined cluster/filter support increased from 68 to 74 in `i2`. The
redshift/richness-adjusted forced-only `i2` relation is imprecise
(`r=0.147`, 95% connected-block interval −0.485–0.425). On the 68 pairs in
common, the paired correlation change from 60% to 20% is +0.003
(−0.070–0.126); adding six pairs changes it by −0.025 (−0.127–0.079).
The coverage floor is post-hoc, masks can be non-random, and CAMIRA locations
are not confirmed galaxies. This sensitivity does not establish a physical
ICL–satellite association or replace the science sample. Keep
`science_measurement=false`. See
[`hsc_debug/icl_satellite_forced_coverage_sensitivity_assessment_2026-11-06.md`](hsc_debug/icl_satellite_forced_coverage_sensitivity_assessment_2026-11-06.md).


### 2026-11-06: radial selection in forced satellite-position support

A follow-up audit of the 20%-coverage forced-position run shows severe
host-centric radial selection at the satellite 50–75 kpc profile. None of the
65 retained profiles is from a candidate inside 100 kpc of its cluster center;
41/65 are at 300–500 kpc. The parent candidate-radius median is 195–200 kpc,
versus 339–388 kpc among retained profiles. Only 18/282 candidate-bearing `i`
clusters and 32/439 `i2` clusters retain any profile. This violates positivity
for a correction to the full parent population and means forced profiles must
remain separate from the detected-satellite sample. The combined `i2`
`r=0.257` result is not a completeness-corrected relation. See
[`hsc_debug/icl_satellite_forced_position_selection_assessment_2026-11-06.md`](hsc_debug/icl_satellite_forced_position_selection_assessment_2026-11-06.md).


### 2026-11-06 위성 스택 구성원 수 민감도

기존 94개 공통지원 쌍에서 50–75 proper-kpc 위성 스택에 실제 포함된 검출 구성원 수가 ICL–위성 부분상관을 설명하는지 검사했다. 개별 구성원 표에서 각 쌍의 수를 다시 세고 저장된 `n_satellites`와 대조했다. 기준 모형은 redshift 및 선형 CAMIRA richness를 통제했다. 확장 모형은 log 구성원 수도 통제했다. 두 필터가 공유하는 표적·대조 5도 연결 블록을 10,000회 재표집하도록 요청했으며 9,989회에서 두 필터의 유한 추정치를 얻었다.

`i2` 기준 부분상관 `r=0.293`은 구성원 수 추가 후 `0.295`로 바뀌었다. 짝지은 차이는 `+0.002`이며 95% 구간은 `[-0.023, 0.044]`이다. 구성원 수와 ICL 대비의 부분상관은 `0.057` `[-0.185, 0.313]`이고, 구성원 수와 위성 밝기 초과의 부분상관은 `−0.032` `[-0.226, 0.136]`이다. `i` 부분상관은 `0.079`에서 `0.088`로 바뀌며 차이는 `+0.009` `[-0.081, 0.192]`다. 모든 구성원 수 연관 구간과 모형 변화 구간은 0을 포함한다. 따라서 측정된 구성원 수가 선택된 표본의 양의 `i2` 점추정치를 설명한다는 증거는 없다. 검출 completeness, source-Poisson, PSF wing, sky-transfer 및 coadd 공분산은 해결되지 않았다. 구성원 수 조정은 정밀도 민감도일 뿐 인과 보정이 아니다. `science_measurement=false`를 유지한다.

재현기 `hsc_debug/analyze_icl_satellite_membercount_adjustment.py`, 평가 `hsc_debug/icl_satellite_membercount_adjustment_assessment_2026-11-06.md`, 결과 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_membercount_adjustment_2026-11-06_02/`에 상세 수치와 입력 해시를 기록했다. 기존 로컬 표만 사용했고 새 HSC 질의, 다운로드 또는 영상 픽셀 읽기는 없었다.


### 2026-11-06: 위성 스택의 은하단 중심거리 구성 민감도

고정된 94개 ICL–위성 공통지원 쌍에서 실제 50–75 proper-kpc 스택 구성원의 host-centric 반지름 중앙값을 산출해 기준 부분상관 모형(redshift와 선형 CAMIRA richness 통제)에 추가했다. `i2` 부분상관은 0.293에서 0.322로 바뀌었고 짝지은 변화는 +0.030 (공유 5도 블록 부트스트랩 95% 구간 −0.019–0.101)이다. `i`에서는 0.079에서 0.029로 바뀌며 변화 −0.050 (−0.330–0.044)이다. 반지름 변수의 부분 연관 및 두 모형 변화의 구간은 모두 0을 포함한다. 측정된 위치 구성은 선택된 표본의 `i2` 점추정치를 설명하지 못하며, 이 민감도 분석은 전체 모집단 또는 인과 보정이 아니다. 검출 completeness, source-Poisson, PSF wing, sky-transfer 및 총 영상 공분산은 여전히 미해결이다. `science_measurement=false`를 유지한다.

재현기 `hsc_debug/analyze_icl_satellite_clustercentric_position_adjustment.py`, 평가 `hsc_debug/icl_satellite_clustercentric_position_adjustment_assessment_2026-11-06.md`, 결과 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_clustercentric_position_adjustment_2026-11-06_01/`에 기록했다. 기존 로컬 표를 사용했으며 새 HSC 질의·다운로드나 픽셀 읽기는 없었다.


### 2026-11-06: S20A 확장 coverage metadata 제출 감사

기승인된 S20A 2,596개 후보의 9점 coverage metadata 계획(23,364 probe, 52 batch)을 lageunha에서 제출했으나 첫 CAS 요청이 HTTP 406으로 거절됐다. 조사 중 mosaic 내부 JOIN에 조건이 빠진 SQL 오류를 발견해 `m.band = 'HSC-I'` 조건을 JOIN에 넣었다. 새 계획의 후보 중심과 probe 파일은 승인본과 바이트 단위로 같고 열 및 데이터 범위도 유지했다. 수정 계획도 첫 batch에서 HTTP 406을 반환했다. 두 체크포인트 모두 CAS job ID 없이 첫 batch가 pending이며 결과 파일이 없다. 오류 원인은 미확정이므로 CAS 작업 상태와 거절 이유 확인 전 재시도하지 않는다. 이미지·source row·shape·photo-z 반출은 없고 `science_measurement=false`다. 상세 감사 `hsc_debug/icl_expansion_coverage_submission_assessment_2026-11-06.md`.


### 2026-11-06: ICL 및 weak-lensing 증거 종합

최신 Astrafek 산출물 13개의 SHA-256을 고정해 한 보고서로 종합했다. ICL–위성 공통지원 94쌍에서 `i2` 부분상관은 raw `r=0.293` (95% 블록 구간 0.039–0.500), 유한-lag 픽셀 잡음 후 `0.306` (0.045–0.516)이다. 지원 확률 가중과 같은 잡음 보정을 결합하면 `r=0.238` (−0.023–0.446; 유효 N 59.9/68)로 구간이 0을 포함한다. 구성원 수/host-centric 위치 조정은 점추정치를 크게 설명하지 못하지만, 주입 전달행렬이 불안정하고 sample selection 및 source-Poisson, PSF, sky-transfer 보정은 해결되지 않았다.

weak-lensing cluster-minus-random stack은 435 렌즈와 2,175 무작위 중심을 사용한다. 세 연결 tract 성분 중 최대 성분 261개 렌즈를 빼면 내측 프로파일이 28.6%, 0.4–0.7 Mpc 구간이 26.8% 변한다. ICL 연속관계는 i2 70개 렌즈에서 다중검정 조정 후 유의하지 않고, 제안된 boost 보정은 적용되지 않았다. 따라서 두 결과 모두 탐색적 상태이며 `science_measurement=false`를 유지한다.

재현기 `hsc_debug/build_icl_wl_current_evidence_report.py`, 종합문서 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_03/icl_wl_current_evidence_report.md`에 보고했다. 새 HSC 질의·다운로드·영상 픽셀 읽기는 없었다. S20A coverage metadata 계획은 별도로 CAS HTTP 406으로 미완료 상태다.

### 2026-11-06: HSC warp 단일 위치 분산 검증 및 공통지원 94쌍 확장 계획

사용자가 승인한 target 6966의 HSC-I `warp` 9개 방문 영상·mask·variance를
55,511,040 bytes (52.92 MiB) 수신했다. 두 방문은 유효 영상/분산 픽셀이 없고
한 방문은 59.5% coverage다. FITS 배열과 WCS grid 정렬은 확인했다. 200–300
proper-kpc 고리의 162,352개 primary 픽셀에서 coadd/predicted variance 중앙값은
1.00748 (q16–q84 1.00372–1.01506)이고, 공통 mask 적용 시 면적 48.83%가 남으며
중앙값은 1.00656이다. 인접 픽셀은 상관되므로 이 분위수는 오차구간이 아니다.
한 위치의 제품 일관성 검증으로만 해석하며 source-Poisson/effective gain,
resampling noise, sky-transfer covariance, PSF wing 또는 전체 ICL 공분산을
검증하지 않았다. `science_measurement=false` 유지.

이 결과를 지원 94개 ICL–위성 target/control 쌍(i 26, i2 68) 전체의 warp
분산 비교로 확장하는 새 계획을 준비했다. 승인된 frame metadata에 표적·대조
188개 위치 모두 존재한다. 요청 목록은 188개 HSC-I 2′×2′ 위치 요청이며,
metadata에서 1,471 visit product(i 379, i2 1,092)를 예상한다. pilot 실측에
비례한 전송 예상은 9,072,971,094 bytes (8.45 GiB)이고 요청별 강제 상한 160 MiB, 전체 12 GiB에서
중지한다. 각 요청에 image/mask/variance를 포함하며 catalog, shape, photo-z,
다른 필터 또는 광역 coadd는 포함하지 않는다. 표적/control matched pair와
공간 block 정보를 결과 해석에 유지한다. plan SHA-256
`bd936a94bcdbff83a607acbfca66616c396db526ada50cc89c9e3e312d37007a`;
plan은 `prepared_not_submitted`, 현재까지 네트워크 요청은 0건이다. 다음 단계는
이 명시 범위의 사용자 승인을 받은 뒤에만 제출한다.

계획기 `hsc_debug/prepare_icl_supported94_warp_noise_plan.py`, 승인-gated
다운로더 `hsc_debug/download_icl_supported94_warp_noise.py`는 기본적으로 plan-only다.
제출 경로는 plan hash를 이중 검증하고 순차 다운로드·체크포인트 재개·개별 및
전체 byte cap을 강제한다. 이번 계획 준비에서는 HSC 조회/다운로드를 하지 않았고
`science_measurement=false`다. 파일 목록, hashes, 결과 및 과학적 한계는
`hsc_debug/icl_warp_pilot_coadd_variance_assessment_2026-11-06.md`와
`HANDOVER.md`에 있다.
계획기 SHA-256 `475cb7a3436415f8f920ab2fe04e9c20ee28f2565373ddd2f598c8e4fc2ca059`;
다운로더 SHA-256 `1fe8720f87ae31c951eadd40d83988402d89971f4f2d37ea1bc08da8c2da8da53`.


## 2026-11-06: supported-94 HSC warp variance audit

The user-approved 188-cutout HSC-I warp acquisition completed with 184 archives (9.26 GiB) and four repeat-404 data gaps. All expected visits were present in returned archives. The offline Astrafek comparison passed 184 centers and retains 88 matched pairs (24 `i`, 64 `i2`). The delivered coadd/predicted visit-variance ratio has center-median 1.00669 at 200–300 proper kpc; strict common-mask support covers a median 81.9% of the full annulus. Paired target-minus-control ratios are consistent with zero in both filters under connected-sky-block bootstrap. This only checks approximate variance-plane consistency; nearest-neighbor mapping is not the exact HSC resampling kernel, and source Poisson, correlated resampling, sky covariance, PSF wings, and total ICL-estimator covariance remain unresolved. No ICL association claim follows. Keep `science_measurement=false`. See `hsc_debug/icl_supported94_warp_noise_assessment_2026-11-06.md` and its pinned outputs.


### 2026-11-06: integrated 94-pair warp variance evidence

The consolidated ICL/weak-lensing report now includes the approved 184-center warp variance analysis as a hash-validated 23rd input. It records 4 repeated-404 gaps, 9.26 GiB downloaded, variance ratio 1.00669 at 200–300 proper kpc, and strict shared-mask median support 81.9%. The target/control block-bootstrap intervals include zero. This remains an approximate variance-product check, not exact HSC resampling covariance or a science detection. Updated report: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_17/icl_wl_current_evidence_report.md`; `science_measurement=false`. The next gate is kernel-aware resampling/covariance, then source-Poisson, PSF-wing, and sky-transfer validation.


### 2026-11-06: strict warp-mask ICL association sensitivity

The strict visit-valid mask sensitivity has been measured on already downloaded warp cutouts. Baseline coadd ICL aperture fluxes reproduce the saved measurements; applying the at-least-three-valid-visits/shared-mask criterion and preserving the 40% annulus floor retains 81 of 90 warp-complete pairs (23 i, 58 i2). On those same pairs, partial-r changes by +0.0020 (`i`) and −0.0048 (`i2`), with spatial-block paired intervals containing zero. Nine additional pairs fall below the coverage floor, so this does not correct selection in the full 94-pair sample. Report as an image-support sensitivity only; keep `science_measurement=false`. Assessment: `hsc_debug/icl_supported94_warp_mask_association_assessment_2026-11-06.md`.


### 2026-11-06 candidate expansion: PDR3 mosaic filter-column correction

Fixed the S20A/S23B candidate-coverage SQL generators to use the HSC PDR3
`mosaic.filter01` column (`filter01 AS band`) instead of `mosaic.band`. This is
the likely cause of the earlier HTTP 406 response, but not confirmed because
the response body was suppressed. Candidate IDs, center and 1-Mpc boundary
probes, geometry, output schema, and metadata-only scope are unchanged.

Offline plan checks passed for the corrected 2,596-candidate S20A plan (23,364
probes, 52 batches; SHA-256
`96ad3b9078beff082d3dce4de1f9bb7608791915ffc7a46a656d9e6cf30e1149`) and
1,333-candidate S23B plan (11,997 probes, 27 batches, upper bound 47,988
metadata rows; SHA-256
`ee6d2b11893dab0e9ceb0c0d2a3c80c2dee27b600683839d2be294441578ef05`). No CAS
query was submitted. The S23B plan still needs approval for this exact hash.
These requests return only patch footprint, depth, and seeing metadata; no
images, source rows, shapes, or photo-z. Follow with actual WCS/mask aperture
support and matched-control checks before estimating any image expansion.
Details: `hsc_debug/icl_expansion_coverage_filter_column_assessment_2026-11-06.md`.
Keep `science_measurement=false`.

## 2026-11-06 update: forced-position support sensitivity

A paired comparison of the detected-only and 20/20%-coverage satellite profiles
now preserves all 224 ICL-supported pairs in its selection ledger. The 20/20
floor retains 50/71 `i` and 110/153 `i2` pairs. Retained pairs have higher
median 200–300 kpc ICL residual than excluded pairs in both filters, while
redshift and richness medians are similar. Therefore the selected support is
associated with the ICL outcome; the analysis cannot be transported to all
ICL-supported systems as if missingness were random. Paired spatial-block
intervals for changes in partial association include zero. Treat the outcome as
a post-hoc support sensitivity only; do not apply a correction or claim a
physical ICL–satellite relation. Preserve `science_measurement=false`. Full
methods and numbers are in
[`hsc_debug/icl_satellite_forced_coverage_both_sides_assessment_2026-11-06.md`](hsc_debug/icl_satellite_forced_coverage_both_sides_assessment_2026-11-06.md).

## Current evidence synthesis: 2026-11-06

The latest hash-pinned Astrafek ICL/weak-lensing report is version 27 at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_27/icl_wl_current_evidence_report.md`
(SHA-256 `af88c214eadad945c5177e701f4d7496f95222b297ad84e0be2ba91f595f8b0f`).
It includes the 20/20 forced-position coverage sensitivity. Outcome-related
selection remains a key limitation, and neither the ICL–satellite relation nor
the weak-lensing profile is promoted to a science measurement.

## 2026-11-06 spatial support audit for forced satellite profiles

The 20/20 coverage rule has selected and excluded ICL-supported pairs in only
9/18 `i` and 21/40 `i2` connected paired-sky blocks. Within mixed blocks, the
equal-block ICL residual contrast remains positive in `i2` (+3.73 nJy arcsec⁻²,
95% block interval +0.55 to +6.93) and is unresolved in `i` (+1.80, −3.26 to
+6.66). This descriptive result is unadjusted for redshift/richness, uses a
post-hoc threshold, and cannot validate transport to blocks without both
selection states. It reinforces the need to model the support mechanism before
interpreting the ICL–satellite relation; it is not a physical detection. See
[`hsc_debug/icl_forced_coverage_spatial_selection_assessment_2026-11-06.md`](hsc_debug/icl_forced_coverage_spatial_selection_assessment_2026-11-06.md).

## 2026-11-06: covariate-adjusted support selection

Within 5-degree blocks that contain both 20/20-selected and excluded pairs, a
block-fixed-effect model controlling redshift and log richness estimates the
selected-minus-excluded ICL residual as +1.47 nJy arcsec⁻² in `i` and +2.13 in
`i2`. The block-bootstrap 95% intervals (−2.10 to +4.79 and −0.64 to +5.91)
include zero. This attenuation means the unadjusted within-block i2 difference
does not survive the current covariate adjustment. The result remains restricted
to 9/21 mixed blocks, with post-hoc support selection and no measurement-error
propagation; it does not supply a selection correction. Details are in
[`hsc_debug/icl_forced_coverage_adjusted_selection_assessment_2026-11-06.md`](hsc_debug/icl_forced_coverage_adjusted_selection_assessment_2026-11-06.md).

## 2026-11-06 radial forced-position sensitivity

The new radial analysis relates the 200–300 proper-kpc ICL residual to mean
light at unmatched CAMIRA positions in the existing 20/20%-coverage profiles.
Across nine supported filter/radius tests with redshift and richness control,
connected 5-degree blocks, and 10,000 shared restricted Webb draws, no
single-step max-T p-value is below 0.05. The i 50–75 kpc nominal association
(p=0.030) has max-T p=0.059 and Holm p=0.274. Since the positions lack unique
SEP identities, this is unmatched-position light, not a confirmed satellite
measurement. Support is outcome-related, inner bins are sparse, and image
measurement/noise/PSF/sky uncertainty is not propagated. Keep
`science_measurement=false`; see
[`hsc_debug/icl_forced_position_radial_relations_assessment_2026-11-06.md`](hsc_debug/icl_forced_position_radial_relations_assessment_2026-11-06.md).

The hash-pinned integrated evidence report is now version 30 with 91 validated
inputs at `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_30/icl_wl_current_evidence_report.md`.
Report SHA-256: `62a41b863e1ee2d2150c8647829626972fb3a99d5b123cd993289cc130a65515`.
The approved S23B metadata submission received HTTP 401 before job creation;
no rows were returned. The analysis status remains exploratory and
`science_measurement=false`.

The same approved metadata scope later completed after fixing the CAS host and
geometry SQL. It returned 13,575 patch-coverage/QA rows for 1,333 candidates
and 11,997 probes in 27 batches. 1,043 candidates have HSC-I mosaic matches at
their centers and all eight 1-Mpc boundary points; product-availability flags
are present at all nine points for these candidates. This does not verify
full-area WCS/masks or controls, and no images were downloaded. The separate
image step still requires explicit approval. Full audit:
[`hsc_debug/icl_s23b_coverage_metadata_assessment_2026-11-06.md`](hsc_debug/icl_s23b_coverage_metadata_assessment_2026-11-06.md).

The current integrated evidence report is version 31 with 124 validated
inputs: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_31/icl_wl_current_evidence_report.md`
(SHA-256 `3d6aebd01cb2727ca21d2617ef2de1fd17a85c526199315510710f167a1acf15`).

## 2026-11-06 S23B 표본 확장 게이트

기존 800개 영상 표본의 외부 타당도와 약한렌즈 공간 공분산 한계를 줄이기 위해,
HSC-I metadata screen을 통과한 S23B 1,043개 후보의 대조 위치 후보를 찾는
metadata-only HSC PDR3 random query plan을 준비했다. 각 표적 중심에서 30′ 안의
random 위치 최대 200개와 ID/좌표 및 I-band quality/depth scalar 29개 필드를
요청하는 범위로 제한했다(최대 208,600행, 약 80 MB, 105 batch). SQL 및 입력은
오프라인 해시 검증을 통과했고 2026-11-06 사용자 승인 후 LagEunha에서 CAS
제출·수신을 시작했다. 실행 상태와 체크포인트는 HANDOVER.md에 기록한다.
이 자료는 대조 후보 선별용이며 images, source detections, shape, photo-z를
포함하지 않는다. 표적/대조 이미지 다운로드에는 별도 승인과 별도 쿼리 계획이
필요하다. 승인·실행 후에도 S23B는 기존 discovery sample의 독립 공간 블록을
크게 늘리지 못할 수 있으므로, 표본 수와 지역 중첩을 다시 검증한다.

계획 생성기: `hsc_debug/prepare_s23b_icl_control_candidate_plan.py`.
고정 계획 SHA-256:
`143a1bd0a3a3c341576ce420f77ba0089fbfe80fdf58fbfd7afe30ea8e3f6714`.

## 2026-11-06 radial ICL profile contrast

The saved exploratory 2.5-Kron/40%-common-coverage stack was audited across all
six annuli from 50 to 500 proper kpc. A fixed contrast subtracts the area-
weighted 300–500 kpc target-minus-control profile from the 200–300 kpc value.
The complete-case support contains 20 `i` pairs in 10 connected sky blocks and
48 `i2` pairs in 30 blocks. The `i` contrast is positive in raw and
outer-plane-corrected profiles. Both `i2` block-bootstrap intervals include
zero. The 300–500 kpc signal can contain real cluster light, and the cut is
coverage-selected. The analysis gives a relative radial measurement, not an
absolute ICL estimate. The measurement-error budget remains incomplete, so the
science gate stays closed. Reproducibility is documented in
`hsc_debug/icl_radial_outer_baseline_contrast_assessment_2026-11-06.md` and
integrated report version 32.

## 2026-11-06 ICL–satellite outer-profile relation

The ICL relative radial contrast was joined to the existing detected-member
satellite-centered 50–75 proper-kpc target-minus-azimuthal-control residual
profile, controlling cluster redshift and log CAMIRA richness. Partial
correlations are −0.140 (raw) and −0.126 (outer-plane corrected) in `i`
(20 clusters/10 connected blocks), and +0.417 and +0.391 in `i2` (47/30).
All four 95% connected-block intervals include zero. One i2 profile is
unsupported and remains missing. This does not resolve an ICL–satellite
association; the detected-member sample is incomplete and the uncertainty
budget omits measurement error, source-Poisson noise, coadd covariance, sky
transfer, and extended PSF wings. Keep `science_measurement=false`. Assessment:
`hsc_debug/icl_satellite_outer_profile_relation_assessment_2026-11-06.md`.
The integrated report is version 33 with 940 independently verified inputs at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_33/icl_wl_current_evidence_report.md`
(SHA-256 `1ef9acc5bff699fed8e1d75131192cc03e88d48832e8ea485918b70c94832742`).

## 2026-11-07 S23B CAMIRA member-catalog intake gate

The approved S23B target manifest has 1,043 star-masked CAMIRA clusters and no
member coordinates (`member_catalog_sha256=null`). A metadata-only listing of
the Drive folder linked by the catalog authors shows the candidate file
`camira_s23b_wide_sm_member.dat.gz` (Drive file ID
`1DEnjGj9P1ilUooEm2U5_A-qXUrxeiJdE`, 10,626,590 compressed bytes). Its `sm`
variant matches the selected `camira_s23b_wide_sm.tbl` cluster catalog. The
file contents and external sharing visibility have not been verified, and it
has not been downloaded. The exact one-file intake plan is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_import_plan_2026-11-07_01/plan.json`
(SHA-256 `724860324ee069e1b07c1e324ca48c14fde9c6278343b4f645044c002e4eeb56`),
generated by `hsc_debug/prepare_s23b_camira_member_import_plan.py`. It excludes
the unmasked member table, LRG FITS catalogs, images, shapes, and photo-z rows;
explicit approval is still required before receiving the file. After approval,
inspect its schema and verify cluster-ID/member-coordinate coverage before
using it in satellite masks or profiles. Until then, the S23B member-light
analysis gate remains closed and `science_measurement=false`.

## 2026-11-07 S23B target pixel-support audit

The approved target-image transfer is complete at 1,043/1,043. The independent
header/WCS/file audit passed all 1,043 frames, with primary filter counts `i`
251 and `i2` 792. The subsequent local mask/area pass found four cutouts with
zero usable pixels: S23B_10647, S23B_11786, S23B_432, and S23B_8226. In all
four, every pixel is flagged both `MP_INTRP` and `MP_NO_DATA`, so the
post-exclusion support is 1,039 targets (`i`: 251, `i2`: 788).

The automated K=2.5 masks leave median annular area fractions of 0.350 at
50–100 kpc and 0.417–0.423 across 200–500 kpc. At 50–100 kpc, 144/1,039
targets have zero unmasked area and 331/1,039 have at least 50%. The exact
descriptive distribution is in
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_mask_support_assessment_2026-11-07_01/summary.json`
(SHA-256 `c140d6e02862b6d7654494a9c243eba5e664fb04234d25d5ae3e18e2ecb27c47`).
Filter-separated K=2.5/3.5/5.0 diagnostic plot:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_mask_support_assessment_2026-11-07_01/radial_unmasked_area_by_filter_v4.png`.
This pass used no CAMIRA member coordinates; 434 centers have a unique
automated central-source candidate, which is not a manual BCG confirmation.
The exact center lies in a base quality mask for 446 frames and in the K=2.5
source mask for 824. These support values are diagnostic only: no coverage
threshold or pair-selection rule has been chosen, and the `i`/`i2` images must
remain separate. Matched-control images, paired common-support coverage,
extinction, sky transfer, PSF wings, and covariance are still outstanding.
Keep `science_measurement=false` and do not interpret the support output as an
ICL stack.

## 2026-11-07 S23B matched-control image transfer checkpoint

The explicitly approved v2 control manifest remains fixed at 1,041 random
positions (SHA-256
`ae9eb406cfc179d9980031bb67f85f8ee71c4580ba59e97713278823db5c06df`). The
transfer reached 627/1,041 before its validator stopped on
`CTRL_S23B_8547`. The 66,349,440-byte FITS response is locally present as a
`.part` file and has the expected 2,335×2,367 HSC-I cutout layout, but its
MASK header maps `MP_NOT_DEBLENDED` to bit 12 and `MP_UNMASKEDNAN` to bit 11;
the other downloaded controls and all targets map those two names to bits 11
and 12, respectively. HSC's pipeline workshop says bit-to-event correspondence
can vary and directs users to the MASK header
(https://hsc.mtk.nao.ac.jp/pipedoc/DataWorkshop2022/_downloads/2682edf0a494c8b005acf4efa51952ae/hsc_kaiseki_20221128.pdf).

The validator now checks the complete known plane-name set and a unique
0–16 bit assignment. Mask preparation and matched-control masking derive the
bit values from each FITS header, so a quality-pixel selection does not depend
on a global bit-number assumption. The transfer is being resumed from its existing checkpoint; no
new position or product was introduced. After it completes, run the approved-v2
control FITS/WCS audit, then generate control mask/area diagnostics. The S23B
stack driver is `hsc_debug/stack_s23b_icl_target_control_2d.py`; it requires
both complete independent image audits and pinned provenance, excludes only
the four zero-support targets already identified, separates `i` from `i2`,
and writes exploratory products with `science_measurement=false`. It has not
been run.

The exact-scope resume recovered the saved `CTRL_S23B_8547` response after
reviewing its alternate header-declared mask map and continued the same v2
transfer. At the latest checkpoint, 734/1,041 are complete and the approved
manifest SHA remains unchanged. The pre-stack sequence is now explicit:
complete the transfer and FITS/WCS audit; prepare all 1,041 control masks;
audit target/control/common warped mask area without reading science flux;
review the radial support and filter counts; only then run the exploratory
flux stack. The stack driver requires the support-audit artifact and excludes
only audited zero-support targets, filter mismatches, and pairs with no common
mask support within 1,000 proper kpc. Per-bin coverage remains reported; no
minimum coverage cut is adopted in advance. No S23B flux stack has been run.

### 2026-11-07 — S23B exploratory ICL target/control stack

The approved 1,041-frame matched-control download is complete and independently
verified. The mask-only paired support audit yields 1,022 filter-matched pairs
(`i`: 241; `i2`: 781). Common support is sparse in the center: only 23–25% of
pairs have support at 0–25 proper kpc and median common area is zero there.
Median common area is about 10–20% from 50 to 1,000 kpc. The exploratory
target-minus-control stack shows a candidate positive excess at 50–300 kpc
that declines toward zero farther out, but current systematics and the
between-pair SEM are insufficient for an ICL measurement. Keep the result at
`science_measurement=false`; complete Galactic-extinction, PSF-wing,
residual-sky, coadd-covariance, spatial-jackknife, and mask/BCG sensitivity
work before scientific interpretation. Full quantitative assessment:
`hsc_debug/s23b_icl_stack_quality_assessment_2026-11-07.md`.

## 2026-11-07 S23B spatial-block uncertainty follow-up

The existing S23B matched target/control stack now has a pair-level spatial
block-bootstrap sensitivity analysis in
`hsc_debug/analyze_s23b_icl_spatial_blocks.py`. It resamples connected sky
cells at 2.5/5/10-degree widths and two grid phases, using only the saved
radial-profile table and pinned manifests. The analysis adds spatial clustering
to the uncertainty estimate but is not a substitute for official HSC region
geometry or pixel-level covariance.

Across all six block layouts, the 50–100 and 100–300 proper-kpc target-minus-
control profiles remain positive in both `i` and `i2`. The 300–800 kpc interval
contains zero throughout. The apparent `i2` offset at 800–1,000 kpc is
layout-sensitive, which reinforces that the outer profile cannot yet be used
as a sky correction. The block-bootstrap intervals still omit residual-sky
transfer, extinction, extended PSF, correlated coadd noise, and centering/mask
uncertainty. Keep `science_measurement=false`.

Results are in
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_spatial_blocks_2026-11-07_02/`.
The CAMIRA member file remains a separate, approval-gated input; do not infer
satellite positions from the cluster-center catalog or substitute the LRG
catalog. Once approved and schema-validated, the next scientific branch is a
member-centered, matched-control satellite stack joined to this ICL profile,
with completeness and injection-recovery carried into the association model.


### Paired inner-minus-outer radial contrast (2026-11-07)

The saved S23B target/control profiles now have a paired contrast analysis:
for each pair, compare the common-support inner radial ΔSB (50–300, 100–300,
or 200–300 proper kpc) against 300–800 proper kpc, then bootstrap connected
sky blocks across six grid layouts. The 50–300 and 100–300 kpc contrasts stay
positive in both HSC-i and HSC-i2 across all layouts. The 200–300 kpc contrast
is smaller and crosses zero in at least one layout for each filter. See the
runner `hsc_debug/analyze_s23b_icl_relative_radial_contrasts.py` and the
assessment in
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_relative_radial_contrasts_2026-11-07_02/`.
The 300–800 kpc annulus may contain genuine cluster light, so this is a radial
contrast, not an absolute ICL or blank-sky estimate. Residual sky, extinction,
PSF wings, pixel covariance, and center/mask review remain gates;
`science_measurement=false`. The CAMIRA member catalog remains approval-gated
and has not been downloaded, so no satellite association is measured yet.


### S23B matched-control image transfer and integrity gate (2026-11-07)

The explicitly approved v2 control-image transfer is now complete: all 1,041
HSC PDR3 `pdr3_wide` HSC-I `coadd/bg` IMAGE/MASK/VARIANCE cutouts were received
(58,271,158,080 bytes; 248 `i`, 793 `i2`). The offline integrity audit checked
checksums, FITS planes, filter headers, WCS, and the one-to-one control/target
mapping against the approved manifest SHA-256
`ae9eb406cfc179d9980031bb67f85f8ee71c4580ba59e97713278823db5c06df`.
The target-image audit is also complete (1,043 frames). Automated masks were
made for 1,039 targets; four failed with too few valid pixels. The mask-only
common-support audit retained 1,022 same-filter pairs and found no pair with
zero common support inside 1,000 kpc; an exploratory 2D stack has been produced
for those pairs. This does not validate sky residuals, extinction, PSF wings,
coadd covariance, or member completeness. Keep `images_ready_for_science=false`
and `science_measurement=false`. At the time of this earlier image-stack
assessment CAMIRA member positions had not yet been imported; the later
approved intake and member-centered analysis below supersede that status.


### S23B ICL contrast versus richness (2026-11-07)

An exploratory association now uses the 1,022 saved target/control pairs. For
each filter separately, the outcome is the 50–300, 100–300, or 200–300 kpc
relative contrast against 300–800 kpc. Equal-pair OLS adjusts for log CAMIRA
richness, redshift, and residual target-control SFD E(B−V) mismatch; connected
RA-Dec blocks are bootstrapped across six grid layouts. The primary 5°
phase-zero slopes are positive but all 95% intervals include zero. For i2,
two wider-annulus estimates become positive only at 10° phase zero and cross
zero under the alternate phase, so the result is block-layout sensitive.
This is not a satellite-profile measurement: richness is only a cluster-scale
proxy. No p-values or detection claim are made. Extinction is not corrected;
`science_measurement=false`. Runner:
`hsc_debug/analyze_s23b_icl_richness_relation.py`; outputs:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_richness_relation_2026-11-07_01/`.
At the time of this richness-only analysis, the member file had not yet been
received; the later intake and direct member-centered analysis below
supersede that status.


### Common-support sensitivity for the S23B richness relation (2026-11-07)

A second exploratory model adds log common target/control pixels in the inner
annulus to the primary richness regression, using the same block bootstrap
draws for both slopes. At 5° phase zero, the 50–300 kpc slope changes by
−0.16 [−5.15, 6.15] in `i` and −0.12 [−1.97, 1.03] in `i2`; all other
paired slope-change intervals also include zero. Log richness and log support
are weakly correlated (|r| ≤ 0.017). This suggests the linear trend is not
explained by a simple linear support correlation, but it does not correct
coverage selection. Keep `science_measurement=false`. See
`hsc_debug/analyze_s23b_icl_richness_coverage_sensitivity.py` and
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_richness_coverage_sensitivity_2026-11-07_01/`.


### S23B nominal Galactic-extinction sensitivity (2026-11-07)

The S23B target/control profiles were recomputed after applying separate
nominal hscPipe SFD corrections to each endpoint (`i=1.641`, `i2=1.617`;
[official hscPipe coefficient table](https://hsc.mtk.nao.ac.jp/pipedoc/pipedoc_7_e/colorterms.html)).
The estimator preserves the exact common-pixel support and uses the existing
connected target/control blocks. At 5° phase zero, the 50–300 minus 300–800
kpc contrast increases by 0.369 nJy arcsec⁻² (95% block interval 0.190–0.516)
in `i` and 0.444 (0.337–0.569) in `i2`, about 5.6–5.8% of the raw contrast.
The correction interval is positive in all six spatial layouts for both
filters. This is a nominal correction sensitivity only: it also rescales
residual-sky terms and does not correct PSF wings, sky transfer, source-Poisson
noise, coadd covariance, or mask/center uncertainty. Keep
`science_measurement=false` and `images_ready_for_science=false`.

Runner: `hsc_debug/analyze_s23b_icl_extinction_sensitivity.py`; assessment:
`hsc_debug/s23b_icl_extinction_sensitivity_assessment_2026-11-07.md`; output:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_extinction_sensitivity_2026-11-07_04/`.
The separate S23B CAMIRA member catalog was subsequently approved, received,
and used for the exploratory member-centered profile analysis documented in
the later 2026-11-07 section below. The existing weak-lensing sample remains
spatially under-supported for calibrated field covariance.
## 2026-11-07: Direct S23B ICL–satellite candidate relation

The separately approved single-file member-catalog intake is complete. The
earlier approval-gated status above records the state before this intake and
is superseded here. Only camira_s23b_wide_sm_member.dat.gz was received
(10,626,590 bytes); no other file or credential was downloaded/saved.

The importer verified the gzip stream and matched all 12,162 member groups
to the pinned cluster catalog, including all 1,043 target clusters. It
produced 45,634 target member positions. Of these, 39,166 fall within the
target-image pixel-center footprints. This audit establishes coordinates
and geometric coverage only; it does not establish usable pixels, member
purity, or completeness. The two trailing headerless fields remain opaque.

The first direct exploratory model uses 1,022 matched target/control pairs
and the common-support relative ICL contrast between 50–300 and 300–800
proper kpc. Its predictor is the unweighted CAMIRA candidate-member surface
density ratio in those annuli with a 0.5-count continuity correction. The
model adjusts for richness, redshift, target/control SFD E(B−V) mismatch,
and minimum common-pixel support; uncertainty uses 5,000 connected sky-block
bootstrap draws per layout.

At 5° phase zero, the adjusted slopes are 5.72 [−0.06, 13.63] for HSC-i
(232 pairs, 26 blocks) and 4.24 [0.99, 7.80] for HSC-i2 (756 pairs, 42
blocks), in nJy arcsec⁻² per dex. The intervals exclude zero in only 3/6
and 4/6 sky-grid layouts, respectively. With at least 10% common pixels
in each annulus, both filters' intervals include zero. The positive point
estimates are suggestive but do not establish a robust association.

The outer reference may contain cluster light; common-pixel support is
sparse; member candidates are unweighted and photometric; and positions
inside a footprint need not be on usable pixels. Residual sky, extinction,
extended PSF wings, complete coadd covariance, centering/masking, and member
purity/completeness remain uncalibrated. Keep science_measurement=false.

Next analysis gates: audit candidate positions against usable masks and
centering; establish the semantics and calibration of the member fields
before any weighting; measure member-centered satellite light with matched
control positions and injection/recovery; then propagate sky, PSF, and
spatially correlated noise uncertainty before considering a science result.
Use the importer and model in hsc_debug/import_s23b_camira_member_catalog.py
and hsc_debug/analyze_s23b_icl_satellite_relation.py. The dated assessment
records the exact input/output hashes and sensitivity results.

### S23B CAMIRA member-field-8 sensitivity (2026-11-07)

A raw-field weighting sensitivity gives a positive adjusted HSC-i2 relation
between candidate-member density contrast and common-support relative ICL
contrast at 50–300 versus 300–800 proper kpc: 5.01 [2.44, 7.62]
nJy arcsec⁻² dex⁻¹ for the 5-degree phase-zero layout. The interval excludes
zero in all six spatial layouts, but crosses zero after requiring at least
10% common pixels in both annuli. HSC-i is unresolved. The member file's
headerless eighth field is numerically probability-like but remains
unverified; it must not be called CAMIRA membership probability or w_mem yet.
This result is exploratory, not a detection. Keep science_measurement=false.

The member-position mask audit and the exploratory source-aware member-profile
stack have since been completed; they do not resolve the unverified member
column schema or calibrate member completeness. Remaining gates are schema
verification, member-preserving injection/recovery, and calibration of sky,
PSF-wing, extinction, source-Poisson, and coadd-covariance effects. Assessment:
`hsc_debug/s23b_icl_satellite_field8_weight_sensitivity_assessment_2026-11-07.md`.

### S23B member-centered outer-light profile and adjusted ICL relation (2026-11-07)

The member-position mask audit has been followed by a source-aware exploratory
profile stack and a spatial-block adjusted relation. The stack covers all 979
eligible target/filter frames (226 HSC-i and 753 HSC-i2), with 16,503 CAMIRA
candidate positions at 50–500 proper kpc and 13,269 retained candidate-by-bin
measurements. SEP-matched candidate sources are preserved; unrelated sources,
bad/bright pixels, and unmatched-position fallback circles are masked.
Target and rotated-control bins require at least 60% coverage and three
usable controls.

The median candidate-minus-rotated-control surface brightness falls from
459.6 to 0.68 nJy arcsec⁻² between the 5–10 and 50–75 kpc bins in HSC-i; in
HSC-i2 it falls from 583.9 to 1.46. The outer-bin q16–q84 cluster ranges are
[−3.77, 5.64] for HSC-i (62 clusters) and [−4.23, 8.94] for HSC-i2 (365
clusters). These are descriptive cluster spreads.

For the 50–75 kpc member-light outcome against the 50–300 minus 300–800 kpc
relative ICL contrast, the adjusted slopes at 5° phase zero are 0.086
[−0.098, 0.334] in HSC-i (62 clusters, 19 connected blocks) and −0.041
[−0.145, 0.090] in HSC-i2 (365 clusters, 43 blocks). None of the six spatial
layout intervals excludes zero in either filter. At 5%, 10%, and 20% ICL
support floors, none of the available intervals excludes zero; the fitted
sign changes with the floor. The association remains unresolved. Keep
`science_measurement=false` and `images_ready_for_science=false`.

Assessment:
`hsc_debug/s23b_camira_member_outer_profile_relation_assessment_2026-11-07.md`.
Full outputs are under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_outer_profiles_full_2026-11-07_01/`
and
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_outer_icl_relation_2026-11-07_01/`.
The next gates are CAMIRA schema verification, member-preserving
injection/recovery, PSF-wing and sky-transfer calibration, coadd covariance,
and center/mask/control sensitivity. No new data was downloaded.

### S23B raw member-field-8 profile weighting sensitivity (2026-11-07)

The member-centered profiles were re-aggregated using the imported member
file's eighth raw field as an **unverified sensitivity weight**. The runner
reproduces the saved unweighted outer profiles exactly and pairs 5,000
connected sky-block bootstrap draws between weighted and unweighted fits.
No new data or image pixels were used. Although the S23B CAMIRA paper
describes a membership-weighted richness, the public catalog README does not
map the downloaded member file's headerless field 8 to that quantity. Keep it
named raw field 8 until the file-specific schema is verified.

The weighted 50–75 kpc candidate-light medians are 1.03 nJy arcsec⁻² for
HSC-i (62 clusters; q16–q84 −3.87 to 6.36) and 1.99 for HSC-i2 (365 clusters;
−4.43 to 10.89). Adjusted outer-light/ICL slopes at the 5° phase-zero layout
are 0.017 [−0.247, 0.334] for i and −0.013 [−0.118, 0.129] for i2. All
weighted-slope intervals include zero in all six spatial layouts. The paired
i2 weighting shift is positive in only 3/6 layouts; its interval crosses zero
in the other three. This does not support an ICL–satellite detection or
validate field 8. Keep `science_measurement=false` and
`images_ready_for_science=false`.

Assessment:
`hsc_debug/s23b_camira_member_outer_raw_field8_sensitivity_assessment_2026-11-07.md`.
Runner: `hsc_debug/analyze_s23b_member_outer_raw_field8_sensitivity.py`.
Output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_outer_field8_sensitivity_2026-11-07_01/`
(summary SHA-256 `7156397a478932a9d1fb133b88a758a814bfc0543708719cd36f2a21c1b14875`).

### S23B CAMIRA member-position mask audit (2026-11-07)

The saved-mask WCS reproduces member positions within 1.08e-7 pixel, but the
combined Kron masks flag 98.4–99.1% of member centers. The mask products lack
segmentation labels, so they cannot tell a candidate's own source mask from a
neighbor's. The BASE_MASK center flag fractions are 48.9% in HSC-i and 41.4%
in HSC-i2; median whole-image BASE_MASK area fractions are 55.2% and 46.5%,
respectively, and are only a coarse non-spatially-matched comparison. Four
images have no mask product due to too few valid pixels; 255 member positions
in them remain inside the geometric footprint but unsampled.

The combined source masks cannot be used unchanged for member-centered outer
light. Next, retain SEP segmentation labels and detections, match each CAMIRA
position to its own detection, preserve that candidate in a per-member mask,
and assess annular support against matched control positions. Keep
science_measurement=false. Full assessment:
`hsc_debug/s23b_icl_satellite_member_mask_audit_assessment_2026-11-07.md`.

### Weak-lensing and S23B member-profile sample identity

The existing 800-lens S19A weak-lensing sample has no positional match within
60 arcsec to the 1,043 S23B CAMIRA targets. The closest pair is 364.846
arcsec apart. This remains true for the 435 lenses with common cluster/random
shape support and the 427 S23B targets represented in the saved member outer
profiles. Numeric IDs collide for 66 pairs, including 29 profile-target IDs,
but no collision is positionally valid; numeric ID equality must not be used
for this join.

Therefore the current S19A and S23B profile products cannot support a
same-cluster three-way weak-lensing–ICL–satellite relation. A local count-only
plan now covers the 427 S23B member-profile clusters in 18 HSC PDR3 CAS
batches, using 3 proper-Mpc cones and returning only one integer count per
cluster (122,198 SQL bytes). It has not been submitted; its SHA-256 is
`b713eebded489c280fd9ce1b5840bc8d7d3e42d029035f8dbd414b72fb31e984`.
Submission and any later source-row export each require separate explicit
approval. After counts, validate shape support and matched random centers.
Keep `science_measurement=false`.

Cross-match assessment:
`hsc_debug/wl_s23b_camira_crossmatch_assessment_2026-11-07.md`.

### CAMIRA member raw-field schema consistency

The approved member file maps 637,146 rows uniquely to the 12,162 S23B Wide
star-mask clusters. Its raw field 8 lies entirely in [0.001000069, 1.0],
matching the published CAMIRA member-catalog threshold `w_mem >= 10^-3`; this
is strong evidence that field 8 is probably `w_mem`. The catalog README does
not document the member-file column order, so preserve it as `field8_raw`
until the specific layout is confirmed. The 50–75 kpc profile weighting
remains an exploratory sensitivity, not a science detection. See
`hsc_debug/s23b_camira_member_raw_field_schema_assessment_2026-11-07.md`.

### S23B ICL–satellite raw-field-8 weighted-density sensitivity (2026-11-07)

The 1,043-target member-position table was used to repeat the direct
ICL–candidate-density relation with sums of raw member field 8 in the
50–300/300–800 proper-kpc annuli. Public CAMIRA literature makes field 8
strongly consistent with `w_mem`, but no source located states the exact
column order of this headerless S23B file. It remains `field8_raw`; the fit is
a weighting sensitivity, not a verified membership-probability measurement.

The same supported sample and paired 5,000-draw connected-block bootstrap were
used for weighted and unweighted fits. In the primary 5-degree phase-zero
layout, HSC-i slopes are 5.13 [−0.80, 13.05] weighted and 5.72 [−0.06, 13.63]
unweighted; the paired shift is −0.59 [−3.26, 1.60]. HSC-i2 slopes are
5.01 [2.44, 7.62] and 4.24 [0.99, 7.80], with paired shift +0.77
[−1.05, 2.46]. The paired shift includes zero in all six layouts for each
filter. At a 10% common-pixel floor, the primary i2 weighted interval crosses
zero (−0.21 to 6.57), and only two of six layouts exclude zero; i remains
unresolved. The full-sample i2 result is therefore not robust to support
selection. Keep `science_measurement=false`.

Assessment:
`hsc_debug/s23b_icl_satellite_field8_weighted_density_assessment_2026-11-07.md`.
Runner:
`hsc_debug/analyze_s23b_icl_satellite_field8_weighted_density.py`.
The run used the already received approved member file; no CAS request,
download, or image-pixel read occurred. The separate 427-target S19A
count-only plan remains unsubmitted pending its own approval.

### S23B PSF-candidate field pilot (2026-11-07)

A CPU-only pilot on 50 already received target images (25 per actual coadd
filter) found 492 provisional compact-source candidates across 49 usable
frames; 44/49 frames contain at least four candidates. Input image hashes
matched the verified receipt. One frame, `S23B_10647`, had too few valid
pixels. The separate local S19A Gaia subset had no positional overlap. The
exact 492-position Gaia DR3 query was later approved and completed: 77
positions matched, and 48 passed a conservative diagnostic star-evidence
screen. From already received target images, 36 wide stamps were available;
four frame-weighted profile stacks show substantial sky-annulus sensitivity
and sparse outer support. They do not establish a control-matched PSF or a
stable extended-wing correction. Member profiles remain exploratory and
`science_measurement=false`.

Assessment:
`hsc_debug/s23b_detection_psf_pilot_assessment_2026-11-07.md`.
Follow-up assessment:
`hsc_debug/s23b_gaia_empirical_psf_wing_assessment_2026-11-07.md`.

#### 801-pixel PSF-wing extension

Using only the approved Gaia-supported target candidates and already received
images, a larger local extraction produced 22 stamps (8 `i`, 14 `i2`) across
13 frames. The >=500-proper-kpc subset contains 4 `i` stars/frames and 10
`i2` stars across 8 frames. Profiles reach 33.6 arcsec and use two background
annuli at 43.7–50.4 and 57.1–65.5 arcsec. At radii above 16 arcsec the `i2`
stack has only 1–2 frames per bin and is consistent with zero; the outer `i`
profile is supported by a single frame. This does not validate the apparent
positive `i2` wing in the smaller 401-pixel stack, since the wider stamp sample
is smaller and not paired. The PSF wing remains inseparable from sky-annulus
transfer, and `science_measurement=false` remains in force.

Detailed results:
`hsc_debug/s23b_gaia_empirical_psf_wing_assessment_2026-11-07.md`.

### 2026-11-07 matched PSF expansion and Gaia completion

The 49-pair pilot could not estimate a paired PSF-wing difference because no
radial bin had more than one valid matched frame pair. The local expansion
selected 100 new matched pairs in each actual HSC-I subfilter, `i` and `i2`,
across five target-seeing strata. The selection excluded the 49 pilot pairs and
17 pairs with different target and control subfilters. The target-side screen
used 7,505 CAMIRA member coordinates inside selected images.

Local CPU preparation succeeded for 199/200 target frames and 200/200
controls. `S23B_11786` had fewer than 100 valid pixels, and three prepared
target frames yielded no compact candidates. The morphology-only candidate
list contains 2,069 target and 2,230 control positions. Counts by filter are
1,016/1,000 target/control in `i` and 1,053/1,230 in `i2`.

The user approved the exact 4,299-position Gaia DR3 plan with SHA-256
`c00eda5f896f8160022a7e4198b3f8ce69803f0e30adc8f560ce04c8beba0e7d`. The
18 hash-pinned batches returned 706 candidate-source rows in a 138,892-byte
aggregate CSV. The response SHA-256 is
`151f7c3f227ab7a16b1a6015563afd46fd53ccf4067d44d0264e27be75f665fb`.
The local audit reproduced every association from the saved Gaia cone-union
responses and verified each batch query, response, and match-table hash. The
706 matched candidates each have one unique Gaia source. No candidate cone
contains multiple Gaia sources, and no Gaia source is reused for another HSC
candidate. The maximum separation is 0.2876 arcsec within the approved
2-arcsec radius.

A conservative diagnostic screen found 448 strict candidates, with 225
in the target sample and 223 in the control sample. The criteria require one
Gaia source inside 2 arcsec, RUWE below 1.4, no duplicate, QSO, or galaxy
flag, significant parallax or proper motion, and combined star probability
at least 0.5. The criteria do not form a calibrated star or PSF classifier.
The remaining counts are 3,593 Gaia non-matches, 255 matches without clean
astrometry, and 3 matches with clean astrometric evidence but without the
combined-star threshold. A Gaia non-match cannot classify an HSC candidate as non-stellar.

The local CPU stamp pass made 168 target and 166 control stamps, each
401×401 pixels, from existing approved HSC images. The target pass blocked
57 edge-crossing candidates. The control pass blocked 56 edge-crossing
candidates and one known variable. Stamps retain raw image values and record
the fractional SExtractor center. The radial measurement did not shift the
stamp pixels. The profile calculation masks detected neighbors and defects,
subtracts one of two local sky estimates from each annulus, and divides by
masked core flux. The sky annuli span 120–150 and 160–190 pixels. Radial bins
extend to 96 pixels, or 16.1 arcsec.

Eight frame-weighted profiles cover target and control samples, both filters,
and all-star or center-distance-at-least-500-kpc subsets. Each profile first
averages stars within a frame. The uncertainty calculation resamples matched
frame pairs 5,000 times. The `i` filter has no radial bin with at least two
valid matched pairs. The `i2` all-star profile has at least two pairs in 43
of 44 bins, with a maximum of 14 pairs. The `i2` outer-500-kpc profile has at
least two pairs in 30 of 44 bins, with a maximum of eight pairs. Control
center distance uses a random sky position rather than a cluster center.

The `i2` all-star frame-weighted means are negative from 32 to 96 pixels on
both target and control sides under both sky annuli. The paired control-minus-
target difference in the 24–32-pixel bin changes from −6.7
[−13.3, 0.0] × 10⁻⁶ with the 120–150-pixel sky annulus to +1.7
[−8.0, +12.0] × 10⁻⁶ with the 160–190-pixel annulus. The quoted intervals
span the 16th to 84th percentiles of the frame-pair bootstrap. The sign
reversal demonstrates strong sky-annulus sensitivity. Two `i2` outer-500-kpc
bins have nominal 95% intervals outside zero under the farther sky annulus,
but each bin has only three matched pairs among 44 profile bins. The small
support and annulus dependence do not establish a target-control PSF-wing
difference.

The expanded sample supports a larger empirical PSF diagnostic, but no
calibrated extended-wing correction. The saved profiles must not correct the
ICL stack. The effective HSC coadd epoch remains unknown, and the Gaia
positions were not propagated to a presumed epoch. The approved 67,071-row
CCD metadata export covers 800 S19A lens centers and 798 controls. A prior
coordinate audit found no match within 60 arcsec between the S19A lens sample
and the S23B CAMIRA centers. The saved MJD rows cannot assign coadd epochs to
the 200 S23B matched pairs. Any new metadata request requires separate exact
approval. Keep `science_measurement=false` and `psf_approved=false`.

Outputs

- Gaia audit runner `hsc_debug/audit_s23b_psf_expansion_gaia_matches.py`
- Gaia audit `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_expansion_2026-11-07_01/gaia_audit_2026-11-07_01/`
- Target and control stamps under `s23b_psf_expansion_2026-11-07_01/{target,control}_wing_stamps_2026-11-07_01/`
- Eight filter/subset profile stacks under `s23b_psf_expansion_2026-11-07_01/{target,control}_wing_stack_{i,i2}_{all,outer500}_2026-11-07_01/`
- Paired comparison `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_expansion_2026-11-07_01/target_control_profile_comparison_2026-11-07_01/`
- Detailed assessment `hsc_debug/s23b_gaia_empirical_psf_wing_assessment_2026-11-07.md`

### 2026-11-07 matched-control PSF-candidate pilot

The already downloaded matched random-sky control images have now passed the
same local morphology-candidate preparation as the target subset. Across 49
matched target/control frames, the controls contain 458 provisional compact
candidates (180 `i`, 278 `i2`) in 48 fields; the targets contain 492 (229 `i`,
263 `i2`). Candidate FWHM-moment medians are similar by filter, but the screen
does not classify stars. Candidate positions have no target/control match
within 2 arcsec (nearest separation 58.12 arcsec), and there are no duplicate
control cones at 2 arcsec.

The separate Gaia DR3 positional cross-match for the 458 control candidates
was subsequently approved and completed. The exact plan used two batches of
250 and 208 positions, 2-arcsec cones, and 16 Gaia source fields; its SHA-256
is `cbf9bb94c73738b06cd5f616ba1c190f4f0e11ec71a5d914333b4b291d0e7daf`.
It returned 58 match rows (rough estimate 72; response SHA-256
`37a8d515b415990a3444bd4654856b0fe10e2d2c2aa2e10c0c517d699eccb4b7`). The
audit classified 38 as strict diagnostic candidates, 20 as matches without
clean astrometry, and 400 as Gaia non-matches. Together with the target-side
audit, this gives 48 strict target candidates and 38 controls; the thresholds
are not a calibrated PSF-star classifier.

The control candidates produced 26 usable 401×401 stamps (11 `i`, 15 `i2`)
from 38 Gaia diagnostic candidates; 12 were blocked by detector edges. Four
control frame stacks used the same radial bins, two local background annuli,
and 5,000 frame-bootstrap draws as the target pilot. In the paired profile
comparison, across both filters, all radii and the center-distance-≥500-kpc
subset, no radial bin had two valid matched frame pairs; the maximum support
was one pair. A paired PSF-wing difference or uncertainty interval is
therefore not estimable. The unpaired stacks remain diagnostic only. The
control center-distance cut is measured from its random sky coordinate, not a
cluster center. Keep `science_measurement=false` and `psf_approved=false`.

Outputs are under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_control_wing_stamps_2026-11-07_01/`,
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_control_wing_stack_{i_all,i_outer500,i2_all,i2_outer500}_2026-11-07_01/`,
and
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_target_control_profile_comparison_2026-11-07_03/`.

Detailed checkpoint:
`hsc_debug/s23b_gaia_empirical_psf_wing_assessment_2026-11-07.md`.

## 2026-11-07 S23B coadd residual-covariance diagnostic

The 1,022-pair exploratory target/control stack (241 `i`, 781 `i2`) was
checked using 64×64-pixel outer-field patches, saved K=2.5 masks, and the HSC
VARIANCE planes. A finite 9×9 lag-covariance sum gives median approximate
variance factors of 1.586/1.580 (target/control) in `i` and 1.572/1.555 in
`i2`; median standard-deviation factors are about 1.25. Zero-lag residual
variance is about 1.14 times the delivered VARIANCE prediction in both
populations. Paired target-minus-control factor medians are only +0.023 (`i`)
and +0.029 (`i2`), with descriptive frame-pair spreads that include zero.
One `i2` target frame lacked enough valid outer patches. These results
indicate a likely shared short-range covariance contribution, but the finite
kernel is not propagated through reprojection or the exact annular estimator
and is not a defensible error-bar correction. Sky, source-Poisson, extended
PSF, and selection uncertainties remain. Keep `science_measurement=false`.
Full method, limitations, and hashes are in
[`s23b_icl_residual_covariance_assessment_2026-11-07.md`](hsc_debug/s23b_icl_residual_covariance_assessment_2026-11-07.md);
output tables and plot are under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_residual_covariance_2026-11-07_01/`.

### Propagation through the saved relative-contrast estimator

The finite short-lag kernel was subsequently propagated through the actual
bilinear metric-grid warp and common-support 50–300, 100–300, or 200–300 kpc
minus 300–800 kpc annular estimator for all 1,022 matched pairs. The
reconstructed means and all six 5-degree block-support counts match their
saved references. Among supported contrasts, the finite-lag conditional
pixel-noise σ is 1.089–1.098 times the delivered-VARIANCE diagonal σ, and is
18–34% of the spatial-block 95% interval half-width. This scale comparison
does not define a complete uncertainty or a correction. The kernel is measured
in outer patches and assumed stationary in cluster annuli; long-range sky,
source-Poisson, PSF, mask-selection, and other systematics remain. Keep the
block intervals as the population-level reference and keep
`science_measurement=false` and `images_ready_for_science=false`.

Assessment:
`hsc_debug/s23b_icl_relative_covariance_propagation_assessment_2026-11-07.md`.

Runner: `hsc_debug/propagate_s23b_icl_relative_contrast_covariance.py`.
Output (outside Git):
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_relative_covariance_propagation_2026-11-07_02/`.
Summary SHA-256:
`9c1ef1740490f7cbc4811a9bd2573d69405f5482d7ce754a09642108507fc94e`.

## 2026-11-07 CAMIRA member outer-profile injection recovery

The source-aware member profile pipeline was replayed and tested with a
seeded-random, circular exponential halo injection around one already-profiled
member in each of 427 target/filter frames (62 `i`, 365 `i2`). Baseline
replays were exact. Across amplitudes 0.5, 1, and 2 nJy arcsec⁻² at 50 kpc,
1,221 of 1,281 frame-amplitude rows had adequate support. Median recovered to
nominal input ratios were 0.838–0.896 in `i` and 0.807–0.819 in `i2`; paired
amplitude differences have cluster-bootstrap intervals including zero.
However, injection changed the segmentation footprint in most frames and the
recovery remains conditional on already-detected members. The halo was not
PSF-convolved, source Poisson noise was omitted, and the rotated controls were
not injected. Treat this as measurement-transfer calibration only; it does
not establish completeness or an ICL–satellite association. Keep
`science_measurement=false` and `images_ready_for_science=false`.

Assessment: `hsc_debug/s23b_camira_member_outer_injection_recovery_assessment_2026-11-07.md`.
Runner: `hsc_debug/run_s23b_member_outer_profile_injection_recovery.py`.
Output (outside Git):
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_outer_injection_recovery_2026-11-07_01/`.
Summary SHA-256:
`e9d4cced59f0b0f8ea3f189989886f28af74a339cc5a9b3421f8bca187defc77`.
No network request, download, commit, or push was made for this analysis.

## 2026-11-07 S23B member-profile position-specific PSF plan

Prepared, but did not submit, an HSC PDR3 `pdr3_wide` centered coadd model-PSF
request for the exact 427 CAMIRA member positions used in injection recovery
(62 `i`, 365 `i2`). None exactly overlaps the three audited existing PSF
banks, and the nearest old position is 69.6869 arcsec away. The estimated
transfer size is 8,608,320 bytes if all files match the largest previously
observed PSF; this is a review cap, and processing must stop if it would be
exceeded. Images, shape rows, and photo-z are excluded. The plan remains
`prepared_not_submitted`, `network_submitted=false`, and
`science_measurement=false`; separate explicit approval is required.

Assessment: `hsc_debug/s23b_member_profile_psf_request_assessment_2026-11-07.md`.
Planner: `hsc_debug/prepare_s23b_member_profile_psf_plan.py`.
Plan (outside Git):
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_profile_psf_plan_2026-11-07_03/`.
Plan SHA-256: `03f974ad24ead32eabd103ef7a6e0d6b52317aa2647023a5f7402003a2e3b25d`.

## 2026-11-07 Existing S19A source-row reuse audit for S23B WL

The complete 800-lens S19A source-row export was audited locally against its
32 batch receipts, byte counts, SHA-256 hashes, pair counts, source-level
field consistency, and original 3 proper-Mpc query cones. Its 3,069,544
lens-source rows reduce to 2,929,303 consistent unique source IDs. Comparing
those rows with the 427 S23B member-profile centers found geometric query-cone
overlap for 43 targets, but existing source rows inside only 21 target
apertures. No aperture is fully covered by one old cone. Estimated old-cone
union area coverage has median 0, 90th percentile 0.00713, and maximum 0.604;
none reaches 95%. The prior export therefore cannot replace target-specific
support counts or a source-row export.

A second local audit added all 9,010,324 rows from the already approved
3,100-center S19A random export, after checking all 124 batch receipts,
schemas, pair uniqueness, and source-field consistency against the lens
export. The combined query-cone union contains source rows in 32/427 target
apertures (21 lens-export, 16 random-export, with 11 random-only); 395 have no
saved rows. Estimated combined geometric coverage has median 0, p90 0.310,
maximum 0.994, and only 2/427 apertures reach 95%. This partial union still
cannot support a common 427-target measurement. Geometry is an upper bound
on reusable catalog coverage, not a depth or completeness assessment.
`science_measurement=false` remains in force.

Combined audit runner:
`hsc_debug/audit_wl_s23b_reuse_existing_s19a_sources.py`.
Output:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_existing_source_union_2026-11-07_01/`.
Summary SHA-256:
`876bed400f93a55c391c89128aa9e9214ed489e664ea0e618ea3df6920f91ed7`.

The already prepared S23B count-only CAS plan remains unsubmitted and needs
its own exact-plan approval. It returns 427 integer source counts in 18
batches; it is a separate data action from the completed Gaia positional
query. A later shape-row export, matched random-center support, calibrated
tomographic n(z), and B-mode selection validation remain separate gates.
Keep science_measurement=false.

Assessment:
hsc_debug/wl_s23b_existing_export_reuse_assessment_2026-11-07.md.
Runner:
hsc_debug/audit_wl_s23b_reuse_existing_s19a_shape_rows.py.
Output (outside Git):
/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_existing_export_overlap_2026-11-07_01/.
Summary SHA-256:
61026a5d82d47303724f576922415c451589f2516f39337adee3abe021c2c54a.
No network query, data transfer, commit, or push was made for this audit.

## 2026-11-07 Integrated S23B ICL–satellite–PSF–WL evidence checkpoint

`hsc_debug/build_s23b_icl_wl_evidence_report.py` consolidates the saved S23B
stack, CAMIRA position and member-profile relations, injection recovery, Gaia
PSF diagnostics, and the existing S19A shape-sample overlap. The report is
hash-pinned to its source summaries and the completed, approved 4,299-position
Gaia plan. Generating it made no network request, read no image pixels, and
downloaded no data.

The 1,022-pair target/control stack has positive 50–300 minus 300–800 kpc
relative contrasts in both filters. The unweighted satellite-density slope is
positive in the primary `i2` fit, but its interval includes zero at a 10%
common-support floor and it is significant in only four of six block layouts.
The alternative field-8 weighting remains sensitivity-only because the
headerless field is not authoritatively identified. The separate adjusted
member-centered outer-light relation has intervals including zero in both
filters. These results do not yet establish an ICL–satellite association.

The existing S19A shape sample has no positional overlap within 60 arcsec with
the S23B targets. Existing lens and random source exports have rows in only
32/427 member-profile apertures, leaving 395 without rows. The prepared
S23B-specific count-only plan remains unsubmitted and requires separate
explicit approval; it returns only one lens identifier and one integer count
per target. Any source-row export or matched-random support query is a separate
data action.

- Report: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_wl_evidence_2026-11-07_01/assessment.md`
- Machine-readable summary SHA-256: `5684f7e788c72e482ec5f6bbff9e205ffae67a124c0f36fac89ee0fd1055d0fe`
- Builder: `hsc_debug/build_s23b_icl_wl_evidence_report.py`
- Pending HSC count-plan SHA-256: `b713eebded489c280fd9ce1b5840bc8d7d3e42d029035f8dbd414b72fb31e984`
- Keep `science_measurement=false`, `images_ready_for_science=false`, and
  `psf_approved=false`.

## 2026-11-07 S23B ICL–satellite restricted wild-cluster inference

A local CPU analysis re-fit the pinned adjusted model for all eight prespecified
filter/support combinations, using 10,000 Webb six-point restricted wild-
cluster bootstrap-t draws and one shared 5-degree connected target/control
block partition (49 blocks over the union). All eight groups had at least six
blocks. Three unadjusted bootstrap p-values were below 0.05; one remained below
0.05 after Holm adjustment (`i2`, both annuli at >=20% common support: slope
4.21 nJy arcsec^-2 dex^-1, raw p=0.0057, Holm p=0.0456). No group passed the
single-step max-T diagnostic under the complete null (smallest p=0.2666).
Support-threshold subsets reuse the same selected images and are not independent
replication; the one Holm result is not treated as robust evidence. The relative
ICL outcome, CAMIRA candidate counts, approximate sky blocks, and unmodeled sky,
PSF, mask, and coadd-covariance terms keep `science_measurement=false`.

- Runner: `hsc_debug/analyze_s23b_icl_satellite_relation_wild.py`
- Assessment: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_satellite_relation_wild_2026-11-07_01/assessment.md`
- Summary SHA-256: `d95f894c6cf54e92d906063fdc920c24e80bdf9ad4194d2ff63dacd8bdd42494`
- Pair-table SHA-256: `98ed9e0022bfa68d5921ead62dce55be7df1cfa169e6a777d66d1e9760f43994`
- No network request, image-pixel read, new data download, commit, or push was made for this analysis.

## 2026-11-09 S23B S19A count-only audit

The exact approved count-only plan (SHA-256
`b713eebded489c280fd9ce1b5840bc8d7d3e42d029035f8dbd414b72fb31e984`) completed
18/18 batches and returned one integer count for each of 427 centers. The audit
verified the batch receipts, result hashes, SQL center IDs, and target-manifest
mapping. It recorded 356,191 summed lens-source candidate pairs; 329/427
centers have zero entries, 82 have at least 1,000, 62 at least 3,000, and 25
at least 5,000. Median count is zero; maximum is 8,959. Descriptive Spearman
rho is 0.060 versus redshift and -0.145 versus richness; no p-values were
computed. A local cross-tab shows that zero counts vary strongly by sky region,
so footprint and common-support selection need resolution before defining a
usable lens sample.

The `weaklensing_hsm_regauss` count is from the released S19A shape-catalog
subset, whose release selection already includes the weak-lensing flag and star
mask. It precedes project-specific source-redshift, tomographic, B-mode, and
radial cuts. These are availability counts, not unique-source totals, final
effective counts, signal-to-noise estimates, or a weak-lensing measurement.
No source rows, coordinate rows, shape values, photo-z, or images were returned;
credentials were not saved. Keep `science_measurement=false`.

- Audit runner: `hsc_debug/audit_s23b_wl_shape_counts.py`
- Assessment: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_shape_count_audit_2026-11-07_01/assessment.md`
- Assessment SHA-256: `f2ce0c5f22c2f6f894386f0c682d08aaf21236d5dd56188cf896b65064e6d672`
- Per-center table SHA-256: `d650564b4870d4af93ce04b9ebd82475a99366aa68c2dfe236e7f070826c0f7f`

This local assessment and plan-preparation step was completed on 2026-11-10.
The 400-control count-only query and any later shape-row submission remain
separate approval gates.

## 2026-11-10 S23B weak-lensing common-support design

The approved S19A count-only audit covers 427 S23B target apertures, with 329
zero-count positions and 82 positions containing at least 1,000 released shape
catalog entries within 3 proper Mpc. The count is a raw availability screen;
it does not include project-specific source-quality, tomographic, B-mode, or
radial cuts and is not a shear measurement.

The provisional common-support candidate contains 80 lenses with 400 unique
HSC random-table control centers (five per lens). The other two targets above
the count screen lack five qualifying saved candidates. Controls satisfy the
audited center-flag and depth/seeing/SFD-EBV matching thresholds, have no
overlap in the previous CAMIRA 2-Mpc screen, and lie more than 6 proper Mpc
from their paired lens. This avoids geometric overlap between paired 3-Mpc
cones but does not guarantee an uncontaminated large-scale environment.
Selected lenses have 340,256 verified raw count pairs (per-lens min/median/max
1,040/4,323/8,959). Only 29 selected apertures contain any rows from the
already-exported S19A lens/random union; none has 95% estimated geometric
coverage.

The exact next count-only plan has 400 centers in 16 batches and returns only
the center ID, paired lens ID, and integer source count. It is validated by the
existing HSC random-count plan checker but remains unsubmitted. Its exact
plan-file SHA-256 is
`fe19ac426272d4005e3ee74a902dfdb52d1fe83ecb2cc5ec901c10b880d28e8f`; submission
requires separate explicit approval. No data transfer is included.

The source-row plan, which is not yet approved for execution, has 80 lens
centers (4 batches) and 400 control centers (16 batches). The expected lens
rows from the approved counts are 340,256; the control total must first be
measured by the count-only query. A prior, different sample supplies only a
planning proxy of 1.05/1.83/2.23 million combined pairs (q16/median/q84), or
0.25/0.55/1.00 GiB uncompressed CSV (low/point/high). Lens and control schemas
contain 22 and 23 columns respectively, including shape/calibration fields,
PSF moments, raw B-mode mask, and HSC Y3 source-z bin. Full object photo-z PDFs
and images are excluded.

Before weak-lensing interpretation, use identical source cuts for lenses and
controls, establish the B-mode field polarity, provide calibrated source n(z)
and lensing efficiency, and evaluate footprint, PSF, control-overlap, and
spatial-covariance support. Keep `science_measurement=false` until those gates
and image/shape quality checks pass. The local planner and complete assessment
are in `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_wl_common_support_plan_2026-11-10_01/`;
runner: `hsc_debug/prepare_s23b_wl_common_support_plan.py`.

### 2026-11-10 approved control-count query and post-count row plan

The separately approved count-only plan
`fe19ac426272d4005e3ee74a902dfdb52d1fe83ecb2cc5ec901c10b880d28e8f` completed
on `lageunha`: all 16 batches and 400 integer-count rows passed center-ID,
count, and receipt-hash validation. The result schema contained only random
center ID, paired lens ID, and integer source count. No shape values, source
coordinates, photo-z, or images were returned; credentials were not saved.

The 400 controls sum to 1,572,453 raw catalog counts (18 zero, 382 positive).
The 80 lenses sum to 340,256, giving 1,912,709 expected lens/control pairs
for a future row export. The per-lens mean of five control counts divided by
the lens count has q16/median/q84 0.690/0.936/1.295; seven lenses have at least
one zero-count control and none has five zeros. Counts precede the project
shape, source-redshift, B-mode, and radial cuts, so they do not establish final
common support or a shear measurement.

The post-count row manifest expects 80 lens centers in four batches and 400
controls in 16 batches. Estimated uncompressed CSV size is 0.46/0.58/0.86 GiB
(low/point/high), using the exact count totals and historical bytes per pair;
these are not observed transfer bytes. The row query remains unsubmitted and
requires separate explicit approval tied to file SHA-256
`7bfb998149f30a50b5ee7a537bb703ed302b7ea269a3296c35204a836e75aabc`.
`science_measurement=false` remains in force.

- Count audit, post-count assessment, and plan: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_wl_common_support_plan_2026-11-10_01/`
- Finalizer: `hsc_debug/finalize_s23b_wl_common_support_plan.py`

### 2026-11-10 approved S23B weak-lensing source-row export

The user-approved row plan (file SHA-256
`7bfb998149f30a50b5ee7a537bb703ed302b7ea269a3296c35204a836e75aabc`) completed
on `lageunha`: 4/4 lens batches and 16/16 matched-control batches. Returned
data contains 340,256 lens-source pairs and 1,572,453 control-source pairs,
1,912,709 center/source pairs total. Actual CSV size is 102,103,795 bytes for
lenses plus 523,062,067 bytes for controls, 625,165,862 bytes (0.582 GiB)
combined. These are pair rows; repeated source IDs across apertures are kept.

All batches matched the approved count-only totals per center and passed
duplicate-pair, schema, finite-number, coordinate-range, binary B-mode-mask,
tomographic-bin, and positive-weight checks. An independent audit recomputed
the SHA-256 and size for every CSV against the job checkpoint. The schemas are
22 lens columns and 23 control columns, including source coordinates, shapes,
calibration terms, magnitude, resolution, PSF moments, raw B-mode mask, and the
HSC Y3 source-z bin. Images and full photo-z PDFs were not included; credentials
were not saved. Results are in
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_wl_common_support_plan_2026-11-10_01/source_rows_post_count_01/`.

This is validated raw released-catalog support, not the final selected source
sample or a lensing measurement. Next, apply identical documented source cuts
to lenses and controls, report aperture-level retained counts and effective
weights, and verify the raw B-mode-mask semantics. Calibrated source n(z),
lensing efficiency, PSF, spatial covariance, and injection/recovery remain
interpretation gates. Keep `science_measurement=false`.

## 2026-11-10 S23B weak-lensing source-row QA

The approved local-only audit re-read the 20 received CSVs (1,912,709
center/source pairs), verified each hash and count-only row total, and applied
the same documented HSC cuts to 80 lens apertures and 400 matched controls.
After shape, magnitude, resolution, z-bin, B-mode, positive fiducial n(z)
efficiency, 0.1–3 proper-Mpc, and catalog-input-validity gates, 239,565 lens
pairs and 1,105,779 control pairs remain. Pair rows include repeated sources
across apertures and are not globally unique objects.

Support is adequate for preparing a calibrated estimator, but not uniform in
every paired set: two lenses have zero retained rows; 24 control centers have
zero rows, affecting nine of the 80 sets. None of the 80 sets has all five
controls at zero. Among lens sets with nonzero lens support, the
q16/median/q84 control-to-lens ratios are
0.690/0.925/1.311 for count, 0.686/0.929/1.311 for HSC weight sum, and
0.690/0.956/1.387 for the `weight*(D_l beta)^2` geometry proxy. The latter is
not the final estimator weight. Shape and weight medians are similar between
roles, and the raw B-mode mask matches the published GAMA09H rectangle exactly,
confirming stored 0 as masked in this export.

No shear profile was computed. The four-bin SACC n(z) is fiducial; its
uncertainty, per-object calibration/responsivity, selection-bias response,
PSF spatial matching, control contamination, repeated-source/field covariance,
and injection/recovery still need explicit treatment. Keep
`science_measurement=false`. Before a profile, prepare and review a
pre-specified calibrated, random-subtracted estimator and covariance plan,
including treatment of zero-support matched sets.

- QA runner: `analysis/icl_wgl_membership/hsc_debug/audit_s23b_wl_common_shape_rows.py`
- Results (outside Git): `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_wl_common_support_plan_2026-11-10_01/quality_audit_2026-11-10_01/`
- Summary SHA-256: `9ef472c1ae0888e187fce05f3449b1614d1da7c99e7a400403b92f8aa7e2dc54`

## 2026-11-11 radial/tomographic support and zero-control sensitivity

Using only the 20 previously approved S23B source-row CSVs and the existing HSC
Y3 SACC n(z), the local follow-up audit resolved support across seven proper-
radius bins and source-z bins 1–4. It added no CAS query, download, image read,
or shear-profile calculation. The 13,440 per-center cells reconcile exactly
with the earlier radial totals: 239,565 lens pairs and 1,105,779 control pairs
after common source and input-validity cuts. Pair counts are center/source
pairs; the same source can recur in multiple apertures.

Support is weakest at 0.1–0.2 proper Mpc: lens apertures with support by source
bin 1/2/3/4 are 65/63/59/52 of 80, while control coverage is 285/303/293/262
of 400. At 2–3 Mpc, coverage is 78/78/78/78 lenses and 376/376/376/376
controls. The detailed radial × z-bin pair counts, shape-weight sums, effective
counts, and per-center cells are in the external audit outputs below.

Among the 78 lenses with nonzero retained lens support, a requirement for at
least 3 positive-support controls retains 74 sets (92.5%), 367 positive control
centers, and 233,064 lens pairs. Requiring 4/5 gives the same sample. Requiring
5/5 retains 71 sets (88.8%), 355 controls, and 224,375 lens pairs. The median
positive-control/lens retained-count ratio is 0.940 for 3/5 or 4/5 and 0.950
for 5/5. These are support sensitivities only; none is designated as the final
selection or estimator policy.

No shear profile or mass fit was calculated. The HSC Y3 n(z) remains fiducial,
and its posterior uncertainty, shear calibration response, PSF matching,
control contamination, repeated-source/field covariance, and injection/recovery
remain unresolved. Keep `science_measurement=false`.

- Runner: `analysis/icl_wgl_membership/hsc_debug/audit_s23b_wl_common_shape_rows.py`
- Results (outside Git): `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_wl_common_support_plan_2026-11-10_01/quality_audit_2026-11-11_01/`
- Summary SHA-256: `eb73aa0365813d880d3edf60b12e4039f2575f8a97eab8868fc918121a34d445`
- Assessment SHA-256: `35369bbc74fa576570860af09d95d29c33c32cf8b1dfb96c30bb5949d5ddcf95`
- Per-center radial/tomographic TSV SHA-256: `f8d604e91c43705cb24152c7e98553efd8838933034a8657bcc6030c2f8c11bed`
- Aggregate radial/tomographic TSV SHA-256: `ae84048ee707ffe531ecff80ede270544a7475846079db7520ee59bd62706911`
- Zero-support policy TSV SHA-256: `ba89948681cb0560408030a4df66a232a2989cd5415c34ca1cb99c59d8bc26da`

### 2026-11-11 S23B matched-sample random-subtracted 약중력 프로파일

승인된 80개 S23B 렌즈와 400개 대조 중심의 HSC S19A shape 행, 고정된 HSC Y3
SACC `n(z)`만 사용해 Astrafek 보정 추정기로 탐색 프로파일을 계산했다. 4/5
양의 지지 대조 기준에서 74개 렌즈 세트와 367개 대조 중심이 남았다. 주 분석은
이전 프로토콜대로 source bin 2–4를 사용했다. 무작위 중심은 양의 지지 표본을
pooled stack했으므로 대조가 5개인 세트가 4개인 세트보다 한 aperture 더 기여한다.
렌즈 세트별 동일 가중 민감도 비교가 필요하다.

5도 equal-area RA/sin(Dec) paired endpoint jackknife에서 tangential
`DeltaSigma`는 0.1–0.2 proper Mpc의 `23.845 ± 3.432`에서 2–3 Mpc의
`0.762 ± 0.400`으로 감소했다 (단위 `10^13 Msun Mpc^-2`). 조건부 zero-profile
진단은 `chi2=221.85`, rank 7, 명목 `p=2.67e-44`이나 검출 확률이 아니다. 5/5
대조 민감도에서 cross-null 명목 `p=0.0341`, 2.5도 민감도에서는 `p=0.3182`로
달라져 null 결과가 대조 기준과 block 폭에 안정적이지 않다. source bin 1을
포함하면 가장 안쪽 tangential 값이 22.565로 바뀐다.

따라서 이 결과는 80-lens pilot이며 질량 또는 ICL–위성 관계 검출로 해석하지
않는다. fiducial `n(z)`/photo-z 희석 오차, contiguous HSC 영역 공분산, random
aperture의 depth/mask/PSF/sky 일치, pooled 대조 가중, reduced shear,
magnification, cluster-member dilution/boost, selection-response 오차가 남았다.
`science_measurement=false`를 유지한다. 우선 기존 자료만으로 렌즈 세트별 동일
가중 random stack과 pooled 기준의 차이를 평가하고, 보유 중인 S19A 전용 `n(z)`
DESI 곡선/범위를 fiducial SACC와 비교한다 (공동 posterior draw는 없어 완전한
불확도 전파는 불가능하다). 이후 footprint 기반 공분산과 영상/대조 일치 관문을
통과한 뒤에만 질량 모델이나 ICL–위성 연결 분석을 고려한다.

- 실행기: `analysis/icl_wgl_membership/hsc_debug/analyze_s23b_wl_random_subtracted_profile.py`
- 평가서: `analysis/icl_wgl_membership/hsc_debug/wl_s23b_random_subtracted_profile_assessment_2026-11-11.md`
- 산출물 (Git 외): `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_random_subtracted_profile_2026-11-11_01/`
- summary SHA-256: `2c899b0090dd94a34f1fd6b0cb111e3d1ada696da2178f98e9fd5b659e699edf`
- 평가서 SHA-256: `d152c22f21a244e9848f453523c0a346120e30f4897b9a427d65776aa804ba2a`
- TSV SHA-256: `3dd1ef2b7f4ac7aa33d6b046ef9b528ae6c8818d1f3992037b6b70aa6e5548b2`
- 그림 SHA-256: `4c3593236d5a46dd504d0f54b3f2396e0a925460c27f3f879bc6fef7bef92ce9`
