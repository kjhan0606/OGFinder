# S23B HSC shape-count 공간·release 선택 감사 (2026-11-11)

## 범위와 요약

앞서 승인된 S16A count-only 결과, 기존 S19A 427중심 count-only 결과, 고정된
S23B 표적 명단을 로컬에서 결합했다. 분석은 1,043개 은하단 중심의 하늘
위치, 적색편이, CAMIRA 풍부도, 두 release의 3 proper-Mpc 원뿔 count만
사용했다. 고정 S19A 52개 렌즈와 ICL 필터별 ID도 해시가 고정된 기존
산출물에서 읽었다. 새 CAS 조회나 네트워크 전송은 없었다.

| 집단 | 중심 수 | S16A count > 0 | 0 count | 양수 count 중앙값 |
|---|---:|---:|---:|---:|
| 전체 S23B 표적 | 1,043 | 177 (17.0%) | 866 (83.0%) | — |
| 고정 S19A 렌즈 교집합 | 52 | 21 (40.4%) | 31 (59.6%) | 4,708 |
| ICL i 필터인 고정 렌즈 | 25 | 19 (76.0%) | 6 | 4,708 |
| ICL i2 필터인 고정 렌즈 | 27 | 2 (7.4%) | 25 | 3,608 |

이미 S19A count를 조회했던 동일 427개 중심에서는 다음과 같이 겹친다.

| 동일 427중심 count 비교 | 중심 수 |
|---|---:|
| 두 release 모두 0보다 큼 | 32 |
| S16A만 0보다 큼 | 3 |
| S19A만 0보다 큼 | 66 |
| 두 release 모두 0 | 326 |
| 두 release 모두 count ≥1,000 | 27 |
| S16A만 count ≥1,000 | 3 |
| S19A만 count ≥1,000 | 55 |

기존 S19A 427개 count 표에서는 98곳이 양수이고 82곳이 count ≥1,000이다.
동일 중심에서 S16A는 양수 35곳, count ≥1,000 30곳이다. 따라서 이
427중심 범위에서 S16A가 더하는 raw-count 후보는 3곳이다. 전체 1,043개
S16A 표적 중 S19A count 감사에 포함되지 않은 616곳에서는 S16A 양수
142곳, count ≥1,000 126곳이다. 이 616곳의 S19A count는 조회되지 않았으므로
S16A 전용 영역으로 간주하지 않는다.

126개 고-count 후보 중 116곳에는 기존 ICL 표적–대조 영상 쌍이 있고, 나머지
10곳에는 이 산출물의 ICL 대조 쌍이 없다. 기존 대조 위치는 116개 모두
고유하며, 필터는 `i` 98곳과 `i2` 18곳이다. 표적–대조 SFD E(B−V) 차이의
절댓값은 중앙값 0.000713, 최댓값 0.003329다. 이 대조 좌표들은 ICL 영상
검토용으로 이미 보유한 것이며, 약중력렌즈용 5개 random-center 집합과
동일하다고 가정하지 않는다. 대조 위치의 S16A source count도 아직 조회하지
않았다.

전체 양수 count는 177개 렌즈 원뿔에서 중앙값 4,622개(범위 4–10,138)다.
765,197은 렌즈 원뿔별 count 합이므로 같은 shape object가 여러 원뿔에
들어갈 수 있다. 전체 표본의 고유 원천 수가 아니다.

## 공간 및 군집 속성에 따른 선택

양수 count 위치는 전체 S23B 표적의 하늘 위치에 무작위로 흩어져 있지 않고
몇 개의 조밀한 RA/Dec 영역에 모여 있다. 0 count가 S16A catalog object가
해당 원뿔 안에서 없다는 뜻은 되지만, count만으로 survey footprint, 개별
마스크, 관측 깊이 중 어느 원인이 지배적인지 분리할 수 없다.

양수 count 집단과 0 count 집단의 CAMIRA 풍부도 중앙값은 각각 17.46과
18.84다. 적색편이 중앙값은 0.4455와 0.4340이다. 전체 1,043개 중심을
풍부도 순위 5분위로 나누면 S16A 양수 비율은 Q1–Q5에서 20.6%, 23.4%,
17.7%, 14.4%, 8.7%로 감소한다. 적색편이 순위 5분위의 양수 비율은
18.7%, 12.9%, 13.9%, 15.9%, 23.6%다. Spearman 순위상관은 count–richness
`rho=-0.123`, count–redshift `rho=0.028`이며 p-value는 계산하거나
사용하지 않았다. 공간상관과 survey footprint를 통제하지 않은 기술
통계이므로 물리적 상관으로 해석하지 않는다.

고정 52개 ICL/WL 렌즈에서는 i와 i2 표본의 S16A 겹침이 특히 다르다.
i2에서 양수인 두 중심은 RA 211.99–212.23°, Dec 0.94–1.16° 안에 모여
있다. 이 표본에서 필터별 차이를 천체물리 효과로 해석할 수 없으며,
하늘 위치 및 영상 표본 선택과 혼동돼 있다. 같은 427개 중심의 release
비교에서도 S16A-only 양수는 3곳뿐이고 S19A-only 양수는 66곳이다. S16A와
S19A의 count를 합산해 coverage 증가로 세거나 shape 행을 합치는 것은
타당하지 않다.

![S16A aperture count와 S23B 표적의 위치·특성 비교](/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_s16a_shape_coverage_selection_audit_2026-11-11_04/s16a_count_coverage_selection.png)

## 연구상 의미와 다음 관문

이 결과는 1,043개 전체를 대상으로 한 S16A 약중력렌즈 스태킹을 뒷받침하지
않는다. 양수인 177개는 공간적으로 선택됐고, count는 source-redshift,
shape-quality, shear calibration, mask, radial-bin, random-control 컷 전의
catalog availability 상한이다. 특히 기존 52개 ICL/WL 표본 중 S16A 양수
count가 있는 곳은 21개뿐이며, i/i2 필터별 지지율이 76.0% 대 7.4%로
불균형하다.

116개 대조 좌표의 S16A count-only 사전 조회 계획을 준비했다. S16A PDR1의
`s16a_wide` 세 테이블에서 각 control 중심과 짝지은 은하단 redshift 기준
3 proper-Mpc 원뿔 안의 distinct shape-object 개수만 반환한다. 5개 SQL batch
(각 최대 25중심), 116개 결과 행, SQL 36,139 bytes, 결과 CSV는 보수적 상한
4,100 bytes다. 계획 상태는 `prepared_not_submitted`; 아직 CAS에 제출하지
않았다. 정확한 계획 SHA-256은
`832f2de18f5d512396cc84673e458932635c1f766ea91cdc04c7b47125738dfa`다.

이 1:1 ICL 대조는 기존 weak-lensing 5 random centers/lens 절차를 대체하지
않으며, 이번 조회도 단순 count-only 진단이다. source-row 다운로드나
스태킹을 승인하지 않는다.
1,000 count 기준 역시 기존 S19A 후보 선별 기준과의 비교용일 뿐, 유효
shape 수 기준은 아니다. 이후 약중력렌즈용 대조군·shape-row 반출은 각각
정확한 범위와 용량을 제시해 별도 승인을 받아야 한다. S16A 행은 S19A와
바로 합칠 수 없으며, release/calibration, source selection, random-control
및 image-support 검증이 별도 관문이다.
이 감사에서는 weak-lensing stack이나 과학 측정을 수행하지 않았고
`science_measurement=false`를 유지한다.

## 재현 정보

- 실행기: [audit_s23b_s16a_shape_count_selection.py](audit_s23b_s16a_shape_count_selection.py)
- 출력 디렉터리 (서버, Git 밖): `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_s16a_shape_coverage_selection_audit_2026-11-11_04/`
- 요약 JSON SHA-256: `c08602e1b45492a92389dabd074e5f854050e8e14c2cd55b4683b5ae6fc0b398`
- 중심별 coverage TSV SHA-256: `a7ff068f8d4ef81352598865ddc922938bcf2d5b4ba9feb8af0ea02b8a58a667`
- S16A/S19A 427중심 비교 TSV SHA-256: `8d67a18316e1823a06d5bdb44edc4fc2dd3db3a93130cc85d2c67e8499a82fdd`
- 고-count 후보의 저장된 ICL 대조 위치 TSV SHA-256: `81c702166871e75570918f770ffad44702cebe9a29920c80d3d73d2263c21706`
- 다음 단계 count-only 계획 (미제출): `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_s16a_icl_control_count_plan_2026-11-11_01/`
- 다음 단계 계획 SHA-256: `832f2de18f5d512396cc84673e458932635c1f766ea91cdc04c7b47125738dfa`
- 계획 생성기: [prepare_s23b_s16a_icl_control_count_plan.py](prepare_s23b_s16a_icl_control_count_plan.py)
- 계획 생성기 SHA-256: `5364f01f606a56354ac89a505e971a9a2d3893822e9205646698a3a37fdecb46`
- quintile TSV SHA-256: `0e0c1f9cdaf178696103232c02d1fa4e7ccf16ac376f3c565eb24a67a4d1a755`
- 고정 52 필터별 TSV SHA-256: `3ed200d6544cfcdde874e421158e3753b8876ae528225cc546203a6451db4c7f`
- 그림 SHA-256: `37f8ba9519044436a218cd120cb5edccc6c39c65c205e7a12b9840e98b2b2683`
- 실행기 SHA-256: `5bd15bb4c18392f6c9f666edde84389945f69511d4b9eb16d8cbc82ea3b2f1a6`
- 입력 count 결과 SHA-256: `0945329119d1100b783211fc128cd0a329b449f5acf6648313a1484ac55a7b12`
- 입력 S19A count 결과 SHA-256: `d650564b4870d4af93ce04b9ebd82475a99366aa68c2dfe236e7f070826c0f7f`
- 입력 S19A count summary SHA-256: `748fea398f144165a25210071a5f6ab1866792df3be231ccfdb75282829935b8`
- 입력 target manifest SHA-256: `e487582ba40cd3e58da67704242db4d7219e381fb08f4b33fed47ce4ff268b44`
- 이 조회 계획 생성 중 새 CAS 요청·source-row 다운로드·영상 읽기 없음; `science_measurement=false`.
