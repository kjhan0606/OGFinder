# 보충 HSC 무작위 중심 후보 고리 계획 — 2026-10-20

## 동기

기존 사용자 승인 800-lens 좌표 내보내기의 첫 두 배치 부분 기하 파일럿은 50개 렌즈 중 19개에서 기하 필터 뒤 후보가 0개였고, 전체 33,514행 가운데 2,129행만 남았다. 첫 50개는 대표 표본이 아니므로 전체 후보 부족을 확정하지는 않지만, 기본 30′ 탐색 고리가 일부 렌즈에 충분하지 않을 가능성을 드러냈다.

## 준비한 보충 범위

기존 조회의 2×3 Mpc–30′ 범위와 겹치지 않도록 고정된 30′–60′(1,800–3,600 arcsec) annulus를 정의했다. 이는 최종 random 중심을 선택하는 계획이 아니라, 각 800개 표적 중심 주변에서 이 추가 고리에 들어오는 HSC PDR3 `pdr3_wide.random`의 적격 후보 수를 렌즈마다 하나씩 세는 count-only 사전 계획이다. 기존과 같은 `adjust_density < 0.01`, `isprimary`, `i_pixelflags_edge`, `i_pixelflags_bright_objectcenter` 조건을 사용한다. 출력은 `lens_id`, `n_random_candidates` 두 열뿐이며 좌표, source 행, shape, photo-z, 영상은 반환하도록 구성하지 않았다.

Astrafek은 count query와 후속 coordinate query가 같은 고리 경계를 사용하도록 `prepare_wl_random_control_counts.py`와 `submit_wl_random_control_coordinates.py`를 확장했다. 기본 2×3 Mpc–30′ 계획의 동작은 그대로 두었다. 보충 계획의 `random_control_count_plan.json` SHA-256은 `546da14b240ed7168f0115225eff01413a0b88c02d249e6f73892f160bd78cdb`이며 결과 디렉터리는 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_point_count_expansion_plan_2026-10-20_01/`다. 계획은 32개 SQL 배치로 구성되고 `submit_wl_random_control_counts.py --check` 검증을 통과했다.

## 승인과 과학적 한계

계획 생성·검증은 로컬에서만 했고 HSC CAS 요청은 제출하지 않았다. 사용자 승인 전에는 실행하지 않는다. 현재 기존 2×3 Mpc–30′ coordinate export는 lageunha PID `1049265`로 실행 중이며 2/32 배치(33,514행), CAS 작업 `1062312` 제출 상태다. 이 기존 export가 완료되기 전 보충 조회를 병렬 제출하지 않는다.

보충 count 결과가 있어도 이는 후보 개수만 나타내며 random 중심 선택을 승인하지 않는다. 이후 좌표 반출은 별도 확인하고, 전체 CAMIRA 조리개 중첩·좌표 중복·image footprint·마스크·깊이·PSF·sky·소광·중심 간 중첩·공간 공분산을 검증해야 한다. Random shape 행 반출과 최종 차감 스택은 각각 정확한 범위와 전송량을 제시해 별도 승인을 받아야 한다. 이 계획은 `science_measurement=false`이며 어떤 weak-lensing 결과도 추가하지 않는다.

단위검사 10개 및 전체 `hsc_debug` 회귀검사 288개와 6개 하위 검사가 통과했다. 새 count-only plan은 계획 검사기로 800개 렌즈·32개 batch·두 출력 열임을 재확인했다.
