# 스케일 ICL 스택, 멤버십, 약중력 질량

Astrafex Desktop에서 돌린 합성 자기일관성 실행이다. 두 은하단은 적색이동, 질량, 픽셀 크기가 다르다. 표면밝기는 스케일 반지름 `u = θ / θ_Δ`의 함수로 그려 두었다. `θ_Δ`는 `M_Δ`와 각지름거리가 정한다. `Δ = 200`이므로 질량은 `M_200c`이다.

멤버십은 한 목록을 세 갈래가 나눠 쓴다.

- ICL 표면밝기 스택은 BCG를 남기고, 위성·배경·별·미분류의 픽셀을 가린다.
- ICL 광량비는 BCG와 위성을 분모의 은하 빛으로 둔다. 배경·별·미분류는 빼다.
- 약중력은 라벨이 배경이고, `Σ_crit`에 쓸 적색이동이 유한하며, 렌즈보다 적색이동 간격 뒤에 있는 열만 남긴다.

분광 적색이동이 유한하면 그 창이 라벨을 정하고, 측광은 그 열을 바꾸지 않는다. 분광 적색이동이 없는 열만 점추정 측광적색이동 또는 `P(z)`로 채운다. 호출자가 준 은하단 창과 배경 창 중 하나에만 들어가면 그 라벨이다. 둘 다이거나 둘 다 아니면 미분류로 남는다. 그래도 미분류인 열은, 분광 위성으로만 맞춘 적색서열 위에 있고 등급 한계보다 밝으면 위성이 된다. 서열 밖의 푸른 은하는 미분류로 남는다. 측광 라벨은 `Σ_crit`의 적색이동을 대신하지 않는다.

약중력 부분은 전단 목록을 평균한 은하단 스택이 아니다. 렌즈 하나의 접선 전단으로 `M_200c`를 맞춘다. 접선 방향으로 돌린 우표 스택도 이 실행에 없다. 이 실행의 확산광 산출물은 BCG를 남긴 평균이므로, 이후 단계의 이름은 BCG+확산광이다.

이후 구현 순서는 `implementation_plan.md`에 있다. Astra 자문(배경은하 영상 스택, 그리고 2026-10-12 현실성·방법·데이터 자문)을 각 단계의 완료 조건으로 넣었다. 관측 과학 주장은 그 계획이 이름 있는 자료의 유효 표본과 계통오차를 확인한 뒤에만 붙인다.

1절 후속인 공통 kpc 스택, 질량 오차 전방모형, 반지름별 참여 표본, 서베이 점검 관문은 `icl/wgl.py`에 있다. 점검 항목이 비어 있으면 과학 측정 플래그는 내려가 있다. 이 합성 실행은 그 목록을 채우지 않는다.

2절의 전단 프로필 평균은 `ogfkit/wlcluster.py`의 `stack_shear_profiles`와 `plugins/wlcluster`의 `--task stack`에 있다. 픽셀 영상 스택은 그 함수가 거절한다. 이 디렉터리의 합성 실행은 렌즈 하나의 접선 전단만 맞추며, 그 결과에는 `science_measurement`가 없다.

3절의 접선 우표 진단은 `ogfkit/wlcluster.py`의 `stack_tangential_stamps`에 있다. 접선 방향으로 돌린 평균이 진단이고, 하늘 좌표 평균은 영점 그림이다. 이 함수는 질량을 맞추지 않는다. 이 디렉터리의 합성 실행에는 우표 스택이 없다.

4절의 사용자 안내 15.5·15.6절은 `ds10-web/docs/user_guide/sections/own_files.tex`에 있고, 방법론 초고는 `methods_draft.md`에 있다. 합성 실행이 기록하는 주장과 기록하지 않는 주장은 그 두 곳에 같은 표로 있다.

5절의 이름 있는 자료 관문은 `icl/named_gate.py`의 `named_data_gate`에 있다. 호출자가 준 행에서 HSC-SSP PDR3 `deepCoadd` global-sky, CAMIRA, S19A의 공통 영역과 분광 제한을 통과한 소스 수를 센다. 이 함수는 릴리스를 읽지 않는다. 세어진 수는 호출자가 준 행의 수이고, 과학 측정 플래그는 내려가 있다. 관측 과학 문장은 아직 없다.

스케일 스택의 메뉴 항목은 ICL의 `Stack BCG+diffuse`이다. 프로그램은 `plugins/icl/icl_stack.py`이고 `icl.wgl.stack_clusters`를 호출한다. 한 프레임 프로파일은 `ds9/library/ds9_icl.py`에 그대로 있다. local-sky 프레임과, 제품 이름이 `deepCoadd`가 아닌 프레임은 거절된다. 산출물 키는 `bcg+diffuse`이고, 서베이 점검표가 비어 있으면 과학 측정 플래그는 내려가 있다.

## 실행

OGFinder 루트에서:

```bash
python analysis/icl_wgl_membership/run_analysis.py
```

스크립트가 같은 디렉터리에 `history.json`을 쓴다. 그 파일의 숫자는 이 스크립트가 정의한 합성 은하단에서 나온 출력이다. 논문의 관측 질량이나 표본 복원율로 옮기지 않는다.

## Data availability

The synthetic self-consistency run for the scaled intracluster-light stack, the membership split, and the weak-lensing mass is in the Astrafex Desktop repository, [https://github.com/kjhan0606/OGFinder](https://github.com/kjhan0606/OGFinder), under `analysis/icl_wgl_membership/`. `run_analysis.py` writes `history.json`. That file records the caller settings and the numbers produced by the run. Cite the commit that contains this directory. `history.json` names the parent commit that was HEAD when the file was written; the analysis files sit on top of that commit until they are added.
