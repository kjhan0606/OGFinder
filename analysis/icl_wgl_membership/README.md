# 스케일 ICL 스택, 멤버십, 약중력 질량

Astrafex Desktop에서 돌린 합성 자기일관성 실행이다. 두 은하단은 적색이동, 질량, 픽셀 크기가 다르다. 표면밝기는 스케일 반지름 `u = θ / θ_Δ`의 함수로 그려 두었다. `θ_Δ`는 `M_Δ`와 각지름거리가 정한다. `Δ = 200`이므로 질량은 `M_200c`이다.

멤버십은 한 목록을 세 갈래가 나눠 쓴다.

- ICL 표면밝기 스택은 BCG를 남기고, 위성·배경·별·미분류의 픽셀을 가린다.
- ICL 광량비는 BCG와 위성을 분모의 은하 빛으로 둔다. 배경·별·미분류는 빼다.
- 약중력은 라벨이 배경이고, `Σ_crit`에 쓸 적색이동이 유한하며, 렌즈보다 적색이동 간격 뒤에 있는 열만 남긴다.

분광 적색이동이 유한하면 그 창이 라벨을 정하고, 측광은 그 열을 바꾸지 않는다. 분광 적색이동이 없는 열만 점추정 측광적색이동 또는 `P(z)`로 채운다. 호출자가 준 은하단 창과 배경 창 중 하나에만 들어가면 그 라벨이다. 둘 다이거나 둘 다 아니면 미분류로 남는다. 그래도 미분류인 열은, 분광 위성으로만 맞춘 적색서열 위에 있고 등급 한계보다 밝으면 위성이 된다. 서열 밖의 푸른 은하는 미분류로 남는다. 측광 라벨은 `Σ_crit`의 적색이동을 대신하지 않는다.

약중력 부분은 전단 목록을 평균한 은하단 스택이 아니다. 렌즈 하나의 접선 전단으로 `M_200c`를 맞춘다. 접선 방향으로 돌린 우표 스택도 이 실행에 없다.

## 실행

OGFinder 루트에서:

```bash
python analysis/icl_wgl_membership/run_analysis.py
```

스크립트가 같은 디렉터리에 `history.json`을 쓴다. 그 파일의 숫자는 이 스크립트가 정의한 합성 은하단에서 나온 출력이다. 논문의 관측 질량이나 표본 복원율로 옮기지 않는다.

## Data availability

The synthetic self-consistency run for the scaled intracluster-light stack, the membership split, and the weak-lensing mass is in the Astrafex Desktop repository, [https://github.com/kjhan0606/OGFinder](https://github.com/kjhan0606/OGFinder), under `analysis/icl_wgl_membership/`. `run_analysis.py` writes `history.json`. That file records the caller settings and the numbers produced by the run. Cite the commit that contains this directory. `history.json` names the parent commit that was HEAD when the file was written; the analysis files sit on top of that commit until they are added.
