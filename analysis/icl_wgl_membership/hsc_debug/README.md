# HSC 디버그 절차

이 디렉터리의 스크립트는 ICL 스케일 스택을 실제 목록으로 준비하는 절차이다. 논문에 함께 제공한다. 실행 결과는 과학 측정이 아니다. `science_measurement` 는 거짓으로 남는다.

작업 디렉터리는 `/gpfs/kjhan/Astrafek/data` 이다. 다른 경로를 쓰려면 `ASTRAFEK_DATA` 를 둔다. OGFinder 위치는 `ASTRAFEX_OGF_ROOT` 이다.

```bash
bash scripts/run_debug.sh
```

순서:

1. `01_fetch_catalog.sh` 는 CAMIRA S20A Wide 별마스크 목록을 받는다. SHA-256 이 기록된 값과 같으면 다시 받지 않는다.
2. `02_select_clusters.py` 는 분광 적색이동 창을 통과한 가장 촘촘한 0.05 적색이동 구간의 은하단을 고르고 `logs/selected_clusters.tsv` 를 쓴다. 6각분보다 가까운 쌍은 richness 가 큰 쪽만 남긴다.
3. `03_mass_gate.py` 는 목록의 질량이 `M_200c` 가 아니라는 점(Oguri et al. 2018, arXiv:1701.00818 은 `M_200m`)을 기록한다. `clusters.json` 을 만들지 않는다.
4. `04_cutout_status.py` 는 `cutouts/` 안의 FITS 유무만 적는다. 영상을 내려받지 않는다. HSC 영상은 로그인 뒤 아카이브에서 가져온다.
5. `05_stack.py` 는 `masses.tsv` 와 `cutouts/frames.tsv` 가 있을 때 메뉴와 같은 명령으로 `icl_stack.py` 를 실행한다. 둘 중 하나가 없으면 스택을 실행하지 않는다.

`masses.tsv` 의 머리글은 `id`, `m200c_msun`, `mass_definition`, `citation` 이다. `mass_definition` 값은 `200c` 이다. 식은 그 인용에서 그대로 옮긴다.

`cutouts/frames.tsv` 의 머리글은 `id`, `filename`, `sky`, `product` 이다. `sky` 는 `global` 또는 `global-sky` 이다. `product` 는 비우거나 `deepCoadd` 이다. 중심은 astropy `world_to_pixel` 의 0-기반 픽셀이고, 픽셀 크기는 헤더의 CD 또는 CDELT 이다.
