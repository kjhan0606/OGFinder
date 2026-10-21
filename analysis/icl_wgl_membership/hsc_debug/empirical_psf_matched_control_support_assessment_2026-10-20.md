# Gaia PSF 별의 매칭 대조장 위치 지원 검사 — 2026-10-20

## 결과

기존 HSC-I 표적 영상에서 Gaia 지지 별 우표 333개(227개 표적 프레임)를 골라,
각 별의 표적 은하단 중심 기준 각거리 오프셋을 이미 수신한 매칭 대조장 중심에
옮길 수 있는지 검사했다. Astrafek은 target/control WCS, 하위 필터, 원본 우표의
pixel 좌표와 SHA-256을 확인했다. 원본 science 영상의 픽셀은 읽지 않았다.

| 검사 항목 | 결과 |
|---|---:|
| 저장된 별 우표 | 333 |
| 고유 표적 프레임 | 227 |
| 같은 HSC 하위 필터의 target/control 쌍 | 319 (217 프레임) |
| `i`/`i2` 필터 불일치로 제외 | 14 (10 프레임) |
| 401×401 control 우표 전체가 영상 안에 위치 | 319/319 |
| HSC 요청·자료 수신 | 없음 |

대조 우표는 별의 절대 RA/Dec를 옮긴 것이 아니라 target cluster center에 대한
east/north 각도 오프셋을 control center에 적용한 뒤 control WCS로 pixel 변환했다.
따라서 표적 별과 대조장 우표는 cluster-centric 상대 위치를 보존한다. 유효 쌍의
가장자리 최소 여유는 x 1.94 pixel, y 3.78 pixel이고, 일부 위치는 영상 경계에
가깝지만 요청한 전체 우표가 들어온다.

## 해석과 다음 단계

이 결과는 매칭 대조장 영상에 같은 각거리의 placebo/background 우표를 추출할
기하·필터 지원이 있음을 보인다. PSF 날개 측정이나 ICL sky 보정은 아니다. 대조장
우표에는 중심 별이 없으며 source mask·잔류 sky·target 쪽 ICL 곡률을 별도로
다뤄야 한다. 다음 분석은 기존 global-sky target/control 픽셀에서 동일한 반경·하늘
추정기를 적용해 대조장 잔차 분포를 만들고, 별 프로필의 국소 sky 민감도를 비교한다.
그 차이를 PSF wing으로 단독 해석하거나 ICL 보정값으로 쓰지 않는다.

HSC FAQ는 coadd PSF 모델을 별도 `deepCoadd_psf` 자료로 설명하고 PSF Picker를
모델 획득 경로로 안내한다. HSC 문서는 확장천체 날개 분석에는 global-sky coadd를
사용하도록 안내한다. 따라서 이번 매칭 대조 검사는 global-sky 영상 잔차를 직접
평가하는 보조 실험이며 모델 PSF를 대체하지 않는다 ([HSC PDR3 FAQ](https://hsc-release.mtk.nao.ac.jp/doc/index.php/faq__pdr3/)).

Astrafek 재현기:
[`audit_empirical_psf_matched_control_support.py`](audit_empirical_psf_matched_control_support.py)
와 focused tests
[`test_audit_empirical_psf_matched_control_support.py`](test_audit_empirical_psf_matched_control_support.py).
테스트 3개 통과. 실행은 lageunha CPU에서 기존 FITS header/WCS만 읽었으며, 실행 명령은
다음과 같다.

```bash
python3 analysis/icl_wgl_membership/hsc_debug/audit_empirical_psf_matched_control_support.py \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/empirical_psf_matched_control_support_2026-10-20_01
```

산출물은
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/empirical_psf_matched_control_support_2026-10-20_01/`에
있다. `summary.json` SHA-256은
`a6d80adbfe8d01e466930dc68d44355287124b1b5ad8d8221a0c025f4cd7346f`, 표 SHA-256은
`fc255bf0704ddfd43e817e6b21e4a66f28c593e6748fda50976be12558f5da15`다. 별도 검산에서
333행, 319개 동일 필터 우표, 14개 필터 불일치, 출력 표 체크섬을 다시 확인했다.
`science_measurement=false`, `psf_approved=false`를 유지한다.
