# OGFinder / Astrafex Desktop 인수인계 (HANDOVER)

작성: 2026-10-08 (Asia/Seoul), 소유자 Juhan Kim <kjhan0606@gmail.com>.
추가: 2026-10-12 (Asia/Seoul). 6장. 1–5장은 2026-10-08 기록이다. 전체 서술은 ds10-web `HANDOVER.md` 11장이다.

**전체 인수인계 문서는 Astrafex Web 저장소의 [`HANDOVER.md`](https://github.com/kjhan0606/ds10-web/blob/master/HANDOVER.md)** (`kjhan0606/ds10-web`, private)이다.
프로젝트 목적, 이름(Astrafex), 아키텍처, 기능별 검증 결과, 작업 규칙, 남은 일은 모두 거기에 있다. 이 문서는 OGFinder에만 해당하는 부분만 다룬다.
두 저장소는 같은 상위 폴더에 나란히 clone한다 (`<WS>/OGFinder`, `<WS>/ds10-web`; 원래 box는 `<WS>=/workspace`).

## 1. 이 저장소는 무엇인가

* SAOImageDS9 8.7 소스 트리 + OGFinder 분석 기능 = **Astrafex Desktop** ("based on the OGFinder tool set"). 라이선스 GPL-3.0 (소유자의 잠정 결정, `README.md` License 절, `docs/licensing.md`, `docs/license_audit.md`). 번들된 third-party 구성요소는 각자 라이선스 유지.
* 저장소/폴더/라이선스 표기/ds9 쪽 GUI 문자열의 "OGFinder" 이름은 상표 검토가 끝날 때까지 그대로 둔다 (ds10-web `docs/naming.md`).
* 주요 디렉터리:

| 경로 | 내용 |
|---|---|
| `ds9/`, `tksao/`, `tcl8.6/`, `tk8.6/`, `tcllib/`, `ast/`, ... | ds9 8.7과 vendored Tcl/Tk/C 라이브러리 (빌드 대상) |
| `ds9/library/*.tcl` | ds9 GUI + OGFinder catalog panel, session recorder, layout (`docs/architecture.md`) |
| `plugins/<name>/` | 분석 plugin (manifest `plugin.json` + CLI `*.py` + Tcl GUI hook + `tests/`): extract, photometry, isophote, multifit, psfex, daophot, completeness, lsbg, icl, moving, trails, stacking, photoz_sed, sedcodes, lensmodel, spectra, xmatch, lightcurves, cluster, batch, repro, report, ai_services, ... (`docs/plugins.md`) |
| `plugins/ds10core/` | **Astrafex Core의 Desktop shell**. `ds10.py`가 ds10-web의 `server/ds10core`를 찾아 **별도 process**로 `python -m ds10core ...` 실행 (import하지 않음 = 라이선스 경계). 코어 탐색: `--core-dir`, `DS10_CORE`, `~/.ds9/ogf_params/ds10core.json`, `../ds10-web/server`, `/workspace/ds10-web/server`. GUI: Measure 탭 chip "Astrafex Core"(More 메뉴). 단계 목록과 검사: `docs/ds10core.md` |
| `ogfkit/` | 공용 Python 라이브러리 (cliexpand, trails, photerr, ...) |
| `moving/`, `icl/`, `lsbg/`, `photo_z/`, `ai_bridge/`, `ai_merge/`, ... | 분석 패키지 (이동천체, ICL, LSBG, photo-z, AI 서비스 연결) |
| `regression/` | 공개 데이터 회귀 검증 세트 (`docs/testing.md` 끝) |
| `scripts/` | `run_all_checks.sh`(모든 검사), `verify_*.tcl/.sh/.py`(GUI 검사), `astrafex`(Core CLI), 재현 스크립트 |
| `docs/` | 기능 문서, `testing.md`, `progress_log.md`(R1–R25 진행 기록), `phase2_web_design.md`(웹 설계), `shots/`(스크린샷) |

## 2. 빌드와 환경

```bash
# Debian/Ubuntu 빌드 의존성 (BUILD.txt)
sudo apt-get install -y automake autoconf libx11-dev zlib1g-dev libxml2-dev libxslt1-dev libxft-dev zip tcl8.6-dev tk8.6-dev \
     libcfitsio-dev xvfb x11-utils xdotool tcl
cd OGFinder
unix/configure && make            # -> bin/ds9   (원래 box: 8 core로 수십 분)
./build_sextract.sh               # -> bin/ds9_sextract (SEP 기반 추출기; batch/repro/report 테스트에 필요)
python3 -m venv --system-site-packages ~/ogf_venv
~/ogf_venv/bin/pip install --extra-index-url https://download.pytorch.org/whl/cpu -r scripts/ogf_venv_freeze.txt
export OGFINDER_PYTHON=~/ogf_venv/bin/python3     # run_all_checks 기본값은 /workspace/ogf_venv/bin/python3
export OGF_TEST_FITS=/path/to/fits                # m51.fits 가 있는 폴더 (기본 /workspace/fits); HUDF 검사는 hudf_f160w.fits 도 사용
```

* `scripts/ogf_venv_freeze.txt` = 원래 box에서 검사가 돌던 venv의 package 목록 (제품 dependency 목록이 아니라 기록). GPL(`rebound`, `assist`), LGPL(`sep`) 주의.
* **빌드가 tracked 파일을 바꾸거나 지운다** (`ast/config.h.in`, `tcllib/modules/*/build/*`, `tkimg/compat/libtiff/build/*` 등). 이런 변경은 **절대 commit하지 않는다**: 항상 `git add <지정 파일>`로만 commit. 특히 `ast/config.h.in`은 commit 금지.
* `m51.fits`는 600x600 2D 영상이다 (출처는 저장소에 기록되어 있지 않음). GUI 검사 일부(`cat_behavior` golden 등)는 원래 box의 그 파일 기준이라, 다른 영상이면 golden 비교가 다를 수 있다.
* `ds10core_*` 검사와 Astrafex Core 단계는 sibling `../ds10-web` checkout이 필요하다 (없으면 SKIP).

## 3. 검사 (`scripts/run_all_checks.sh`)

```bash
flock /tmp/astrafex_heavy.lock scripts/run_all_checks.sh            # 전체 (long 검사 포함)
flock /tmp/astrafex_heavy.lock scripts/run_all_checks.sh --quick    # long/network 검사(session_replay, link_bench, moving_session, regression_data) 제외
scripts/run_all_checks.sh --list                                     # 검사 이름
scripts/run_all_checks.sh --only geometry,mouse                      # 일부만
```

* GUI 검사는 `--display`(기본 `:77`)에 X가 없으면 Xvfb(1400x1000x24)를 직접 띄우고, 끝날 때(정상/에러/INT/TERM/HUP) PID로 정리한다. Xvfb는 fd>2를 닫고 띄워서 `flock`의 lock fd를 물려받지 않는다.
* 로그: `$OGF_CHECK_OUT` (기본 `/tmp/ogf_checks.<pid>/NAME.log`). 끝나면 지운다.
* 검사 표와 각 검사의 필요조건: `docs/testing.md`. live AI agent 테스트는 `OGF_LIVE=1`일 때만.
* `link_bench`는 `/workspace/work/inj*.json.pkl` 주입 세트, `moving_session`은 캐시된 BB89 노출(`~/.ds9/moving_cache`)이나 네트워크, `regression_data`는 `$OGF_DATA_CACHE`(기본 `~/.cache/ogfinder_regression`)에 공개 데이터를 받는다. 없으면 SKIP (실패 아님).
* 실행 중에는 저장소 파일을 수정하지 않는다.

### 3.1 최신 결과

인수인계 시점의 `run_all_checks.sh` 전체 표는 이 저장소에 없다. 검사 로그는 끝나면 지우고, `docs/progress_log.md` R24도 "full-run result is in the final report of this session (log kept at scripts run output, not in the repo)"라고 적는다.

저장소에 숫자가 있는 마지막 전체 실행은 `docs/progress_log.md`의 "Final full run"이다.

* HEAD `5d3271b7e` (2026-10-02): 0 failures, 35 PASS, 1 SKIP (`agent_real`, agent CLI 없음).
* 시간: newplugins 708 s, session_replay 183 s, link_bench 191 s, moving_session 629 s.
* 출력 파일은 원래 box의 `/workspace/work/run_all_checks_full.txt`이다. 저장소 밖이라 여기에는 없다.

그 뒤에 `ds10core`, `aperphot_gui`, `icl_tests`, `nightlink_tests`, trails stack이 들어왔으므로, 35 PASS는 인수인계 시점의 검사 목록과 같지 않다. 통과로 적힌 `geometry`는 R1 (HEAD `ac4fff4be`, table 181/769/154)까지다.

### 3.2 알려진 문제

* **geometry 검사의 Tcl `ALPHA_J2000` 오류**: `scripts/verify_geometry.tcl`이 추출 결과 행에서 `dict get $r ALPHA_J2000`을 하는데, 환경(추출 결과/WCS, 테스트 영상)에 따라 그 key가 없으면 Tcl 오류로 실패한다. 코드 회귀가 아니라 환경 문제로 분류해 왔다. 실패 로그 자체는 저장소에 없다. 스크립트는 `CatalogPanelExtract` 뒤 카탈로그의 3·4·5·7번째 행에 `ALPHA_J2000`과 `DELTA_J2000`이 있다고 가정한다 (`scripts/verify_geometry.tcl` 43–46행). 추출 결과에 그 열이 없으면 `dict get`이 Tcl 오류로 끝난다. 사용 영상은 `$OGF_TEST_FITS/m51.fits`(기본 `/workspace/fits/m51.fits`)이다.
* **`scripts/verify_ai_gui.py`, `scripts/verify_agent_gui.py`가 자기가 띄운 Xvfb를 남긴다**: 해당 display에 Xvfb가 없으면 `subprocess.Popen(['Xvfb', display, ...])`으로 띄우고 끝날 때 정리하지 않는다 (fd도 닫지 않아 lock fd를 물려받을 수 있음). `run_all_checks.sh` 안에서는 먼저 Xvfb를 띄우므로 문제 없고, **단독 실행할 때** 남는다. 끝난 뒤 `ps -ef | grep Xvfb`로 PID 확인 후 `kill <PID>`. 고칠 일: `run_all_checks.sh`처럼 띄운 Xvfb의 PID를 기억해 종료 시 kill, fd>2 닫기.
* trails: SDSS 후보 29개의 육안 확인은 2026-10-08에 끝났다 (spike 16, 별 5, 연속 streak 없음 8, 위성 trail로 확인된 것 없음). 표는 `plugins/trails/validation/trails_validation.md`. 실제 satellite-trail truth sample은 2026-10-08에 ACS/WFC flc 8장, 라벨 9개로 기록했다 (recall 6/9, precision 6/12). `jc8m32j5q`에서 이전에 적힌 3개는 그 노출에서 검출기가 낸 선분이다. 이동천체 linker의 검색은 two-body이고, 연결된 그룹의 N-body 보정은 별도 프로세스다 (`--nl-nbody`, `docs/moving_objects.md`). 이 머신에서는 ASSIST 없이 CODES로 입자 하나만 측정했다. `--nl-ranging`은 그 그룹의 궤도요소 분포다. 1시간 호는 `iod.ranging`(a 16/50/84% = 0.78/1.77/5.18 AU), 21일 호는 국소 Laplace 샘플(a 2.6941/2.6961/2.6985 AU, 스케일된 공분산으로 다시 잰 값)이다. OpenOrb 실행 파일은 쓰지 않았다. (1862) Apollo 관측 12줄의 예제는 `examples/apollo_1862_ranging/`이다. 자세한 남은 일은 ds10-web `HANDOVER.md` 8장.

## 4. 최근 commit (이번 인수인계 시점)

아래 목록은 2026-10-08, master `ad8761063` 기록이다. 2026-10-12의 HEAD `4dee078b467d6f1dac1aaef8ef2f80a55f577222`와 그 사이 24개 commit은 6장이다.

master `ad8761063` (2026-10-08, 이 문서) 기준 최근 commit. 한 줄은 `git log` 제목이다. 전체 이력은 약 5573 commit이고, 기능 단위 기록은 `docs/progress_log.md`다.

* `16798dd79` (2026-10-08) trails stack을 실제 HST dither (jc8m32010)에서 돌리며 WCS 왜곡, CR 제거, 스케일을 고침.
* `fa0c7cfe1` (2026-10-07) `run_all_checks`와 moving-session 검사가 띄운 Xvfb를 모든 종료 시그널에서 PID로 정리. lock fd를 물려받지 않게 fd>2를 닫고 띄움.
* `aaf0114a3`, `1fa6218ba`, `c6c5dc819` (2026-10-07) ICL: 공유 마스크와 밴드별 영점, hot+cold 마스크, μ가 정의된 픽셀만으로 f_ICL. Abell S1063 숫자를 설명서에 반영.
* `8f9e9bd28`, `885162e78` (2026-10-07) Astrafex Core 별 구경 측광 단계와 `aperphot_gui` 검사, 데스크톱 스크린샷.
* `2932ee850`, `9f6f71461`, `0b821cfa4` (2026-10-07) 하늘 카탈로그, MAST, 사용자 플러그인 단계와 `verify_ds10core` 검사.
* `f54dda631` (2026-10-06) PSF 측광 / CMD 단계.
* `60373f050`, `0c9d0fb8a`, `d28f5b6a3` (2026-10-05–06) LSBG finder와 `--classify`를 데스크톱에 연결.

## 5. 규칙 요약 (전체는 ds10-web `HANDOVER.md` 7장)

파일 지정 commit만; `ast/config.h.in` commit 금지; 무거운 작업은 `flock /tmp/astrafex_heavy.lock`으로 하나씩 (RAM 15 GB, swap 없음); process는 PID로 kill (`pkill -f` 금지); push는 사용자가 지시할 때만; git identity `Juhan Kim <kjhan0606@gmail.com>`; 기능 추가 시 `docs/`와 ds10-web 사용자 설명서도 갱신.

## 6. 2026-10-12 추가

이 머신의 작업 폴더는 `/home/kjhan/BACKUP/Astrafek`이다. HEAD는 `origin/master`와 같은 `4dee078b467d6f1dac1aaef8ef2f80a55f577222`이다. 스케일 ICL, 멤버 맞춤, HSC 디버그가 멈춘 지점, 커밋에서 빼 둘 파일은 ds10-web `HANDOVER.md` 11장에 있다. 여기에는 이 저장소에서 이어서 볼 것만 적는다.

### 6.1 `ad8761063` 다음 commit

`git log --oneline ad8761063..HEAD` 24개. 한 줄은 제목이다.

* `4dee078b4` Add the CAMIRA member match and the masked kpc stack step.
* `893b6e3af` Keep the full spectroscopic slice for the HSC stack list.
* `7956b45b6` Save the HSC stack debug procedure as scripts for the paper.
* `7e57dcd68` Wire the scaled BCG+diffuse stack into the ICL menu.
* `c8d3e44a0` Count lens sources at the named-data gate before any observed sentence.
* `0b81eb9bc` Draft the methods note for the scaled light stack and the lensing mass.
* `dd26ab3f0` Add a tangential postage-stamp diagnostic beside the shear catalogue.
* `5c72e470f` Average tangential shear profiles before the next NFW fit.
* `dd4a87e9c` Name the scaled stack BCG+diffuse and keep a kpc comparison beside it.
* `e441eb7f5` Record the synthetic ICL and weak-lensing self-consistency run.
* `e7cfadc2a` Label satellites and background galaxies for ICL masks and weak lensing.
* `1bbecb041` Scale cluster images by apparent size and M_200c before an ICL stack.
* `3bbaa6dd2` Add cluster weak lensing from metacalibration through an NFW mass.
* `9a4adcb8f` Add learned spectrum steps for a temperature letter, galaxy neighbours, and a user template match.
* `01e6610a3` Add a supernova type from one host redshift and the line depths.
* `8b2325e68` Add a synthetic spectrum example for preprocessing and the science step.
* `0282c122c` Classify one spectrum from the lines that actually fit.
* `92ef9a1f0` Draw the default galaxy stamp with an in-tree Sérsic, shear, and Gaussian PSF.
* `b7c00b6e0` Use in-tree models for the standalone measurements, orbits, and SED fits.
* `7b3b757cb` Record orbit ranging for linked groups and stop the GUI checks leaving Xvfb up.
* `e2f307324` Refine accepted nightlink groups in an N-body child process.
* `1d3edf8b6` Record the ACS satellite-trail truth sample (recall 6/9, precision 6/12)
* `2984e9971` Record the visual inspection of the 29 SDSS trail candidates
* `8984e74f7` HANDOVER.md: fill check-result, geometry, and git-log placeholders from the repo

### 6.2 이 저장소에서 이어서 볼 경로

* `analysis/icl_wgl_membership/`: 합성 자기일관성 (`run_analysis.py`, `history.json`), `implementation_plan.md`, `methods_draft.md`. `history.json`의 `parent_commit`은 `9a4adcb8f421838115acbff9c8612671d4fff723`이다. 그 파일은 고치지 않고, 스크립트는 다시 돌리지 않는다.
* `analysis/icl_wgl_membership/hsc_debug/`: 논문에 함께 제공하는 디버그 스크립트 01–08과 `run_debug.sh`. 실행 데이터는 `/gpfs/kjhan/Astrafek/data`에만 있다. 카탈로그, FITS, 로그는 커밋하지 않는다.
* 메뉴에 연결된 스택은 `plugins/icl/icl_stack.py`의 `Stack BCG+diffuse`이다. 산출물 키 `bcg+diffuse`는 BCG를 남긴 평균이다. 가린 kpc 평균의 이름은 `bcg-and-members-masked`이고 `hsc_debug/08_mask_stack.py`가 쓴다. 2026-10-12 실행은 컷아웃 HTTP 401에서 멈추었고, `cutouts/frames.tsv`가 없어 평균 파일은 쓰이지 않았다. `science_measurement`는 거짓이다.
* 접선 우표 진단은 `ogfkit/wlcluster.py`의 `stack_tangential_stamps`이다. 메뉴에는 연결되어 있지 않다.

### 6.3 작업 트리에 남은 파일

2026-10-12에 커밋하지 않고 남겨 둔 파일이다. 인수인계 커밋에 넣지 않는다.

* 수정: `ds9/library/ds9_psf_deconv.py`, `plugins/galaxy_model/galaxy_model.tcl`, `plugins/star_psf/plugin.json`, `plugins/star_psf/star_psf.tcl`, `psf_deconv/psf/sim_psf.py`
* 미추적: `plugins/galaxy_model/galaxy_fit.py`, `plugins/galaxy_model/tests/`, `psf_deconv/tests/test_diffraction_psf.py`

나선 은하는 exponential disc(preset `exp`)로 둔다. bulge+disk는 연결하지 않는다. `ast/config.h.in`과 빌드 산출물은 커밋하지 않는다.
