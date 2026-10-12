# OGFinder / Astrafex Desktop 인수인계 (HANDOVER)

작성: 2026-10-08 (Asia/Seoul), 소유자 Juhan Kim <kjhan0606@gmail.com>.
추가: 2026-10-12 (Asia/Seoul). 6장. 1–5장은 2026-10-08 기록이다. 전체 서술은 ds10-web `HANDOVER.md` 11장이다.
추가: 2026-10-13 (Asia/Seoul). 7장. ICL–위성은하 외곽구조 장기 연구주제와 도구 설계.
추가: 2026-10-13 (Asia/Seoul). 8장. 분석 목적 기반 범용 파이프라인 요구사항.

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

### 6.4 2026-10-12 승인된 HSC kpc 시범과 인증 재확인

앞 절의 401은 최초 실행 기록이다. 사용자가 새 `hsc_debug/probe_login.py`를
터미널에서 실행한 뒤 HTTP 200을 받았고, 이어 `--global-sky`로 공식 매뉴얼과
`coadd/bg` 한 장의 모든 HDU 헤더를 받았다. 매뉴얼은 이를 `{rerun}/deepCoadd`
의 coadd 및 공격적 배경 차감 전 중간 자료로 설명한다. 헤더는 HSC-I, inherited
primary `FLUXMAG0=63095734448.01944` (AB 영점 27), HDU 1 IMAGE, 2 MASK, 3 VARIANCE를
확인한다. 이에 이 한 프레임의 제품 출처와 영점 매핑은 통과했다. 판정은 실행
데이터의 `logs/cutout_global_review.json`에 남겼다. 영상 픽셀은 probe가 폐기했으므로
대량 다운로드와 `cutouts/frames.tsv` 작성, 관측 스택은 아직 실행하지 않았다.
`MP_*` 헤더 값은 bit 위치이고 실제 제외 비트 값은 `2**n`이다. 검출 마스크에서
`DETECTED`만으로 BCG 픽셀을 제거하지 않는다.

승인된 첫 목표는 질량을 쓰지 않는 HSC-I BCG+확산광 kpc 프로파일이다.
`analysis/icl_wgl_membership/pilot_README.md`에 입력 스키마와 명령이 있고,
`hsc_pilot.py`의 prepare/inspect/measure, `icl/pilot.py`, `pilot_job.sh`가 새 경로다.
BCG를 남기고 모든 다른 검출을 가린다. 개별 원본 고리 프로파일의 동일 은하단
가중 평균이 주 추정량이며, 면적 보존 영상은 별도 대조다. matched random과
영역별 jackknife 경로가 있으나 실제 random·PSF·sky 전달·주입/회수 검증은
미완료이다. 모든 시범 산출물의 `science_measurement`는 false이다.

실행 데이터 아래 `pilot_2026-10-12/pilot_manifest.json`에 richness 순위로 12개를
준비했다. 영상·교정·품질 입력은 비어 있으며 `inspection_initial/inspection.json`
검사에서 12개 모두 멈춘다. 관측 스택은 실행하지 않았다. Slurm 작업 420654의
회귀 검사는 35 passed (8.94 s); 로그는 `logs/hsc-icl-pilot-420654.out`이다.
`mass_interpretation`은 입력 공분산을 기록할 뿐이며 `covariance_recorded`와
`covariance_propagated=false`로 구별한다. 자세한 전체 갱신은 Web 인수인계 11.9다.
이 추가 작업은 아직 커밋·푸시하지 않았다.

## 7. 2026-10-13 장기 연구주제

사용자는 **ICL과 위성은하 외곽구조의 관계를 공동 스태킹으로 연구**하고 장기
주제로 유지하도록 지시했다. 기준 문서는
[과학 연구계획](analysis/icl_wgl_membership/science_research_plan.md)이며,
[분석 도구 설계](analysis/icl_wgl_membership/research_software_design.md)와
[실행 계획](analysis/icl_wgl_membership/implementation_plan.md)을 함께 읽는다.

은하단 중심의 BCG+확산광과 위성 중심의 외곽광을 같은 은하단 표본에서 측정한다.
ICL의 조작적 정의, 공동 주입/회수, 광량 분배가 만드는 가짜 상관, 은하단 단위
공분산을 검증한 뒤 관계를 해석한다. 약중력은 후속 확장이다. 기존 213개 영상의
검사와 총 800개 후보 준비는 실행 계획에 있다. 실제 전송 상태는 실행 데이터의
상태 파일로 확인하며, 위 6장의 최초 401·헤더 시범 기록을 현재 상태로 쓰지 않는다.

다음 우선 작업은 추가 측광·PSF 입력 감사와 30–50개 예비 표본 설계다.
CPU 분석은 `lageunha`, 대용량 데이터는 `/gpfs/kjhan/Astrafek/data`에 둔다.
장기 연구주제는 등록됐으나 공동 분석 코드는 설계 단계이고 관측 과학 결과는
미완료다. 인증정보를 문서·설정에 남기거나 기존 다운로드 방식을 새로 작성하지 않는다.

## 8. 2026-10-13 목적 기반 파이프라인 도구

사용자는 “구상성단 별 측광”, “LSBG” 등 분석 목적에 따라 파이프라인을 자동
구성하고, 그 흐름을 따라 최종 결과를 얻는 도구를 요청했다. 요구사항·구조·
예시·개발 순서는 [science_workflow_designer.md](docs/science_workflow_designer.md)에
등록했다. `astrispex`는 사용자 표현을 따른 임시 명칭이며 제품명 변경은 미확정이다.

기존 plugin/CLI/headless, batch, session, repro를 기반으로 목적·자료 감사,
과학 recipe, 입출력 계약을 검증하는 DAG planner와 단계별 실행을 설계한다.
없는 기능을 자동 실행 가능으로 표시하지 않고, 자료 부족·상한·검증 실패도
결과에 기록한다. ICL–위성 공동 연구는 이 도구의 장기 과학 recipe다.
현재는 설계 문서만 등록했으며 planner·신규 UI·공동 분석 엔진은 미구현이다.

### 8.1 같은 날 추가 구현: Desktop / Web 순서도와 저장

위 문서 작성 이후 두 저장소에 11개 분야의 순서도, 단계 도구 연결, 검토 메모,
중간 저장·복원을 구현했다. 사용법은 [science_pipelines.md](docs/science_pipelines.md).
Desktop은 `plugins/workflow_designer`, `ogfkit/science_workflow.py`이며
Measure → Workflows에서 연다. Web은 별도 구현으로 동일한 선언적 recipe 데이터만
공유하고 Desktop 코드를 import하지 않는다. 모든 레이아웃에서 Ctrl-K → Science
pipelines로 열 수 있다.

Desktop은 완료 산출물과 카탈로그를 새 스냅샷에 복사하고 SHA-256을 검증한다.
원본 영상은 경로·크기·수정시간 참조이므로 별도 보존한다. 자동 기록되지 않은
결과 파일은 사용자가 연결한다. 복원은 작업 사본을 만들어 저장본을 보존한다.
Web은 계정 세션의 저장 시점까지 완료 산출물을 ZIP으로 내려받고, 서버의 세션
저장 지점으로 복원한다. 누락·변조 파일은 저장 성공으로 처리하지 않는다.

Desktop 저장/복원 6건, Web API 3건과 브라우저 선택·메모·ZIP·특정 저장 지점
복원 1건이 통과했다. 타입 검사·Web 빌드·Desktop manifest 검사도 통과했다.
검사는 CPU로 lageunha에서 수행했다. Desktop Tk 화면의 실제 클릭 검사는
아직 하지 않았다. 자연어/입력 감사 기반 DAG 자동 생성, 미연결 분석 adapter,
모든 분야의 최종 과학 검증은 후속 개발이며 현재 UI의 `reviewed`는 수동 검토다.
이 작업은 커밋·푸시하지 않았다.

### 2026-10-13 사용자 UI 수정: 본창 파이프라인

사용자는 팝업 대신 본창의 분석 메뉴를 클릭하면 왼쪽 위 PIPELINE 단계 목록이
해당 분야에 맞게 바뀌도록 지시했다. Web은 11개 분야를 분석 메뉴에 직접 넣고,
선택 단계의 설명·메모·저장을 중앙 분석 영역으로 옮겼다. 일반 도구는 기본
단계 목록으로 연결한다. Desktop도 별도 toplevel을 제거하고 기존 분석 패널에
순서도와 선택기를 넣었다. Web의 이전 캡처 4장은 같은 파일 경로에 본창 캡처로
갱신했다. 저장/복원 동작은 유지하며 과학 검증 완료와 수동 검토를 구별한다.

### 2026-10-13 추가 구현: 단계 실행과 작업 연결

Web 중앙 단계 영역에 직접 실행과 작업 ID·산출물 목록을 추가했다. 단계별
`runs` 참조는 계정·기록 세션·분석기 일치를 서버가 검사하며 실행 상태는 실제
작업 기록에서 확인한다. 이전 저장 지점은 당시 완료 산출물만 포함한다. 완료
기록이 아직 없는 저장 지점은 완전한 결과 저장으로 표시하지 않는다. 서버에
남은 작업은 재개 시 기존 SSE로 관찰하며 새 작업을 중복 제출하지 않는다.
Desktop은 순서도에서 직접 실행한 도구의 recorder 항목을 단계에 연결하고
write trace/idle callback으로 상태를 갱신한다. 실행 중 스냅샷을 새 프로그램에서
열면 재검토로 바뀌며 계산 프로세스가 복원됐다고 표시하지 않는다. 앞 단계
재실행·입력/설정 변경은 후속 결과를 재검토 상태로 바꾼다. 수동 검토와 과학
검증은 계속 구별한다. 사용법은 `docs/science_pipelines.md`에 추가했다.

이번 추가 구현 검사는 Desktop 7건, Web API 6건, 브라우저 3건 통과다. 실제
합성 영상 검출·산출물 표시·새로고침 재개·설정 변경 재검토를 포함한다. 타입
검사·빌드도 통과했고 CPU 검사는 lageunha에서 수행했다. 본창 캡처 4장은
갱신했다. Desktop 실행 상태는 Tcl로 검사했으며 실제 Tk 클릭 검사는 아직
하지 않았다. 커밋·푸시와 관측 과학 분석 실행은 하지 않았다.

### 2026-10-13 미연결 분석단계 추가 연결

영상 입력 검사와 결과 저장을 실제 연결했다. 원본 공유 코어에 동일 은하단 가중 kpc 프로파일 stack을 추가했다. Web에는 계정 소유 TSV import 및 통합 PSF/LSBG·구경 측광·PSF stamp·검출 주입 회수 adapter를 연결했다. 통합 실행과 개별 검토를 구별한다. Desktop 명시적 산출물도 recorder에 등록한다. 계약과 남은 전용 분기는 docs/science_pipelines.md 마지막 절 참조. 합성자료 검사: Desktop 9, 서버 15, 브라우저 4개 통과. 실제 HSC 계산/다운로드는 수행하지 않았으며 science_validated=false다.

### 2026-10-13 다음 구현: 1D 분광 준비

공유 코어 spectrum_prepare의 단일 1D TSV 입력 검사·오차 가중 Chebyshev 연속광을 양쪽에 연결했다. 웹 계정 소유 import와 Desktop catalog 로드, 산출물·입력 표·공분산 저장/재개를 포함한다. 파장계와 flux 단위는 명시 필수이며 AGN/host/Fe II 분해를 대신하지 않는다. 웹의 미구현 후속 1D 과학 adapter는 버튼을 비활성화한다. Desktop workflow의 close 이름 충돌도 ::close로 수정했다. 입력 계약과 남은 단계는 docs/science_pipelines.md 마지막 절 참조. 실제 관측 데이터 계산은 수행하지 않았다.

1D 분광 최종 검증: Desktop 11, 서버 25, 브라우저 5개(총 41) 통과. Web typecheck/build 통과; Desktop 39 plugins/289 steps/873 params manifest 검사 문제 0. 검증은 lageunha CPU에서 수행했다. 연속광 출력에도 원 MASK를 유지하며 각 단계의 입력 snapshot 파일은 서로 다른 이름으로 보존한다.


### 2026-10-13 후속: 동정된 고립선/적색편이 후보

두 제품의 항성/은하 선 측정 및 은하/AGN 적색편이 노드를 original shared-core `spectrum_lines.py`에 연결했습니다. 사용자 rest 파장·초기 z·반창 폭·방출/흡수 입력, 국소 연속광+Gaussian 공동 적합, covariance/모델/원본/측정표 보존. 여러 accepted 선이 일치할 때만 조건부 z 후보를 보고하며 자동 선 동정·전역 검색 및 AGN 분해는 미구현입니다. 상세 사용법·판정 범위는 docs/science_pipelines.md 마지막 절에 있습니다. 실제 HSC 데이터 분석·다운로드, 커밋·푸시는 수행하지 않았습니다.

추가로 Web `sessionController.ts`의 영상 없는 workflow 복원에서 활성 분석 표를 생략하던 문제를 수정했습니다. 저장→새로고침→후속 선 측정→다시 ZIP 저장 흐름을 브라우저 회귀 검증에 포함합니다.

검증(lageunha, CPU/BLAS 1 thread): 서버 기존 25개 및 새 선 측정 15개, Desktop 12개, 실제 브라우저 5개 통과(총 57개). Web typecheck/build 통과. Desktop manifest 39 plugins/291 steps/878 params, 0 problems. 변경 파일 diff whitespace 검사 통과. 합성 자료 검증이며 실관측 과학 검증은 수행하지 않았습니다.


### 2026-10-13 후속: 초기 z 없는 피크 탐색/적색편이 후보

항성·은하·AGN 분광에 search 노드를 추가(11 recipes/110 nodes). Web ds10_spectrum.search / Desktop ds10core.spectrum-search, original NumPy/SciPy shared core. 사용자 기준선 목록, 구간별 median detrending/피크 대조, 다중 대응과 ambiguity 보존; 후보표·피크표·JSON·입력 원본을 기록/체크포인트/ZIP에 포함. 과학 승인·설정 자동 변경을 하지 않으며 후보를 검토해 다음 정밀 측정 설정에 반영합니다. AGN 분해는 아직 미구현입니다. docs/science_pipelines.md 마지막 절에 계약·한계가 있습니다. 실제 관측 데이터 계산·커밋·푸시는 하지 않았습니다.

검증(lageunha, CPU/BLAS 1 thread): 서버 52개·Desktop 13개·브라우저 5개, 총 70개 통과. Web typecheck/build 통과. Desktop manifest 39 plugins/292 steps/882 params, 0 problems. 초기 z 없는 합성 방출/흡수선 회수, 모호한 대응 보존, 잡음/마스크/스파이크, 실제 owned job 및 저장/재개 검증.


### 2026-10-13 후속: AGN 방출선 공동 적합

AGN lines를 Web ds10_spectrum.joint / Desktop ds10core.spectrum-joint에 실제 연결했습니다. 공통 z에 묶인 여러 방출 Gaussian+국소 선형 연속광, 사용자 observed sigma bounds/fit-range, 결정적 3-start, 차원 정규화 퇴화 검사, 전체 parameter/성분 flux covariance, 판정·잔차·MASK·원본 저장. 이름 broad/narrow의 물리적 해석이나 host/Fe II 템플릿 분해를 대신하지 않습니다. docs/science_pipelines.md 마지막 절에 계약과 한계가 있습니다. 실제 관측 계산·커밋·푸시 없음.

검증(lageunha, CPU/BLAS 1 thread): 서버 66개·Desktop 14개·브라우저 6개, 총 86개 통과. Web typecheck/build 및 Desktop manifest(39 plugins/293 steps/885 params, 0 problems) 통과. 합성 겹침선/동일 중심 좁은+넓은 성분/작은 물리 flux 단위/마스크/잘못된 모형과 계약을 확인했습니다. 실제 owned job, Desktop 체크포인트, Web ZIP 저장→새로고침→동일 단계 재실행도 통과했습니다.


### 2026-10-13 후속: 사용자 host/Fe II 템플릿 분해

AGN decompose를 두 제품에 실행 연결했습니다. 같은 관측 TSV의 TEMPLATE_HOST/FEII, 출처·frame/grid/redshift/resolution template-note 필수. 고정 F_lambda index/pivot + nonnegative AGN/host/FeII 세기, 퇴화·경계 검사와 조건부 공분산. joint continuum-model=templates에서 template 세기와 선을 원본 flux에 함께 재적합하여 전체/선 flux 공분산에 반영합니다. 분해표·원본·JSON은 저장/재개합니다. 템플릿 모양/폭/AGN 지수 자체 탐색·블랙홀 질량은 남아 있습니다. 실제 관측 계산·다운로드·커밋·푸시 없음. 상세는 docs/science_pipelines.md 마지막 절.

검증(lageunha, CPU/BLAS 1 thread): 서버 75개·Desktop 15개·브라우저 7개, 총 97개 통과. Web typecheck/build, Desktop manifest(39 plugins/294 steps/889 params, 0 problems) 통과. 합성 계수/선 회수·작은 물리 flux 단위·퇴화/경계·입력 계약, owned job·체크포인트·ZIP 저장→새로고침→template joint 재실행 검증.


### 2026-10-13 후속: 공동 적합 선비 진단

AGN diagnostic 노드를 두 버전에 연결했습니다. 선행 lines의 joint.json을 자동 입력하며 활성 스펙트럼 유지. 성분 합계 ratio/log10와 전체 flux covariance를 전파한 오차·비 사이 covariance, S/N/품질/행렬 검증, source report snapshot/SHA 기록. Desktop 기록·복원 artifact 경로 연결, Web bound owned job/artifact 연결. 낮은 S/N/flagged 성분에 숫자나 상한을 만들지 않습니다. 자동 분류·소광·광도·블랙홀 질량은 남아 있습니다. 실제 관측 계산·커밋·푸시 없음. 상세는 docs/science_pipelines.md 마지막 절.

검증(lageunha, CPU/BLAS 1 thread): 서버 87개·Desktop 17개·브라우저 7개, 총 111개 통과. Web typecheck/build 및 Desktop manifest(39 plugins/295 steps/892 params, 0 problems) 통과. 선비·log 및 교차 covariance의 해석식/작은 물리 단위·잘못된 행렬/flags/SNR/소유권 검사, Desktop 복원 파일 연결·최신 실패 차단, Web 활성 스펙트럼 유지·ZIP 저장→새로고침→진단 재실행 검증.


### 2026-10-13 후속: 분류·소광·선 광도·블랙홀 질량 통합 해석

남아 있던 네 항목을 AGN inference 노드로 묶어 두 제품에 연결했습니다(11개 목적/111개 단계). 선행 lines/joint.json 자동 연결, Web은 popup 없이 단계 내 JSON 편집기, Desktop은 spectrum-inference-config. 명시한 좁은선 BPT와 모호 경계, 외부 또는 Balmer E(B-V)/사용자 k, 명시한 광도거리·cgs 변환의 선 광도, 확인된 넓은 Hα와 기기 Gaussian 폭의 Greene-Ho 2005 식6 조건부 질량을 계산합니다. 전체 공동 적합 covariance 및 독립 보조입력 오차를 전파하고 측정오차와 사용자 채택 scatter는 분리합니다. 생략은 not_configured, 부적절/부족한 입력은 blocked+이유이며 job 성공은 보고서 생성 성공입니다. science_validated=false/requires_review=true. 네 산출물과 원본 적합·설정은 체크포인트/ZIP 저장·복원에 포함됩니다. 상세 계약·합성 예시·1차 논문 링크는 docs/science_pipelines.md 마지막 절.

검증(lageunha, CPU/BLAS 1 thread): 서버 100개·Desktop 18개·실제 브라우저 8개, 총 126개 통과. Web typecheck/build, Desktop manifest(39 plugins/296 steps/894 params, 0 problems) 통과. 네 결과 값/교차 covariance, 질량오차의 독립 유한차분 검증, 경계 모호성/음의 소광/오류 가정/해상도 차단, owned API·Desktop subprocess와 체크포인트·Web inline 입력→네 결과→ZIP 저장→새로고침→재실행 검증. 실제 관측 계산·HSC 데이터 변경·커밋·푸시 없음. 원본 코드/허용 NumPy·SciPy만 사용, 추가 dependency/외부 템플릿 없음. 확정 AGN 판정/자동 선 식별, cosmology 자동선택, 소광법칙·host/Fe II 모양/virial factor 등 체계오차 posterior는 여전히 범위 밖입니다.


### 2026-10-13 통합 해석 추가 검증과 결함 수정

정상 회귀 126개 재통과 후 독립 합성 주입(z guess=0/.35, truth=guess+.0008, flux scale=1/1e-17)의 참값 회수 및 외부/Balmer 소광 질량오차 유한차분 검증을 추가했습니다. 운동학 보고서 변조 세 경우에서 질량을 차단하지 않는 결함을 발견하여 finite/양의 기준파장·유효 z·fitted shift와 common z의 일치 및 finite uncertainty 검사로 수정했습니다. 최종 서버 108개, Desktop 18개, 실제 브라우저 8개, 총 134개 통과(lageunha CPU/BLAS 1 thread). Web typecheck/build 및 Desktop manifest 0 problems. 범위·허용오차·결함 재현/수정·검증 한계는 docs/validation_spectrum_2026-10-13.md에 기록했습니다. 실제 관측 승인 검증·HSC 데이터 변경·커밋·푸시 없음.

커밋 준비 추가 검증: 커밋 대상 파일만 추출한 lageunha 체크아웃에서도 서버 126개·Desktop 18개·실제 브라우저 8개, 총 152개 및 Web typecheck/build 통과. 파이프라인 직접 의존 morphology/PSF stamp/counts 포함, 기존 런처·요금제·persona·회절 PSF 변경 제외. 상세와 기존 zero-noise 테스트 RuntimeWarning 3개는 검증 보고서 마지막 절.
