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

인증 조회 후 확인: catalogue_status는 161개 중 patch_qa 1개만 완료했고,
56,994행의 파일 해시는 일치했다. extinction_000(CAS ID 1059397)의 `.part`는
CAS CSV 첫 열의 `# ` 접두사 때문에 검사에서 중단됐다. 파서를 수정해 실제 파일
10개 대상/11행을 읽고 회귀 10건을 통과했다. 기존 `.part`는 완료로 가정하지 않고
동일 명령 재실행 시 기존 CAS ID에서 다시 받는다. patch_qa는 건너뛴다.
아직 소광 전체·random 목록은 완료되지 않았다. 상세 ancillary_README.md.

다음 준비에서 `validation_sample_2026-10-13_01`에 면적 조건 충족 40개와 마스크
stress 10개의 층화 예비 표본을 만들었다. 원본 800개는 유지한다. 활성 소광/대조장은
여전히 미입력이며 관측 측정을 실행하지 않았다. 표본·경로 독립 검증 및 선택 회귀
2건 통과. 사용자 지시로 숨김 암호를 실행당 1회 입력하는
`analysis/icl_wgl_membership/hsc_debug/fetch_ancillary_catalogues.py`를 추가했다.
161개 순차 CAS 쿼리의 결과/작업 ID/해시를 저장·재개하며 기존 HSC 영상 다운로더는
변경하지 않았다. 모의 CAS 회귀 5건 통과, 실제 인증 조회는 사용자 실행 대기다.
실행 명령과 결과 경로는 hsc_debug/ancillary_README.md에 있다.

후속 자료 상태: HSC-I 800개 다운로드 완료, 실제 파일 누락·크기 불일치 0,
42.24 GiB. `sample_800_2026-10-13/quality_2026-10-13`의 전체 CPU QA에서
구조·WCS 오류 0이다. 기존 마스크/고리별 면적 50% 기준으로 0–500 kpc 모든 고리를
통과하는 동일 표본 165개, 고리별 360–433개다. 800개를 모두 전 반경에 쓸 수 있다고
가정하지 않는다. `ancillary_2026-10-13_01`에는 SFD/i·i2 소광 후보 800행과 출처,
독립 조회 검증, 공식 random/patch_qa/a_i 조회 준비 파일이 있다. 소광 후보는 아직
활성 보정에 적용하지 않았다. 대조장 좌표·영상은 0이며 현재 환경에 HSC 인증이
남아 있지 않아 공식 CAS 조회는 미실행이다. 암호 재요청 없음.
`analysis/icl_wgl_membership/hsc_debug/ancillary_README.md` 및 실행 계획의
"후속: 다운로드 완료와 소광·대조장 준비"를 현재 기준으로 읽는다.

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

### 2026-10-13 사용자 AI 도구 Import 사전검사

Tools 메뉴에 AI submission kit 내보내기와 local Import review를 추가했다.
`tools/user_tool_import.py --kit <새 ZIP>` / `--review <제출폴더>` CLI도 제공한다.
Python 3.10 이상 필요, Web 서버나 sibling checkout 없이 동작한다. 검증기는
Web first-party `tools/validate_user_tool.py`의 사본이며 제출 코드를 실행하지
않는다. `tools/user_tool_kit/`에 규약·예제·런타임 정보를 포함했다. 내보내는
kit의 target은 standalone이다. 실행 승인/자동 설치/분석 단계 등록은 별도다.

검증: `python3 regression/test_user_tool_import.py` 4개(CLI, kit/덮어쓰기 거부,
검증기 동일성, Tcl 메뉴 명령 stub 실행). Desktop 바이너리를 재빌드하지 않았으므로
기존 실행본의 메뉴 사용은 Tcl 변경 반영 후 가능하다. 상세 `docs/user_tool_import.md`.
Web 측에서는 해당 기능을 열면 로컬 모드로 자동 전환한다. 커밋·푸시 없음.


### 2026-10-13 사용자 도구의 실행 위치 분리 (최신 지시)

이전 절의 전체 로컬 모드 자동 전환 설명을 대체한다. Web 사용자 도구 Import는
독립된 브라우저 Worker에서만 검사하며 전역 서버/로컬 모드, localStorage 모드,
일반 분석의 실행 위치를 변경하지 않는다. 모드 변경 이벤트를 발생시키지 않고
사용자의 전체 모드 변경에도 검사 Worker를 취소하지 않는다. 후속 등록·실행도
도구별 browser 실행 위치로 구현해야 한다. 현재는 사전검사까지이며 사용자 도구
자동 등록·실행 연결은 여전히 미구현이다. 규약과 standalone 동봉 사본도 갱신했다.

### 2026-10-13 개인 도구 native 실행·분석 단계 편입

앞선 사용자 도구 자동 등록/실행 미구현 설명을 대체한다. Tools의 native Python
등록으로 고정 manifest 해시 버전과 plugin.json을 생성하고, Tools의 selected
pipeline 메뉴로 현재 공식 분석 단계 뒤에 사용자 단계를 넣는다. 입력·매개변수
버튼에서 inputs-json / params-json을 설정한다. 각 실행은 별도 작업 폴더에서
해시를 재검사하고 실행하며 로그·보고서·선언 산출물을 기록한다.
시간 제한과 SIGTERM 시 자식 프로세스 그룹을 정리한다. 네이티브 Python은
사용자 권한으로 실행되며 OS 보안 샌드박스는 아니다. 자동 의존성 설치는 없다.

연구 세션에는 사용자 단계 정의, 원본 코드/manifest, 결과와 실제 실행 입력
사본을 보관한다. 복원은 코드 검증/재등록과 입력 경로 재연결만 수행하며 실행하지
않는다. 등록 버전이 바뀌면 기존 노드는 실행을 막는다. 기존 공식 recipe 파일은
수정하지 않고 사용자 세션 안에서만 순서도를 확장한다.
사용법: docs/user_tool_execution.md. 검사: regression/test_personal_pipeline.py,
regression/test_user_tool_import.py, plugins/workflow_designer/test_workflow.py의
27건이 lageunha CPU에서 통과했다. Web 관련 실제 CSP E2E 10건도 통과했다.
Desktop 바이너리는 재빌드하지 않았으며 Tcl/CLI 소스 수준에서 검증했다.
실제 과학 자료 계산, 커밋·푸시는 하지 않았다.


### 2026-10-13 HSC 50장 검출·마스크·PSF 입력 준비 완료

lageunha CPU에서 검증 표본 50/50장 처리, 오류 0. 검출 156,984개, 임시 PSF 후보
569개를 저장했다. 면적 집합 40장과 마스크 스트레스 10장은 분리 유지한다.
BCG 위치 1 arcsec 내 고유 검출 11장, 검출 없음 39장이므로 후자는 BCG 보존 미확정으로
차단했다. 원본 품질 마스크의 BCG 중심 가림은 18장이다. PSF 후보 4개 미만은
3806(1), 6940(3), 3700(2)이며 대체 Gaussian을 생성하지 않았다.
Kron 2.5/3.5/5.0 마스크, 검출/멤버 표, 원본 PSF 우표와 방사 진단, 전체/BCG 확대
이미지를 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/detection_psf_2026-10-13_01`
에 저장했다. 합성 회귀 3건과 독립 FITS/표/WCS/마스크 산출물 검사 50건 통과.
`review_queue.tsv`에 후속 검토를 기록했다. 과학 승인 마스크·PSF가 아니며
`science_measurement=false`, `mask_ready_for_science=false`다. 원본 영상 변경,
소광 적용, 대조장 차감, 관측 스태킹은 하지 않았다.
상세: `analysis/icl_wgl_membership/hsc_debug/detection_psf_README.md`.
동시 진행 중인 사용자 CAS 수집은 이번 확인 시 running/random_008이었다.
앞선 인증 조회 미실행 기록은 당시 상태이며 현재 상태로 해석하지 않는다.
커밋·푸시는 하지 않았다.


### 2026-10-13 후속: BCG 정점 대조와 멤버 마스크 연결 보완

lageunha CPU에서 50장 검출을 재현하고 deblend contrast 0.005/0.001의 원본·필터
정점과 중심 분할 번호를 대조했다. 두 설정에서 안정된 위치 후보는 35장으로,
기존 centroid 후보 11장을 유지하고 24장을 추가했다. 남은 15장은 중심 품질 마스크
8장과 정점 조건 불충족 7장이다. 고유 BCG 멤버가 미검출 위성으로 재마스킹되는
연결을 보완했으며 실제 이웃 마스크와 원본 품질 마스크는 보존했다.

산출물은 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/bcg_review_2026-10-13_01`
및 `bcg_masks_2026-10-13_01`이다. 후자가 멤버 처리까지 보완한 검토용 마스크다.
Kron 2.5/3.5/5.0에서도 중심이 가려진 영상은 31/32/35장이므로 후보 35장을
사용 가능 과학 표본으로 해석하지 않는다. 모든 과학 승인 상태는 false로 유지했다.
PSF 후보 569개의 외부 대조 좌표를 준비했으나 별 분류 확인은 아직 수행하지 않았다.
회귀 8건과 독립 저장 산출물 검사(마스크 50장, PSF 좌표 569개) 통과.
기존 결과·원본·활성 manifest를 보존했고 실제 스태킹·새 인증 요청·커밋·푸시는 없다.
상세는 `analysis/icl_wgl_membership/hsc_debug/bcg_review_README.md`.
후속은 중심부 이웃 분리/마스크 유효면적 검토, 외부 별 분류 및 넓은 PSF 날개,
주입 회수와 공식 random에 기반한 대조장 검증이다.


### 2026-10-13 후속: 마스크 면적/분리 민감도와 Gaia PSF 후보 대조

lageunha CPU에서 50장의 0–500 kpc, 25 kpc 고리·8방향 면적을 계산했다.
품질 마스크와 분리 2종×Kron 3종의 7,000행, 기존 마스크 200개의 독립 재계수
검사 통과. 신규 합성 테스트 5건 통과. 기준 contrast .005/Kron 2.5에서 BCG 후보와
모든 고리 유효면적≥50%를 동시에 만족한 것은 5724·619·5814 세 장, contrast .001에서는
619·5814 두 장이다. Kron 3.5/5.0 및 스트레스 집합은 통과 0이다. 과학 승인은 아니다.
5724는 BCG 정점/모멘트가 같아도 이웃 검출 차이로 최대 면적 차이 약 72.7%p 발생.
마스크를 임의로 축소하거나 중심을 강제로 비우지 않았으며, 이웃 분리/광량 모형과
주입 회수 검증을 다음 단계로 유지한다.

공개 CDS Gaia DR3 위치 대조 후 ESA source_id 조회로 외부 정보를 확보했다.
569개 후보 중 대응 73개, 보수적 별 근거 충족 32개(23개 영상), 근거 불충분 41개,
미대응 496개다. 입력 좌표/ID·응답 해시·천구 거리 독립 검사를 통과했다.
PSF 승인과 HSC epoch 전파는 미완료다. 지연된 ESA 공간 결합은 ABORTED로 정리했다.
HSC 별도 CAS 조회 161개도 모두 완료: 소광 800개 대상 1,212행, random 후보
800개 대상 160,000행, patch QA 56,994행. 완료 파일 감사 기록을 저장했다.

결과는 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/` 아래
`mask_coverage_2026-10-13_02`, `mask_case_5724_2026-10-13_01`,
`psf_cds_2026-10-13_01/evidence`다. 상세와 그림·검증 범위는
`analysis/icl_wgl_membership/hsc_debug/mask_psf_review_README.md`.
원본·활성 manifest 보존. 과학 스태킹·대조 영상 다운로드·재인증·커밋·푸시 없음.


### 2026-10-13 후속: 이웃 외곽·확산광의 주입/회수 예비 실험

lageunha CPU에서 실제 HSC 품질 마스크/분산/멤버 위치 6개 조건을 사용한 합성
432장면×5방법=2160개 회수 기록을 완료했다. science 픽셀은 합성이고 Gaussian
PSF·고정 기하/형태·독립 잡음을 가정했다. 실제 관측 및 upstream 하늘 전달 검증이 아니다.
결과는 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/neighbor_injection_2026-10-13_01`.
단위 검사 5건, 입력/격자 감사, 독립 열 정규화 SVD 재계산 48건 통과.

면적 집합 5개 조건의 이웃 중심+외곽 공동 모형에서 확산광 상대 오차 중앙값은
같은 모형 대조 −3.82%, 외곽 크기 오차 +11.38%, 합성 PSF 날개 누락 +3.39%,
하늘 곡률 누락 −51.65%다. 이 값은 관측 편향 보정값이 아니다. 스트레스 6940은
따로 보고했고 외곽광 형식 오차가 최대 주입량의 약 33배여서 과학 회수로 승인하지 않는다.
마스크 방식은 조건별 108장면 중 41–45개만 성공했다. 알려진 위치를 사용하는 모형
방식과 성공 사례만 비교해 우월성을 주장하지 않는다. 완전한 진폭 격자만 전달계수 계산.
독립 주입 성분의 오차에도 반상관(면적 집합 약 −0.281)이 생겨 공동 공분산 처리가 필요하다.

외곽 포함 모형은 후속 후보지만 하늘 잔차/외곽 형태/경험 PSF/상관잡음 검증은 남아 있다.
science_measurement=false, observational_validation_passed=false 유지. 원본·활성
manifest 보존, 실제 스태킹·추가 인증·커밋·푸시 없음. 상세:
`analysis/icl_wgl_membership/hsc_debug/neighbor_injection_README.md`.


### 2026-10-13 후속: 외곽 크기·하늘 제약 공동 검증

lageunha CPU/1 thread에서 이전 합성 실험과 같은 입력·잡음으로 864장면,
5방법=4320개 회수 기록, 크기별 적합 25,056건을 완료했다(약 19.2초).
4–16 kpc 공통 외곽 크기 격자, 평면/이차 하늘, 독립 합성 하늘 제약 및
편향된 제약을 비교했다. 실제 대조장 측정이나 경험적 PSF를 사용한 실험은 아니다.

면적 집합 5개 조건에서 외곽 12 kpc+하늘 곡률 동반 시 확산광 상대 오차 중앙값은
기존 고정 크기/평면 −34.42%, 자유 크기/이차 하늘 −4.63%, 합성 외부 제약 −4.02%다.
각각의 절대 상대 오차 중앙값은 34.42%, 17.00%, 9.21%이므로 장면별 정확도 4%를
주장하지 않는다. 편향된 제약은 +7.78%, 모형에 없는 국소 하늘 잔차는 합성 제약에도
+27.96%가 남았다. PSF 날개 오차도 미해결이다. 6940의 720개 방법 기록은
수치적으로 성공해도 가려진 성분의 정보 부족으로 별도 차단 표시했다.

이전 기준 432건 재현(최대 차이 1.93e-10), 독립 SVD 검산 60건(최대 3.23e-9),
신규 6건 포함 테스트 11건과 해시/격자/실패/profile 구간 감사 통과.
결과: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/halo_sky_validation_2026-10-13_01`.
상세: `analysis/icl_wgl_membership/hsc_debug/halo_sky_validation_README.md`.
크기 모호성을 반영한 구간도 보정된 신뢰구간은 아니다. 실제 대조장 하늘 공간 변동,
경험적 PSF/상관잡음과 추가 독립 주입 실험을 다음 조건으로 유지한다.
science_measurement=false, observational_validation_passed=false. 원본/활성 manifest
보존, 관측 스태킹·추가 인증·다운로드·커밋·푸시는 하지 않았다.


### 2026-10-14 후속: 잠정 대조장 좌표와 실제 외곽 잔차

lageunha CPU/1 thread에서 공식 random 160,000행(고유 ID 109,402개), patch QA,
기존 SFD 지도와 전체 보유 CAMIRA 7,319개를 사용했다. 은하단 배제 반경
0.5/1/2 Mpc에서 각각 대상 787개의 비중첩 잠정 대조 좌표를 배정했다. 기본 1 Mpc
대비 좌표 교체는 316/505개다. 13개 대상은 5각분 안의 유효 QA 대용점이 없어 남겼다.
대상 patch는 근방 random의 QA 대용값이며 정확한 patch 식별은 미확인이다.
전체 대조 영상 footprint/마스크/필터/PSF 검증과 대조 영상 다운로드는 아직 없다.
`request_plan.json`은 잠정 좌표 계획이고 과학용 또는 실행용 manifest가 아니다.

기존 검증 영상 50장(면적 집합 40, 스트레스 10)의 600–900 kpc 외곽을
Kron 2.5/3.5/5.0으로 검사했다. 150개 진단, 타일 26,906개, 상관 기록 1,200개다.
면적 집합의 타일 평면 잔차 MAD/명목 평균 오차 중앙값은 2.79/2.58/2.41,
스트레스 집합은 3.15/2.85/2.57이다. 면적 집합 외곽 유효면적은 44.7/37.3/26.0%다.
1픽셀 상관은 약 .054–.055. 은하단 빛/PSF/미검출 소스/마스크 선택도 포함되므로
하늘 보정값, 순수 잡음 커널 또는 최종 오차 배율로 편입하지 않는다.

독립 검산에서 기존 거리 격자 보간과 Astropy의 약 .002–.005% 차이로 고리 경계
픽셀 불일치를 확인했다. 이번 외곽 진단을 FlatLambdaCDM으로 통일한 `_02`가 최종이다.
첫 결과와 공용 코어/과거 주입 결과는 보존했다. 신규 테스트 8건, 입력 해시 190건,
반경별 787개 대조 좌표의 전체 CAMIRA/상호 겹침 검사, 원본 타일 독립 재계산
54건이 통과했고 타일 평균 최대 차이는 0이다.

자료 상위 경로 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13` 아래:
`control_candidates_2026-10-14_01`, `outer_sky_diagnostics_2026-10-14_02`,
`control_sky_audit_2026-10-14_02`. 상세는
`analysis/icl_wgl_membership/hsc_debug/control_sky_README.md`.
원본/활성 manifest 보존, science_measurement=false. 다음은 정확한 patch 연결,
대조장 전체 영상 검증, 경험적 PSF/배경 공분산 및 같은 타일의 마스크 비교다.
추가 인증·새 영상 다운로드·관측 스태킹·커밋·푸시는 하지 않았다.


### 2026-10-14 후속: 공식 patch 연결과 경험적 PSF 진단

lageunha CPU/1 thread에서 공개 tract·patch 다각형 91,449개를 사용해 근방 QA
대용값을 교체했다. 대상 800개 중 단일 유효 patch 607개, 복수 patch 192개,
미연결 1개다. 복수 patch 모두의 화질 조건을 만족하도록 해 비중첩 대조 좌표
798개를 선정했다. 1497은 공개 경계 미연결, 6470은 겹친 네 patch와 동시에
맞는 후보가 없어 남겼다. 실제 primary tract 및 전체 cutout 품질은 미확인이다.

`control_download_manifest.json`은 원래 성공한 `download_pilot.py`를 그대로
사용하는 실행 입력이다. matched_control 역할·대응 대상 ID·자체 좌표·BCG 없음과
미검증 소광을 보존한다. 예상 42.10 GiB이고 실제 새 대조 영상은 0개다.
현재 프로세스에 HSC 인증이 남아 있지 않아 전송하지 않았으며 인증을 재요청하지 않았다.

Gaia 별 근거 32개/23영상에 401×401 native 우표와 마스크 2종×배경 2종을
검사했다. 우표 20개 확보, 12개 경계 실패. 128개 방법 기록 중 진단 가능 30개,
배경 면적 부족 49개, 영상 경계 48개, 중심부 마스크 1개다. 일부 조건 통과 별
11개, 네 조건 모두 통과 별 4개다. 바깥 10.75–16.13초각 고리 면적 통과 기록
25개 중 23개가 배경 차감 후 0 이하다. 중심 광량 안정성을 PSF 날개 승인으로
대신하지 않았고 음수·실패 사례를 모두 보존했다.

신규 테스트 8건, 입력 해시 80건, 접평면 다각형 재검산 2,840건, 798개 전체
CAMIRA/대조장 겹침 검사, 원본 우표 20개·방사 구간 330개·마스크 표본점
3,240건의 독립 재계산을 통과했다. 원본 픽셀 평균 최대 차이는 0이다.

자료는 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13` 아래
`control_patches_2026-10-14_02`, `empirical_psf_2026-10-14_01`,
`patch_psf_audit_2026-10-14_02`. 상세는
`analysis/icl_wgl_membership/hsc_debug/patch_psf_README.md`.
빈 QA 중심을 처리하지 못한 최초 patch 실행과 마스크 HDU 참조 키를 수정하기 전
감사 `_01`은 보존했다. 원본·활성 manifest 보존, science_measurement=false,
psf_approved=false. 다음은 실제 대조 영상의 배경/공분산과 PSF 입력·epoch 검증이다.


### 2026-10-14 후속: 공개 확장 PSF 실제 확보·전달 시험

lageunha CPU에서 Garate-Nuñez et al. (2024)의 공개 i-band PSF 대칭/비대칭 FITS를
실제 내려받았다. Zenodo doi:10.5281/zenodo.14915578, CC-BY-4.0. 두 파일 합계
19,618,560 bytes(18.71 MiB), 배포 MD5와 SHA-256을 확인했다. 4001×4001 배열의
중심·유한성·단위 합을 검사했고 원본 음수도 보존했다. 외부 R 코드 실행 없음.

동일한 마스크와 배경 연산으로 기존 관측 성공 30조건×참조 2종, 방사 구간
660개를 비교했다. 참조 날개의 10.75–16.13초각 평균에서 배경 추정이 제거하는
비율은 10.32–28.39%(중앙값 26.54%)다. 남는 참조 신호는 명목 고리 평균 오차의
.210–.709배다. 관측 23/25개 방법의 0 이하 값은 그대로 두었다. 이 자료로
개별 PSF 날개를 승인하거나 공개 평균을 800장 전체의 참 PSF로 고정하지 않는다.

원본에 참조 신호를 더한 1,980개 차분 기록이 선형 연산과 일치했다.
분산·마스크 고정, photon noise 추가 없음이므로 HSC upstream sky나 재검출
완전성 검증은 아니다. 신규 6개·기존 4개 테스트, 입력 해시 31건, 독립 방사
비교 660개/신호 증분 1,980개 검산 통과. 최대 독립 증분 오차 5.68e-14다.

기존 800장 위치/크기 선별 후 관련 24개 헤더에서 별 좌표 33조합을 WCS로
확인했다. 경계 때문에 잘린 별 12개를 대체할 보유 영상은 없다. 날짜/MJD 키는
0개이고 PSF archive ID는 있어도 모델 본체는 없다. 최초 inventory `_01`의
상대 영상 파일명 해석을 고쳐 `_02`에 실행했으며 첫 폴더를 보존했다.

자료 상위 경로 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13` 아래
`extended_psf_reference_2026-10-14_01`, `extended_psf_comparison_2026-10-14_01`,
`psf_support_2026-10-14_02`, `extended_psf_audit_2026-10-14_01`.
설명은 `analysis/icl_wgl_membership/hsc_debug/extended_psf_README.md`.
대조장 798개 입력은 준비 상태이고 실제 대조 영상은 인증 세션이 없어 0개다.
science_measurement=false, psf_approved=false, reference_approved=false 유지.


### 2026-10-14 사용자 다운로드 승인과 실제 전송 재개

사용자는 새 자료 다운로드 시 대상·수량·예상 용량을 제시하고 명시적 승인을
받도록 지시했다. 이후 **HSC 대조장 남은 794장, 예상 약 41.9 GiB의 다운로드를
명시적으로 승인**했다. 승인한 동일 묶음의 전송/재개에 이 승인을 유지하며,
다른 자료나 수량 확대는 별도 승인을 받는다. 암호 자체는 이 문서에 기록하지 않는다.

이전에 대화로 제공된 인증정보를 실제로 사용한 결과 인증에 성공했고 대조장
4장의 FITS·헤더 검증을 완료했다(240,894,720 bytes). 새 승인 조건을 확인하는
동안 전송을 중단했다가, 명시적 승인 후 기존 `download_pilot.py --resume`으로
lageunha CPU에서 남은 794장을 재개했다. 중단된 다섯 번째 응답은 보존하고 재수신한다.

입력: `sample_800_2026-10-13/control_patches_2026-10-14_02/control_download_manifest.json`.
실행: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/control_download_01`.
실시간 상태는 `download_status.json`, 전체 완료 후에는 `download_manifest.json`이다.
앞 절의 "인증 세션이 없어 0개"는 재개 전 기록이며 현재 상태는 이 실행 파일로 확인한다.
인증정보는 실행 중인 프로세스의 숨김 입력으로만 전달했고 파일에 저장하지 않았다.

### 2026-10-14 대조장 798장 다운로드 완료 확인

기존 downloader가 798/798장 완료 후 정상 종료했다. 2026-10-14 17:23 KST에
승인 입력 SHA-256, 입력/상태/출력 manifest의 전체 고유 ID, 실제 파일 존재·크기,
control 역할·대상 ID·좌표·각척도 보존을 대조했다. 누락·크기 불일치·역할 불일치 0건.
실제 합계는 45,254,626,560 bytes(42.1467 GiB)다. 다운로드 당시 FITS 구조·WCS
검증 기록은 798건이며 헤더 필터는 i 312개/i2 486개다.

완료 자료는 `sample_800_2026-10-13/control_download_01`의 `images/`,
`download_manifest.json`, `download_status.json`이다. 확인 보고서는
`confirmation_20261014T082353Z.json`. 이 확인에서 45 GB 전체 해시를 다시 읽지는
않았으며 다운로드 중 기록한 SHA-256을 보존한다. 전체 픽셀·품질 마스크·유효면적·
대조장 배경/공분산과 실제 passband에 따른 소광 검증은 다음 단계다.
science_measurement=false, credentials_saved=false 유지.

### 2026-10-14 대조장 전체 품질검사와 대응 외곽 배경 비교

승인된 다운로드의 후속 분석을 lageunha CPU에서 완료했다. 추가 다운로드·GPU
사용은 없다. 기존 `hsc_quality` 픽셀 검사와 SEP/외곽 타일 함수를 재사용하고,
대조장 중심에는 BCG 보존 예외를 두지 않는 검증 스크립트를 추가했다.
상세 방법·수치·실행 명령은 `analysis/icl_wgl_membership/hsc_debug/control_images_README.md`다.

`sample_800_2026-10-13/control_quality_2026-10-14_01`:
798/798장의 45,254,626,560 bytes 전체 SHA-256을 다시 읽어 다운로드 기록과
대조했고 구조·WCS 오류 0, 비유한 science 픽셀을 가진 영상 0이다.
26장의 일부 분산 픽셀 1,205개를 제외했다. random 중심 BRIGHT_OBJECT는 0개다.
전체 품질 유효면적 중앙값 47.44%, 밝은별 마스크 면적 중앙값 42.13%다.
품질 마스크만 적용한 모든 0–500 kpc 고리의 50% 유효면적 문턱을 통과한
동일 대조장 표본은 234개이며, 고리별 수는 716개에서 338개까지 달라진다.
이는 천체 확장 마스크 적용 전 기술 진단이며 최종 과학 표본 선정이 아니다.

`control_outer_2026-10-14_01`:
기존 면적 집합 40개·마스크 스트레스 10개의 은하단과 대응 대조장 50쌍을
600–900 kpc에서 비교했다. Kron 2.5/3.5/5.0의 150개 방법 행,
대조장 타일 25,229개, 간격/방향별 상관 1,200행을 완료했다. 50쌍 모두
헤더 i/i2가 일치하고 각 방법에서 양쪽 유효 타일 10개 이상이다.
Kron 3.5의 타일 잔차 MAD/명목 평균 오차 중앙값은 면적 집합 은하단/대조장
2.577/2.687, 스트레스 집합 2.854/2.826이다. 유효면적 중앙값은 각각
37.31%/34.18%, 35.68%/28.03%다. 이는 순수 잡음 배율이나 ICL 신호가 아니다.
마스크 크기별 살아남는 타일도 달라지므로 같은 유효 영역에서 후속 비교가 필요하다.

`control_image_audit_2026-10-14_01/audit.json`:
입력 해시 108건·코드 해시 5건, 반경별 면적 1,000건, 마스크 150개,
WCS 표본점 600개, 원본 타일 450건, 쌍별 광량 환산 150건 검산이 통과했다.
타일 평균 최대 차이 0. 새 테스트 4개와 기존 회귀 7개 총 11개도 통과했다.
비교 그림은 `control_outer_2026-10-14_01/paired_outer_diagnostics.png`다.

다음은 이 50쌍에서 같은 유효 영역의 마스크 민감도, 실제 대조장 신호 주입·회수,
배경 공분산을 검증하는 단계다. PSF/소광/유효 passband 검증도 남아 있다.
원본·기존 결과·활성 science manifest는 보존하며 science_measurement=false,
controls_ready_for_science=false, sky_correction_approved=false를 유지한다.

### 2026-10-15 대조장 비교 확대 재검증

사용자의 재검증 요청에 따라 lageunha CPU에서 기존 테스트 11개와 원본 감사
전체를 다시 실행해 통과했다. 감사 결과는
`sample_800_2026-10-13/control_image_audit_2026-10-15_01/audit.json`이다.
같은 폴더의 `summary_recheck.json`에서 품질 집계와 그룹별 통계 60개도 확인했다.
234개 유효면적 조건 통과, 26장의 분산 문제 픽셀 1,205개 등 보고 수치는 그대로다.

기존 감사의 일부 타일 표본 검사에 더해, 별도 코드
`analysis/icl_wgl_membership/hsc_debug/recheck_control_statistics.py`로
대조장 798장 전체 45,254,626,560 bytes 해시와 양쪽 영상 50쌍의 모든 통계를
재계산했다. 기존 타일 통계 함수를 호출하지 않고 가중 정규방정식과 영상 전체의
평균 차감 픽셀쌍 공분산을 사용했다. 타일 경계 밖의 쌍은 명시적으로 제외한다.

`control_full_recheck_2026-10-15_01/verification.json`:
타일 52,135개·평면 적합 300건·상관 2,400건이 통과했다.
타일 평균·중앙값·명목 평균 오차 차이는 0이고, 잔차/명목 오차 비율의 최대
절대 차이는 1.42e-14, 상관계수는 1.11e-16으로 반올림 수준이다.
별도 계산식 테스트 4개를 포함해 총 15개 테스트가 통과했다.
원본·기존 결과·과학용 활성 상태는 유지했고 추가 다운로드는 없다.
이번 단계는 수치 재현성 검증이며 신호 주입·회수와 과학용 공분산 검증은 다음 단계다.

### 2026-10-15 실제 영상 주입·회수와 스태킹 전 공분산 검증

사용자 지시에 따라 lageunha CPU에서 대조장 798개 전체를 재검출하고
Kron 2.5/3.5/5.0 마스크, 50–500 kpc 여섯 구간 프로필을 추출했다.
기존 50개 검증 대조장에 원형 지수 h=100 kpc와 타원 h=250 kpc/q=0.65,
100 kpc의 관측 SB 1/3/10 nJy/arcsec²를 주입해 총 300장을 재검출했다.
이 단계는 코애드 주입이며 추가 photon noise·HSC 이전 처리·PSF 날개 전달을
검증한 것은 아니다. 기존 은하단 50개의 마스크 기반 프로필도 별도 추출했다.
추출/주입 시간 약 1,254초, 오류 0, 회수 행 10,800개, 집계 432개다.

기준/방법/결과는 `analysis/icl_wgl_membership/hsc_debug/prestack_validation_README.md`다.
Kron 3.5에서 면적 조건 집합 40개의 50–300 kpc는 모형별 유효 13–25개로,
원본 고리 추출 24개 집계의 회수 편향 중앙값이 −3.31%~+3.39%였다.
300–500 kpc의 유효 영상은 6–8개라 미승인이다. 마스크 스트레스 집합은
구간별 0–4개여서 전체가 표본 부족이다. 전체 고리 참값과 열린 면적 참값을
구분해 보조 회수비도 저장했으며 원래 기하학적 고리 기준 판정은 유지했다.

외곽 평면을 빼면 넓은 모형의 200–300 kpc에서 A=1/3 주입에 각각
−16.44%/−10.02% 편향이 생겼다. 평면 차감은 민감도 시험으로만 남기고
기존 원본 고리/대조장 기반 분석에 새 과학용 보정으로 편입하지 않았다.

여섯 구간 모두 면적 50%를 충족한 대조장은 Kron 2.5/3.5/5.0에서 140/40/0개다.
Kron 3.5 대조장 공분산은 40개·28개 5° 하늘 셀, rank=6, 조건수 약 84.1이다.
독립 bootstrap 대각은 S/N 대비 최대 4.72% 차이, 공간 hold-out의 1.96σ 포함률
80%로 정의한 탐색적 기술 기준을 통과했다. 공간 bootstrap 표준오차는 독립
영상 가정보다 약 7–17% 크다. 원점/격자 크기·필터·각척도 구간 민감도도 저장했다.
이는 대조장/잔차 성분의 조건부 공분산이며 순수 잡음 배율이나 전체 과학 오차가 아니다.

기존 은하단–대조장 50쌍 중 양쪽 모두 여섯 구간을 지원하는 것은 Kron 3.5에서
0쌍이다. 구간별로도 8/6/8/7/4/2쌍이다. Kron 2.5의 공통 3쌍은 rank=2여서
6차원 전체 스태킹 공분산에 부족하다. 대조장 공분산을 전체 오차로 대체하지 않는다.
**원본 추출의 지원 집합에서 50–300 kpc 회수는 확인했지만 전체 과학 스태킹은 보류**다.
다음은 양쪽에 공통 면적을 갖는 짝 표본과 마스크/반경 정의 검증이다.

실행 폴더는 `sample_800_2026-10-13` 아래
`prestack_extraction_2026-10-15_01`, `prestack_clusters_2026-10-15_01`,
`prestack_validation_2026-10-15_01`, `prestack_audit_2026-10-15_01`이다.
검토 그림은 `prestack_validation_2026-10-15_01/prestack_validation_review.png`다.

새 테스트 6개(독립 잡음 500회/공간 공유 오프셋 포함)가 통과했다.
독립 감사도 입력·산출물 해시 1,896건, 기존 마스크 150개, 픽셀 프로필 값 1,512개,
공통 면적 증분 1,296개·고정 면적 모형 1,296개, 공분산 84개 설정을 통과했다.
원본 재계산 최대 차이는 2.51e-13 nJy/arcsec² 미만이다.
감사 폴더에 해시 이름의 코드 사본을 보존했다. 추가 다운로드·GPU 사용은 없고,
원본·기존 결과·science_measurement=false·sky_correction_approved=false를 유지했다.

### 2026-10-16 matched pair 공통 면적 민감도

기존 은하단·대조장 50쌍의 저장된 프로필에서 Kron 2.5/3.5/5.0, 면적 바닥
25/35/40/50%, 50–150/200/300/500 kpc와 개별 고리를 조합해 120개 공분산 지원
설정을 계산했다. 신호 주입 회수 민감도 72개 그룹도 재집계했다. 코드는
`analysis/icl_wgl_membership/hsc_debug/evaluate_prestack_pair_support.py`이고,
상세 설명은 `analysis/icl_wgl_membership/hsc_debug/prestack_validation_README.md`다.

결과는 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_pair_support_2026-10-16_02/`다.
원래 Kron 3.5·50% 기준에서는 양쪽이 50–300 kpc 네 구간을 모두 만족하는 쌍이
1개, 50–500 kpc는 0개였다. 탐색 조건 Kron 2.5·40%는 50–300 kpc 24쌍
(모두 area-qualified), 21개 공동 5° 블록, rank 4를 남겼다. bootstrap 대각
차이 6.8% 이내, 공간 hold-out 포함률 100%로 기술 공분산 진단은 통과했다.
공동 블록 표준오차는 [11.30, 2.80, 1.34, 0.892] nJy/arcsec²이며,
전체 800개 표본의 오차나 순수 픽셀 잡음은 아니다.

같은 Kron 2.5·40% 주입 민감도에서 각 6개 모형이 28개 필드의 네 구간을 모두
지원했고 중앙 회수 편향 절대값은 최대 5.90%다. 이것은 coadd smooth-signal
주입의 회수이며 실제 미검출/잔류 천체 외곽빛은 검증하지 않았다. Kron 2.5/40%
조건은 표본 수를 확보하는 후보일 뿐 과학용 측정 정의로 채택하지 않았다.
기존 50개는 전체 800개 대표 표본이 아니고 PSF·상류 sky·i/i2 소광/응답·위성–ICL
광량 누출이 미승인이다. `science_measurement=false`,
`sky_correction_approved=false`를 유지한다. 전체 matched sample의 잔류 천체광과
면적 민감도 평가 후 본 분석 반경/마스크를 고정한다. 다운로드·GPU 사용 없음.

### 2026-10-16 전체 표본 자동 마스크·paired 면적 지원 계산 완료

전체 `download_01/download_manifest.json` 800개와 798개 matched control 연결을
확인했다. 기존 BCG 검토 마스크는 50개에 한정되어 있어 새 스크립트
`analysis/icl_wgl_membership/hsc_debug/prepare_full_cluster_mask_support.py`로
800개 source-mask 후보와 50–500 kpc 구간별 면적을 계산한다. 이 작업은 밝기
프로필을 만들지 않고 자동 중심 검출 후보를 BCG 확정으로 취급하지 않는다.

Slurm job 424800은 800/800 성공·오류 0·46:03 경과로 완료했다. 자동 중심 후보
258개, 중심 품질 마스크 391개, Kron 2.5/3.5/5.0에서 가려진 중심 671/696/725개다.
후보 중 품질/천체 마스크가 중심에 걸리지 않은 수는 129/104/75개다. 이는 자동
위치 후보이며 BCG 검증이나 측광용 표본 승인이 아니다. cluster-only 50–500 kpc
지원 수는 Kron 2.5에서 면적 25/35/40/50% 기준 331/186/131/40개,
Kron 3.5에서 171/50/27/4개다. paired 면적 집계 job 424865도 완료했다.
전체 798개 matched control 중 2개 대상은 대조장이 없어 분석 연결에서 빠졌다.
탐색 조건 Kron 2.5·40%, 50–300 kpc는 공통 면적 138쌍(동일 필터 쌍 137; 대상
필터 i 44/i2 94), 공동 5° 공간 블록 44개로 표본 수 스크린을 통과했다.
같은 마스크/바닥에서 50–500 kpc는 69쌍(동일 필터 68), 블록 37개로 수 스크린을 통과했지만,
그 반경의 주입 회수와 광도 공분산 검증은 완료되지 않았다. Kron 2.5·50%의
50–300 kpc는 41쌍/30블록으로 스크린을 통과하나 50–500 kpc는 11쌍/10블록으로
실패했다. 기준 Kron 3.5·50%는 50–300 kpc 4쌍, 50–500 kpc 1쌍뿐이다.
이 표는 면적 지원만 계산하며 flux covariance, 광도 프로필, 영상 스택을 만들지
않았다. 자동 중심 후보도 BCG 판정이 아니다. 모든 조건에서
`flux_covariance_available=false`, `science_measurement=false`다.
결과는 `prestack_full_pair_area_support_2026-10-16_01`에 있다. 작업 스크립트는
`run_full_cluster_mask_support.sh`, `run_full_pair_area_support.sh`이며 로그는
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_mask_support_logs_2026-10-16_01/`다.
출력은 `prestack_cluster_masks_2026-10-16_01`과
`prestack_full_pair_area_support_2026-10-16_01`이다. mask summary에 입력/출력 해시
800개를 보존했다. 두 작업 모두 GPU 코드를 사용하지 않고 과학용 마스크/측정을
승인하지 않는다.


### 2026-10-20 후속: 전체 Gaia PSF 후보 진단과 weak-lensing 좌표 반출 상태

기존 Gaia DR3 별 근거 후보 480개에 대해 `vari_summary`와
`vari_classifier_result`를 조회했다. 두 표에서 일치 source 행은 0개였으며,
이는 비변광 판정이 아니라 Gaia DR3 변광 처리 기록이 없다는 뜻이다. 전체
HSC 진단은 288개 프레임, 333개 401×401 우표, 1,920개 방법 조합을 기록했다.
137개 후보가 적어도 한 조합의 유효 프로필을 냈고 42개는 네 조합 모두에서
중심/배경 분석 조건을 통과했다. 64–96픽셀 외곽 고리에서는 기하 지원을 통과한
방법 프로필 274개 중 118개 양수·156개 0 이하였고, 후보별 묶음 113개의 지원
설정 부호는 전체 양수 38, 전체 0 이하 58, 혼합 17이었다. 외곽 PSF 날개는
검출/승인되지 않았다. Gaia 변광 분석 자료 부재, 정확한 epoch별 Gaia↔HSC 후보
식별, 영상별 포화·고립도 및 프레임 전달 검증이 남아 `science_measurement=false`,
`psf_approved=false`다.

Gaia DR3 `ref_epoch=2016.0` proper motion을 구면 접평면에서 ±10년까지 전파하는
민감도 분석도 수행했다. HSC 후보 좌표를 고정한 시나리오에서 480개 모두 기존
2″ 연결반경 안에 남았고 최대 분리는 1.572″였다. HSC coadd의 실제 epoch를
대입하지 않은 시나리오이므로 실제 전파나 천체/PSF 검증이 아니다. 결과와
재현기는 `empirical_psf_gaia_epoch_sensitivity_2026-10-20_01/` 및
`summarize_psf_gaia_epoch_sensitivity.py`다.

재현 스크립트와 수치는 `analysis/icl_wgl_membership/hsc_debug/`의
`diagnose_empirical_psf.py`, `query_psf_gaia_variability.py`,
`summarize_empirical_psf_star_level.py`, `extended_psf_README.md`, 통합 연구평가에
있다. 원자료 기반 진단과 Gaia 조회 산출물은 Git 밖의 `/gpfs/kjhan/Astrafek/data/`
아래에 보존한다.

Weak-lensing 대조장용 좌표 export 프로세스 PID 1049265는 2026-10-20 확인 시
`lageunha`에서 살아 있었다. checkpoint에서 SQL batch 000과 001은 완료,
002는 CAS job 1062312에 제출된 상태다. 프로세스를 재시작하지 않는다. Slurm
제어기 DNS SRV 조회가 실패해 `squeue`와 새 `sbatch` 작업을 사용할 수 없다.

### 2026-10-20 후속: 승인 shape 행과 렌즈 단계 재확인

현재 GPFS의 authoritative shape-row run은 실패한 `_02` 시도가 아니라
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_shape_rows_export_2026-10-18_03/`이다.
체크포인트와 요약을 다시 읽어 32/32 batch, 3,069,544 lens–source pair,
909,833,598 bytes, 22개 승인 열, `credentials_saved=false`, 전체 photo-z PDF
미포함, `science_measurement=false`를 확인했다. 각 배치의 SHA-256과 렌즈별
count 대조는 체크포인트 및 `row_summary.json`에 기록되어 있다. 재조회나 추가
반출은 하지 않았다.

이 승인 행과 로컬 HSC Y3 SACC `n(z)`, 좌표로 검증한 raw B-mode mask 극성을 사용한
cluster-only 약한렌즈 스택은 2,153,004 lens–source pair와 609개 렌즈를 포함한다.
52-cell 연결 sky jackknife 및 source-bin 일관성 진단이 Astrafek 스크립트로
저장되어 있다. 0.1–0.7 Mpc source-bin 간 차이가 남고 random 중심 차감 및 대응
random-field 공분산은 아직 없으므로 해당 스택 역시 `science_measurement=false`다.
수치와 해시는
`analysis/icl_wgl_membership/hsc_debug/wl_tomographic_consistency_assessment_2026-10-20.md`
및 `wl_source_reuse_covariance_validation_2026-10-20_02` 검증 산출물을 참조한다.

random-control 좌표 내보내기 PID 1049265는 이번 확인에서도 `lageunha`에서 살아
있었다(상태 `hrtimer_nanosleep`). 좌표 체크포인트는 batch 000·001 완료
(33,514/780,191행), batch 002는 CAS job 1062312 제출 상태로 그대로다. 이 CAS
job을 재제출하거나 살아 있는 프로세스를 재시작하지 않는다. 현재 output은 좌표
내보내기뿐이며 random source shape 행은 포함하지 않는다. 다음 실행 단계는
기존 좌표 체크포인트가 32/32가 된 뒤 행 수와 해시를 검증하고 전체 CAMIRA 조리개
중첩·중복좌표·렌즈별 지원을 다시 screening하는 것이다. 이후 random shape count와
shape-row 반출은 각각 실제 필터링된 중심 수·행 수·전송량을 확정해 별도 승인을
받은 뒤 진행한다.

2026-10-20 21:44 KST의 read-only SSH 재확인에서 PID 1049265는 34시간 37분째
살아 있었고 `hrtimer_nanosleep` 상태였다. 표준입출력은 원격 `/dev/pts/0`에 연결돼
있으며, 실행 중 메모리에 올라간 구버전 submitter는 CAS 상태를 확인한 뒤 60초
주기로 대기한다. 저장된 `coordinate_status.json`은 여전히 2/32 배치이며 최상위
`HTTPError` marker가 남아 있다. 성공 poll 시각을 저장하지 않는 구버전이므로
서버 측 작업 상태는 체크포인트로 판별할 수 없다. 프로세스와 제출 배치를 유지하고
중복 제출하지 않는다.

### 2026-10-20 후속: 전체 Gaia PSF 스택 및 승인 DESI 적색이동 보정

확장 Gaia ledger의 480 후보 중 hash-verified stamp 333개(227 프레임)를
기존 800개 HSC-I 컷아웃으로부터 스택했다. 이 표본에서도 외곽 64–96픽셀은
배경 annulus/source-mask 선택에 따라 부호가 달라져 확장 PSF 날개를 측정하지
못했다. 표본 수 증가로 두 설정은 0과 양립하고 한 설정은 양의 평균을 보였으나,
나머지 두 설정은 음수다. 국소 sky와 PSF 날개를 분리할 수 없어 `psf_approved=false`,
`science_measurement=false`다. 권위 산출물은
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/empirical_psf_frame_weighted_full_gaia_2026-10-20_02/`
이며 summary SHA-256은
`0be06f1a827823bcd04ce23db4113dca1bd368e1a2b7205d0bab0ed76911a7ba`다.

사용자가 DESI DR1/DR2 HSC PDR3 공개 보정 자료 4개만 다운로드하도록 승인했다.
고정 commit `abc70af29910e49fce17a63ef3083160d515fb6c`의 두 scale-cut NPZ와
JSON sidecar만 수신했고 총 549,818 bytes, 파일별 URL/크기/SHA-256은
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/desi_y3_hsc_nz_abc70af2_2026-10-20/manifest.json`에
있다. 4개 bin의 1,000점 `n(z)` mean/median와 점별 lower/upper, magnification
corrected 곡선을 확인했다. joint posterior 샘플은 포함하지 않는다. 계산 재현기는
`audit_wl_desi_nz_distribution_sensitivity.py`; SACC 기준과 두 scale-cut mean
곡선의 동일 shape rows 및 5도 sky jackknife, 점별 envelope stack을 포함한다.
random-lens shape rows가 없어 아직 cluster-only 진단이다. 결합 PSF/DESI Slurm job
430248은 69분 계산 후 pointwise envelope 행에 없는 jackknife-error 열을 요구하는
출력 조합 버그로 종료됐다. 승인된 기존 SACC 기준 summary가 3,069,544행 계획,
각 배치 해시, SACC 해시, 5도 52셀, 반경 구간, 우주론과 정확히 일치함을 확인해
재사용하도록 고쳤다. 반복 렌즈 적색편이의 `Sigma_crit^{-1}` 적분도 lens-z/bin별
한 번만 계산하도록 최적화했다. 회귀검사 356개와 subtest 6개가 통과했다. DESI만
계산하는 Slurm job 430481은 37분 41초 만에 완료했고 독립 validator도 통과했다.
기존 SACC stack과 비교해 DESI central curve는 7개 반경 median `DeltaSigma`를
source bin 1에서 +12.5–14.8%, bins 2–4에서 약 −1–3.5% 바꿨다. 네 central
model의 52셀/56×56 covariance는 rank 51이고 상대 비대칭은 4.1e-17 이하이다.
0.1–0.7 Mpc의 내부 토모그래피 불일치는 낮아지지만 계속 남았다(세 radial raw
4-bin p 범위 약 4.4e-8–1.2e-3). random lens subtraction이 없으므로 모두
cluster-only diagnostic이며 `science_measurement=false`다. 정량표와 해시는
[`wl_desi_nz_distribution_sensitivity_assessment_2026-10-20.md`](analysis/icl_wgl_membership/hsc_debug/wl_desi_nz_distribution_sensitivity_assessment_2026-10-20.md)
및 GPFS의 `spatial_jackknife_2026-10-20_01/summary.json`, `validation.json`에 있다.

### 2026-10-20 후속: 초기 20개 별 우표 PSF 파일럿 (확장 표본으로 대체)

초기 파일럿에서 승인된 기존 HSC 별 우표 20개(15개 프레임)를 Astrafek 재현 스크립트
`analysis/icl_wgl_membership/hsc_debug/stack_empirical_psf_star_profiles.py`로
프로필 스택했다. Kron 2.5/3.5와 120–150/160–190픽셀 배경 설정 모두 48–96픽셀
외곽 net profile이 음수였다. 따라서 현 국소 sky estimator와 별 표본으로는
양의 확장 PSF 날개를 측정하거나 보정할 수 없다. 이는 날개가 없다는 검출이
아니며 배경·ICL 혼입과 epoch 문제도 남는다. 결과는
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/empirical_psf_frame_weighted_stack_2026-10-20_02/`
및 `analysis/icl_wgl_membership/hsc_debug/empirical_psf_frame_weighted_stack_assessment_2026-10-20.md`에
해시와 함께 기록했다. 이 파일럿은 333개 stamp/227개 프레임의 확장 표본 평가에
의해 해석상 대체됐다. 확장 결과와 해시는 위 절에 기록했다. 두 결과 모두
`science_measurement=false`, `psf_approved=false` 상태를 유지한다.

### 2026-10-20 후속: Gaia PSF 별 우표의 matched-control WCS 지원

확장 Gaia ledger의 333개 별 우표(227개 표적 프레임) 각각에 대해 표적 은하단 중심
기준 east/north 각거리 오프셋을 기존 매칭 대조장 중심에 옮겨 WCS 대응을 확인했다.
319개 후보(217개 표적 프레임)는 같은 하위 필터 `i` 또는 `i2`였고 401×401 전체
control 우표가 영상 안에 들어왔다. 14개 후보(10개 프레임)는 target/control 하위
필터가 달라 제외했다. 검사기는 영상 헤더와 기존 star stamp만 읽고 coadd science
픽셀·새 HSC 자료는 읽지 않았다. Astrafek 재현기
`audit_empirical_psf_matched_control_support.py`와 테스트
`test_audit_empirical_psf_matched_control_support.py` 3개가 있다. 산출물은
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/empirical_psf_matched_control_support_2026-10-20_01/`에
있으며 summary SHA-256은 `a6d80adbfe8d01e466930dc68d44355287124b1b5ad8d8221a0c025f4cd7346f`,
행 표 SHA-256은 `fc255bf0704ddfd43e817e6b21e4a66f28c593e6748fda50976be12558f5da15`다.
이는 다음 대조장 잔차 profile의 입력 지원만 확인한다. PSF wing 보정 및
`science_measurement`는 승인되지 않았다. 상세는
`analysis/icl_wgl_membership/hsc_debug/empirical_psf_matched_control_support_assessment_2026-10-20.md`를
참조한다.

후속 paired profile은 기존 표적·대조 영상만 사용해 같은 필터 319개 별 위치(217 프레임)를 세 mask/coverage 설정으로 평가했다. outer-ring 중심 mask 조건을 조정한 뒤에도 최대 지원은 i2 64–96픽셀 구간의 12 프레임이며, sky annulus에 민감하다. control 위치는 별 없는 placebo이고 target에는 ICL 곡률이 남아 PSF wing을 분리할 수 없다. `_03` summary SHA-256은 `7478c905a4e46c0c82f3fe569e4ee1b4cfb53280c9acb635267bc4188a98cd5b`다. 방법·상태별 제외·모든 실행 해시는 `analysis/icl_wgl_membership/hsc_debug/empirical_psf_matched_control_profiles_assessment_2026-10-20.md`에 있다. `psf_approved=false`, `science_measurement=false`.

2026-10-20 22:04 KST read-only 원격 재확인: lageunha PID 1049265는 `Ss+`로 살아 있으며 누적 경과는 34시간 55분이다. 승인된 `wl_random_point_coordinate_plan_2026-10-19_01` 체크포인트는 `running`, 32개 중 2개 완료(33,514/780,191 예상 좌표), batch `002`는 CAS job `1062312` 제출 상태다. 기존 `HTTPError` 표시는 여전히 요청/시각이 특정되지 않아 terminal 실패로 판단하지 않았다. 재제출·프로세스 재시작은 하지 않았다. 이어질 전체 기하 스크리너는 승인된 좌표 32/32와 계획 해시 검증 뒤에만 실행한다.

전체 좌표 기하 스크리너 회귀검사를 보강해 target과 CAMIRA 이웃의 각각의 적색편이에서 계산한 proper-3-Mpc 각반경 합으로 교차/비교차 지점을 검증하고, chunk 크기 1과 10에서 결과 동일성을 확인했다. 좌표 체크포인트 감사기와 함께 focused test 7개가 통과했으며, 실제 전체 좌표 자료가 도착하기 전에는 기하 스크리닝을 실행하지 않는다.

Gaia↔HSC 후보 연결의 epoch 범위를 점검하려고 기존 333개 별 우표의 227개 target coadd FITS header를 조사했다. `DATE-OBS`/`DATE-*`/`MJD-*` 키는 0/227이다. 단, PSF 우표 중심은 Gaia가 아닌 HSC 검출의 `X_IMAGE/Y_IMAGE`이므로 이 누락이 영상 중심정렬을 막지 않는다. ±2/5/10년 Gaia 이동 결과는 실제 epoch 확인이 아닌 cross-match 민감도이며, 480개 모두 2초각 내에 남았다. 명확화한 Astrafek header-only 검사기는 `analysis/icl_wgl_membership/hsc_debug/audit_empirical_psf_hsc_coadd_epochs.py`, 테스트 3개, 상세·체크섬은 `analysis/icl_wgl_membership/hsc_debug/empirical_psf_hsc_coadd_epoch_audit_2026-10-20.md`에 있다. `_02` 요약 SHA-256은 `fb5ff0c7b87d1d17bdeb42ef824980b3a50eb7c7f0a487d7449b1fd9fe88d93c`이며 픽셀과 HSC 서비스는 읽지 않았다. `psf_approved=false`, `science_measurement=false`를 유지한다.

같은 333개 별 우표에 대해 HSC 포화/`BRIGHT_OBJECT` 비트, 검출 이웃, Kron source-mask 지지율을 추가 점검했다. 중심 12픽셀 내 포화는 0/333, 밝은천체 비트는 3/333이고, 12픽셀 이웃 검출은 없다. 중심 source-mask 개방률 99% 기준을 통과한 수는 Kron 2.5에서 290개/201프레임, Kron 3.5에서 242개/179프레임이다. 하지만 48픽셀 내 이웃은 317개, 96픽셀 내 이웃은 333개다. 중심 적격 subset에서도 48–64 및 64–96픽셀 개방률 중앙값은 Kron 2.5에서 0.718/0.694, Kron 3.5에서 0.605/0.579로 외곽 날개 측정에 강한 마스크·혼잡도 통제가 필요하다. 검사기 `analysis/icl_wgl_membership/hsc_debug/audit_empirical_psf_star_stamp_quality.py`, 관련 테스트 3개 통과. 결과는 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/empirical_psf_star_stamp_quality_audit_2026-10-20_02/`에 있고 summary SHA-256은 `bdb176c6ede927369e0c9b37bbc9ec88b0e1ba607bce89dfbab18795fe925bc0`, 행 표 SHA-256은 `ca2b908c717358ca0c3bb581da87ad717af6206b40ef0c6333d462a80ed48187`이다. 중심부 후보 진단이며 `psf_approved=false`, `science_measurement=false`를 유지한다. 해석과 한계는 `analysis/icl_wgl_membership/hsc_debug/empirical_psf_star_stamp_quality_assessment_2026-10-20.md`에 있다.

### 2026-10-20 후속: weak-lensing random 후보 count 전량 감사

이미 승인·수신된 count-only 표의 32개 CSV 해시와 800개 렌즈별 row를 좌표 SQL 계획에 대조했다. 30 arcmin 외곽 및 각 렌즈의 두 배 proper-3-Mpc 제외반경 안에서 HSC quality-filtered random 표본은 총 780,191개였다. 렌즈별 count는 364–1,443(중앙값 1,000)이고 전 렌즈가 계획한 200개를 넘는다. 유효 고리 면적으로 나눈 점밀도 중앙값은 0.566 arcmin⁻²다. raw count와 렌즈 적색편이의 고정표본 Spearman `rho=0.595`는 면적 정규화 뒤 `rho=0.080`으로 줄어, count의 z 의존 대부분이 각도 제외면적과 연결됨을 시사한다. 이는 중첩 제거 전 후보 가용성으로 shape 지원이나 최종 random pool을 보장하지 않는다. 전체 결과·입력 해시는 `analysis/icl_wgl_membership/hsc_debug/wl_random_coordinate_candidate_count_assessment_2026-10-20.md` 및 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_point_candidate_count_audit_2026-10-20_02/`에 있다. 검사기 `audit_wl_random_coordinate_candidate_counts.py`, 테스트 2개 통과. CAS 추가 요청·좌표 반출은 없었다.

현재 read-only lageunha 확인에서 PID 1049265는 여전히 살아 있음(`Ss+`, 경과 1일 11시간 32분)이며 좌표 체크포인트는 `running`, 완료 2/32, 33,514/780,191행, batch `random_coordinates_002` CAS job `1062312` 제출 상태다. 최상위 `HTTPError` 표시는 남아 있고 checkpoint mtime은 바뀌지 않았다. 서버 job의 terminal 상태를 확인하지 못했으므로 살아 있는 프로세스를 재시작하거나 batch를 재제출하지 않는다. 남은 전체 기하 스크리닝은 기존 export가 끝나고 32개 batch의 row count·SHA 검증을 통과한 뒤 수행한다.

### 2026-10-20 후속: 경험 PSF radial coverage와 ICL annulus

기존 333개 Gaia 별 우표/227 프레임의 표적 적색편이와 HSC 픽셀 scale을 사용해 경험 PSF 스택의 radial extent를 proper kpc로 환산했다. 실제 stack은 96픽셀까지이며 물리반경 중앙값 88.4 kpc (q10–q90 76.8–101.4 kpc)다. ICL 분석의 200–300 kpc annulus에는 227개 frame 중 하나도 닿지 않는다. 전체 우표의 circular-inscribed radius도 중앙값 184.1 kpc여서 200 kpc까지 닿는 프레임은 52개, 300 kpc는 0개다. 반면 96픽셀 endpoint는 75 kpc 이상을 덮는 stamp 316/333개로 satellite 외곽 분석과 일부 기하 중첩은 있다. 이것은 측정 PSF가 아니라 coverage audit다. 재현기 `analysis/icl_wgl_membership/hsc_debug/audit_empirical_psf_physical_radius_coverage.py`, 테스트 2개 통과. output `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/empirical_psf_physical_radius_audit_2026-10-20_01/`; summary SHA-256 `6d88022648930be6fd6f5f6785d66271cb399b1f6c4107dd7305453fead5d872`, table SHA-256 `b1b4cfee369a20939dbf3b6b040ff4a08b7efb4239e155829f65fc1bd99a1f32`. 상세 범위와 입력 provenance는 `analysis/icl_wgl_membership/hsc_debug/empirical_psf_physical_radius_coverage_assessment_2026-10-20.md`에 있다. PSF 보정과 `science_measurement`는 승인되지 않았다.

### 2026-10-20 후속: full-frame outer-field PSF-wing 파일럿

기존 333개 hash-verified Gaia 별 우표 중 18개가 target/control 양쪽에서 300 kpc 기하 지원, 동일 하위필터, 800 kpc 외곽장, 중심 Kron 2.5 개방률 조건을 통과했다. 이미 승인·다운로드된 800 target 및 798 control HSC 영상만 읽었다. Slurm job 430734는 18쌍의 25 kpc profile을 처리해 중복 없는 864행과 24개 frame-equal paired stack bin을 만들었다.

독립 검증에서 profile arithmetic·SHA-256이 일치했다. 하지만 radial bin별 paired support는 최대 4쌍이며 250–300 kpc에는 없고 200–225 kpc는 1쌍이다. `_01`에는 control profile 중복 버그가 있었으므로 사용하지 않는다. 수정된 `_02`도 PSF 보정 또는 ICL 과학 측정으로 쓸 수 없다: `science_measurement=false`, `psf_approved=false`. 핵심 산출물은 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/empirical_psf_outer_field_wings_2026-10-20_02/`이며 summary SHA-256 `c8085c8a7459a34ab03e551de6d2d37cccf637a80156bff4f70ba61b0dc46369`, validation SHA-256 `c70a144aae21922044d14a230404c8633a7d9be755580b976d7752a0fad860f6`이다. 재현기·실패 원인과 후속 기준은 `analysis/icl_wgl_membership/hsc_debug/empirical_psf_outer_field_wings_assessment_2026-10-20.md`에 있다.

### 2026-10-20 후속: random coordinate radial-yield

이미 승인된 첫 2/32 coordinate batch(50개 계획 순서 렌즈, 33,514행)만 사용해 inner 2×proper-3-Mpc 경계에서 거리에 따른 geometry-clean 후보수를 나눴다. 30′ outer limit은 각 렌즈에 7.14–8.76′(중앙값 8.10′)의 annulus 폭만 남긴다. 0–5′에서 1,199/19,920, 남은 shell에서 930/13,594가 geometry-clean이며 총 2,129(6.35%)다. 기존 50개는 비대표 표본이다. 30–60′ 추가 고리의 이상적 spherical area는 현재 검색 고리의 6.02–7.15배(중앙값 6.42배)이지만 HSC coverage·depth를 반영한 것은 아니다. 보충 count-only CAS 계획은 아직 제출하지 않았다. output `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_coordinate_radial_yield_audit_2026-10-20_03/`, summary SHA-256 `cd5a30122a2beb2ad875dec3f59bbca35286ee0862cbc1473d021cbb240402ad`. Astrafek 재현기와 해석은 `analysis/icl_wgl_membership/hsc_debug/wl_random_coordinate_radial_yield_assessment_2026-10-20.md`에 있다. `science_ready=false`, `science_measurement=false`.

### 2026-10-20 후속: random-coordinate export live audit

이번 read-only lageunha 확인에서 PID 1049265는 살아 있었다(`Ss+`, elapsed 약 1일 11시간 59분). 기존 체크포인트는 `running`; `random_coordinates_000`·`001` complete, `002`는 CAS job 1062312 submitted이며 top-level `HTTPError` marker가 남아 있다. `audit_wl_random_coordinate_checkpoint.py`로 저장된 CSV 2개를 count plan별 33,514/33,514행 및 SHA-256과 대조해 검증했다. 새 CAS 요청·중복 제출은 하지 않았고, 해당 batch의 실제 CAS terminal state는 여전히 모른다. checkpoint audit `_11`은 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_coordinate_checkpoint_audit_2026-10-20_11/`에 기록했다.

### 2026-10-20 후속: ICL×위성 반경 교차 지원

기존 800개 HSC full-sample 프로필의 여섯 ICL 고리(50–500 kpc)와 일곱 위성 고리(5–75 kpc)의 exact cluster/filter 교집합 84개를 작성·해시 검증했다. 가장 큰 지원은 ICL 100–150 kpc × 위성 5–10 kpc로 i 78개, i2 150개지만 위성 내부 구조를 뜻한다. 관심 외곽 조합인 ICL 200–300 kpc × 위성 50–75 kpc는 i 26개, i2 68개이며, 400–500 kpc × 50–75 kpc는 21/57개다. 최대 지원 조합으로 과학 질문을 바꾸지 않았다. 이번 단계는 지원 감사뿐이며 상관·유의성·선택 보정을 하지 않았다. 이 표본도 미검증 PSF 날개·sky 전달 때문에 `science_measurement=false`다. 재현기 `analysis/icl_wgl_membership/hsc_debug/audit_icl_satellite_radial_support_matrix.py`, 테스트 3개 통과; 보고서 `analysis/icl_wgl_membership/hsc_debug/icl_satellite_radial_support_matrix_assessment_2026-10-20.md`. 결과 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_radial_support_matrix_2026-10-20_01/`, summary SHA-256 `7d00bc9fad2e9066cc152c2485c3f5c301637594ac6acd37212cea1264204126`, 84행 table SHA-256 `1ee4f8e4a9c682296c64781d82cb74a31267acc1952b4814b643aa3b779be32d`.

### 2026-10-20 후속: 전체 ICL×위성 반경 가족 검정

상기 84개 반경/필터 조합의 exact 프로필 교집합과 추출된 대조장 좌표를 사용했다. 420개 고유 표적/대조 cluster-filter 쌍을 공통 5도 공간블록에 배정했고, 각 셀의 `satellite ΔSB ~ ICL ΔSB + z + richness` 부분상관을 CR1 공간블록 t 통계와 함께 계산했다. 84개 셀 전부 10개 행·8개 블록 최소기준을 통과했다. 모든 셀에 공통 Webb 6점 wild-cluster bootstrap 5,000회를 적용해 Holm 및 max-T 완전귀무 가족조정을 수행했다.

미보정 `p<0.05` 8개가 있었으나 Holm 조정 최소 p는 0.588, max-T 조정 최소 p는 0.156으로 84개 중 어느 것도 보정 후 유의하지 않았다. 후보 ICL 200–300 kpc × 위성 50–75 kpc의 i2는 `r=0.293`, 68쌍/25블록, `p_raw=0.0288`이지만 `p_Holm=1.000`, `p_maxT=0.987`; i는 `r=0.079`, 26쌍/9블록, `p_raw=0.621`이다. 표본별 지원 차이와 필터 불일치가 남고 PSF 날개·sky 전달·완전도 미해결이어서 탐색 결과이며 `science_measurement=false`다. 상세·한계는 `analysis/icl_wgl_membership/hsc_debug/icl_satellite_radial_matrix_familywise_assessment_2026-10-20.md` 및 재현기 `analyze_icl_satellite_radial_matrix_familywise.py`에 기록했다. 회귀검사 6개와 독립 계수·p조정·해시 검증을 통과했다. 결과 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_radial_matrix_familywise_2026-10-20_01/`; summary SHA-256 `9e4bf2e8369c2c7b3de025a51eafc74538b560b9c42a9cf2fc68de49e1a8ccfb`, association TSV SHA-256 `99ae3afd7a8845de9fb668cb415dcc2b7f8edded60bba550e2d69b5765014b0d`, max-T null TSV SHA-256 `2fcfe2a763d33a1ee12255c35dc77e5623dfdcf5121f4283597629d3368f4c19`.

### 2026-10-20 후속: CAMIRA 위성 가중을 포함한 168개 가족 민감도

기존 hash-verified member-weighted 프로필의 CAMIRA 열 8 `weight`를 사용해 84개 반경/필터 equal-satellite 분석과 CAMIRA-weighted 분석을 같은 168개 검정 가족으로 결합했다. 두 estimator는 모든 셀에서 동일한 cluster/filter 행과 공간블록을 유지했다. 각 클러스터 안에서만 CAMIRA weight를 적용하고 cluster는 동일 가중했다. 5,000회 공통 Webb wild-cluster bootstrap의 max-T와 Holm 보정 후 유의한 셀은 0개였다(17개 원시 `p<0.05`, 최소 Holm `1.000`, 최소 max-T `0.255`). ICL 200–300 × 위성 50–75 kpc의 i2는 equal `r=0.293`, weighted `r=0.241` (둘 다 68쌍/25블록); 원시 p는 각각 `0.030`과 `0.055`지만 168개 보정 후 둘 다 `p_Holm=1.000`, `p_maxT≥0.999`다. CAMIRA weight는 독립 불편 확률이 아니며 효과 차이에 별도 추론을 하지 않았다. PSF·sky·완전도 문제가 남아 `science_measurement=false`다. 재현기, 독립 validator, 결과·제한은 `analysis/icl_wgl_membership/hsc_debug/analyze_icl_satellite_weighted_radial_matrix_familywise.py`, `validate_icl_satellite_weighted_radial_matrix_familywise.py`, `icl_satellite_weighted_radial_matrix_familywise_assessment_2026-10-20.md`에 있다. 테스트 8개 및 168행/5,000 null·Holm/max-T·입력/출력 해시 독립 검증 통과. Slurm 430809는 4분 16초, max RSS 57,896 KB, stderr 0으로 종료했다. 결과 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_weighted_radial_matrix_familywise_2026-10-20_01/`, summary SHA-256 `0a337df80fda6fb92f178d251e09347f27012f5cc53f088fd967311b31ef8068`, 결과표 SHA-256 `6e6c37a44823ab0d6ac5f1795af87f4d47cc3dda1f65b0872ea9e622828142e7`, max-T null SHA-256 `04386baebe37577b40b96d2c859f9202e5cdc73b5eaf01b7ed6c1cc6b53e56a5`.

같은 200–300 × 50–75 kpc endpoint에서 equal 대비 CAMIRA-weighted 상관 변화량을 같은 cluster/filter 지원과 연결 5도 공간블록으로 10,000회 paired bootstrap했다. i: Δr=−0.003, 95% 구간 [−0.225, 0.092] (26개/9블록); i2: Δr=−0.052, [−0.139, 0.061] (68개/25블록). 두 구간 모두 0을 포함하며, 선택된 탐색 endpoint에 대한 민감도이지 독립 유의성 검정이 아니다. 분석 및 입력과 bootstrap 재계산 검증은 `analysis/icl_wgl_membership/hsc_debug/analyze_icl_satellite_weighting_contrast.py`, `icl_satellite_weighted_radial_matrix_familywise_assessment_2026-10-20.md`에 기록했다. 결과 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_weighting_contrast_2026-10-20_01/`; summary SHA-256 `b1a07b7cd46c8a2dc8fd562663a3297ba1d4f46c3a0b25df9ab35f9b56ffc699`, 결과표 SHA-256 `992bf7e7883b3b21babf88d6ba73c0d54e52d1e2ed8b1a7ff2ed07b294d0a012`. `science_measurement=false` 유지.

2026-10-20 cluster-only weak-lensing cross-shear null을 반경/source-bin별로 분해했다. 전체 28성분 χ²=74.75/28 (p=3.93×10⁻⁶)와 일치하지만 어떤 단일 반경 가족도 Holm 후 유의하지 않다(최소 p=0.102); source-bin family도 모두 비유의(최소 p=0.310)이며, 28개 개별 성분 family의 최소 Holm p=0.966이다. 현 요약통계에서 null 오염의 위치를 특정하지 못했다. 재현기·테스트·결과 SHA는 `analysis/icl_wgl_membership/hsc_debug/wl_tomographic_consistency_assessment_2026-10-20.md`에 있고, JSON은 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_cross_shear_null_decomposition_2026-10-20_01/`에 있다. 추가 CAS 조회 없음. `science_measurement=false` 유지.

2026-10-21 PSF support sweep은 outer 600 kpc 기준 91개 star/control 쌍, 78프레임을 찾았다. 첫 `_01` stack에서 `i`/`i2` 필터 혼합 오류를 발견해 이를 제외하고 필터별 `_02` 48셀 집계만 사용한다. 독립 검증기가 프로필 4,368행과 10,000회 frame-bootstrap을 재현했다. `i` 29쌍/25프레임, `i2` 62쌍/53프레임이다. 엄격 지지 기준에서 `i`의 50–75 kpc 표본은 0, `i2`는 4–6프레임이고, 200–300 kpc `i2` 지원은 0이다. mask 지지 완화에도 bin당 독립 프레임 최대 1–6개뿐이라 PSF 보정은 불가하다. 보고서 `analysis/icl_wgl_membership/hsc_debug/empirical_psf_outer_frame_bootstrap_assessment_2026-10-21.md`; 권위 산출물 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/empirical_psf_outer_field_wings_frame_bootstrap_2026-10-20_02/`, summary SHA `556fe0e828a5f2b0ff1600c4e4f2449147e9b383a115d3757eb4b5ea18faa388`, validation SHA `844337eb4fa89f72648ec1d660f57cc89208df2fa3ff9b44f15b406a3116da6f`. `psf_approved=false`, `science_measurement=false`.

### 2026-10-21 재연결 후 read-only 확인

2026-10-21 01:59 KST `lageunha` read-only 확인에서 승인된 random-center coordinate export PID 1049265가 여전히 살아 있었다(`Ss+`, elapsed 1일 14:48:55). 체크포인트는 `running`, mtime `2026-10-19 11:07:53 +0900`, 크기 936 bytes다. `random_coordinates_000`/`001`은 완료(총 33,514행), `002`는 CAS job `1062312` 제출 상태이며 최상위 `HTTPError` marker가 남아 있다. 이후 batch는 제출하지 않았고 PID를 재시작/종료하지 않았다. 완료된 32개 batch의 row/hash 검증 뒤에만 전체 기하 스크리닝을 진행한다.

같은 확인에서 target/control manifest(800/798개)를 다시 읽어 1,598개 고유 중심 중 승인 PSF 코어가 있는 위치 56개, 미포함 위치 1,542개라는 기존 오프라인 계획을 확인했다. 계획 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_psf_picker_target_control_centers_plan_2026-10-21_01/psf_picker_plan.json`은 `prepared_not_submitted`; `pdr3_wide`, `coadd`, `HSC-I` 중심 모델 PSF 1,542개만 요청하며 image/shape/photo-z는 요청하지 않는다. 좌표 목록 SHA-256은 `78bb8b426ad3aa99df68d757b8981e0893517c9e246667aa15d0c891d071d99c`, mapping SHA-256은 `4c2f060f350b65fc93dbb0a8b7c12421517b9be54ddba09c866b08a401dc0212`다. 이전 PSF 파일 크기 분포로 계산한 예상량은 평균 27.6 MiB, 중앙값/관측 최대 기준 29.6 MiB였다. 사용자가 2026-10-21에 이 1,542개 범위만 승인했다. HSC PDR3 문서상 i/i2 visit은 하나의 HSC-I coadd로 합쳐지므로 요청 band는 HSC-I로 유지하되, 분석 진단의 i/i2 표본은 계속 분리한다.

승인 후 기존과 동일한 공식 HSC PDR3 PSF Picker `downloadPsf.py`를 재사용했다. 내려받은 원본은 SHA-256 `17cee33a86de87a77f3f800049b3fdb056f2e0c50b96c2b8cf7e80c1f62f2298`로 앞선 성공 실행과 일치한다. 요청은 좌표목록 방식으로 `rerun=pdr3_wide`, `type=coadd`, `filter=HSC-I`, centered stamp만 받았으며 정상 종료(exit 0)했다. 데이터는 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_psf_picker_target_control_centers_download_2026-10-21_01/`에 보존했다. 승인된 1,542개 FITS가 좌표·매핑 해시와 일대일 대응하고 누락·중복·추가 파일은 0개다. 총 29,269,440 bytes이며 크기는 17,280–20,160 bytes, 형상은 43×43 911개, 43×41 242개, 41×43 266개, 41×41 123개다. 모든 픽셀은 유한하고 전체 합은 0.9999999999999986–1.0000000000000018, 양의 광량 중심의 최대 편차는 0.2761 pixel이다. 음수 픽셀 699,893개를 포함한 원본 부호를 보존해야 한다. 파일별 SHA-256 및 좌표·FITS 검증 manifest는 같은 폴더의 `psf_picker_audit.json`; SHA-256 `b6079bd47f5dd0a716a5841fc7c1cca01f96a857a4cbd7a4d8d9132b4856e6b4`다. 감사기는 `audit_hsc_psf_picker_target_control_batch.py`로 저장하고 3개 테스트를 통과시켰다. 이 스크립트로 실제 1,542개 전체를 다시 감사해 기존 manifest와 바이트 단위로 같은 해시를 얻었다. 재현 실행은 기존 파일을 덮어쓰지 않도록 `psf_picker_audit_reproducible.json`에 기록했으며 SHA-256도 동일하다. 기존 승인 PSF 56개와 합쳐 800 target + 798 matched-control 중심 1,598곳 모두에 모델 PSF 코어가 확보됐다. 그러나 이는 중심 코어 입력의 완전성만 확인한 것으로 PSF 날개·영상별 하늘 전달·주입회수 성능은 검증하지 않았다. `psf_approved=false`, `science_measurement=false`를 유지한다. 영상·shape·photo-z는 받지 않았고 자격정보도 저장하지 않았다. 이 다운로드는 실행 중인 random-coordinate CAS 작업을 건드리지 않았다.

공개 확장 HSC-i PSF 대칭·비대칭 FITS 2개(각 9,809,280 bytes, 합계 18.71 MiB)는 이미 사용자 승인 후 수신되어 Zenodo 체크섬과 독립 재검증을 마쳤다. 직전 대화에서 이 자료가 아직 미수신인 것처럼 말한 것은 부정확했다. 기존 위치별 PSF-picker 744 코어도 이미 수신됐다. 다만 두 종류 모두 PSF/하늘 전달을 과학적으로 승인하지 않았으며 `psf_approved=false`, `science_measurement=false`를 유지한다. 세부 결과는 `analysis/icl_wgl_membership/hsc_debug/extended_psf_public_reference_assessment_2026-10-21.md`에 있다.

### 2026-10-21 전체 target/control 위치별 PSF 코어 및 ICL 의존성

신규 재현 감사기 `audit_hsc_psf_picker_target_control_batch.py`를 저장하고 테스트했다. 기존 감사 manifest와 공식 도구 해시 `17cee33a…f2298`를 확인한 1,542개 수신 PSF에서 이 스크립트를 다시 실행해 바이트 단위로 동일한 audit SHA를 재현했다. 감사 테스트와 새 target/control 연결 테스트를 합쳐 9개가 통과했다. 신규 1,542개와 기존 승인 위치 PSF 중 target/control 중심과 겹치는 56개를 합쳐 1,598개 중심을 `analyze_hsc_target_control_psf_cores.py`로 분석했다. 중심 table은 1,598행, paired table은 798행이며 800 target 중 control이 없는 것은 2개다. target-control 장축 moment FWHM 차이 중앙값은 −0.0070 arcsec, q16–q84는 −0.0572–+0.0392 arcsec, 중앙 절대차는 0.0319 arcsec다. 해당 요약 SHA는 `cb9244c39e5608932434eec8eb0e7aca0a50ebac585cd769a19f996937843bb0`이다.

이 코어 차이가 기존 ICL 200–300 kpc 초과광과 함께 변하는지 공통 면적 40% 통과 137쌍(44 i, 93 i2)에서 추가 진단했다. 적색편이·풍부도·공통 피복률을 조정한 partial r은 i `0.154`, i2 `0.038`; 중앙 마스크 flag까지 넣으면 i `0.016`, i2 `0.036`이다. 5,000회 연결 5° 공간블록 bootstrap에서 여섯 모델의 모든 2.5–97.5 percentile 범위가 0을 포함한다. 독립 validator가 6개 점추정, 총 30,000개 재표집 상관, leave-one-block-out 범위를 재현했다. 보고서는 `analysis/icl_wgl_membership/hsc_debug/icl_psf_core_mismatch_full_sample_assessment_2026-10-21.md`; 분석·독립검증 코드는 `analyze_full_sample_icl_psf_core_dependence.py`, `validate_full_sample_icl_psf_core_dependence.py`다. 결과 요약 SHA-256은 `f67019b8cb6f3c096b45e8640d489022db19f2b1d29f7b4db1ab89441936c34e`, joined table SHA-256은 `d6e0f0b669e3eb2c0f7a57c37669260da4ad17216daaf66fda2c60af7223fcf5`다. Slurm controller DNS 실패로 이 137행 단일 CPU 작업은 로그인 노드에서 약 6초에 끝냈다. 코어가 target-control 중앙값에서 비슷해 보인다는 결과는 확장 PSF 날개나 하늘 전달 검증을 뜻하지 않으므로 `psf_approved=false`, `science_measurement=false`다.

weak-lensing random-control-coordinate export는 재확인 당시 PID 1049265가 살아 있었으나 체크포인트 mtime은 계속 `2026-10-19 11:07:53 KST`였다. CAS job `1062310`, `1062311` 두 개는 33,514행으로 완료되어 있고 `1062312`는 체크포인트에서 `submitted`로 남아 있다. Astrafek에 추가한 read-only `check_wl_cas_job_status.py`로 CAS에 딱 한 번 상태를 조회한 결과 `1062312`는 실제 `running`이었다. 이 상태 조회는 `query_submitted=false`, `rows_downloaded=false`이며 자격정보를 저장하지 않았다. `audit_wl_random_coordinate_checkpoint.py`도 다시 실행해 2/32 배치, 33,514/780,191행과 파일 해시를 확인했다. 새 audit `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_coordinate_checkpoint_audit_2026-10-21_01/checkpoint_audit.json` SHA-256은 `d5072ee7cecd4ec4798892dc6d3be4acf1d916738e0ff65564666234a3ec8212`, checkpoint SHA-256은 `84cccc6443ec8d2fc8bee1611ab6b9b57a3aabdc90978fcdfc558151d2bb92b2`다. 활성 submitter를 중단·재시작하지 않았고 같은 CAS 쿼리를 재제출하지 않았다. 저장된 완료분과 결합하기 전에는 기존 작업의 terminal status 및 체크포인트 갱신을 확인해야 한다. 상태조회기 테스트 2개 통과.

### 2026-10-21 확장 피복 계획 및 CAS 좌표 작업 재확인

lageunha에서 승인된 `submit_icl_expansion_coverage_metadata.py --check`를 실행했다. 기존 저장 plan의 체크섬·SQL 구조·열 제한이 모두 통과했고, 후보 2,596개·probe 23,364개·배치 52개를 확인했다. 이 계획은 PDR3 `mosaic`/`skymap`/`patch_qa` 메타데이터만 반환하며 image, object/source rows, shape, 개별 photo-z 요청 플래그가 모두 false다. 이 검사는 네트워크 요청을 하지 않았고 아직 `prepared_not_submitted` 상태다.

같은 read-only 확인에서 weak-lensing 좌표 export PID 1049265가 `lageunha`에서 계속 살아 있었다(`Ss+`, `hrtime`, elapsed 1일 15:56). `cas_coordinates/coordinate_status.json`은 2026-10-19 11:07 KST 이후 갱신되지 않았으며 32개 중 2개 완료(33,514/780,191행), CAS job 1062312 제출 상태, `error_type=HTTPError`를 기록한다. 현재 `squeue`는 Slurm controller DNS SRV 조회 오류로 상태를 가져오지 못했다. 이 확인에서는 HSC CAS job 상태를 새로 조회하지 않았고 좌표 worker를 종료·재시작하거나 새 CAS 쿼리를 제출하지 않았다. 다음 단계는 기존 worker/job의 권위 있는 terminal 상태를 확인한 뒤, 이미 승인된 피복 메타데이터 52개 batch를 기존 HSC 세션에서 순차 실행하는 것이다.

### 2026-10-21 전체 표적·대조 PSF 코어/공개 확장 PSF 비교

새 Astrafek 분석기 `analyze_public_psf_target_control_overlap.py`로 800개 표적과 798개 대조장의 모델 PSF 코어 1,598개를 공개 HSC-I 확장 PSF와 1.5–3.3 arcsec에서 비교했다. target filter는 기존 검증 스택, control filter는 해당 대조장의 실제 `coadd_header_filter`에서 가져왔다. 동일 filter 쌍은 i 302개, i2 481개(총 783개); 다른 filter 15쌍은 paired 통계에서 제외됐다. 2개 target에는 control이 없다.

전 창 진폭 scale 중앙값(i target/control, i2 target/control)은 0.713/0.740, 0.787/0.797이다. 중앙 scale은 inner 1.5–2.4에서 outer 2.4–3.3 arcsec로 각각 0.731→0.520, 0.754→0.560, 0.803→0.555, 0.817→0.572로 감소한다. full-window 모양 잔차 RMS 중앙값은 11.7–12.8%다. same-filter target-minus-control scale 차이 중앙값은 i −0.0239 (302쌍, q16–q84 −0.1294…0.0644), i2 −0.0033 (481쌍, −0.1181…0.0988)이다. 이는 위치별 분포/민감도이지 보정계수나 추론 검정이 아니다. 한 진폭으로 public wing을 접합할 수 있다는 증거가 없으며 200–300 kpc wing 및 sky 전달은 계속 미검증이다.

초기 `_01` 실행은 control에 target filter를 상속하여 paired 결과에 사용하지 않는다. 실제 control filter를 조인한 `_02`만 권위 출력이다. 독립 FITS 재계산기가 1,598개 전체의 네 fit window, 동일-filter 783쌍, 15 mismatch를 재현했다. 7개 focused test와 syntax/whitespace 검사 통과. 분석기 `analysis/icl_wgl_membership/hsc_debug/analyze_public_psf_target_control_overlap.py`, 독립 검증기 `validate_public_psf_target_control_overlap.py`, 상세 보고서 `extended_psf_target_control_overlap_assessment_2026-10-21.md`. GPFS 결과 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/extended_psf_target_control_overlap_2026-10-21_02/`; summary SHA-256 `8ba338cf139ab469d1c21a10f4401ea3d6a86e3dd402ea07bc49ca512be2e6d8`, position TSV `ffe0a557d47694f73156b8f084d22f5977d658bb1b6932fe1cb5400dba22c62a`, paired TSV `b6276e52678abc91499ce56a65acbaec0ceaa055599aa56fbecb2a6be1686556`, validation JSON `add0ce56323fe32564147f6ab9ec91d05d8ee3d2810f62c69494ac1b38b607ec`. `psf_approved=false`, `science_measurement=false` 유지. 후속은 기존 image mask에서 splice radius·진폭별 점광원/diffuse injection recovery다.

### 2026-10-21 `i` 필터 공개 PSF 주입회수

Astrafek `extended_psf_injection.py`로 동일 필터 PSF 지원과 K2.5/40% ICL 지원을 만족하는 13개 `i` 프레임에서 HSC-I 대칭·비대칭 public PSF를 각각 주입 민감도에 사용했다. 두 run 모두 13/13, error 0, 각각 약 68초. 독립 `validate_extended_psf_injection_batch.py`가 저장 응답을 재계산했고 최대 산술 오차는 0이다. 기존 public PSF는 새 승인 위치의 동일 byte file(SHA `8897af…2d85`)임을 확인했다.

대칭 모델 회수 중앙값/프레임 범위: smooth diffuse→ICL 200–300 kpc `0.215` (`−2.734…0.271`), 위성 외곽광→ICL `0.353` (`0.024…1.269`), diffuse→위성 50–75 kpc `−0.039` (`−0.915…4.324`), 위성 외곽광→위성 `0.294` (`−3.024…1.895`). 13개 표본은 선택된 지원 집합이며 quantile은 모집단 CI가 아니다. 비대칭−대칭 paired 모델 차이는 네 응답에서 최대 절댓값 `4.77×10⁻⁷`; 동일한 공개 별 표본을 쓰므로 독립 검증은 아니다. 음수와 넓은 분산 때문에 scalar recovery correction을 적용하지 않는다. 실제 science coadd는 수정하지 않았고 source photon noise/coadd covariance는 추가하지 않았다. `science_measurement=false`, `psf_approved=false`.

결과 보고서 `analysis/icl_wgl_membership/hsc_debug/extended_psf_i_injection_recovery_assessment_2026-10-21.md`. symmetric run `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_extended_psf_i_sensitivity_2026-10-21_01/` summary SHA `e9853c6626cace7b630d6db17e3d9a8587a26007d726b884801608f2719be228`, validation SHA `aa1c1b3fd80f044d76f66a25ca2da19517b90c2e511ef48d3d2c3033b43df2d1`. asymmetric run `icl_extended_psf_i_asymmetric_sensitivity_2026-10-21_01/` summary SHA `dd3f40c56c55840a6a268eb9192cfaf7acede0183a138951df978bf63f882527`, validation SHA `b62255417a19e069bbb53879b002bde3932a40d10a98cf2f3f8cba80db382aa8`. Filter-aware paired comparison `variant_comparison_filter_aware.json` SHA `9281dca86f14aaa4ce7e63537b021e2ce861e1c3663f5b4a17b7ff86f5fba575`; earlier `variant_comparison.json` is retained but superseded. Comparator `compare_extended_psf_injection_variants.py`를 i/i2 필터 선택형으로 일반화; focused test 8개 통과. 다음은 여러 splice 반경·진폭의 real-mask 주입과 더 넓은 local wing 검증이다.


### 2026-10-21 고정 마스크 주입 민감도 검증 및 Slurm 위치

`icl_extended_psf_i_fixed_masks_2026-10-21_01/`의 13개 프레임을 독립 validator로 다시 계산했다. 13/13 완료, 오류 0, 저장 산술 최대 오차 0이다. 고정 baseline mask에서 diffuse→ICL 회수 중앙값은 0.22495 (q16–q84 0.22175–0.22717), diffuse→위성 50–75 kpc는 0.00062 (−0.01161–0.01260)다. 위성-envelope→ICL은 0.494 (0.165–0.713), 위성-envelope→위성은 0.278 (−0.174–0.509)으로 여전히 분산이 크다. 재탐지−고정마스크 diffuse→ICL 차이는 중앙값 −0.0115지만 최솟값 −2.955이며 diffuse→위성은 범위 −0.900–4.326이다. 따라서 일부 프레임의 극단값은 재탐지/마스크 선택에 민감하다. 이는 합성 주입 민감도일 뿐 PSF 보정이나 과학 측정이 아니다. 보고서 `extended_psf_i_injection_recovery_assessment_2026-10-21.md`에 추가했다. 결과 SHA-256: summary `a732fbc5fdaf41932f5685a329143f8d794f543885b5a5c0917ea321c3c264c5`, validation `2bcd49fa584c18e2415be56ed4fa75fa4894f2762b839ec9cdf689b2147c0b45`. 관련 focused pytest 11개 통과.

운영 위치를 확인했다. Slurm 작업의 조회·제출은 syntax에서 한다. 단, 1062312는 Slurm job ID가 아니라 HSC CAS job ID이므로 `squeue`/`sacct` 조회 결과는 이 CAS 요청의 상태 증거가 아니다. lageunha 좌표 exporter PID 1049265는 계속 살아 있으며 checkpoint는 2/32 배치, 33,514/780,191행, 마지막 갱신은 2026-10-19 11:07 KST다. CAS 상태를 별도 인증 조회하지 않았고, 프로세스를 종료·재시작하거나 요청을 재제출하지 않았다.


### 2026-10-21 공개 PSF 주입의 재탐지 영향도

Astrafek `audit_extended_psf_injection_reselection_influence.py`와 독립 focused tests 2개를 추가했다. 13개 `i` 프레임의 입력 파일 해시, summary-to-frame 응답값, redetection−fixed-mask 산술을 대조하고 프레임별 지원 변화와 leave-one-frame-out 중앙값을 기록한다. 가장 큰 재탐지 outlier는 cluster 4982다. baseline에서 matched satellite 4개/control 16개였으나 diffuse-only·joint injection에서 2개/8개로 감소했고 ICL 200–300 kpc coverage는 0.5395에서 0.4982로 떨어졌다. 이 프레임의 diffuse→ICL 재탐지−고정마스크 응답차는 −2.955, diffuse→위성은 +4.326이다. cluster 4982 제외 시 diffuse→ICL 재탐지 median은 0.2150→0.2181, 고정마스크 median은 0.22495→0.22522로 유지되지만 satellite-envelope→satellite median은 0.294→0.218로 변한다. 선택된 13프레임의 기술적 leave-one-out일 뿐 모집단 추론이 아니다. 보고서 `extended_psf_i_injection_recovery_assessment_2026-10-21.md`; 산출물 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_extended_psf_i_reselection_influence_2026-10-21_01/reselection_influence_audit.json` SHA-256 `09dc009995b8dc8bbbc29cff0b392812f47b5351e399c57c943a26d8a56cee29`, 도구 SHA-256 `6b9f2a1e9f9f97dcbc1c27a21e1383524a42aed8a023c5a8b5fe086c8129ccad`. `science_measurement=false`, `psf_approved=false` 유지.

기존 random-control 좌표 exporter PID 1049265를 60초 간격 재확인했다. lageunha에서 여전히 `hrtimer_nanosleep`이며 checkpoint는 2026-10-19 11:07:53 KST 그대로(2/32 배치, 33,514/780,191행; batch 002 CAS 1062312 `submitted`)다. CAS 재조회·재제출·프로세스 재시작은 하지 않았다. 1062312는 Slurm ID가 아니라 HSC CAS ID이며 Slurm 작업은 syntax에서 조회·제출한다.

### 2026-10-21 위치별 PSF 코어/공개 날개 하이브리드 구성기

`position_psf.py`의 `build_core_wing_hybrid_psf`는 signed local core와 public extended wing을 smooth radial taper로 결합해 단위 정규화한다. 중앙 코어·음수 픽셀 보존 및 범위 검사를 단위시험으로 확인했다. 이 커널은 synthetic-light transfer 민감도 시험에 적용됐으나 observed coadd를 수정하지 않았고 새 HSC 요청도 없다. 결과는 `psf_approved=false`, `science_measurement=false`다. 재현 및 제한은 `analysis/icl_wgl_membership/hsc_debug/core_wing_hybrid_psf_design_2026-10-21.md`에 있다.

### 2026-10-21 i2 확장 PSF 프록시/위치별 코어 주입 민감도 비교

이미 완료된 HSC-I public symmetric/asymmetric proxy 주입을 위치별 PSF core 주입의 동일 44개 `i2` 프레임과 paired 비교했다. 기존 배치 각각 44/44 프레임, 오류 0이고 저장 응답을 독립 재계산한 최대 산술 차이도 0이다. 프레임·지원·주입 규격을 맞춘 재현 감사기 `audit_extended_psf_core_vs_wide.py`에 `i2` 선택 경로를 추가했고 focused tests 3개가 통과했다. 산출물 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_extended_psf_i2_proxy_variant_comparison_2026-10-21_01/core_vs_wide_audit.json`, SHA-256 `ba1287493d5d3db4d8050f601ac982a659e1ca8588c024c90d8a7814c30e0f97`; 감사기 SHA-256 `8d08096c9d4babab32aa61fc6bef34d710267e65873c7cc1296ac3c36009f625`.

Diffuse→ICL 200–300 kpc 회수율 중앙값은 local core `0.2245`, public proxy `0.2235`, paired 차이 `−0.00036` (q16–q84 `−0.00092–0.00093`)였다. Satellite-envelope→ICL 중앙값은 `0.4904`와 `0.4849`, paired 차이 `0.00846` (`−0.00773–0.02847`)였다. 반면 satellite 50–75 kpc 회수 분산은 두 모델에서 매우 컸다. Satellite-envelope→satellite paired 차이의 q16–q84는 `−0.191–0.081`, 전체 범위는 `−3.551–0.377`이다. 대칭·비대칭 public PSF 모델 차이는 최대 `9.5×10⁻⁷`로 수치 정밀도 수준이며, 두 모델은 같은 별 표본을 써 독립 검증이 아니다. HSC-I public PSF는 i2 전용 보정이 아니므로 프록시 민감도에 한정하고 회수율 보정이나 과학 측정으로 승격하지 않는다. 보고서 `analysis/icl_wgl_membership/hsc_debug/extended_psf_i2_proxy_core_comparison_assessment_2026-10-21.md`; `science_measurement=false`, `psf_approved=false` 유지. 다음 ICL gate는 기존에 설계한 multi-radius/amplitude local-core/public-wing hybrid 점광원·diffuse real-mask 주입 파일럿이다.

### 2026-10-21 위치별 하이브리드 위성 PSF 경로 및 i-band 파일럿 실행

하이브리드 주입에서 위성 envelope에 은하단 중심 PSF 코어를 재사용하던 것을 수정했다. 확산 ICL에는 중심 위치 코어+public wing hybrid를, 각 위성 envelope에는 해당 위성 위치 코어+같은 wing을 적용한다. 결과 파일과 독립 validator는 각 위성 코어의 filename/SHA, 정규화, wing parameters, 승인 플래그를 확인한다. `test_position_psf.py`, `test_extended_psf_injection.py`, `test_validate_extended_psf_injection_batch.py`, 새 grid 테스트를 포함한 focused test 22개 통과했다.

미리 고정한 12개 조합은 splice 1.8/2.4/3.0 arcsec, blend 0.6 arcsec, wing amplitude 0.50/0.65/0.80/1.00이다. cluster 1056의 10개 위성 위치를 이용한 grid job `432792`는 syntax `syn04`에서 13m58s 만에 정상 완료됐다(1 CPU, 32 GB, GPU 미요청; peak RSS 2,866,960 KiB). 12/12 설정과 독립 response 산술 검증이 통과했고 최대 오차는 0이다. Diffuse→ICL 응답은 grid에서 0.22111–0.22128로 안정적이었지만 위성 외곽광→위성은 mask 재탐지 시 0.09297–0.16836, 고정 mask에서 0.31063–0.32470이었다. 단일 cluster라 모집단 범위는 아니다. 재현 보고서 `core_wing_hybrid_i_pilot_assessment_2026-10-21.md`; `grid_audit.json` SHA-256 `18606dcd49e8b918186eda6158e2209f2755a76b767806756e3571d7b4990ee1`. 다음은 13개 supported i 프레임의 frame-to-frame 민감도다. 이 grid에는 point source model profile이 포함되지 않았고 별도 noiseless point-source mask diagnostic만 cluster 1056에서 수행됐다(아래 절). 모든 과학·PSF 승인 플래그는 false다.


### 2026-10-21 추가: position-specific hybrid 단일 프레임/점광원 진단 완료

`lageunha`에서 사전 정의한 단일 프레임 cluster 1056의 새 고유 출력 경로로 splice 2.4 arcsec, blend 0.6 arcsec, wing amplitude 0.65를 실행했다. 앞서 기록한 syntax Slurm array(job 432792)와는 별도 실행이다. 이번에는 Slurm 명령을 실행하지 않았다. 1/1 프레임, 10개 satellite 위치 코어, elapsed 12.12초, 오류 0으로 완료했다. 독립 validator는 모든 저장 응답을 재계산해 최대 오차 0, `validated_exploratory_transfer_only`를 반환했다.

고정 mask 응답은 diffuse→ICL 0.230718, 위성 envelope→ICL 0.506852, diffuse→satellite −0.012496, satellite envelope→satellite 0.314188이다. 동일 cluster의 기존 public PSF fixed-mask 값과 차이는 각각 약 +0.000105, −0.000554, −0.000093, −0.012785다. 별도 noiseless unit point-source 진단에서 50–100 kpc annulus의 mask coverage/마스크 내 flux fraction은 46.3%/0.001996, 200–300 kpc는 55.4%/0.000834다. 이는 단일 프레임 민감도 진단이며 Poisson noise, coadd covariance, sky transfer, 관측원 제거를 반영하지 않는다.

결과 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_hybrid_psf_i_pilot_2026-10-21_r2p4_a0p65_perpos_point_01/`; `summary.json` SHA-256 `8a2fa26dd197aae9f29fafaaaac7e91afeab0bcc9c1d2299b2a73704a367d16e`, `cluster_1056.json` `ddb650883c7e664fb9c026d28a83768c195d7b5ccaa867b6b5c58d8048666711`, `validation.json` `d547fba7bfbae4869f48c6946a316190c860afe4dd711b301c2205e55d03c568`. focused tests 23개 통과; `py_compile`, array script `bash -n`, `git diff --check` 통과. 설계와 다음 grid 단계는 `analysis/icl_wgl_membership/hsc_debug/core_wing_hybrid_psf_design_2026-10-21.md`. 해당 결과도 `psf_approved=false`, `science_measurement=false`다.

2026-10-22 `validate_extended_psf_injection_batch.py`로 위 단일 hybrid+point-profile batch를 독립 재검증했다. 출력 `validation_repro_2026-10-22.json`은 기존 검증과 동일 SHA-256 `d547fba7bfbae4869f48c6946a316190c860afe4dd711b301c2205e55d03c568`, 1/1 프레임, 최대 응답 산술 오차 0이다.

### 2026-10-22 supported i 프레임 전체 hybrid 민감도 grid 실행 중

2026-10-21 cluster 1056 단일 프레임에서 12개 PSF 모델 grid가 끝난 뒤, 13개 supported `i` 프레임으로 동일한 frozen grid를 확대했다. 프레임별 결과에는 redetection/fixed-mask diffuse·satellite 응답과 별도의 noiseless unit point-source profile/mask transmission이 포함된다. 재현 가능한 `--all-supported` driver, updated point-profile grid summarizer, CPU-only Slurm wrapper, 범위 문서 `core_wing_hybrid_i_supported_plan_2026-10-22.md`를 추가했다. 관련 grid/summarizer/validator focused test 14개 통과.

Slurm job `432803`은 `syntax`의 `a10`/`syn04`에서 2:07:02 만에 정상 완료했다(1 CPU, 32 GB, GPU 미요청, peak RSS 4,363,184 KiB). 12개 설정×13 프레임=156개 frame record, 오류 0, 각 저장 응답의 독립 재계산 최대 오차 0이다. 점광원 profile은 별도 noiseless unit source 모델을 saved masks에 통과시킨 것으로 source detection/noise 회수가 아니며 독립 FITS 재구성은 하지 않았다.

중앙 `r2.4_a0.65` 모델에서 fixed-mask diffuse→ICL median 0.2251 (q16–q84 0.2219–0.2272), redetection median 0.2141이며 12개 모델 medians의 envelope는 각각 0.2250–0.2252 및 0.2141–0.2150이다. Satellite envelope→ICL redetection−fixed paired median은 −0.1022이고 PSF grid에 따른 median 범위는 −0.1045…−0.0825다. Satellite envelope→satellite paired 차이는 중앙 모델 median −0.0536이나 q16–q84 −0.3158…+0.2065, 전체 범위 −3.0894…+2.1677 (cluster 430 및 4982 극단)으로 단일 보정값을 쓸 수 없다. Point source 유효 flux fraction은 중앙 모델 median 50–100 kpc 0.00259, 200–300 kpc 0.000957다. 결과는 과학 측정·PSF 승인 또는 ICL–위성 상관의 증거가 아니다.

독립 Astrafek 요약은 `analysis/icl_wgl_membership/hsc_debug/core_wing_hybrid_i_supported_assessment_2026-10-22.md`. v2 audit은 12 configuration 전체가 동일한 13개 고유 cluster ID를 썼는지도 확인한다. `grid_audit_v2.json` SHA-256 `c0545e9f36c1b0d47589744abfcf21264461e8bf467f6f63fd379436ee35ba2c`, checkpoint SHA-256 `91f5aed2dbff64fd5ae07290c30b8f77b69fde15535c881d1822352ba25747b6`. 최초 `grid_audit.json`도 보존했다. 출력 `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_core_wing_hybrid_i_supported13_2026-10-22_01/`. 다음 단계는 matched-control sky/PSF-wing leakage와 pixel/noise response를 관련 물리반경까지 검증하는 것이며, 현재 응답 중앙값을 보정값으로 사용하지 않는다. `science_measurement=false`, `psf_approved=false` 유지.
## 2026-10-21 ICL 주입 응답의 마스크 재선택 분리

### 고정 마스크와 재검출 비교

이전 2026-10-19 위치별 HSC PSF 코어 민감도 시험을 이어, PSF를 합성 은하단 영상에
합성곱한 뒤 확산 성분과 위성 성분의 주입 응답을 다시 계산했다. 기존 방식은 각
주입 시나리오마다 source detection을 다시 수행하므로, 주입된 빛의 전달과 검출
마스크·위성 표본의 변경 효과가 한 응답에 섞여 있었다.

이를 분리하기 위해 로컬 작업 트리에서
`analysis/icl_wgl_membership/hsc_debug/extended_psf_injection.py`에 기준 검출 마스크를
고정한 분기와 `calculate_transfer_response`를 추가했다. 프레임 JSON schema v2에는
`measurements_by_scenario_fixed_masks`, `fixed_mask_transfer_response`,
`redetection_minus_fixed_mask_response`를 저장한다. 독립 검사기
`validate_extended_psf_injection_batch.py`는 v1/v2를 읽고 두 분기 및 차이를
재계산한다. 이 코드 변경은 이 인수인계 커밋에 포함하지 않는다. 재개 시 현재 로컬
작업 트리 상태와 테스트를 확인하고 코드 변경은 별도로 검토한다.

### 13개 HSC-I 프레임의 탐색 결과

lageunha CPU 실행은 `icl_extended_psf_i_fixed_masks_2026-10-21_01`에서 13/13 프레임,
오류 0, 72.627초로 완료했다. 독립 검사 결과는
`validated_exploratory_transfer_only`, 13개 검증, 저장 응답과 재계산 응답의
최대 차이 0이다. 요약 SHA-256은
`a732fbc5fdaf41932f5685a329143f8d794f543885b5a5c0917ea321c3c264c5`,
검증 SHA-256은
`2bcd49fa584c18e2415be56ed4fa75fa4894f2762b839ec9cdf689b2147c0b45`이다.
산출물은 실행 데이터의
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_extended_psf_i_fixed_masks_2026-10-21_01`
아래에 있으며 Git에 넣지 않는다.

재검출 마스크 결과에서 확산광→ICL 응답 중앙값은 0.215였으나 범위는
−2.734–0.271이었다. 기준 마스크를 고정하면 같은 응답 중앙값 0.22495,
16–84백분위 0.22175–0.22717, 전체 범위 0.22130–0.23061로 좁아졌다.
고정 마스크 확산광→위성 측정 응답 중앙값은 0.000624
(16–84백분위 −0.0116–0.0126)였다. 이는 재검출 결과의 극단값 상당 부분이
영상 내 광량 전달 자체보다 재검출·측정 대상 변경에서 왔음을 보여준다.

주요 사례인 cluster 4982에서 확산광→ICL은 재검출 −2.734, 고정 마스크 0.2213,
차이 −2.955였고, 확산광→위성 지표는 각각 4.324, −0.00287, 4.326이었다.
재검출 중 위성 표본은 4개에서 2개, ICL 면적 지원율은 0.539에서 0.498로 변했다.
cluster 411은 위성 측정에 usable satellite가 1개뿐이므로 극단 응답을 일반화하지
않는다. 위성→위성 자체 응답도 고정 마스크 중앙값 0.278이지만 표본 간 산포가
커서 보정값으로 사용하지 않는다.

이 결과는 선택된 13개 프레임의 모델 PSF 합성곱 주입 민감도일 뿐이다. 국소 PSF
보정, 확장 날개, coadd 상관 잡음과 주입 광자잡음, upstream sky 전달, 전체 표본
대표성은 검증하지 않았다. 따라서 PSF 보정값·과학 측정·ICL–위성 상관의 증거로
사용하지 말고 `science_measurement=false` 및 PSF 승인 보류를 유지한다.
새 고정 마스크 테스트를 추가하려던 마지막 편집 이후 전체 검사는 다시 실행하지
못했으므로, 앞서 통과한 9개 테스트 결과를 최종 코드 상태의 검사로 간주하지 않는다.

### 작업환경·HSC 반출 상태와 다음 순서

2026-10-21 세션에서는 `exec_command`가 Astrafek 경로, `/tmp`, `/` 및 명시적
`/bin/sh` 실행에서도 셸 프로세스를 시작하지 못하고
`CreateProcess: No such file or directory`를 반환했다. 따라서 이 세션에서
공유 지침, 로컬 `git status`, Slurm 상태 및 실행 데이터 경로를 다시 읽지 못했다.
이는 `/tmp`의 파일 권한 거부를 뜻하지 않는다. 도구가 셸을 띄우기 전에 실패한
것이며 실제 `/tmp` 권한 모드는 확인되지 않았다.

사용자는 추가 HSC PDR3 Wide HSC-I PSF 코어 1,542개 반출을 명시적으로 승인했다.
이 작업 세션에서는 셸·HSC 전송을 시작하지 못했으므로 추가 요청이나 다운로드를
완료했다고 기록하지 않는다. 재개 시 기존 manifest, 완료 상태, 파일 크기·해시를
먼저 확인하고 이미 수신된 파일을 중복 요청하지 않는다. 암호를 다시 요청하거나
문서·로그에 저장하지 않는다.

Weak-lensing 대조장용 좌표 export 프로세스 PID 1049265는 2026-10-20 확인 시
`lageunha`에서 살아 있었다. checkpoint에서 SQL batch 000과 001은 완료,
002는 CAS job 1062312에 제출된 상태다. 이 상태는 마지막 확인 시점의 것이므로
재개할 때 같은 PID와 Slurm 작업 ID의 현재 상태를 먼저 확인한다. 살아 있는 작업을
재시작하거나 중복 제출하지 않는다.

다음 순서:
1. 작업환경 복구 후 `SHARED_CONTEXT.md`, 저장소 로컬 지침, branch/status와 변경 파일을 확인한다.
2. 위 고정 마스크 테스트 추가분을 확인하고 관련 테스트·validator를 재실행한다.
3. 기존 프로세스·checkpoint 및 PSF manifest를 확인한다. 승인 범위를 넘는 반출은 하지 않는다.
4. ICL 및 위성 공동 분석의 전체 표본 마스크·광량 누출·공분산·주입 회수를 계속 검증한다.
5. 관측 weak-lensing 스택과 최종 과학 해석은 별도 검증을 마친 뒤 진행한다.
   이 인수인계 커밋에는 `HANDOVER.md`만 포함한다.

### 2026-10-22 후속: weak-lensing random-control 좌표 export 중간 감사

`syntax`에서 `lageunha`를 read-only 확인했을 때 이전 export PID 1049265는 더 이상
존재하지 않았다. 기존 HSC CAS checkpoint는 여전히 `running`으로 표시되지만, 승인된
32개 배치 중 `random_coordinates_000`과 `_001`만 완료되어 있었다. 이 두 CSV를
Astrafek 재현 감사기 `audit_wl_random_coordinate_partial_export.py`로 렌즈별
count-only 결과, 행 수, checkpoint SHA-256에 대조했다. 15,879행 및 17,635행, 총
33,514행이 모두 통과했고 기존 파일 해시도 일치한다. 승인 계획의 전체 예상은
780,191개 coordinate row다.

배치 `_002`는 CAS job `1062312`로 `submitted` 기록이며 결과 CSV는 없다. 해당 CAS
작업의 terminal 상태는 조회되지 않았으므로 재제출하거나 export를 재시작하지 않았다.
나머지 29개 배치는 아직 제출되지 않았다. 오프라인 감사 결과는
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_point_coordinate_plan_2026-10-19_01/cas_coordinates/coordinate_partial_audit_2026-10-22.json`
(SHA-256 `29c24052e09eca4102dcd2837ebdcba701087208767b82ebbe55384218b4d327`)이다.
감사기는 HSC에 연결하지 않으며 제출·다운로드·과학 측정을 하지 않는다는 것을
3개 단위시험으로 검증했다. 다음 export 단계 전에 기존 job `1062312`의 CAS 상태를
확인해야 한다. Random-lens subtraction 및 spatial-block covariance가 없으므로 weak
lensing 결과는 계속 `science_measurement=false`다.
