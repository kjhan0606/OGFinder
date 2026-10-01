Measured on the final build (`bin/ds9`), display :77 (Xvfb), `-geometry 1300x950`.  Screenshots were only checked through
measured geometry, pixel counts and logs, not looked at.

| check | baseline (before) | final build |
|---|---|---|
| full replay harness (`scripts/verify_session_replay.sh`) | 78 checks, 0 failed | 78 checks (T1 21, T2 30, T3 12, T4 15), 0 failed |
| step signatures of the GUI sessions (hudf 12 steps, m51 15 steps) | - | identical to the baseline (`diff` empty) |
| `ai_bridge` pytest / `moving` pytest | 49 / 12 passed | 49 / 12 passed |
| `scripts/verify_ai_gui` | 25/25 | 25/25 |
| ICL export smoke | OK | OK |
| leaf menu entries of the old 11 menus reachable | 145 | 145 of 145 (`Show Results Table` -> kind filter) |
| catalogue panel width | 620 px | 559 px |
| table y / h / info h | 181 / 769 / 154 | 181 / 769 / 154 in every state tested (6 tabs, 5 table kinds, single/tile, detached/reattached) |
| image canvas | 675x772 | 736x745 (wider panel gone, 24 px tab strip added) |
| detach / reattach | 1328x709 -> 738x709 + 585x709 -> 1328x709 | 1300x950 -> 736x950 + 559x950 -> 1300x950 |
| `layout.tcl` lines | 12958 | 12198 |
