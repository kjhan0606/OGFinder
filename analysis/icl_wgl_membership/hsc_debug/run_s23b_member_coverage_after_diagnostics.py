#!/usr/bin/env python3
"""Run local coverage sensitivity only after the exact live diagnostic run exits.

A Linux pidfd or /proc start-time identity tracks the runner's actual lifetime. A successful terminal
status and completed matching historical audit are required before analysis.
No query, download, credentials or science-image reads are performed here.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import select
import socket
import subprocess
import sys
import time


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()



def runner_identity(pid):
    try:
        fields=(Path('/proc')/str(pid)/'stat').read_text().rsplit(')',1)[1].split()
    except FileNotFoundError:
        return None
    if fields[0] == 'Z':
        return None
    return fields[19] # kernel start-time ticks; PID reuse must not extend a wait


def wait_runner(pid, identity, handle=None):
    if handle is not None:
        poll=select.poll()
        poll.register(handle,select.POLLIN)
        poll.poll()
    else:
        while runner_identity(pid) == identity:
            time.sleep(5)


def run(root, runner_pid, replicas):
    if not socket.gethostname().lower().startswith('lageunha'):
        raise ValueError('this CPU follow-up must run on lageunha')
    identity=runner_identity(runner_pid)
    if identity is None:
        raise ValueError('requested diagnostic runner is no longer live')
    opener=getattr(os,'pidfd_open',None)
    handle=opener(runner_pid) if opener is not None else None
    command=(Path('/proc')/str(runner_pid)/'cmdline').read_bytes().split(b'\0')
    if os.fsencode(root) not in command or not any(b'run_s23b_member_bin_diagnostics.sh' in part for part in command):
        if handle is not None:
            os.close(handle)
        raise ValueError('runner PID does not identify the requested diagnostic run')
    work=root/'coverage_followup'
    work.mkdir(exist_ok=False)
    status=work/'status'
    status.write_text('WAITING_FOR_DIAGNOSTIC_RUNNER\n')
    here=Path(__file__).resolve().parent
    scripts=[here/'analyze_s23b_member_coverage_sensitivity.py',
             here/'plot_s23b_member_coverage_sensitivity.py',
             here/'analyze_s23b_member_outer_icl_relation.py',
             here/'analyze_s23b_icl_spatial_blocks.py']
    pinned={str(p):sha(p) for p in scripts}
    plan=dict(schema='astrafek-member-coverage-followup-v1',run_root=str(root),
        runner_pid=runner_pid,runner_start_ticks=identity,
        wait_method='pidfd' if handle is not None else 'proc_start_time',created_utc=datetime.now(timezone.utc).isoformat(),
        script_sha256=pinned,bootstrap_replicas=replicas,
        science_measurement=False,network_used=False)
    (work/'plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    try:
        wait_runner(runner_pid,identity,handle)
        if handle is not None:
            os.close(handle)
            handle=None
        if (root/'status').read_text().strip() != 'COMPLETE':
            raise ValueError('diagnostic runner exited without successful terminal status')
        audit=json.loads((root/'audit/summary.json').read_text())
        measurement=root/'measurement'
        if (audit.get('status')!='complete_local_reconstruction_and_reference_audit'
                or audit.get('n_frames')!=979
                or audit.get('historical_profile_rows_identical') is not True
                or audit.get('measurement_summary_sha256')!=sha(measurement/'summary.json')):
            raise ValueError('full historical audit did not pass')
        if any(sha(Path(p))!=digest for p,digest in pinned.items()):
            raise ValueError('analysis sources changed after scheduling; review before executing')
        status.write_text('RUNNING_LOCAL_SENSITIVITY\n')
        env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
        with (work/'run.log').open('x') as log:
            subprocess.run([sys.executable,str(scripts[0]),'--run-root',str(root),
                '--output',str(work/'analysis'),'--replicas',str(replicas)],
                env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
            subprocess.run([sys.executable,str(scripts[1]),'--analysis',str(work/'analysis'),
                '--output',str(work/'figures')],env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
        status.write_text('COMPLETE\n')
    except BaseException as error:
        status.write_text('FAILED: '+type(error).__name__+': '+str(error)+'\n')
        raise
    finally:
        if handle is not None:
            os.close(handle)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-root',type=Path,required=True)
    parser.add_argument('--runner-pid',type=int,required=True)
    parser.add_argument('--replicas',type=int,default=5000)
    args=parser.parse_args()
    if args.replicas<100 or args.runner_pid<1:
        parser.error('positive runner PID and >=100 bootstrap replicas required')
    run(args.run_root.resolve(),args.runner_pid,args.replicas)
