"""Small labelled sample of held-out (BB89, never trained on) detections for the real/bogus regression test:
    python make_rb_fixture.py OUT.json inj1.json.pkl inj2.json.pkl ...     (plain sets; 60 real + 240 bogus with SNR >= 6, fixed seed)"""
import json
import os
import pickle
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..'))
from moving import realbogus as RB  # noqa: E402

out, sets = sys.argv[1], sys.argv[2:]
real, bogus = [], []
for p in sets:
    dets = pickle.load(open(p, 'rb'))[0]
    for d in dets:
        if d['sign'] > 0 and d['snr'] >= 6 and 'real' in d:
            (real if d['real'] else bogus).append([float(v) for v in RB.feature_vector(d)])
rng = np.random.default_rng(0)
ir = rng.choice(len(real), 60, replace=False)
ib = rng.choice(len(bogus), 240, replace=False)
json.dump(dict(features=list(RB.FEATURES), real=[real[i] for i in sorted(ir)], bogus=[bogus[i] for i in sorted(ib)], note='held-out BB89 injection detections'), open(out, 'w'), separators=(',', ':'))
print(len(real), len(bogus), os.path.getsize(out))
