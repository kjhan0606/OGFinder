#!/usr/bin/env python3
"""xshot.py DISPLAY OUT.png [X Y W H] - grab an X display (or a region) with PIL; used by the GUI checks for documentation screenshots."""
import sys
from PIL import ImageGrab
disp, out = sys.argv[1], sys.argv[2]
box = tuple(int(v) for v in sys.argv[3:7]) if len(sys.argv) >= 7 else None
img = ImageGrab.grab(xdisplay=disp, bbox=None if box is None else (box[0], box[1], box[0] + box[2], box[1] + box[3]))
img.save(out)
print(out, img.size)
