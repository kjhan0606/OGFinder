set -e
P=/workspace/ogf_venv/bin/python
D=/workspace/OGFinder/ds9/library/ds9_moving.py
W=/workspace/work/wd_real
$P $D --mode align --workdir $W --files /home/box/.ds9/mast_cache/j8pu38c7q_flc.fits /home/box/.ds9/mast_cache/j8pu38caq_flc.fits /home/box/.ds9/mast_cache/j8pu38ceq_flc.fits /home/box/.ds9/mast_cache/j8pu38ciq_flc.fits 
$P $D --mode difference --workdir $W --files /home/box/.ds9/mast_cache/j8pu38c7q_flc.fits /home/box/.ds9/mast_cache/j8pu38caq_flc.fits /home/box/.ds9/mast_cache/j8pu38ceq_flc.fits /home/box/.ds9/mast_cache/j8pu38ciq_flc.fits  --ra 150.1375 --dec 2.3610 --half-pix 700 --snr 5
$P $D --mode link --workdir $W --files /home/box/.ds9/mast_cache/j8pu38c7q_flc.fits /home/box/.ds9/mast_cache/j8pu38caq_flc.fits /home/box/.ds9/mast_cache/j8pu38ceq_flc.fits /home/box/.ds9/mast_cache/j8pu38ciq_flc.fits  --snr 8 --tol 1.0 --max-per-exposure 900 --max-tracklets 400
F="/home/box/.ds9/mast_cache/j8pu38c7q_flc.fits /home/box/.ds9/mast_cache/j8pu38caq_flc.fits /home/box/.ds9/mast_cache/j8pu38ceq_flc.fits /home/box/.ds9/mast_cache/j8pu38ciq_flc.fits"
$P $D --mode identify --workdir $W --radius 6 --tol 1.5 --max-tracklets 15 --files $F
