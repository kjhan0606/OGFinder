#!/bin/bash
# Stand-in for the python interpreter in scripts/verify_cli_templates.tcl (OGFINDER_PYTHON points here).  It records its argv in $FAKE_LOG and
# prints what the real driver would: a TSV with NUMBER + the driver's result columns (canned values, same formula as verify_cat_behavior.tcl).
echo "$@" >> "${FAKE_LOG:-/dev/null}"
cat=""; prev=""
for a in "$@"; do [ "$prev" = --catalog ] && cat="$a"; prev="$a"; done
case "$(basename "$1")" in
  ds9_sersic.py) cols="SERSIC_N SERSIC_RE SERSIC_IE SERSIC_ELLIP SERSIC_THETA SERSIC_CHI2";;
  ds9_morphometry.py) cols="CONC ASYM GINI M20 R_PETRO";;
  ds9_bulge_disk.py) cols="BT_RATIO BULGE_RE BULGE_MAG DISK_RS DISK_MAG BD_CHI2 BD_FLAG";;
  ds9_psf_phot.py) cols="FLUX_PSF FLUXERR_PSF MAG_PSF MAGERR_PSF CHI2_PSF X_PSF Y_PSF";;
  ds9_crowded_phot.py) cols="FLUX_CROWD FLUXERR_CROWD MAG_CROWD X_CROWD Y_CROWD N_NEIGHBORS";;
  ds9_multiband.py) cols="MB_MAG_G MB_MAG_R";;
  ds9_crossmatch.py) cols="MATCH_DIST MATCH_ID";;
  ds9_photo_z.py) cols="PHOTO_Z PHOTO_Z_ERR PHOTO_Z_Q68 PHOTO_Z_OUTLIER";;
  ds9_sed_fit.py) cols="LOG_MASS LOG_MASS_ERR LOG_AGE LOG_AGE_ERR LOG_Z AV SFR";;
  # drivers without --catalog: canned stdout
  ds9_segmap.py) printf 'OK 7 /workspace/fits/m51.fits\n'; exit 0;;
  ds9_completeness.py|ds9_dual_extract.py) printf 'NUMBER\tMAG_AUTO\tFRAC\n1\t21.5\t0.99\n2\t22.5\t0.95\n'; exit 0;;
  ds9_psf_deconv.py) exit 0;;
  *) cols="";;
esac
[ -n "$cat" ] || exit 0
awk -v cols="$cols" 'BEGIN{FS="\t"; n=split(cols,c," ")} NR==1{printf "NUMBER"; for(i=1;i<=n;i++) printf "\t%s",c[i]; printf "\n"; next} NF{printf "%s",$1; for(i=1;i<=n;i++) printf "\t%.3f",($1%7)+0.25*i; printf "\n"}' "$cat"
