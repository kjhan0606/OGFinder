#!/usr/bin/env python3
"""
OGFinder catalog edit helper (standard library only).

Re-implements, string for string, the catalog manipulations done by the
Tcl catalog panel in ds9/library/layout.tcl and ogf_bands.tcl so that a GUI
session can be replayed from the shell:

  CatalogPanelAddColumnsFromTSV  -> add-columns
  OGFBandsRenameCols             -> forced-rename
  OGFBandsAddColors              -> band-colors
  CatalogPanelTrimApply          -> apply (op "trim")
  CatalogPanelSort               -> apply (op "sort")
  CatalogPanelMergeSources       -> apply (op "merge")
  CatalogPanelDeleteSelected     -> apply (op "delete")
  CatalogPanelMorphParseResults/MorphAddColumns -> morph-columns
  (ICL "BCG" default = brightest MAG_AUTO)       -> brightest
  CatalogPanelLoadTSV            -> set  (needs >= 2 lines, else keeps old catalog)

A "catalog" is the exact text held in catpanel(alldata): TSV, rows separated by
"\\n", no trailing newline.  All subcommands read/write such files.

Usage:
  ds9_catalog_edit.py set          --catalog C --input RAW [--strip-newline]
  ds9_catalog_edit.py add-columns  --catalog C --result R --columns A,B [--strip-comments]
  ds9_catalog_edit.py forced-rename --input RAW --band B --output OUT
  ds9_catalog_edit.py band-colors  --catalog C --order A,B,C
  ds9_catalog_edit.py apply        --catalog C --edits edits.json [--lenient]
  ds9_catalog_edit.py morph-columns --catalog C --result R
  ds9_catalog_edit.py brightest    --catalog C
  ds9_catalog_edit.py count        --catalog C
"""
import argparse
import functools
import json
import math
import os
import re
import sys

_DBL = re.compile(r'^[+-]?(?:\d+\.?\d*(?:[eE][+-]?\d+)?|\.\d+(?:[eE][+-]?\d+)?|inf(?:inity)?|nan)$',
                  re.I)


def is_double(s):
    """Tcl: string is double -strict"""
    return bool(_DBL.match(s.strip())) if s.strip() != '' else False


def tcl_trim(s):
    return s.strip(' \t\n\r\x0b\x0c')


def lindex(seq, i):
    if i < 0 or i >= len(seq):
        return ''
    return seq[i]


def read_text(path):
    with open(path, 'r', newline='', encoding='utf-8') as f:
        return f.read()


def write_text(path, text):
    tmp = path + '.tmp'
    with open(tmp, 'w', newline='', encoding='utf-8') as f:
        f.write(text)
    os.replace(tmp, path)


def nonblank(line):
    return tcl_trim(line) != ''


def header_index(headers, name):
    for i, h in enumerate(headers):
        if tcl_trim(h) == name:
            return i
    return -1


# ------------------------------------------------------------ add columns
def add_columns(alldata, result_data, col_names):
    rlines = result_data.split('\n')
    rheaders = rlines[0].split('\t')
    r_num_col = header_index(rheaders, 'NUMBER')
    if r_num_col < 0:
        return alldata
    rdata = {}
    for line in rlines[1:]:
        if not nonblank(line):
            continue
        fields = line.split('\t')
        num = tcl_trim(lindex(fields, r_num_col))
        vals = []
        for cn in col_names:
            cidx = header_index(rheaders, cn)
            if cidx >= 0 and cidx < len(fields):
                vals.append(tcl_trim(fields[cidx]))
            else:
                vals.append('')
        rdata[num] = vals
    lines = alldata.split('\n')
    headers = lines[0].split('\t')
    num_col = header_index(headers, 'NUMBER')
    if num_col < 0:
        return alldata
    existing = any(header_index(headers, cn) >= 0 for cn in col_names)
    out = []
    if existing:
        col_indices = [header_index(headers, cn) for cn in col_names]
        out.append(lines[0])
        for line in lines[1:]:
            if not nonblank(line):
                continue
            fields = line.split('\t')
            sn = tcl_trim(lindex(fields, num_col))
            if sn in rdata:
                vals = rdata[sn]
                for v in range(len(col_names)):
                    ci = col_indices[v]
                    if ci >= 0:
                        if ci >= len(fields):
                            raise IndexError('list index out of range (row too short)')
                        fields[ci] = vals[v]
            out.append('\n' + '\t'.join(fields))
    else:
        head = lines[0] + ''.join('\t' + cn for cn in col_names)
        out.append(head)
        for line in lines[1:]:
            if not nonblank(line):
                continue
            fields = line.split('\t')
            sn = tcl_trim(lindex(fields, num_col))
            row = '\n' + line
            if sn in rdata:
                row += ''.join('\t' + v for v in rdata[sn])
            else:
                row += '\t' * len(col_names)
            out.append(row)
    return ''.join(out)


def strip_star_finder(raw):
    """CatalogPanelStarFinderParse: header_line + result lines (drops '#' lines)"""
    header_line = ''
    res = []
    for line in raw.split('\n'):
        if line.startswith('#'):
            continue
        if line.startswith('NUMBER'):
            header_line = line
            continue
        if tcl_trim(line) == '':
            continue
        res.append(line)
    return '\n'.join([header_line] + res), len(res)


# ----------------------------------------------------------- forced bands
def forced_rename(tsv, b):
    lines = tsv.split('\n')
    h = lines[0].split('\t')
    ci_n = h.index('NUMBER') if 'NUMBER' in h else -1
    ci_m = h.index('MAG_AUTO_' + b) if ('MAG_AUTO_' + b) in h else -1
    ci_e = h.index('MAGERR_AUTO_' + b) if ('MAGERR_AUTO_' + b) in h else -1
    out = 'NUMBER\tMAG_%s\tMAGERR_%s' % (b, b)
    for line in lines[1:]:
        if not nonblank(line):
            continue
        f = line.split('\t')
        out += '\n%s\t%s\t%s' % (lindex(f, ci_n), lindex(f, ci_m), lindex(f, ci_e))
    return out


def band_colors(alldata, order):
    if len(order) < 2:
        return alldata
    lines = alldata.split('\n')
    h = lines[0].split('\t')
    pairs = []
    for i in range(len(order) - 1):
        a, b = order[i], order[i + 1]
        ia = h.index('MAG_' + a) if ('MAG_' + a) in h else -1
        ib = h.index('MAG_' + b) if ('MAG_' + b) in h else -1
        if ia >= 0 and ib >= 0:
            pairs.append((a, b, ia, ib))
    if not pairs:
        return alldata
    ni = h.index('NUMBER') if 'NUMBER' in h else -1
    res = 'NUMBER'
    names = []
    for p in pairs:
        res += '\t%s-%s' % (p[0], p[1])
        names.append('%s-%s' % (p[0], p[1]))
    for line in lines[1:]:
        if not nonblank(line):
            continue
        f = line.split('\t')
        res += '\n' + lindex(f, ni)
        for a, b, ia, ib in pairs:
            ma, mb = lindex(f, ia), lindex(f, ib)
            if is_double(ma) and is_double(mb) and float(ma) < 90 and float(mb) < 90:
                res += '\t%.3f' % (float(ma) - float(mb))
            else:
                res += '\t99.000'
    return add_columns(alldata, res, names)


# ------------------------------------------------------------- morphology
def morph_columns(alldata, raw):
    """CatalogPanelMorphParseResults + CatalogPanelMorphAddColumns"""
    morph = {}
    order = []
    for line in raw.split('\n'):
        if line.startswith('#') or line.startswith('NUMBER') or tcl_trim(line) == '':
            continue
        f = line.split('\t')
        if len(f) < 11:
            continue
        morph[f[0]] = (f[1], f[2], f[3])
        order.append(f[0])
    if not order:
        return alldata
    lines = alldata.split('\n')
    if len(lines) < 2:
        return alldata
    headers = lines[0].split('\t')
    col_num = header_index(headers, 'NUMBER')
    if col_num < 0:
        return alldata
    has = header_index(headers, 'MORPH_TYPE') >= 0
    out = []
    if has:
        c_t, c_d, c_c = (header_index(headers, k) for k in ('MORPH_TYPE', 'MORPH_DESC', 'MORPH_CONF'))
        out.append(lines[0])
        for line in lines[1:]:
            if not nonblank(line):
                continue
            f = line.split('\t')
            sn = tcl_trim(lindex(f, col_num))
            if sn in morph:
                info = morph[sn]
                if c_t >= 0: f[c_t] = info[0]
                if c_d >= 0: f[c_d] = info[1]
                if c_c >= 0: f[c_c] = info[2]
            out.append('\n' + '\t'.join(f))
    else:
        out.append(lines[0] + '\tMORPH_TYPE\tMORPH_DESC\tMORPH_CONF')
        for line in lines[1:]:
            if not nonblank(line):
                continue
            f = line.split('\t')
            sn = tcl_trim(lindex(f, col_num))
            mt = md = mc = ''
            if sn in morph:
                mt, md, mc = morph[sn]
            out.append('\n' + line + '\t' + mt + '\t' + md + '\t' + mc)
    return ''.join(out)


# ----------------------------------------------------------------- edits
def op_trim(alldata, conditions, lenient, log):
    lines = alldata.split('\n')
    header = lines[0]
    headers = header.split('\t')
    conds = []
    for c in conditions:
        idx = header_index(headers, c['col'])
        if idx < 0:
            msg = 'trim: column %s not in catalog' % c['col']
            if lenient:
                log('WARNING: ' + msg + ' (condition skipped)')
                continue
            raise KeyError(msg)
        has_min = c.get('min') not in (None, '')
        has_max = c.get('max') not in (None, '')
        conds.append((idx, has_min, float(c['min']) if has_min else 0.0,
                      has_max, float(c['max']) if has_max else 0.0))
    if not conds:
        return alldata            # "Trim cleared - showing all sources"
    filtered = [header]
    for line in lines[1:]:
        if not nonblank(line):
            continue
        fields = line.split('\t')
        ok = True
        for idx, has_min, mn, has_max, mx in conds:
            val = tcl_trim(lindex(fields, idx))
            if not is_double(val):
                ok = False
                break
            v = float(val)
            if has_min and v < mn:
                ok = False
                break
            if has_max and v > mx:
                ok = False
                break
        if ok:
            filtered.append(line)
    new = '\n'.join(filtered)
    # CatalogPanelLoadTSV refuses data with fewer than two lines ("No sources detected")
    if len(new.split('\n')) < 2:
        log('WARNING: trim left no sources; catalog unchanged (same as the GUI)')
        return alldata
    return new


def op_sort(alldata, col, direction, lenient, log):
    lines = alldata.split('\n')
    header = lines[0]
    headers = header.split('\t')
    colidx = header_index(headers, col)
    if colidx < 0:
        msg = 'sort: column %s not in catalog' % col
        if lenient:
            log('WARNING: ' + msg + ' (skipped)')
            return alldata
        raise KeyError(msg)
    rows = [l for l in lines[1:] if nonblank(l)]
    isnum = True
    for r in rows:
        val = tcl_trim(lindex(r.split('\t'), colidx))
        if val != '':
            if not is_double(val):
                isnum = False
            break

    def cmp_num(a, b):
        va = tcl_trim(lindex(a.split('\t'), colidx))
        vb = tcl_trim(lindex(b.split('\t'), colidx))
        va = float(va) if is_double(va) else 0.0
        vb = float(vb) if is_double(vb) else 0.0
        return -1 if va < vb else (1 if va > vb else 0)

    def cmp_str(a, b):
        va = tcl_trim(lindex(a.split('\t'), colidx))
        vb = tcl_trim(lindex(b.split('\t'), colidx))
        return -1 if va < vb else (1 if va > vb else 0)

    cmp = cmp_num if isnum else cmp_str
    if direction == 'descending':
        key = functools.cmp_to_key(lambda a, b: cmp(b, a))
    else:
        key = functools.cmp_to_key(cmp)
    srt = sorted(rows, key=key)
    return '\n'.join([header] + srt)


def _fnum(s, default=0):
    return float(s) if is_double(s) else default


def op_merge(alldata, nums, mag_zp, lenient, log):
    nums = [str(n) for n in nums]
    lines = alldata.split('\n')
    header = lines[0]
    headers = header.split('\t')
    ncols = len(headers)
    idx = {k: -1 for k in ('NUMBER', 'X_IMAGE', 'Y_IMAGE', 'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE',
                           'ISO_RADIUS', 'FLUX_AUTO', 'MAG_AUTO', 'NPIX_ISO')}
    for i, h in enumerate(headers):
        h = tcl_trim(h)
        if h in idx:
            idx[h] = i                       # last occurrence wins, like the Tcl switch loop
    i_num, i_x, i_y, i_a, i_b, i_th = (idx[k] for k in ('NUMBER', 'X_IMAGE', 'Y_IMAGE', 'A_IMAGE',
                                                       'B_IMAGE', 'THETA_IMAGE'))
    i_ir, i_flux, i_mag, i_npix = idx['ISO_RADIUS'], idx['FLUX_AUTO'], idx['MAG_AUTO'], idx['NPIX_ISO']
    merge_rows, other_rows = [], []
    max_number = 0
    brightest_idx, brightest_flux = -1, -1e30
    for line in lines[1:]:
        if not nonblank(line):
            continue
        fields = line.split('\t')
        num_val = tcl_trim(lindex(fields, i_num))
        if re.match(r'^[+-]?\d+$', num_val) and int(num_val) > max_number:
            max_number = int(num_val)
        if num_val in nums:
            merge_rows.append(fields)
            if i_flux >= 0:
                fv = tcl_trim(lindex(fields, i_flux))
                if is_double(fv) and float(fv) > brightest_flux:
                    brightest_flux = float(fv)
                    brightest_idx = len(merge_rows) - 1
        else:
            other_rows.append(line)
    if len(merge_rows) < 2:
        msg = 'merge: sources %s not found in catalog' % ','.join(nums)
        if lenient:
            log('WARNING: ' + msg + ' (skipped)')
            return alldata
        raise KeyError(msg)
    if brightest_idx < 0:
        brightest_idx = 0
    new_num = max_number + 1

    def getflux(row):
        flux = 1.0
        if i_flux >= 0:
            fv = tcl_trim(lindex(row, i_flux))
            if is_double(fv) and float(fv) > 0:
                flux = float(fv)
        return flux

    total_flux = 0.0
    wx = wy = 0.0
    total_npix = 0
    for row in merge_rows:
        flux = getflux(row)
        x = _fnum(tcl_trim(lindex(row, i_x)))
        y = _fnum(tcl_trim(lindex(row, i_y)))
        total_flux = total_flux + flux
        wx = wx + x * flux
        wy = wy + y * flux
        if i_npix >= 0:
            nv = tcl_trim(lindex(row, i_npix))
            if re.match(r'^[+-]?\d+$', nv):
                total_npix += int(nv)
    if total_flux <= 0:
        total_flux = 1.0
    new_x = wx / total_flux
    new_y = wy / total_flux
    Ixx = Iyy = Ixy = 0.0
    PI = 3.14159265358979
    for row in merge_rows:
        flux = getflux(row)
        x = _fnum(tcl_trim(lindex(row, i_x)))
        y = _fnum(tcl_trim(lindex(row, i_y)))
        ak = bk = tk = 0.0
        if i_a >= 0:
            ak = _fnum(tcl_trim(lindex(row, i_a)), 0.0)
        if i_b >= 0:
            bk = _fnum(tcl_trim(lindex(row, i_b)), 0.0)
        if i_th >= 0:
            tk = _fnum(tcl_trim(lindex(row, i_th)), 0.0)
        rad = tk * PI / 180.0
        cosT, sinT = math.cos(rad), math.sin(rad)
        a2, b2 = ak * ak, bk * bk
        ixx_k = a2 * cosT * cosT + b2 * sinT * sinT
        iyy_k = a2 * sinT * sinT + b2 * cosT * cosT
        ixy_k = (a2 - b2) * sinT * cosT
        dx, dy = x - new_x, y - new_y
        Ixx = Ixx + flux * (ixx_k + dx * dx)
        Iyy = Iyy + flux * (iyy_k + dy * dy)
        Ixy = Ixy + flux * (ixy_k + dx * dy)
    Ixx /= total_flux
    Iyy /= total_flux
    Ixy /= total_flux
    trace = Ixx + Iyy
    disc = math.sqrt(abs((Ixx - Iyy) * (Ixx - Iyy) + 4.0 * Ixy * Ixy))
    lam1 = (trace + disc) / 2.0
    lam2 = (trace - disc) / 2.0
    if lam1 < 0:
        lam1 = 0.0
    if lam2 < 0:
        lam2 = 0.0
    new_a, new_b = math.sqrt(lam1), math.sqrt(lam2)
    new_theta = 0.5 * math.atan2(2.0 * Ixy, Ixx - Iyy) * 180.0 / PI
    zp = float(mag_zp) if is_double(str(mag_zp)) else 25.0
    tf = total_flux
    new_mag = -2.5 * math.log10(tf) + zp
    new_ir = 5.0
    if total_npix > 0 and new_a > 0 and new_b > 0:
        ratio = new_b / new_a
        if ratio <= 0:
            ratio = 1.0
        new_ir = math.sqrt(total_npix / (3.14159265 * ratio))
    base_row = merge_rows[brightest_idx]
    nf = []
    for c in range(ncols):
        val = tcl_trim(lindex(base_row, c))
        if c == i_num: val = str(new_num)
        if c == i_x: val = '%.4f' % new_x
        if c == i_y: val = '%.4f' % new_y
        if c == i_flux and i_flux >= 0: val = '%.6g' % total_flux
        if c == i_mag and i_mag >= 0: val = '%.4f' % new_mag
        if c == i_npix and i_npix >= 0: val = str(total_npix)
        if c == i_ir and i_ir >= 0: val = '%.4f' % new_ir
        if c == i_a and i_a >= 0: val = '%.4f' % new_a
        if c == i_b and i_b >= 0: val = '%.4f' % new_b
        if c == i_th and i_th >= 0: val = '%.4f' % new_theta
        nf.append(val)
    return '\n'.join([header] + other_rows + ['\t'.join(nf)])


def op_delete(alldata, number, lenient, log):
    lines = alldata.split('\n')
    header = lines[0]
    headers = header.split('\t')
    num_idx = header_index(headers, 'NUMBER')
    new_lines = [header]
    deleted = False
    for line in lines[1:]:
        if line == '':
            continue
        fields = line.split('\t')
        this_num = ''
        if num_idx >= 0 and len(fields) > num_idx:
            this_num = tcl_trim(fields[num_idx])
        if this_num == tcl_trim(str(number)):
            deleted = True
        else:
            new_lines.append(line)
    if not deleted:
        msg = 'delete: source %s not found' % number
        if lenient:
            log('WARNING: ' + msg)
            return alldata
        raise KeyError(msg)
    return '\n'.join(new_lines)


def apply_edits(alldata, ops, lenient, log):
    for op in ops:
        k = op['op']
        if k == 'trim':
            alldata = op_trim(alldata, op['conditions'], lenient, log)
        elif k == 'sort':
            alldata = op_sort(alldata, op['col'], op['dir'], lenient, log)
        elif k == 'merge':
            alldata = op_merge(alldata, op['nums'], op.get('mag_zp', '25.0'), lenient, log)
        elif k == 'delete':
            alldata = op_delete(alldata, op['number'], lenient, log)
        else:
            raise ValueError('unknown edit op %r' % k)
    return alldata


def brightest(alldata):
    """The default the ICL 'Set BCG Center' dialog proposes: smallest MAG_AUTO (< 99).
    Returns (number, x_image, y_image) as strings, or None."""
    lines = alldata.split('\n')
    headers = lines[0].split('\t')
    ncol, xc, yc, mc = (header_index(headers, k) for k in ('NUMBER', 'X_IMAGE', 'Y_IMAGE', 'MAG_AUTO'))
    if min(ncol, xc, yc) < 0:
        return None
    best, best_mag = None, 99.0
    for line in lines[1:]:
        row = line.split('\t')
        if len(row) <= xc:
            continue
        if mc >= 0 and len(row) > mc:
            m = tcl_trim(row[mc])
            if is_double(m) and float(m) < best_mag:
                best_mag = float(m)
                best = (tcl_trim(row[ncol]), tcl_trim(row[xc]), tcl_trim(row[yc]))
    return best


def count_rows(alldata):
    return sum(1 for l in alldata.split('\n')[1:] if nonblank(l))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('cmd', choices=['set', 'add-columns', 'forced-rename', 'band-colors', 'apply',
                                    'morph-columns', 'brightest', 'count'])
    ap.add_argument('--catalog')
    ap.add_argument('--input')
    ap.add_argument('--result')
    ap.add_argument('--columns')
    ap.add_argument('--band')
    ap.add_argument('--output')
    ap.add_argument('--order')
    ap.add_argument('--edits')
    ap.add_argument('--strip-newline', action='store_true',
                    help='drop one trailing newline (what Tcl exec does)')
    ap.add_argument('--strip-comments', action='store_true')
    ap.add_argument('--force', action='store_true',
                    help='set: replace the catalog even when the result has < 2 lines (LSBG steps)')
    ap.add_argument('--lenient', action='store_true')
    a = ap.parse_args()

    def log(m):
        print(m, file=sys.stderr)

    if a.cmd == 'set':
        raw = read_text(a.input)
        if a.strip_newline and raw.endswith('\n'):
            raw = raw[:-1]
        if len(raw.split('\n')) < 2 and not a.force:
            log('WARNING: fewer than 2 lines in result; catalog not replaced (GUI: "No sources detected")')
            print('#CATALOG ROWS=0 REPLACED=0')
            return 0
        write_text(a.catalog, raw)
        print('#CATALOG ROWS=%d REPLACED=1' % count_rows(raw))
    elif a.cmd == 'add-columns':
        res = read_text(a.result)
        if a.strip_newline and res.endswith('\n'):
            res = res[:-1]
        if a.strip_comments:
            res, nres = strip_star_finder(res)
            if nres == 0:
                log('WARNING: no classified sources in result; catalog unchanged (same as the GUI)')
                print('#CATALOG ROWS=%d' % count_rows(read_text(a.catalog)))
                return 0
        cols = a.columns.split(',')
        new = add_columns(read_text(a.catalog), res, cols)
        write_text(a.catalog, new)
        print('#CATALOG ROWS=%d COLS=%d' % (count_rows(new), len(new.split('\n')[0].split('\t'))))
    elif a.cmd == 'forced-rename':
        raw = read_text(a.input)
        if raw.endswith('\n'):
            raw = raw[:-1]
        write_text(a.output, forced_rename(raw, a.band))
    elif a.cmd == 'band-colors':
        new = band_colors(read_text(a.catalog), a.order.split(','))
        write_text(a.catalog, new)
        print('#CATALOG ROWS=%d COLS=%d' % (count_rows(new), len(new.split('\n')[0].split('\t'))))
    elif a.cmd == 'apply':
        with open(a.edits) as f:
            ed = json.load(f)
        new = apply_edits(read_text(a.catalog), ed['ops'], a.lenient, log)
        write_text(a.catalog, new)
        print('#CATALOG ROWS=%d COLS=%d' % (count_rows(new), len(new.split('\n')[0].split('\t'))))
    elif a.cmd == 'morph-columns':
        res = read_text(a.result)
        if res.endswith('\n'):
            res = res[:-1]
        new = morph_columns(read_text(a.catalog), res)
        write_text(a.catalog, new)
        print('#CATALOG ROWS=%d' % count_rows(new))
    elif a.cmd == 'brightest':
        b = brightest(read_text(a.catalog))
        if b is None:
            print('NONE')
            return 2
        print('%s\t%s\t%s' % b)
    elif a.cmd == 'count':
        print(count_rows(read_text(a.catalog)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
