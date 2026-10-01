#!/usr/bin/env python3
"""Rewrite direct catpanel(...) accesses of a plugin's Tcl files into ::ogf::cat accessor calls (item 5 of the decoupling work).

usage: cat_migrate.py FILE.tcl [...]        rewrites in place and prints the residual catpanel references
Rules (all keys are kept verbatim, so storage and behaviour are the same):
  ![info exists catpanel(alldata)] || $catpanel(alldata) eq {}   ->  ![::ogf::cat::has]
  $catpanel(alldata)  -> [::ogf::cat::tsv]       $catpanel(sel,nums) -> [::ogf::cat::selection]
  $catpanel(K)        -> [::ogf::cat::get K]
  set catpanel(K) W   -> ::ogf::cat::set K W     (W = one Tcl word)
  append|lappend catpanel(K) args -> ::ogf::cat::append|lappend K args
  info exists catpanel(K) -> ::ogf::cat::exists K
  array unset catpanel P  -> ::ogf::cat::unset_glob P
  unset catpanel(K)       -> ::ogf::cat::unset K
and drops 'catpanel' from `global` lines of procs that no longer mention it.
"""
import re, sys

def parse_word(s, i):
    """return end index of the Tcl word starting at s[i]"""
    n = len(s)
    if s[i] == '"':
        j = i + 1
        while j < n:
            if s[j] == '\\': j += 2; continue
            if s[j] == '[':
                j = skip_brackets(s, j); continue
            if s[j] == '"': return j + 1
            j += 1
        return n
    if s[i] == '{':
        d = 0; j = i
        while j < n:
            if s[j] == '\\': j += 2; continue
            if s[j] == '{': d += 1
            elif s[j] == '}':
                d -= 1
                if d == 0: return j + 1
            j += 1
        return n
    j = i
    while j < n and s[j] not in ' \t\n;':
        if s[j] == '\\': j += 2; continue
        if s[j] == '[':
            j = skip_brackets(s, j); continue
        if s[j] in '}]': break   # end of the enclosing braces / brackets (a bare word never contains them here)
        j += 1
    return j

def skip_brackets(s, j):
    d = 0; n = len(s)
    while j < n:
        if s[j] == '\\': j += 2; continue
        if s[j] == '[': d += 1
        elif s[j] == ']':
            d -= 1
            if d == 0: return j + 1
        j += 1
    return n

KEY = r'catpanel\(([^()\s]+)\)'   # keys never contain spaces or parentheses in this code base

def rewrite(s):
    s = re.sub(r'!\[info exists (?:::)?catpanel\(alldata\)\]\s*\|\|\s*\$(?:::)?catpanel\(alldata\)\s+eq\s+\{\}', '![::ogf::cat::has]', s)
    s = re.sub(r'\[info exists (?:::)?catpanel\(alldata\)\]\s*&&\s*\$(?:::)?catpanel\(alldata\)\s+ne\s+\{\}', '[::ogf::cat::has]', s)
    # set / append / lappend with one-word (set, append) or rest-of-command (lappend) values
    out = []; i = 0
    pat = re.compile(r'(?<![\w:$])(set|append|lappend)\s+(?:::)?' + KEY + r'(?=[ \t])')
    while True:
        m = pat.search(s, i)
        if not m:
            out.append(s[i:]); break
        out.append(s[i:m.start()])
        cmd, key = m.group(1), m.group(2)
        j = m.end()
        while j < len(s) and s[j] in ' \t': j += 1
        if cmd == 'set':
            e = parse_word(s, j); val = s[j:e]
            out.append('::ogf::cat::set %s %s' % (key, val)); i = e
        else:
            # rest of the command: words up to newline / ; / unmatched }
            e = j
            while e < len(s) and s[e] not in '\n;}':
                if s[e] in ' \t': e += 1; continue
                e = parse_word(s, e)
            out.append('::ogf::cat::%s %s %s' % (cmd, key, s[j:e].rstrip())); i = e
    s = ''.join(out)
    s = re.sub(r'\binfo exists (?:::)?' + KEY, r'::ogf::cat::exists \1', s)
    s = re.sub(r'\barray unset (?:::)?catpanel (\S+?)(?=[\s;}\]]|$)', r'::ogf::cat::unset_glob \1', s)
    s = re.sub(r'(?<![\w:])unset (?:-nocomplain )?(?:::)?' + KEY, r'::ogf::cat::unset \1', s)
    s = re.sub(r'\$(?:::)?catpanel\(alldata\)', '[::ogf::cat::tsv]', s)
    s = re.sub(r'\$(?:::)?catpanel\(sel,nums\)', '[::ogf::cat::selection]', s)
    s = re.sub(r'\$(?:::)?' + KEY, r'[::ogf::cat::get \1]', s)
    return s

def drop_globals(s):
    lines = s.split('\n')
    # proc extents: start at 'proc ' at col 0, end at the line before the next 'proc ' / EOF
    starts = [i for i, l in enumerate(lines) if l.startswith('proc ')] + [len(lines)]
    for a, b in zip(starts, starts[1:]):
        body = '\n'.join(l for l in lines[a:b] if not re.match(r'\s*global\b', l))
        has = re.search(r'\bcatpanel\b', body) is not None
        if has: continue
        for k in range(a, b):
            if re.match(r'\s*global\b', lines[k]):
                w = lines[k].split()
                w2 = [x for x in w if x != 'catpanel']
                if w2 != w:
                    lines[k] = None if len(w2) == 1 else re.match(r'\s*', lines[k]).group(0) + ' '.join(w2)
    return '\n'.join(l for l in lines if l is not None)

if __name__ == '__main__':
    for f in sys.argv[1:]:
        t = open(f, encoding='utf-8', errors='surrogateescape').read()
        t2 = drop_globals(rewrite(t))
        open(f, 'w', encoding='utf-8', errors='surrogateescape').write(t2)
        res = [(n + 1, l.strip()) for n, l in enumerate(t2.split('\n')) if re.search(r'\bcatpanel\b', l)]
        print('%s: %d residual catpanel lines' % (f, len(res)))
        for n, l in res: print('   %d: %s' % (n, l[:140]))
