#!/usr/bin/env python3
"""Documentation checks for the dkge package.

The extractor is the contract; the rules are regexes over what it yields.
Tiers: GATE fails the build, REPORT prints a number, WARN prints for a human.
"""
import re, sys, glob, os

def extract(path):
    L = open(path, errors="ignore").read().split("\n"); i = 0
    if L and L[0].strip() == "---":
        i = 1
        while i < len(L) and L[i].strip() != "---": i += 1
        i += 1
    out, incode = [], False
    while i < len(L):
        t = L[i]
        if t.startswith("```"):
            incode = not incode; i += 1; continue
        if not incode and not t.lstrip().startswith(("|", "<")):
            out.append((i + 1, t))
        i += 1
    return out

def content(rows):
    return [(n, t) for n, t in rows
            if t.strip() and not t.lstrip().startswith(("#", ">"))]

def dequote(t):
    t = re.sub(r"`[^`]{0,200}`", " ", t)
    return re.sub(r'"[^"]{0,300}"', " ", t)

C1 = re.compile(r", not [a-z]|is not (a|an|the|evidence|proof|intended)"
                r"|does not (imply|mean|establish|prove|constitute)"
                r"|not a (test|substitute|benchmark|guarantee|proof|claim)"
                r"|rather than (a|an|the)", re.I)
# lecturer markers only: construction verbs (create/add/omit/apply/achieve) are
# example-constructor voice, which the style rule preserves.
C4 = re.compile(r"\bwe (will|begin|now|can now|need to|might want|proceed)"
                r"|\blet us\b|\ballow(s)? us to\b|\bour (exploration|purposes|design|primary)\b"
                r"|build understanding|(First|Next|Then|Finally|Now),? we", re.I)
C8 = re.compile(r"\b(localis\w*|summaris\w*|normalis\w*|behaviour\w*|colour\w*|customis\w*"
                r"|modelling|labelled|centred|visualis\w*|unnormalis\w*)\b", re.I)
NAV = re.compile(r"where to go next|next steps?|which page should you read next", re.I)

def check_file(path):
    rows = extract(path); cont = content(rows)
    r = {}
    r["C1"] = sum(len(C1.findall(dequote(t))) for _, t in cont)
    r["C4"] = sum(len(C4.findall(dequote(t))) for _, t in cont)
    r["C8"] = sum(len(C8.findall(dequote(t))) for _, t in cont)
    r["C9"] = sum(1 for _, t in cont if len(t) > 90)
    first = cont[0][1].strip() if cont else ""
    heads = [t for _, t in rows if t.lstrip().startswith("##")]
    firstline = next((t for _, t in rows if t.strip()), "")
    r["C2"] = int(bool(re.match(r"^This (vignette|page|guide|article)", first))
                  or first.startswith("**See also:")
                  or firstline.lstrip().startswith("#"))
    # C7 tests the closing section's BEHAVIOUR, not its heading text:
    # fail when the text from the last ## to EOF carries no onward link.
    idx = [i for i, (n, x) in enumerate(rows) if x.lstrip().startswith("##")]
    tail = rows[idx[-1]:] if idx else []
    has_link = any(re.search(r'vignette\("|\]\([a-z0-9-]+\.html', x) for _, x in tail)
    r["C7"] = int(bool(idx) and not has_link)
    # C12: the closing section navigates with relative .html links, never with
    # vignette() calls, which render as <code> with no anchor.
    # C12 INVERTED after the pkgdown build test: pkgdown autolinks vignette("x")
    # into a real anchor, so vignette() is correct in the navigational slot and a
    # raw .html link is the non-idiomatic form. Only list items are the slot.
    r["C12"] = sum(len(re.findall(r'\]\([a-z0-9-]+\.html\)', x))
                   for _, x in tail if x.lstrip().startswith(("- ", "* ")))
    src = open(path, errors="ignore").read().split("\n")
    incode, comments = False, 0
    for t in src:
        if t.startswith("```"): incode = not incode; continue
        if incode and t.lstrip().startswith("#") and re.search(r"\byou\b|\byour\b|in a real analysis", t, re.I):
            comments += 1
    r["C6"] = comments
    # prose lines only: the extractor already removed chunk bodies and fences
    prose_ln = {n for n, t in cont}
    figs = [i for i, t in enumerate(src, 1) if t.startswith("```{r") and "fig." in t]
    bare = 0
    for i in figs:
        if "fig.cap" in src[i - 1]: continue
        j = i + 1
        while j <= len(src) and not src[j - 1].startswith("```"): j += 1
        pb = any(n in prose_ln for n in range(max(1, i - 6), i))
        pa = any(n in prose_ln for n in range(j + 1, min(len(src), j + 7)))
        if not pb and not pa: bare += 1
    r["C5"] = bare
    return r

GATE, REPORT, WARN = ["C2", "C4", "C5", "C12"], ["C1", "C7", "C9"], ["C6", "C8"]

def main(paths):
    tot, per = {}, {}
    for p in sorted(paths):
        r = check_file(p); per[p] = r
        for k, v in r.items(): tot[k] = tot.get(k, 0) + v
    print(f"{'file':36s}" + "".join(f"{k:>6s}" for k in ["C1","C2","C4","C5","C6","C7","C8","C9","C12"]))
    for p, r in per.items():
        print(f"{os.path.basename(p):36s}" + "".join(f"{r[k]:6d}" for k in ["C1","C2","C4","C5","C6","C7","C8","C9","C12"]))
    print(f"{'TOTAL':36s}" + "".join(f"{tot[k]:6d}" for k in ["C1","C2","C4","C5","C6","C7","C8","C9","C12"]))
    print("\nGATE  ", {k: tot[k] for k in GATE})
    print("REPORT", {k: tot[k] for k in REPORT})
    print("WARN  ", {k: tot[k] for k in WARN})
    return 1 if any(tot[k] for k in GATE) else 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:] or glob.glob("vignettes/*.Rmd")))


def check_bytes(path):
    """C13: byte-level properties the line-based checks cannot see.

    Every other rule in this module operates on extracted, line-split text, so
    anything living in the bytes is invisible to all of them by construction.
    Found the hard way: a rewrap helper silently ate trailing newlines in six
    files and no check noticed.
    """
    b = open(path, "rb").read()
    out = []
    if b and not b.endswith(b"\n"): out.append("no-trailing-newline")
    if b.startswith(b"\xef\xbb\xbf"): out.append("BOM")
    if b"\r\n" in b: out.append("CRLF")
    ws = sum(1 for l in b.split(b"\n") if l.rstrip(b" \t") != l)
    if ws: out.append(f"trailing-whitespace:{ws}")
    nb = b.count(" ".encode())
    if nb: out.append(f"nbsp:{nb}")
    return out


def lost_identifiers(removed_lines, current_text):
    """Content-loss audit for a re-voicing edit.

    Scope the match to a single line: a whole-file regex for backticked spans
    runs across newlines and reports nonsense, which is how this was first
    written and why it is now documented.
    """
    ids = set()
    for line in removed_lines:
        ids.update(re.findall(r"`([^`\n]{1,60})`", line))
    return sorted(i for i in ids if i not in current_text)
