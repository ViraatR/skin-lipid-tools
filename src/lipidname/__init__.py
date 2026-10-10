"""Parse lipid annotations from MS-DIAL, LIPID MAPS shorthand, RefMet and Skyline.

Everything downstream (enzyme-activity proxies in insilico_pace, RT models in
lipid_rt_qc) works from the structured `Lipid` record, so name handling lives in
one place.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

# MS-DIAL ceramide subclass codes -> (sphingoid base type, acyl hydroxylation)
#   base: S = sphingosine (d18:1), DS = dihydrosphingosine (d18:0),
#         P = phytosphingosine (t18:0), H = 6-hydroxysphingosine
#   acyl: N = non-hydroxy, A = alpha-OH, B = beta-OH, O = omega-OH,
#         EO = esterified omega-OH (acylceramide)
CER_SUBCLASSES = {
    "NS", "NDS", "NP", "NH", "AS", "ADS", "AP", "AH", "BS", "BDS",
    "OS", "ODS", "EOS", "EODS", "EOP", "EOH",
}

# Aliases seen in MS-DIAL / RefMet / vendor exports -> canonical class
CLASS_ALIASES = {
    "LYSOPC": "LPC", "LYSOPE": "LPE", "LYSOPS": "LPS", "LYSOPI": "LPI",
    "LYSOPG": "LPG", "LYSOPA": "LPA",
    "PLASMENYL-PE": "PE-P", "PLASMENYL-PC": "PC-P",
    "PLASMANYL-PE": "PE-O", "PLASMANYL-PC": "PC-O",
    "ETHERPE": "PE-O", "ETHERPC": "PC-O",
    "ACYLCARNITINE": "CAR", "CHOLESTERYL ESTER": "CE", "CHE": "CE",
    "HEXCER": "HexCer", "GLCCER": "HexCer", "GALCER": "HexCer",
    "HEX2CER": "Hex2Cer", "LACCER": "Hex2Cer",
    "TAG": "TG", "DAG": "DG", "MAG": "MG",
    "FFA": "FA", "SPH": "SPB", "SO": "SPB", "SA": "SPB",
}

CHAIN_RE = re.compile(r"(?:(?P<prefix>[dtmOP])-?)?(?P<c>\d{1,2}):(?P<db>\d{1,2})(?:;(?:O(?P<o>\d?)|(?P<o2>\d)O))?")
SUM_RE = re.compile(r"^(?P<c>\d{1,3}):(?P<db>\d{1,2})(?:;(?:O(?P<o>\d?)|(?P<o2>\d)O))?")


@dataclass
class Lipid:
    raw: str
    lipid_class: str            # canonical, e.g. "Cer_NS", "PC", "PE-P", "TG"
    total_c: int | None = None
    total_db: int | None = None
    total_oh: int = 0           # extra oxygens beyond the class backbone
    chains: list[tuple[int, int]] = field(default_factory=list)  # acyl/alkyl chains
    base: tuple[int, int] | None = None   # sphingoid base (C, DB) for sphingolipids
    adduct: str | None = None
    labelled: bool = False      # isotope-labelled internal standard, e.g. '(d7)'

    @property
    def is_sphingolipid(self) -> bool:
        return self.lipid_class.startswith(("Cer", "SM", "HexCer", "Hex2Cer", "SPB"))

    @property
    def acyl_c(self) -> int | None:
        """Carbon count of the N-acyl chain for sphingolipids (total minus base)."""
        if self.base and self.total_c is not None:
            return self.total_c - self.base[0]
        return None

    @property
    def ecn(self) -> float | None:
        """Equivalent carbon number, the classic reverse-phase RT predictor."""
        if self.total_c is None or self.total_db is None:
            return None
        return self.total_c - 2 * self.total_db

    @property
    def key(self) -> str:
        """Species-level name, e.g. 'Cer_NS 42:1'."""
        if self.total_c is None:
            return self.lipid_class
        return f"{self.lipid_class} {self.total_c}:{self.total_db}"


def _canon_class(token: str) -> str:
    t = token.strip()
    up = t.upper()
    if up in CLASS_ALIASES:
        return CLASS_ALIASES[up]
    m = re.fullmatch(r"(Cer|HexCer|GlcCer|GalCer|Hex2Cer)\[(\w+)\]", t, flags=re.I)
    if m:
        head = "Cer" if m.group(1).upper() == "CER" else ("Hex2Cer" if m.group(1).upper() == "HEX2CER" else "HexCer")
        return f"{head}_{m.group(2).upper()}"
    return t


def _oh(m: re.Match) -> int:
    if m.group("o2"):
        return int(m.group("o2"))
    if m.group("o") is not None:
        return int(m.group("o") or 1)
    return 0


def _parse_chains(text: str) -> tuple[list[tuple[int, int]], tuple[int, int] | None, int]:
    """Parse '(d18:1/24:0)' or '18:1;O2/24:0' style chain lists."""
    chains, base, oh = [], None, 0
    parts = re.split(r"[/_]", text.strip("() "))
    for i, p in enumerate(parts):
        # omega-esterified acylceramide chain like '32:0-O-18:2'
        for j, sub in enumerate(re.split(r"-O-", p)):
            m = CHAIN_RE.search(sub)
            if not m:
                continue
            c, db = int(m.group("c")), int(m.group("db"))
            pre = m.group("prefix")
            o = _oh(m)
            if i == 0 and j == 0 and (pre in ("d", "t", "m") or ";O" in sub or "O2" in sub or "O3" in sub):
                base = (c, db)
                o_base = {"m": 1, "d": 2, "t": 3}.get(pre, o)
                oh += max(o_base - 2, 0)  # d-bases carry the 2 backbone OH
            else:
                chains.append((c, db))
                oh += o
    return chains, base, oh


def parse(name: str) -> Lipid | None:
    """Parse one annotation. Returns None for non-lipids / unparseable names."""
    if not isinstance(name, str) or not name.strip():
        return None
    raw = name.strip()
    s = raw.rstrip("*").strip()
    adduct = None

    # MS-DIAL 5: 'Cer 42:1;O2|Cer 18:1;O2/24:0' (sum level | molecular species)
    if "|" in s:
        sum_part, detail = (x.strip() for x in s.split("|", 1))
        l = parse(detail) or parse(sum_part)
        if l is not None:
            l.raw = raw
        return l

    # Isotope-labelled standards: 'DG 15:0_18:1(d7)', 'PC 15:0-18:1(d7)'
    labelled = bool(re.search(r"\(d\d+\)|-d\d+\b|\[D\d+\]", s))
    s = re.sub(r"\(d\d+\)|\[D\d+\]", "", s)
    # Avanti-style hyphenated chains 'PC 15:0-18:1' mean sn-1/sn-2
    s = re.sub(r"(\d+:\d+)-(?=\d+:\d+)", r"\1/", s)

    # Estolide / FAHFA-esterified glycerolipids: 'TG 18:1_18:2_18:1;O(FA 16:0)'
    # ...and acyl-sphingomyelin 'SM 34:1;O2(FA 18:0)' (MS-DIAL 'ASM')
    estolide = re.search(r";?O?\(FA (\d+):(\d+)\)", s)
    if estolide and not s.startswith(("Cer", "HexCer")):
        s = s[:estolide.start()] + s[estolide.end():]
        l = parse(s)
        if l is None:
            return None
        ec, ed = int(estolide.group(1)), int(estolide.group(2))
        l.raw, l.lipid_class, l.labelled = raw, f"{l.lipid_class}_EST", labelled
        if not l.chains or sum(c for c, _ in l.chains) + (l.base[0] if l.base and l.chains else 0) == l.total_c:
            l.total_c, l.total_db = l.total_c + ec, l.total_db + ed
        l.chains = l.chains + [(ec, ed)]
        l.total_oh += 1
        return l

    # MS-DIAL 5 hydroxy / acylceramide notation -> forms the chain parser knows
    #   'Cer 18:1;O2/24:0;(2OH)' alpha-OH, ';(3OH)' beta-OH, ';(FA 18:2)' omega-esterified
    s = re.sub(r";\(2OH\)", ";O", s)
    s = re.sub(r";\(3OH\)", ";OB", s)
    s = re.sub(r";\(FA (\d+:\d+)\)", r"-O-\1", s)

    # MS-DIAL 4: 'Cer[NS] 50:3; Cer[NS](d18:1/32:2); [M+Hac-H]-'
    segs = [x.strip() for x in s.split(";")] if "; " in s else [s]
    if len(segs) > 1:
        if segs[-1].startswith("["):
            adduct = segs[-1]
            segs = segs[:-1]
        # keep the most detailed segment (the one with chains), fall back to the sum
        s = next((x for x in reversed(segs) if "(" in x or "/" in x or "_" in x), segs[0])
        s_sum = segs[0]
    else:
        s_sum = s

    # '20:1 Cholesteryl ester' (number first)
    m = re.fullmatch(r"(\d+:\d+)\s+(Cholesteryl ester|CE)", s_sum, flags=re.I)
    if m:
        c, db = map(int, m.group(1).split(":"))
        return Lipid(raw, "CE", c, db, chains=[(c, db)], adduct=adduct)

    # Class token, then either '(chains)', ' chains' or ' sum'
    m = re.match(r"^\s*([A-Za-z][\w\-]*(?:\[[A-Za-z]+\])?|Cholesteryl ester)\s*(.*)$", s)
    if not m:
        return None
    cls_tok, rest = m.group(1), m.group(2).strip()

    # 'PE(P-18:0/22:6)' or 'PE P-18:0/22:6' -> plasmalogen / ether classes
    ether = re.match(r"^\(?\s*([OP])-", rest)
    lipid_class = _canon_class(cls_tok)
    if ether and lipid_class in ("PC", "PE", "PS", "PG", "PI", "LPC", "LPE", "TG", "DG", "MG"):
        lipid_class = f"{lipid_class}-{ether.group(1)}"
        rest = rest.replace(f"{ether.group(1)}-", "", 1)

    # MS-DIAL 5 / LIPID MAPS: 'Cer 18:1;O2/24:0' -> infer ceramide subclass
    if lipid_class in ("Cer", "HexCer"):
        lipid_class = lipid_class + _infer_cer_subclass(rest)[3:]

    chains, base, oh = _parse_chains(rest) if re.search(r"\d+:\d+.*[/_]", rest) else ([], None, 0)

    # Totals: prefer the sum segment, else add up chains
    total_c = total_db = None
    msum = re.search(r"(?<![/_\d])(\d{1,3}):(\d{1,2})(?:;(?:O(\d?)|(\d)O))?", s_sum.split("(")[0]) \
        or re.fullmatch(r"\(?(\d{1,3}):(\d{1,2})\)?", rest)
    if msum and not chains:
        total_c, total_db = int(msum.group(1)), int(msum.group(2))
    if chains or base:
        allc = chains + ([base] if base else [])
        total_c = sum(c for c, _ in allc)
        total_db = sum(d for _, d in allc)
    if total_c is None:
        return None

    if lipid_class.startswith(("Cer", "SM", "HexCer", "Hex2Cer")) and base is None:
        # Species-level only: assume d18:x base for acyl-length math
        base = (18, 1 if not lipid_class.endswith(("DS", "P")) else 0)

    if lipid_class in ("TG", "DG", "MG", "FA", "CAR", "LPC", "LPE", "LPS", "LPI", "LPG", "LPA") and not chains:
        if lipid_class in ("FA", "CAR", "LPC", "LPE", "LPS", "LPI", "LPG", "LPA", "MG"):
            chains = [(total_c, total_db)]

    if not lipid_class.startswith(("Cer", "SM", "HexCer", "Hex2Cer", "SPB")):
        mo = re.search(r"\d+:\d+;O(\d?)(?![\w(])", s_sum) or re.search(r"\d+:\d+;O(\d?)\s*$", s)
        if mo:
            oh = max(oh, int(mo.group(1) or 1))

    return Lipid(raw, lipid_class, total_c, total_db, oh, chains, base, adduct, labelled)


def _infer_cer_subclass(rest: str) -> str:
    """'18:1;O2/24:0' -> Cer_NS ; '18:0;O3/24:0' -> Cer_NP ; '18:1;O2/24:0;O' -> Cer_AS ;
    '18:1;O2/30:0-O-18:2' -> Cer_EOS ; species level '42:1;O2' -> Cer_NS, '42:1;O3' -> Cer_O3
    (ambiguous: NP, AS or ADS cannot be told apart without the chains)."""
    m = re.match(r"\(?\s*(\d+):(\d+);O(\d?)\s*/\s*(\d+):(\d+)(;OB?\d?)?(-O-)?", rest)
    if not m:
        m2 = re.match(r"\(?\s*(\d+):(\d+);O(\d?)\)?\s*$", rest)
        if not m2:
            return "Cer"
        db, o = int(m2.group(2)), int(m2.group(3) or 1)
        return "Cer_O3" if o >= 3 else ("Cer_NS" if db >= 1 else "Cer_NDS")
    db, o = int(m.group(2)), int(m.group(3) or 1)
    base = "P" if o == 3 else ("S" if db >= 1 else "DS")
    if m.group(7):
        return f"Cer_EO{base}"
    acyl = "N"
    if m.group(6):
        acyl = "B" if "B" in m.group(6) else "A"
    return f"Cer_{acyl}{base}"
