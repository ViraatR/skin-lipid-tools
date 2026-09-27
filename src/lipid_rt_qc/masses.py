"""Species-level formulas and adduct m/z for common lipid classes.

Each class is a linear rule in total chain carbons T and double bonds D:
C = T + c0, H = 2T - 2D + h0, plus fixed N/O/P. The rules were checked against
LIPID MAPS reference formulas (e.g. PC 34:1 = C42H82NO8P, Cer 34:1;O2 =
C34H67NO3, SM 34:1;O2 = C39H79N2O6P, TG 52:2 = C55H102O6).
"""
from __future__ import annotations

from dataclasses import dataclass

from lipidname import Lipid

MONO = {"C": 12.0, "H": 1.00782503207, "N": 14.0030740048, "O": 15.99491461956, "P": 30.97376163}
PROTON = 1.007276466812

# class -> (c0, h0, N, O, P, db_offset). db_offset=1 for plasmalogens written P-x:y
# (the vinyl-ether double bond is not counted in y). Sphingolipid T/D include the base.
CLASS_RULES: dict[str, tuple[int, int, int, int, int, int]] = {
    "PC": (8, 16, 1, 8, 1, 0), "PE": (5, 10, 1, 8, 1, 0), "PS": (6, 10, 1, 10, 1, 0),
    "PG": (6, 11, 0, 10, 1, 0), "BMP": (6, 11, 0, 10, 1, 0), "PI": (9, 15, 0, 13, 1, 0),
    "PA": (3, 5, 0, 8, 1, 0),
    "PC-O": (8, 18, 1, 7, 1, 0), "PE-O": (5, 12, 1, 7, 1, 0),
    "PC-P": (8, 18, 1, 7, 1, 1), "PE-P": (5, 12, 1, 7, 1, 1),
    "LPC": (8, 18, 1, 7, 1, 0), "LPE": (5, 12, 1, 7, 1, 0), "LPS": (6, 12, 1, 9, 1, 0),
    "LPG": (6, 13, 0, 9, 1, 0), "LPI": (9, 17, 0, 12, 1, 0), "LPA": (3, 7, 0, 7, 1, 0),
    "TG": (3, 2, 0, 6, 0, 0), "DG": (3, 4, 0, 5, 0, 0), "MG": (3, 6, 0, 4, 0, 0),
    "CE": (27, 44, 0, 2, 0, 0), "FA": (0, 0, 0, 2, 0, 0), "CAR": (7, 13, 1, 4, 0, 0),
    "SM": (5, 13, 2, 6, 1, 0), "HexCer_NS": (6, 11, 1, 8, 0, 0), "HexCer_NDS": (6, 11, 1, 8, 0, 0),
    "Cer_NS": (0, 1, 1, 3, 0, 0), "Cer_NDS": (0, 1, 1, 3, 0, 0),
    "Cer_NP": (0, 1, 1, 4, 0, 0), "Cer_NH": (0, 1, 1, 4, 0, 0),
    "Cer_AS": (0, 1, 1, 4, 0, 0), "Cer_ADS": (0, 1, 1, 4, 0, 0), "Cer_BS": (0, 1, 1, 4, 0, 0),
    "Cer_BDS": (0, 1, 1, 4, 0, 0), "Cer_OS": (0, 1, 1, 4, 0, 0),
    "Cer_O3": (0, 1, 1, 4, 0, 0),  # species-level ;O3 (NP / AS / ADS unresolved)
    "Cer_AP": (0, 1, 1, 5, 0, 0), "Cer_AH": (0, 1, 1, 5, 0, 0),
    "Cer_EOS": (0, -1, 1, 5, 0, 0), "Cer_EODS": (0, -1, 1, 5, 0, 0),
}

ADDUCTS = {
    "[M+H]+": PROTON, "[M]+": PROTON, "[M+Na]+": 22.989218, "[M+NH4]+": 18.033823,
    "[M+H-H2O]+": PROTON - 18.010565, "[M+K]+": 38.963158,
    "[M-H]-": -PROTON, "[M+HCOO]-": 44.998201, "[M+FA-H]-": 44.998201,
    "[M+CH3COO]-": 59.013851, "[M+Hac-H]-": 59.013851, "[M+Cl]-": 34.969402,
}
# Plausible composition windows (total C, max DB) so isobar search ignores chemistry
# that does not occur in mammalian lipidomes (e.g. a ceramide with 11 double bonds).
def _window(cls: str) -> tuple[int, int, int]:
    if cls.startswith("Cer_EO"):
        return 50, 90, 6
    if cls.startswith(("Cer", "SM", "HexCer")):
        return 28, 72, 4
    if cls == "TG":
        return 36, 66, 14
    if cls in ("LPC", "LPE", "LPS", "LPG", "LPI", "LPA", "FA", "CAR", "CE", "MG"):
        return 10, 30, 6
    return 24, 48, 12  # diacyl glycerophospholipids and DG


POSITIVE = {a for a in ADDUCTS if a.endswith("+")}
NEGATIVE = {a for a in ADDUCTS if a.endswith("-")}


@dataclass(frozen=True)
class Formula:
    C: int
    H: int
    N: int
    O: int
    P: int

    @property
    def mass(self) -> float:
        return self.C * MONO["C"] + self.H * MONO["H"] + self.N * MONO["N"] + self.O * MONO["O"] + self.P * MONO["P"]

    def __str__(self) -> str:
        s = f"C{self.C}H{self.H}"
        for el in ("N", "O", "P"):
            n = getattr(self, el)
            s += "" if n == 0 else (el if n == 1 else f"{el}{n}")
        return s


def formula(lipid_class: str, total_c: int, total_db: int) -> Formula | None:
    r = CLASS_RULES.get(lipid_class)
    if r is None:
        return None
    c0, h0, n, o, p, dbo = r
    d = total_db + dbo
    return Formula(total_c + c0, 2 * total_c - 2 * d + h0, n, o, p)


def lipid_formula(l: Lipid) -> Formula | None:
    return formula(l.lipid_class, l.total_c, l.total_db)


def mz(lipid_class: str, total_c: int, total_db: int, adduct: str) -> float | None:
    f = formula(lipid_class, total_c, total_db)
    if f is None or adduct not in ADDUCTS:
        return None
    return f.mass + ADDUCTS[adduct]


def enumerate_isobars(target_mz: float, adducts: set[str], ppm: float = 10.0,
                      classes: list[str] | None = None):
    """All plausible (class, C, DB, adduct, ppm error) within `ppm` of target m/z."""
    out = []
    tol = target_mz * ppm / 1e6
    for cls in classes or CLASS_RULES:
        c0, h0, n, o, p, dbo = CLASS_RULES[cls]
        c_lo, c_hi, db_max = _window(cls)
        for add in adducts:
            base = ADDUCTS[add] + n * MONO["N"] + o * MONO["O"] + p * MONO["P"]
            for d in range(db_max + 1):
                # mass is linear in T: T*(C + 2H) + const, so solve for T directly
                k = MONO["C"] + 2 * MONO["H"]
                const = base + c0 * MONO["C"] + (h0 - 2 * (d + dbo)) * MONO["H"]
                t = round((target_mz - const) / k)
                if not (c_lo <= t <= c_hi):
                    continue
                m = t * k + const
                if abs(m - target_mz) <= tol:
                    out.append((cls, t, d, add, (m - target_mz) / target_mz * 1e6))
    return out
