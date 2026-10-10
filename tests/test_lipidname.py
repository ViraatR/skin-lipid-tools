import pytest

from lipidname import parse


@pytest.mark.parametrize("name, cls, c, db, acyl", [
    ("Cer[NS] 50:3; Cer[NS](d18:1/32:2); [M+Hac-H]-", "Cer_NS", 50, 3, 32),
    ("Cer[NP] 35:0; Cer[NP](t18:0/17:0); [M+Hac-H]-", "Cer_NP", 35, 0, 17),
    ("Cer 18:1;O2/24:0", "Cer_NS", 42, 1, 24),
    ("Cer 18:0;O3/17:0", "Cer_NP", 35, 0, 17),
    ("SM 34:1;O2", "SM", 34, 1, 16),
    ("GlcCer[NS] 42:1; GlcCer[NS](d18:1/24:0); [M+Hac-H]-", "HexCer_NS", 42, 1, 24),
])
def test_sphingolipids(name, cls, c, db, acyl):
    l = parse(name)
    assert (l.lipid_class, l.total_c, l.total_db, l.acyl_c) == (cls, c, db, acyl)


@pytest.mark.parametrize("name, cls, c, db, chains", [
    ("PC 35:3; PC(17:1/18:2); [M+H]+", "PC", 35, 3, [(17, 1), (18, 2)]),
    ("PC 16:0_18:1", "PC", 34, 1, [(16, 0), (18, 1)]),
    ("Plasmenyl-PE 40:6; PE(P-18:0/22:6); [M-H]-", "PE-P", 40, 6, [(18, 0), (22, 6)]),
    ("LysoPC 20:1; [M+Hac-H]-", "LPC", 20, 1, [(20, 1)]),
    ("20:1 Cholesteryl ester; [M+NH4]+", "CE", 20, 1, [(20, 1)]),
    ("TG 56:8; TG(16:1/18:1/22:6); [M+NH4]+", "TG", 56, 8, [(16, 1), (18, 1), (22, 6)]),
    ("FA(23:0)", "FA", 23, 0, [(23, 0)]),
])
def test_glycerolipids(name, cls, c, db, chains):
    l = parse(name)
    assert (l.lipid_class, l.total_c, l.total_db, l.chains) == (cls, c, db, chains)


def test_acylceramide_and_adduct():
    l = parse("Cer[EOS] 67:3; Cer[EOS](d17:1/32:0-O-18:2); [M+Hac-H]-")
    assert l.lipid_class == "Cer_EOS" and l.base == (17, 1) and l.total_c == 67
    assert l.adduct == "[M+Hac-H]-"


@pytest.mark.parametrize("name, cls, c, db, chains", [
    ("DG 15:0_18:1(d7)", "DG", 33, 1, [(15, 0), (18, 1)]),
    ("PC 15:0-18:1(d7)", "PC", 33, 1, [(15, 0), (18, 1)]),
    ("PE 15:0-18:1-d7", "PE", 33, 1, [(15, 0), (18, 1)]),
])
def test_labelled_standards(name, cls, c, db, chains):
    l = parse(name)
    assert (l.lipid_class, l.total_c, l.total_db, l.chains, l.labelled) == (cls, c, db, chains, True)


def test_msdial5_sum_and_species():
    l = parse("Cer 42:1;O2|Cer 18:1;O2/24:0")
    assert (l.lipid_class, l.total_c, l.base, l.acyl_c) == ("Cer_NS", 42, (18, 1), 24)
    assert l.raw == "Cer 42:1;O2|Cer 18:1;O2/24:0"


def test_ecn_and_key():
    l = parse("PC 34:1")
    assert l.ecn == 32 and l.key == "PC 34:1"


@pytest.mark.parametrize("junk", ["", None, "Hypoxanthine", "OHAFA(LA/C36:0)"])
def test_non_lipids(junk):
    assert parse(junk) is None
