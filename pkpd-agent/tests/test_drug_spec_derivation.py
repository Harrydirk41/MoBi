"""The drug-PD re-wire spec is now DERIVED from the paper model's own reactions instead of a
hand-written, drug-specific string (the last piece of scaffolding removed). The derivation runs in
MATLAB (i_deriveDrugSpec in sb_transplant_immune.m); this test is a faithful Python port of that
exact rule, checked against the REAL MTX rate laws captured from sb_drug_mechanism('MTX'). It proves
the derived spec reproduces the old hardcoded spec - so removing the scaffolding changed nothing but
where the knowledge comes from (the model, not the author)."""

import re
import unittest


def derive_drug_spec(reactions, drug="MTX"):
    """Port of i_deriveDrugSpec: reactions is a list of (product, rate, is_influx). Extract each
    '(1 +/- <param-with-drug>)' factor; the '*' wildcard is the most common secretion factor, and a
    cytokine gets an explicit key only if it carries a factor DIFFERENT from the wildcard."""
    rex = re.compile(r"\(\s*1\s*[-+]\s*\w*" + drug + r"\w*\s*\)")
    sec_by_cyt, all_sec, influx = {}, [], []
    for product, rate, is_influx in reactions:
        m = rex.search(rate)
        if not m:
            continue
        f = m.group(0).replace(" ", "")
        if is_influx:
            influx.append(f)
        elif product:
            sec_by_cyt.setdefault(product, []).append(f)
            all_sec.append(f)

    def mode(xs):
        return max(set(xs), key=xs.count) if xs else ""

    star = mode(all_sec)
    parts = ([f"*={star}"] if star else [])
    for cyt, fs in sec_by_cyt.items():
        non_star = [x for x in dict.fromkeys(fs) if x != star]
        if non_star:
            parts.append(f"{cyt}={non_star[0]}")
    return ";".join(parts), mode(influx)


# the REAL rate laws from sb_drug_mechanism('MTX') on the Vantage RA model (abbreviated to the factor)
ANTI = "*(1-Anti_CytSec_MTX)"
PRO = "*(1+Pro_CytSec_MTX)"
INFLUX = "*(1-Anti_CellInflux_MTX)"
PAPER_MTX_REACTIONS = [
    ("VEGF", "F_VEGF*kp_VEGF_byEndo*Endothelial" + ANTI, False),
    ("VEGF", "F_VEGF*kp_VEGF_byFLS*FLS" + ANTI, False),
    ("TNFa", "F_TNFa*kp_TNFa_byFLS*FLS" + ANTI, False),
    ("TNFa", "F_TNFa*kp_TNFa_byTh1*Th1" + ANTI, False),
    ("TNFa", "F_TNFa*kp_TNFa_byMacro*Macrophages" + ANTI, False),
    ("IL6", "F_IL6*Kp_IL6Sec_byFLS*FLS" + ANTI, False),
    ("IL6", "F_IL6*kp_IL6_byMacro*Macrophages" + ANTI, False),
    ("IL1b", "F_IL1b*kp_IL1b_byMacro*Macrophages" + ANTI, False),
    ("IFNg", "F_IFNg*kp_IFNg_byCTL*CTL" + ANTI, False),
    ("BAFF", "F_BAFF*kp_BAFF_byMacro*Macrophages" + ANTI, False),
    ("MCP1", "F_MCP1*kp_MCP1Sec_byMacro*Macrophages" + ANTI, False),
    ("MIP3", "F_MIP3*kp_MIP3_byFLS*FLS" + ANTI, False),
    ("RANTES", "F_RANTES*kp_RANTES_byFLS*FLS" + ANTI, False),
    ("IL17", "F_IL17*kp_IL17_byTh17*Th17" + ANTI, False),
    ("GMCSF", "F_GMCSF*kp_GMCSFSec_byFLS*FLS" + ANTI, False),
    # the two regulatory cytokines the drug BOOSTS from some sources (Pro) - the distinguishing factor
    ("IL10", "F_IL10*kp_IL10_byTreg*Treg" + PRO, False),
    ("IL10", "F_IL10*kp_IL10_byMacro*Macrophages" + ANTI, False),
    ("IL10", "F_IL10*kp_IL10_byBCells*BCells" + ANTI, False),
    ("TGFb", "F_TGFb*kp_TGFb_byMacro*Macrophages" + ANTI, False),
    ("TGFb", "F_TGFb*kp_TGFb_byTreg*Treg" + PRO, False),
    # influx reactions (classified by kIn_ in the real code; here by the is_influx flag)
    ("Th1", "kIn_Th1_Baseline*Hill_LeukoInflux_byCAM" + INFLUX, True),
    ("Treg", "kIn_Treg_Baseline*Hill_LeukoInflux_byCAM" + INFLUX, True),
    ("CTL", "kIn_CTL_Baseline*Hill_LeukoInflux_byCAM" + INFLUX, True),
    ("Macrophages", "kIn_Macrophage_Baseline*Hill_LeukoInflux_byCAM" + INFLUX, True),
]


class TestDeriveDrugSpec(unittest.TestCase):
    def test_derived_spec_matches_the_old_hardcoded_spec(self):
        sec, influx = derive_drug_spec(PAPER_MTX_REACTIONS, "MTX")
        # wildcard is the anti-secretion factor; only IL10/TGFb carry the distinguishing pro factor
        self.assertEqual(influx, "(1-Anti_CellInflux_MTX)")
        parts = sec.split(";")
        self.assertEqual(parts[0], "*=(1-Anti_CytSec_MTX)")
        self.assertIn("IL10=(1+Pro_CytSec_MTX)", parts)
        self.assertIn("TGFb=(1+Pro_CytSec_MTX)", parts)
        # no inflammatory cytokine gets an explicit key (they all match the wildcard)
        self.assertEqual(len([p for p in parts if p and not p.startswith("*=")]), 2)

    def test_no_drug_factors_yields_empty_spec(self):
        clean = [("IL6", "F_IL6*kp_IL6_byMacro*Macrophages", False),
                 ("Th1", "kIn_Th1_Baseline*Hill", True)]
        sec, influx = derive_drug_spec(clean, "MTX")
        self.assertEqual(sec, "")
        self.assertEqual(influx, "")


if __name__ == "__main__":
    unittest.main()
