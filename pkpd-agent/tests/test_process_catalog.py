"""The addable-process catalog: single-compound mechanisms vs documented DDI.

Locks that the action space expanded, that every entry carries a validation
flag + provenance, and that DDI/interaction mechanisms are documented but NOT
offered as single-compound additions (they need a multi-compound setup).
"""

import glob
import json
import os
import unittest

from pkpd_agent.engines import osp_catalog as C

_LIB = os.path.join(os.path.dirname(__file__), "..", "..", "OSP-PBPK-Model-Library")


class TestProcessCatalog(unittest.TestCase):
    def test_single_compound_types_expanded(self):
        # was 8; expanded with kinetic variants + biliary clearance + Hill /
        # tubular secretion / hepatic-t1/2 forms from the PK-Sim v12 docs
        self.assertGreaterEqual(len(C.PROCESS_TYPES), 18)

    def test_no_first_order_active_transport(self):
        # PK-Sim's general membrane transport is MM/Hill only - the spurious
        # linear active-transport type was removed. (Renal tubular secretion is a
        # separate, legitimately first-order renal process and is NOT this.)
        self.assertNotIn("active_transport_first_order", C.PROCESS_TYPES)
        internals = {s.get("internal_name") for s in C.PROCESS_TYPES.values()}
        self.assertNotIn("ActiveTransportSpecific_FirstOrder", internals)

    def test_docs_derived_forms_present_and_flagged(self):
        # Hill, tubular secretion, in-vitro hepatic-t1/2 added from the docs;
        # their inferred InternalNames must be flagged unverified.
        for k in ("metabolization_hill", "active_transport_hill",
                  "tubular_secretion_first_order", "tubular_secretion_mm",
                  "hepatic_clearance_hepatocytes_thalf"):
            self.assertIn(k, C.PROCESS_TYPES)
            self.assertFalse(C.PROCESS_TYPES[k].get("internal_name_verified", False), k)
            self.assertIn("provenance", C.PROCESS_TYPES[k])

    def test_third_permeability_method(self):
        from pkpd_agent.engines.snapshot_edit import PERMEABILITY_METHODS
        self.assertEqual(len(PERMEABILITY_METHODS), 3)
        self.assertTrue(any("normalized" in m for m in PERMEABILITY_METHODS))

    def test_pka_and_compound_type_described(self):
        self.assertIn("ionised", C.describe_parameter("pKa")["description"])
        self.assertEqual(C.param_tier("pKa"), "constant")     # never fitted
        self.assertTrue(C.describe_parameter("Compound type").get("description"))

    def test_recombinant_cyp_types_verified(self):
        # rCYP450_MM is sildenafil's exact metabolism type - must be present AND
        # its InternalName verified from real snapshots (not inferred).
        by = {s["internal_name"]: s for s in C.PROCESS_TYPES.values()}
        self.assertIn("rCYP450_MM", by)
        self.assertIn("rCYP450_FirstOrder", by)
        self.assertTrue(by["rCYP450_MM"]["internal_name_verified"])

    @unittest.skipUnless(os.path.isdir(_LIB), "OSP library not present")
    def test_every_library_process_is_catalogued(self):
        # ground-truth coverage: every process InternalName that appears in any
        # real library snapshot is in our catalog (DDI mechanisms live in the
        # separate interaction catalog).
        known = {s.get("internal_name") for s in C.PROCESS_TYPES.values()}
        ddi = {s.get("internal_name") for s in C.INTERACTION_PROCESS_TYPES.values()}
        missing = set()
        for f in glob.glob(os.path.join(_LIB, "*", "json", "*.json")):
            try:
                d = json.load(open(f, encoding="utf-8"))
            except Exception:
                continue
            for comp in d.get("Compounds") or []:
                for p in comp.get("Processes") or []:
                    inm = p.get("InternalName")
                    if inm and inm not in known and inm not in ddi:
                        missing.add(inm)
        self.assertEqual(missing, set(), f"uncatalogued library processes: {missing}")

    def test_every_type_has_validation_and_provenance(self):
        for k, spec in C.PROCESS_TYPES.items():
            self.assertIn("validated", spec, k)
            self.assertTrue(spec.get("provenance") or k in
                            ("metabolization_first_order", "glomerular_filtration",
                             "metabolization_mm", "active_transport_mm",
                             "specific_binding", "liver_clearance",
                             "kidney_clearance"), k)

    def test_addable_carries_provenance(self):
        rows = C.addable_process_types([{"molecule": "CYP3A4", "type": "Enzyme"}])
        for r in rows:
            self.assertIn("provenance", r)
            self.assertIn("validated", r)

    def test_ddi_documented_but_not_addable(self):
        addable = {r["type"] for r in
                   C.addable_process_types([{"molecule": "CYP3A4", "type": "Enzyme"}])}
        ddi = {r["type"] for r in C.interaction_process_types()}
        self.assertTrue(ddi)                                  # DDI catalog exists
        self.assertFalse(addable & ddi)                      # never addable as single-compound
        self.assertIn("competitive_inhibition", ddi)
        self.assertIn("induction", ddi)

    def test_ddi_internal_names_flagged(self):
        by = {r["type"]: r for r in C.interaction_process_types()}
        # InternalNames confirmed against the library DDI snapshots
        self.assertTrue(by["competitive_inhibition"]["internal_name_verified"])
        self.assertTrue(by["induction"]["internal_name_verified"])
        self.assertTrue(by["mechanism_based_inhibition"]["internal_name_verified"])
        # not seen in the library set -> honestly still unverified
        self.assertFalse(by["uncompetitive_inhibition"]["internal_name_verified"])
        # validated (runs end-to-end through PKSim.CLI) is still False for all DDI:
        # it needs the multi-compound build, not yet exercised here
        self.assertFalse(any(r["validated"] for r in C.interaction_process_types()))

    def test_mbi_uses_k_kinact_half(self):
        by = {r["type"]: r for r in C.interaction_process_types()}
        pnames = [p["name"] for p in by["mechanism_based_inhibition"]["parameters"]]
        self.assertIn("K_kinact_half", pnames)   # ground truth, not "Ki"
        self.assertIn("kinact", pnames)


if __name__ == "__main__":
    unittest.main()
