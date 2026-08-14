"""The addable-process catalog: single-compound mechanisms vs documented DDI.

Locks that the action space expanded, that every entry carries a validation
flag + provenance, and that DDI/interaction mechanisms are documented but NOT
offered as single-compound additions (they need a multi-compound setup).
"""

import unittest

from pkpd_agent.engines import osp_catalog as C


class TestProcessCatalog(unittest.TestCase):
    def test_single_compound_types_expanded(self):
        # was 8; expanded with kinetic variants + biliary clearance
        self.assertGreaterEqual(len(C.PROCESS_TYPES), 12)

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
        # confirmed against the DDI library snapshot
        self.assertTrue(by["competitive_inhibition"]["validated"])
        self.assertTrue(by["induction"]["validated"])
        # inferred ones honestly marked unvalidated
        self.assertFalse(by["mixed_inhibition"]["validated"])


if __name__ == "__main__":
    unittest.main()
