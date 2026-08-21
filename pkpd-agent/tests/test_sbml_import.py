"""SBML import: MathML->infix and full parse feeding QSPModel (no MATLAB, no LLM)."""

import os
import tempfile
import unittest
import xml.etree.ElementTree as ET

from pkpd_agent.engines.sbml_import import mathml_to_infix, sbml_to_network
from pkpd_agent.engines.qsp_model import QSPModel, infer_spec

_M = "{http://www.w3.org/1998/Math/MathML}"

_SBML = """<?xml version="1.0"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version1/core" level="3" version="1">
 <model id="m">
  <listOfSpecies>
   <species id="FLS"/><species id="TNFa"/><species id="IL6"/>
   <species id="Macrophages"/><species id="Endothelial"/>
   <species id="DAS28_CRP"/><species id="MTXDrug_Central"/>
  </listOfSpecies>
  <listOfParameters>
   <parameter id="kd_FLS_Baseline" value="0.1" units="dimensionless"/>
   <parameter id="Pro_IL6Sec_byMacro_effect" value="1" units="dimensionless"/>
  </listOfParameters>
  <listOfRules>
   <assignmentRule variable="Pro_FLSProlif_effect">
    <math xmlns="http://www.w3.org/1998/Math/MathML">
     <apply><min/><cn>10</cn>
      <apply><plus/>
       <apply><ci>MM</ci><ci>TNFa</ci><ci>a</ci></apply>
       <apply><ci>MM</ci><ci>IL6</ci><ci>a</ci></apply>
      </apply></apply></math>
   </assignmentRule>
   <assignmentRule variable="DAS28_CRP">
    <math xmlns="http://www.w3.org/1998/Math/MathML">
     <apply><plus/><ci>FLS</ci><ci>Macrophages</ci><ci>Endothelial</ci></apply></math>
   </assignmentRule>
  </listOfRules>
 </model>
</sbml>"""


def _write(txt):
    fd, p = tempfile.mkstemp(suffix=".xml")
    os.close(fd)
    with open(p, "w") as fh:
        fh.write(txt)
    return p


class TestMathML(unittest.TestCase):
    def _m(self, inner):
        return ET.fromstring(f'<math xmlns="http://www.w3.org/1998/Math/MathML">{inner}</math>')

    def test_binops(self):
        self.assertEqual(mathml_to_infix(self._m(
            "<apply><plus/><ci>a</ci><ci>b</ci></apply>")), "(a + b)")
        self.assertEqual(mathml_to_infix(self._m(
            "<apply><power/><ci>x</ci><cn>2</cn></apply>")), "(x ^ 2)")

    def test_function_application(self):
        self.assertEqual(mathml_to_infix(self._m(
            "<apply><ci>MM</ci><ci>TNFa</ci><ci>k</ci></apply>")), "MM(TNFa,k)")


class TestParse(unittest.TestCase):
    def test_shapes(self):
        p = _write(_SBML)
        try:
            d = sbml_to_network(p)
            self.assertEqual(len(d["species"]), 7)
            self.assertEqual(len(d["rules"]), 2)
            self.assertTrue(any("MM(TNFa" in r["rule"] for r in d["rules"]))
        finally:
            os.remove(p)

    def test_feeds_qspmodel(self):
        p = _write(_SBML)
        try:
            d = sbml_to_network(p)
            m = QSPModel(d, infer_spec(d, "sbml"))
            self.assertIn("FLS", m.nodes)
            self.assertNotIn("DAS28_CRP", m.nodes)
            self.assertNotIn("MTXDrug_Central", m.nodes)
            sset = {(e.source, e.sign, e.target) for e in m.edges}
            self.assertIn(("TNFa", 1, "FLS"), sset)
            self.assertIn(("Macrophages", 1, "IL6"), sset)   # from the _by param name
            self.assertEqual(set(m.readout_drivers),
                             {"FLS", "Macrophages", "Endothelial"})
        finally:
            os.remove(p)


if __name__ == "__main__":
    unittest.main()
