import pandas as pd

from igem_warsaw_rnaseq import conditions

XML = """<?xml version="1.0"?>
<SAMPLE_SET>
  <SAMPLE accession="SAMN1"><TITLE>Microbe sample</TITLE>
    <SAMPLE_ATTRIBUTES>
      <SAMPLE_ATTRIBUTE><TAG>strain</TAG><VALUE>P1</VALUE></SAMPLE_ATTRIBUTE>
      <SAMPLE_ATTRIBUTE><TAG>treatment</TAG><VALUE>light</VALUE></SAMPLE_ATTRIBUTE>
    </SAMPLE_ATTRIBUTES>
  </SAMPLE>
  <SAMPLE accession="SAMN2"></SAMPLE>
</SAMPLE_SET>"""


def test_parse_sample_xml() -> None:
    parsed = conditions.parse_sample_xml(XML)
    assert parsed["SAMN1"] == ("Microbe sample", "strain=P1; treatment=light")
    assert parsed["SAMN2"] == ("", "")


def test_raw_table_groups_runs_per_experiment() -> None:
    meta = pd.DataFrame(
        {
            "run_accession": ["r1", "r2", "r3"],
            "in_ena": [True, True, False],
            "source": ["ena", "ena", "ena"],
            "experiment_accession": ["e1", "e1", ""],
            "sample_accession": ["SAMN1", "SAMN1", ""],
            "study_accession": ["PRJ1", "PRJ1", ""],
            "study_title": ["study", "study", ""],
            "experiment_title": ["exp", "exp", ""],
            "sample_title": ["s", "s", ""],
            "sample_description": ["d", "d", ""],
            "sample_alias": ["a", "a", ""],
            "library_name": ["l", "l", ""],
        }
    )
    table = conditions.raw_table(meta, conditions.parse_sample_xml(XML))
    assert list(table.columns) == list(conditions.RAW_COLUMNS)
    assert len(table) == 1
    assert table.loc[0, "runs"] == "r1;r2"
    assert table.loc[0, "biosample_attributes"] == "strain=P1; treatment=light"
