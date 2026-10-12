from pathlib import Path

import numpy as np
import pandas as pd

from igem_warsaw_rnaseq.normalize import plots


def _summary(n: int) -> pd.DataFrame:
    rng = np.random.default_rng(0)
    p05 = rng.uniform(0, 0.5, n)
    return pd.DataFrame(
        {
            "n_samples": 20,
            "median_pct": p05 + 0.2,
            "p05_pct": p05,
            "p95_pct": p05 + rng.uniform(0, 0.5, n),
            "cons_frac": rng.uniform(0, 1, n),
            "induction_amplitude": rng.uniform(0, 0.5, n),
            "class": rng.choice(["constitutive", "inducible", "variable", "silent"], n),
        },
        index=pd.Index([f"g{i}" for i in range(n)], name="gene_id"),
    )


def test_gene_map_writes_png(tmp_path: Path) -> None:
    out = tmp_path / "map.png"
    plots.gene_map(_summary(50), out, "title", "gene", 20, highlight="g3")
    assert out.stat().st_size > 1000


def test_gene_figure_with_and_without_orthogroup(tmp_path: Path) -> None:
    summary = _summary(30)
    strength = pd.Series(np.random.default_rng(1).uniform(0, 1, 20))
    row = summary.loc["g1"]
    out = tmp_path / "gene.png"
    plots.gene_figure(
        "g1", "ref", "kinase", row, strength, summary, None, None, None, out
    )
    assert out.stat().st_size > 1000
    plots.gene_figure(
        "g1", "ref", "", row, strength, summary, "OG1", row, strength, out
    )
    assert out.stat().st_size > 1000


def test_gene_descriptions_parses_gffread_table(tmp_path: Path, monkeypatch) -> None:
    table = "rna-1\tkinase\tKIN1\nrna-2\t.\t.\nrna-3\tphosphatase\t.\n"
    monkeypatch.setattr(
        plots.subprocess,
        "run",
        lambda *a, **k: type("R", (), {"stdout": table})(),
    )
    desc = plots.gene_descriptions(
        tmp_path / "x.gff", {"rna-1": "g1", "rna-2": "g2", "rna-3": "g3"}
    )
    assert desc == {"g1": "kinase (KIN1)", "g2": "", "g3": "phosphatase"}
