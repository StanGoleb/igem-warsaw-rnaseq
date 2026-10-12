from pathlib import Path

import pytest

from igem_warsaw_rnaseq.config import ConfigError, load_config


def _write(tmp_path: Path, body: str) -> Path:
    cfg = tmp_path / "cfg.yaml"
    cfg.write_text(body)
    return cfg


def test_accessions_file_resolved_relative_to_config(tmp_path: Path) -> None:
    (tmp_path / "runs.txt").write_text("SRR1\n")
    cfg = load_config(
        _write(
            tmp_path, "qc:\n  output_dir: out\n  runs:\n    accessions_file: runs.txt\n"
        )
    )
    assert cfg.qc.runs.accessions_file == tmp_path / "runs.txt"
    assert cfg.qc.output_dir == tmp_path / "out"
    assert cfg.qc.allowed_library_strategies == {"RNA-Seq"}
    assert cfg.qc.gates.min_pct_q30 is None
    assert cfg.qc.rrna.enabled is False


def test_taxon_mode(tmp_path: Path) -> None:
    cfg = load_config(_write(tmp_path, "qc:\n  runs:\n    taxon_id: 63577\n"))
    assert cfg.qc.runs.taxon_id == 63577
    assert cfg.qc.runs.accessions_file is None


@pytest.mark.parametrize(
    ("runs", "message"),
    [
        ("{}", "set accessions_file"),
        ("{accessions_file: runs.txt, taxon_id: 1}", "alternatives"),
    ],
)
def test_run_source_rules(tmp_path: Path, runs: str, message: str) -> None:
    (tmp_path / "runs.txt").write_text("SRR1\n")
    with pytest.raises(ConfigError, match=message):
        load_config(_write(tmp_path, f"qc:\n  runs: {runs}\n"))


def test_rrna_enabled_requires_references(tmp_path: Path) -> None:
    with pytest.raises(ConfigError, match="references"):
        load_config(
            _write(tmp_path, "qc:\n  runs: {taxon_id: 1}\n  rrna: {enabled: true}\n")
        )


def test_gate_values(tmp_path: Path) -> None:
    cfg = load_config(
        _write(
            tmp_path,
            "qc:\n  runs: {taxon_id: 1}\n  gates: {min_pct_q30: 80, max_pct_adapter: 10.5}\n",
        )
    )
    assert cfg.qc.gates.min_pct_q30 == 80.0
    assert cfg.qc.gates.max_pct_adapter == 10.5


def test_bad_threads(tmp_path: Path) -> None:
    with pytest.raises(ConfigError, match="threads"):
        load_config(
            _write(tmp_path, "qc:\n  runs: {taxon_id: 1}\n  fastp: {threads: 0}\n")
        )
