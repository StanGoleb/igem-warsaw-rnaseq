"""Pipeline configuration, loaded and validated from a YAML file."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml


class ConfigError(ValueError):
    """Raised when the configuration file is invalid."""


# ============================================================
# Stage 1: read QC
# ============================================================


@dataclass(frozen=True)
class RunSelection:
    accessions_file: Path | None
    taxon_id: int | None
    local_file: Path | None


@dataclass(frozen=True)
class FastpConfig:
    threads: int
    extra_args: tuple[str, ...]


@dataclass(frozen=True)
class RrnaConfig:
    enabled: bool
    references: tuple[Path, ...]
    threads: int


@dataclass(frozen=True)
class QcGates:
    min_pct_q30: float | None
    max_pct_adapter: float | None


@dataclass(frozen=True)
class QcConfig:
    output_dir: Path
    runs: RunSelection
    allowed_library_strategies: frozenset[str]
    excluded_platforms: frozenset[str]
    fastp: FastpConfig
    rrna: RrnaConfig
    gates: QcGates


# ============================================================
# Stage 2: references
# ============================================================


@dataclass(frozen=True)
class ReferenceSpec:
    name: str
    annotation_source: str
    ncbi_accession: str | None
    genome: Path | None
    annotation: Path | None
    tax_id: int | None
    anchor_genes: tuple[str, ...]


@dataclass(frozen=True)
class SalmonIndexConfig:
    kmer: int
    decoy: bool
    keep_duplicates: bool
    threads: int


@dataclass(frozen=True)
class BuscoConfig:
    enabled: bool
    lineage: str | None
    threads: int


@dataclass(frozen=True)
class RefConfig:
    output_dir: Path
    promoter_bp: int
    salmon: SalmonIndexConfig
    busco: BuscoConfig


# ============================================================
# Stage 4: orthology and phylogeny
# ============================================================


@dataclass(frozen=True)
class OrthoConfig:
    output_dir: Path
    threads: int
    gene_tree_method: str
    tree_model: str
    tree_bootstrap: int
    conservation_seed: int


# ============================================================
# Stage 3: quantification
# ============================================================


@dataclass(frozen=True)
class QuantConfig:
    output_dir: Path
    routing_overrides: Path | None
    threads: int
    bootstraps: int
    detected_min_reads: int
    min_pct_mapped: float | None
    crossmap_min_median_tpm: float | None
    crossmap_min_median_boot_cv: float | None


# ============================================================
# Stage 5: matrices
# ============================================================


@dataclass(frozen=True)
class MatrixConfig:
    output_dir: Path


# ============================================================
# Stage 6: normalisation
# ============================================================


@dataclass(frozen=True)
class AnchorConfig:
    n_auto: int
    min_mean_log2tpm: float
    candidate_top_n: int


@dataclass(frozen=True)
class ClassConfig:
    window: float
    min_cons_frac: float
    silent_min_median_pct: float
    inducible_min_amplitude: float
    top_pct: float


@dataclass(frozen=True)
class CombatConfig:
    enabled: bool
    min_batch_size: int


@dataclass(frozen=True)
class NormalizeConfig:
    output_dir: Path
    conditions_file: Path | None
    anchors: AnchorConfig
    classes: ClassConfig
    combat: CombatConfig


def _parse_normalize(section: dict[str, Any], base: Path) -> NormalizeConfig:
    anchors = section.get("anchors") or {}
    classes = section.get("classes") or {}
    combat = section.get("combat") or {}
    conditions = _path(section.get("conditions_file"), base)
    if conditions is not None and not conditions.is_file():
        raise ConfigError(f"normalize.conditions_file: not a file: {conditions}")

    def frac(key: str, default: float) -> float:
        value = _optional_float(classes.get(key, default), f"normalize.classes.{key}")
        if value is None or not 0.0 <= value <= 1.0:
            raise ConfigError(f"normalize.classes.{key}: expected a number in [0, 1]")
        return value

    return NormalizeConfig(
        output_dir=_path(section.get("output_dir", "data/normalized"), base),
        conditions_file=conditions,
        anchors=AnchorConfig(
            n_auto=_positive_int(anchors.get("n_auto", 3), "normalize.anchors.n_auto"),
            min_mean_log2tpm=float(anchors.get("min_mean_log2tpm", 4.0)),
            candidate_top_n=_positive_int(
                anchors.get("candidate_top_n", 500), "normalize.anchors.candidate_top_n"
            ),
        ),
        classes=ClassConfig(
            window=frac("window", 0.10),
            min_cons_frac=frac("min_cons_frac", 0.95),
            silent_min_median_pct=frac("silent_min_median_pct", 0.75),
            inducible_min_amplitude=frac("inducible_min_amplitude", 0.30),
            top_pct=frac("top_pct", 0.10),
        ),
        combat=CombatConfig(
            enabled=bool(combat.get("enabled", True)),
            min_batch_size=_positive_int(
                combat.get("min_batch_size", 2), "normalize.combat.min_batch_size"
            ),
        ),
    )


@dataclass(frozen=True)
class Config:
    qc: QcConfig | None
    ref: RefConfig | None
    quant: QuantConfig | None
    ortho: OrthoConfig | None
    matrix: MatrixConfig | None
    normalize: NormalizeConfig | None
    final: FinalConfig | None
    references: dict[str, ReferenceSpec]


def load_config(path: Path) -> Config:
    path = path.resolve()
    raw = yaml.safe_load(path.read_text())
    if not isinstance(raw, dict):
        raise ConfigError(f"{path}: expected a YAML mapping")
    base = path.parent
    return Config(
        qc=_parse_qc(raw["qc"], base) if raw.get("qc") is not None else None,
        ref=_parse_ref(raw["ref"], base) if raw.get("ref") is not None else None,
        quant=_parse_quant(raw["quant"], base)
        if raw.get("quant") is not None
        else None,
        ortho=_parse_ortho(raw["ortho"], base)
        if raw.get("ortho") is not None
        else None,
        matrix=MatrixConfig(
            output_dir=_path(
                (raw["matrix"] or {}).get("output_dir", "data/matrix"), base
            )
        )
        if raw.get("matrix") is not None
        else None,
        normalize=_parse_normalize(raw["normalize"], base)
        if raw.get("normalize") is not None
        else None,
        final=_parse_final(raw["final"], base)
        if raw.get("final") is not None
        else None,
        references=_parse_references(raw.get("references") or {}, base),
    )


# ============================================================
# Stage 7: final dataset
# ============================================================


@dataclass(frozen=True)
class SplitConfig:
    train: float
    val: float
    test: float
    seed: int


@dataclass(frozen=True)
class FinalConfig:
    output_dir: Path
    flank_bp: int
    split: SplitConfig


def _parse_final(section: dict[str, Any], base: Path) -> FinalConfig:
    split = section.get("split") or {}
    fractions = {
        k: _optional_float(split.get(k, d), f"final.split.{k}")
        for k, d in (("train", 0.8), ("val", 0.1), ("test", 0.1))
    }
    if (
        any(v is None or v < 0 for v in fractions.values())
        or abs(sum(fractions.values()) - 1.0) > 1e-6
    ):
        raise ConfigError("final.split: train, val and test must be >= 0 and sum to 1")
    return FinalConfig(
        output_dir=_path(section.get("output_dir", "data/final"), base),
        flank_bp=_positive_int(section.get("flank_bp", 1000), "final.flank_bp"),
        split=SplitConfig(
            train=fractions["train"],
            val=fractions["val"],
            test=fractions["test"],
            seed=_non_negative_int(split.get("seed", 1), "final.split.seed"),
        ),
    )


def _parse_ortho(section: dict[str, Any], base: Path) -> OrthoConfig:
    tree = section.get("species_tree") or {}
    method = str(section.get("gene_tree_method", "dendroblast"))
    if method not in ("dendroblast", "msa"):
        raise ConfigError(
            f"ortho.gene_tree_method: expected dendroblast or msa, got {method!r}"
        )
    return OrthoConfig(
        output_dir=_path(section.get("output_dir", "data/ortho"), base),
        threads=_positive_int(section.get("threads", 8), "ortho.threads"),
        gene_tree_method=method,
        tree_model=str(tree.get("model", "MFP")),
        tree_bootstrap=_non_negative_int(
            tree.get("bootstrap", 1000), "ortho.species_tree.bootstrap"
        ),
        conservation_seed=_non_negative_int(
            (section.get("conservation") or {}).get("seed", 1),
            "ortho.conservation.seed",
        ),
    )


def _parse_quant(section: dict[str, Any], base: Path) -> QuantConfig:
    salmon = section.get("salmon") or {}
    gates = section.get("gates") or {}
    crossmap = section.get("crossmap_flag") or {}
    overrides = _path(section.get("routing_overrides"), base)
    if overrides is not None and not overrides.is_file():
        raise ConfigError(f"quant.routing_overrides: not a file: {overrides}")
    return QuantConfig(
        output_dir=_path(section.get("output_dir", "data/quant"), base),
        routing_overrides=overrides,
        threads=_positive_int(salmon.get("threads", 4), "quant.salmon.threads"),
        bootstraps=_non_negative_int(
            salmon.get("bootstraps", 30), "quant.salmon.bootstraps"
        ),
        detected_min_reads=_positive_int(
            section.get("detected_min_reads", 10), "quant.detected_min_reads"
        ),
        min_pct_mapped=_optional_float(
            gates.get("min_pct_mapped", 5.0), "quant.gates.min_pct_mapped"
        ),
        crossmap_min_median_tpm=_optional_float(
            crossmap.get("min_median_tpm"), "quant.crossmap_flag.min_median_tpm"
        ),
        crossmap_min_median_boot_cv=_optional_float(
            crossmap.get("min_median_boot_cv"),
            "quant.crossmap_flag.min_median_boot_cv",
        ),
    )


def _parse_qc(section: dict[str, Any], base: Path) -> QcConfig:
    runs = section.get("runs") or {}
    accessions_file = _path(runs.get("accessions_file"), base)
    taxon_id = runs.get("taxon_id")
    local_file = _path(runs.get("local_file"), base)
    if accessions_file is not None and taxon_id is not None:
        raise ConfigError("qc.runs: accessions_file and taxon_id are alternatives")
    if accessions_file is None and taxon_id is None and local_file is None:
        raise ConfigError("qc.runs: set accessions_file, taxon_id and/or local_file")
    for label, p in (("accessions_file", accessions_file), ("local_file", local_file)):
        if p is not None and not p.is_file():
            raise ConfigError(f"qc.runs.{label}: not a file: {p}")

    fastp = section.get("fastp") or {}
    rrna = section.get("rrna") or {}
    gates = section.get("gates") or {}

    references = tuple(_path(p, base) for p in rrna.get("references") or [])
    enabled = bool(rrna.get("enabled", False))
    if enabled and not references:
        raise ConfigError("qc.rrna: enabled but no references given")
    for ref in references:
        if not ref.is_file():
            raise ConfigError(f"qc.rrna.references: not a file: {ref}")

    return QcConfig(
        output_dir=_path(section.get("output_dir", "data/qc"), base),
        runs=RunSelection(
            accessions_file=accessions_file, taxon_id=taxon_id, local_file=local_file
        ),
        allowed_library_strategies=frozenset(
            section.get("allowed_library_strategies") or ["RNA-Seq"]
        ),
        excluded_platforms=frozenset(
            section.get("excluded_platforms", ["ABI_SOLID"]) or []
        ),
        fastp=FastpConfig(
            threads=_positive_int(fastp.get("threads", 4), "qc.fastp.threads"),
            extra_args=tuple(str(a) for a in fastp.get("extra_args") or []),
        ),
        rrna=RrnaConfig(
            enabled=enabled,
            references=references,
            threads=_positive_int(rrna.get("threads", 4), "qc.rrna.threads"),
        ),
        gates=QcGates(
            min_pct_q30=_optional_float(
                gates.get("min_pct_q30"), "qc.gates.min_pct_q30"
            ),
            max_pct_adapter=_optional_float(
                gates.get("max_pct_adapter"), "qc.gates.max_pct_adapter"
            ),
        ),
    )


def _parse_ref(section: dict[str, Any], base: Path) -> RefConfig:
    salmon = section.get("salmon") or {}
    busco = section.get("busco") or {}
    enabled = bool(busco.get("enabled", False))
    lineage = busco.get("lineage")
    if enabled and not lineage:
        raise ConfigError("ref.busco: enabled but no lineage given")
    return RefConfig(
        output_dir=_path(section.get("output_dir", "data/ref"), base),
        promoter_bp=_positive_int(section.get("promoter_bp", 800), "ref.promoter_bp"),
        salmon=SalmonIndexConfig(
            kmer=_positive_int(salmon.get("kmer", 31), "ref.salmon.kmer"),
            decoy=bool(salmon.get("decoy", True)),
            keep_duplicates=bool(salmon.get("keep_duplicates", True)),
            threads=_positive_int(salmon.get("threads", 8), "ref.salmon.threads"),
        ),
        busco=BuscoConfig(
            enabled=enabled,
            lineage=str(lineage) if lineage else None,
            threads=_positive_int(busco.get("threads", 8), "ref.busco.threads"),
        ),
    )


def _parse_references(section: dict[str, Any], base: Path) -> dict[str, ReferenceSpec]:
    specs: dict[str, ReferenceSpec] = {}
    for name, entry in section.items():
        entry = entry or {}
        accession = entry.get("ncbi_accession")
        genome = _path(entry.get("genome"), base)
        annotation = _path(entry.get("annotation"), base)
        has_files = genome is not None or annotation is not None
        if (accession is None) == (not has_files):
            raise ConfigError(
                f"references.{name}: set either ncbi_accession or genome+annotation"
            )
        if has_files:
            for label, p in (("genome", genome), ("annotation", annotation)):
                if p is None:
                    raise ConfigError(f"references.{name}: {label} is required")
                if not p.is_file():
                    raise ConfigError(f"references.{name}.{label}: not a file: {p}")
        tax_id = entry.get("tax_id")
        if tax_id is not None:
            tax_id = _positive_int(tax_id, f"references.{name}.tax_id")
        specs[str(name)] = ReferenceSpec(
            name=str(name),
            annotation_source=str(entry.get("annotation_source", "")),
            ncbi_accession=str(accession) if accession else None,
            genome=genome,
            annotation=annotation,
            tax_id=tax_id,
            anchor_genes=tuple(str(g) for g in entry.get("anchor_genes") or []),
        )
    seen: dict[int, str] = {}
    for spec in specs.values():
        if spec.tax_id is None:
            continue
        if spec.tax_id in seen:
            raise ConfigError(
                f"references: {seen[spec.tax_id]} and {spec.name} share tax_id "
                f"{spec.tax_id}; routing would be ambiguous"
            )
        seen[spec.tax_id] = spec.name
    return specs


def _path(value: Any, base: Path) -> Path | None:
    if value is None:
        return None
    p = Path(str(value)).expanduser()
    return p if p.is_absolute() else (base / p).resolve()


def _positive_int(value: Any, name: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < 1:
        raise ConfigError(f"{name}: expected a positive integer, got {value!r}")
    return value


def _non_negative_int(value: Any, name: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise ConfigError(f"{name}: expected a non-negative integer, got {value!r}")
    return value


def _optional_float(value: Any, name: str) -> float | None:
    if value is None:
        return None
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigError(f"{name}: expected a number or null, got {value!r}")
    return float(value)
