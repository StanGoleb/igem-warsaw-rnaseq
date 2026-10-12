"""One representative gene per (orthogroup, reference)."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import pandas as pd

from igem_warsaw_rnaseq.ortho import newick

SINGLE_COPY = "single_copy"
GENE_TREE = "gene_tree"
MAX_MEAN_TPM = "max_mean_tpm"
LONGEST_PROTEIN = "longest_protein"

COLUMNS: tuple[str, ...] = (
    "orthogroup",
    "reference",
    "gene_id",
    "criterion",
    "n_paralogs_collapsed",
)


def leaf_owner(leaf: str, references: Sequence[str]) -> tuple[str, str] | None:
    """Split an OrthoFinder leaf ``<reference>_<gene>`` on the longest matching name."""
    for ref in sorted(references, key=len, reverse=True):
        if leaf.startswith(f"{ref}_"):
            return ref, leaf[len(ref) + 1 :]
    return None


def choose_by_tree(
    tree_newick: str,
    reference: str,
    candidates: Sequence[str],
    references: Sequence[str],
) -> str | None:
    """Candidate with the smallest mean patristic distance to other references' genes."""
    tree = newick.parse(tree_newick)
    owners = {leaf: leaf_owner(leaf, references) for leaf in newick.leaves(tree)}
    others = [leaf for leaf, o in owners.items() if o and o[0] != reference]
    wanted = {f"{reference}_{g}": g for g in candidates}
    present = [leaf for leaf in wanted if leaf in owners]
    if not others or len(present) < 2:
        return None
    distances = newick.patristic_distances(tree)
    scores = {wanted[leaf]: distances.loc[leaf, others].mean() for leaf in present}
    return min(scores, key=scores.get)


def choose(
    reference: str,
    genes: Sequence[str],
    tree_newick: str | None,
    mean_tpm: Mapping[str, float] | None,
    protein_length: Mapping[str, int],
    references: Sequence[str],
) -> tuple[str, str]:
    """(gene, criterion) following single_copy > gene_tree > max_mean_tpm > longest."""
    if len(genes) == 1:
        return genes[0], SINGLE_COPY
    if tree_newick:
        picked = choose_by_tree(tree_newick, reference, genes, references)
        if picked is not None:
            return picked, GENE_TREE
    if mean_tpm is not None and any(g in mean_tpm for g in genes):
        return max(genes, key=lambda g: mean_tpm.get(g, 0.0)), MAX_MEAN_TPM
    return max(genes, key=lambda g: protein_length.get(g, 0)), LONGEST_PROTEIN


def table(rows: list[dict[str, object]]) -> pd.DataFrame:
    return pd.DataFrame(rows, columns=list(COLUMNS))
