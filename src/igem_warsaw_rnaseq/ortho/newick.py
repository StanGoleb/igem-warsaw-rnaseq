"""Minimal Newick reader and patristic (branch-length) distances between leaves."""

from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd


@dataclass
class Node:
    name: str = ""
    length: float = 0.0
    children: list[Node] = field(default_factory=list)


def parse(text: str) -> Node:
    tokens = _tokenize(text.strip().rstrip(";"))
    pos = 0

    def node() -> Node:
        nonlocal pos
        n = Node()
        if tokens[pos] == "(":
            pos += 1
            n.children.append(node())
            while tokens[pos] == ",":
                pos += 1
                n.children.append(node())
            if tokens[pos] != ")":
                raise ValueError("newick: expected ')'")
            pos += 1
        if pos < len(tokens) and tokens[pos] not in "(),:":
            n.name = tokens[pos]
            pos += 1
        if pos < len(tokens) and tokens[pos] == ":":
            pos += 1
            n.length = float(tokens[pos])
            pos += 1
        return n

    root = node()
    if pos != len(tokens):
        raise ValueError("newick: trailing tokens")
    return root


def _tokenize(text: str) -> list[str]:
    tokens: list[str] = []
    buf = ""
    for ch in text:
        if ch in "(),:":
            if buf:
                tokens.append(buf)
                buf = ""
            tokens.append(ch)
        elif not ch.isspace():
            buf += ch
    if buf:
        tokens.append(buf)
    return tokens


def leaves(root: Node) -> list[str]:
    return [n.name for n in _walk(root) if not n.children]


def patristic_distances(root: Node) -> pd.DataFrame:
    """Leaf x leaf matrix of summed branch lengths along the connecting path."""
    paths: dict[str, list[Node]] = {}

    def descend(n: Node, trail: list[Node]) -> None:
        trail = [*trail, n]
        if not n.children:
            paths[n.name] = trail
        for c in n.children:
            descend(c, trail)

    descend(root, [])
    names = list(paths)
    matrix = pd.DataFrame(0.0, index=names, columns=names)
    for a in names:
        for b in names:
            if a >= b:
                continue
            shared = 0
            for x, y in zip(paths[a], paths[b], strict=False):
                if x is not y:
                    break
                shared += 1
            d = sum(n.length for n in paths[a][shared:]) + sum(
                n.length for n in paths[b][shared:]
            )
            matrix.loc[a, b] = matrix.loc[b, a] = d
    return matrix


def _walk(n: Node):
    yield n
    for c in n.children:
        yield from _walk(c)
