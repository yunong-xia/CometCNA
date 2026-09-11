#!/usr/bin/env python3
"""Plot a sampled lineage tree with unary lineage nodes suppressed."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_lineages(path):
    parent = {}
    sampled = set()
    with open(path, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample_id", "id", "ancestor"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Lineage file must contain sample_id, id, and ancestor columns")
        for row in reader:
            sample_id = int(row["sample_id"])
            cell_id = int(row["id"])
            ancestor = int(row["ancestor"])
            sampled.add(sample_id)
            previous = parent.setdefault(cell_id, ancestor)
            if previous != ancestor:
                raise ValueError(f"Cell {cell_id} has conflicting ancestors")
    if not sampled:
        raise ValueError("Lineage file contains no samples")
    return parent, sampled


def reduce_tree(parent, sampled):
    children = defaultdict(set)
    for cell_id, ancestor in parent.items():
        if ancestor != 0:
            children[ancestor].add(cell_id)

    # First restrict the tree to ancestors of sampled leaves.
    active = set(sampled)
    for sample_id in sampled:
        current = sample_id
        seen = set()
        while current != 0 and current not in seen:
            seen.add(current)
            active.add(current)
            if current not in parent:
                raise ValueError(f"Missing ancestry for cell {current}")
            current = parent[current]

    # Keep sampled leaves, root, and nodes with at least two active children.
    keep = {1} | sampled
    keep.update(node for node in active if len(children.get(node, set()) & active) >= 2)

    # Walk each retained node's descendants until the next retained node.
    reduced_children = defaultdict(set)
    for start in sorted(keep):
        for child in children.get(start, set()) & active:
            current = child
            while current not in keep:
                next_nodes = children.get(current, ())
                if len(next_nodes) != 1:
                    break
                current = next(iter(next_nodes))
            if current in keep:
                reduced_children[start].add(current)

    if not (keep & sampled):
        raise ValueError("No sampled leaves found in the lineage tree")
    return keep, reduced_children


def layout_tree(keep, children, sampled):
    leaves = sorted((node for node in keep if node in sampled), key=int)
    y = {node: float(index) for index, node in enumerate(leaves)}
    x = {1: 0.0}

    def visit(node, depth):
        x[node] = depth
        child_ids = sorted(children.get(node, ()))
        for child in child_ids:
            visit(child, depth + 1.0)
        if node not in y:
            y[node] = sum(y[child] for child in child_ids) / len(child_ids)

    visit(1, 0.0)
    return x, y, leaves


def plot_tree(parent, sampled, output, dpi=160):
    keep, children = reduce_tree(parent, sampled)
    x, y, leaves = layout_tree(keep, children, sampled)
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    height = max(8.0, len(leaves) * 0.045)
    width = max(12.0, max(x.values(), default=1.0) * 1.25 + 2.0)
    fig, ax = plt.subplots(figsize=(width, height), constrained_layout=True)

    # Orthogonal branches make the retained branching structure easy to follow.
    for node, child_ids in children.items():
        if not child_ids:
            continue
        ordered = sorted(child_ids, key=lambda child: y[child])
        ax.plot([x[node], x[node]], [y[ordered[0]], y[ordered[-1]]],
                color="#555555", linewidth=0.35, zorder=1)
        for child in ordered:
            ax.plot([x[node], x[child]], [y[child], y[child]],
                    color="#555555", linewidth=0.35, zorder=1)

    branch_nodes = [node for node in keep if len(children.get(node, ())) >= 2]
    ax.scatter([x[node] for node in branch_nodes], [y[node] for node in branch_nodes],
               s=7, color="#c23b22", zorder=3, label="branching node")
    ax.scatter([x[node] for node in leaves], [y[node] for node in leaves],
               s=3, color="#1f5a99", zorder=3, label="sampled cell")
    for node in branch_nodes:
        ax.annotate(str(node), (x[node], y[node]), xytext=(3, 2),
                    textcoords="offset points", fontsize=5, color="#8b2415")
    ax.annotate("MRCA 1", (x[1], y[1]), xytext=(4, 4),
                textcoords="offset points", fontsize=7, fontweight="bold")
    ax.set_xlabel("Retained branching depth")
    ax.set_ylabel("Sampled cells (ordered by cell ID)")
    ax.set_title(f"Reduced lineage tree: {len(leaves):,} sampled cells, "
                 f"{len(branch_nodes):,} branching nodes")
    ax.set_yticks([])
    ax.legend(loc="upper right", markerscale=2, fontsize=7)
    ax.spines[["top", "right", "left"]].set_visible(False)
    fig.savefig(output, dpi=dpi)
    if output.suffix.lower() != ".svg":
        fig.savefig(output.with_suffix(".svg"))
    plt.close(fig)
    return len(keep), len(branch_nodes), len(leaves)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("lineages", help="Grouped lineage TSV")
    parser.add_argument("-o", "--output", required=True, help="PNG or SVG output path")
    args = parser.parse_args()
    try:
        parent, sampled = read_lineages(args.lineages)
        nodes, branches, leaves = plot_tree(parent, sampled, args.output)
    except (OSError, ValueError, RecursionError) as exc:
        parser.error(str(exc))
    print(f"nodes={nodes} branching_nodes={branches} sampled_leaves={leaves}")


if __name__ == "__main__":
    main()
