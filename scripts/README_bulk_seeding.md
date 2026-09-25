# Compare a bulk MRCA with seeding cells

`compare_bulk_seeding.py` uses only the Python standard library. Inputs can be
TSV or TSV.GZ and need not be sorted. It stores the combined parent map in memory.

```bash
python3 /path/to/CometCNA/scripts/compare_bulk_seeding.py \
  --population population.tsv.gz \
  --dead-cells dead_cells.tsv.gz \
  --bulk-ids bulk_spatial_k100.tsv \
  --seeding-ids sampled_seeding_cells.tsv.gz \
  --cna cna.tsv.gz \
  --outdir bulk_seeding_comparison
```

Use the spatial sampler's selected-cell table for `--bulk-ids`, not its lineage
table. Bulk IDs must be extant in the population. To select fewer seeding cells,
supply a TSV with an `id` header and one seeding cell ID per row. The script
traces these IDs as given; it does not validate that they were actually seeded.
Dead-cell records are needed because seeding removes cells from the population.
Conflicting parents and incomplete ancestry cause errors.

Outputs (existing files with these names in the output directory are overwritten):

- `bulk_mrca.tsv`: bulk size and MRCA ID.
- `bulk_mrca_lineage.tsv`: bulk MRCA to founder, index 0 at the bulk MRCA.
- `seeding_lineages.tsv`: each seeding cell to founder, grouped by sample_id.
- `pair_mrcas.tsv`: common ancestor of the bulk MRCA and each seeding cell.
- `branches.tsv`: cells strictly after each pair's common ancestor, ordered
  from ancestor toward endpoint, with primary/seeding branch labels.
- `branch_cna_events.tsv`: one row per recorded CNA event per comparison,
  retaining multiple events in the same cell. A branch without events has no
  rows here; it is still represented in pair_mrcas.tsv.

The primary endpoint is the bulk MRCA; subsequent changes in individual bulk
cells are excluded. The seeding endpoint is the seeding cell; subsequent
metastatic growth is not reconstructed. Events assigned to the pair's common
ancestor are excluded because both branches inherit them. Branch membership is
based on ancestry, not merely event time or cell-ID ranges. Event rows describe
historical events, not net copy-number differences; later events can reverse
earlier changes. Events shared across different comparisons appear once per
comparison. A branch can be empty when its endpoint is the common ancestor.
