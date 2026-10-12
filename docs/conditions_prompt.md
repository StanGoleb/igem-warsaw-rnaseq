# Condition annotation with a language model

The pipeline exports, for every sample, all text that ENA and BioSample hold
about it (`igem-warsaw-rnaseq conditions export` writes `conditions_raw.tsv`
and an empty `conditions_template.tsv` into `qc.output_dir`). Assigning a
condition category and a normalised description to each sample is curation.
It can be done by hand or with a language model; the prompt below produces the
file the pipeline reads back (`normalize.conditions_file`).

## Procedure

1. Run `igem-warsaw-rnaseq conditions export --config my.yaml`.
2. Open `conditions_raw.tsv`. For large sets, split it into batches of at most
   100 rows per model call so that nothing is truncated.
3. Paste the prompt below, then the batch (the header line plus the rows), into
   the model. Copy the returned TSV into `conditions.tsv`.
4. Review `other` and every row where the model reported low confidence.
   Replicates of one experiment must receive identical `condition` and
   `description`.
5. Set `normalize.conditions_file: conditions.tsv` and run stage 6 and 7.

The category set below is the one used for the Trichoderma dataset that this
pipeline was developed on. Replace it with categories suited to your data; the
pipeline does not depend on specific category names.

## Prompt

```
You annotate RNA-seq samples. Each input row describes one sample (one ENA
experiment) with the text fields available in public archives: study title,
experiment title, sample title, sample description, sample alias, library name,
BioSample title and BioSample attributes (tag=value pairs separated by "; ").
Rows that share a bioproject come from the same study; read them together,
because the contrast between samples of one study is usually only visible when
they are compared.

For every row return exactly one output row with these tab-separated columns:

sample_id   copied unchanged from the input
condition   exactly one category from the list below
description one sentence, plain English, in the form
            "<genotype>, <treatment or growth condition>, <time point if any>",
            without replicate numbers, without the organism name, without
            instrument or library details; identical for replicates of the
            same experiment
confidence  high, medium or low

Categories:
control              untreated, wild-type reference, standard growth medium;
                     the reference against which other samples of the study
                     are compared
carbon source        growth on a defined or alternative carbon source
                     (cellulose, wheat bran, glucose limitation, no carbon)
biotic antagonism    confrontation with another fungus, nematode, bacterium
                     or other organism, or exposure to its extracts
plant association    growth with or on a plant, root colonisation, plant
                     exudates
light response       light or darkness treatment, light pulses, circadian
                     time courses
mechanical injury    wounding or mechanical damage of the mycelium
other                anything else (stress, chemicals, developmental stages,
                     temperature, pH); state the actual condition in the
                     description

Rules:
- Use "control" only for the reference samples of a study, not for every
  wild-type sample; a wild-type sample under a treatment belongs to the
  treatment category.
- A mutant strain under standard conditions is "control" with the genotype in
  the description; a mutant under a treatment takes the treatment category.
- When two categories apply, choose the one the study was designed to test and
  mention the other in the description.
- Do not invent information. If the texts do not allow a decision, use "other",
  write what is known in the description and set confidence to "low".
- Output only the TSV, with a header line, no commentary.

Input:
```

Paste the batch of `conditions_raw.tsv` rows after the word `Input:`.

## Output file

`conditions.tsv` is tab-separated with the columns `sample_id`, `condition`,
`description`; the `confidence` column and any other columns are carried into
`samples_metadata.tsv` unchanged. `condition` is used as the ComBat covariate
and for `rank_max_category` in `gene_summary.parquet`.
