#!/usr/bin/env bash
# Runs every stage on the yeast example. Downloads ~0.7 GB of reads and two genomes.
set -euo pipefail
cd "$(dirname "$0")"
cfg=config.yaml

igem-warsaw-rnaseq qc metadata --config $cfg
igem-warsaw-rnaseq qc run --config $cfg --all
igem-warsaw-rnaseq qc collect --config $cfg
igem-warsaw-rnaseq ref prepare --config $cfg --all
igem-warsaw-rnaseq ref promoters --config $cfg
igem-warsaw-rnaseq quant route --config $cfg
igem-warsaw-rnaseq quant run --config $cfg --all
igem-warsaw-rnaseq quant collect --config $cfg
igem-warsaw-rnaseq ortho run --config $cfg
igem-warsaw-rnaseq ortho tree --config $cfg
igem-warsaw-rnaseq ortho conservation --config $cfg
igem-warsaw-rnaseq matrix build --config $cfg
igem-warsaw-rnaseq normalize run --config $cfg
igem-warsaw-rnaseq final build --config $cfg
