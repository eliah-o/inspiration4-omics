# Inspiration4 direct RNA transcriptomics and methylation analysis

A collection of scripts, Jupyter notebooks, assets, and LaTeX files for analysis
and processing of direct RNA data, and for its visualization.



### 01-process-fast5/01_01-basecall.sh

Usage:

```
./01-process-fast5/01_01-basecall.sh \
    path/to/input-fast5-directory \
    path/to/output-temporary-blow5-directory \
    path/to/output-temporary-fastq-gz-directory \
    path/to/output-merged.blow5 \
    path/to/output-merged.fastq.gz
```

Runs, per-sample:

* guppy_basecaller: convert fast5 files to fastq.gz
* slow5tools: convert fast5 files to blow5
* f5c: index blow5/fastq.gz files (for eventalign later on)

The second and the third argument specify temporary directories that the script
removes after completing other steps.



### 02-qc/02_01-pycoQC.sh

Usage:

```
./02-qc/02_01-pycoQC.sh \
    path/to/input-guppy-summary.txt \
    path/to/output-pycoqc.json
```

Runs pycoQC per-sample.



### 02-qc/02_02-multiqc.sh

Usage:

```
./02-qc/02_02-multiqc.sh path/to/directory/with/all/input-pycoqc-jsons
```

Is run once on all samples; combines pycoQC JSONs into one report.



### 03-SARTools-pipeline/03_01-align.sh

Usage:

```
./03-SARTools-pipeline/03_01-align.sh \
    path/to/input-genomic-reference.fasta \
    path/to/input-genomic-reference.mmi \
    path/to/input-merged.fastq.gz \
    path/to/output-genomic.bam
```

Runs, per-sample: minimap2 and samtools to generate a genomic alignment file.
If the mmi version of the reference is not present, creates it automatically to
speed up the iterations.



### 03-SARTools-pipeline/03_02-featureCounts.sh

Usage:

```
./03-SARTools-pipeline/03_02-featureCounts.sh \
    path/to/input-index-file \
    path/to/input-annotation-file \
    path/to/input-genomic.bam \
    path/to/output-featureCounts-file
```

Runs, per-sample: featureCounts; generates a counts matrix.



### 03-SARTools-pipeline/03_03-SARTools.sh

Usage: ...



### 04-pipeline-transcriptome-de/04_01-pipeline-transcriptome-de.sh

Usage: ...



### 05-m6anet-and-methylKit/05_01-align.sh

Usage:

```
./05-m6anet-and-methylKit/05_01-align.sh \
    path/to/input-transcriptomic-reference.fasta \
    path/to/input-transcriptomic-reference.mmi \
    path/to/input-merged.fastq.gz \
    path/to/output-transcriptomic.bam
```

Runs, per-sample: minimap2 and samtools to generate a transcriptomic alignment
file.
If the mmi version of the reference is not present, creates it automatically to
speed up the iterations.



### 05-m6anet-and-methylKit/05_02-eventalign.sh

Usage:

```
./05-m6anet-and-methylKit/05_02-eventalign.sh \
    path/to/input-transcriptomic.bam \
    path/to/input-transcriptomic-reference.fasta \
    path/to/input-guppy-summary.txt \
    path/to/input-merged.fastq.gz \
    path/to/input-merged.blow5 \
    path/to/output-eventalign.tsv
```

Runs, per-sample: f5c eventalign.



### 05-m6anet-and-methylKit/05_03-m6anet.sh

Usage:

```
./05-m6anet-and-methylKit/05_03-m6anet.sh
    path/to/input-eventalign.tsv \
    path/to/output-m6anet-directory \
    path/to/output-m6anet.tsv.gz
```

Runs, per-sample: m6anet. Generates a gzip-compressed TSV file with per-site
methylation probabilities.



### 05-m6anet-and-methylKit/05_04-methylKit.sh

Usage:

```
./05-m6anet-and-methylKit/05_04-methylKit.sh \
    COMPARISON_NAME \
    path/to/output-methylKit.csv
```

Runs methylKit.  
NOTE: this combines multiple samples per comparison.

Comparison names that were used internally differ from the profile identifiers
in the paper; additionally, direction of comparisons (i.e. which is considered
a control, and which is considered treatment) were originally inverted for
all "postVrec" (I4-RP) comparisons:

Profile name | Internal COMPARISON_NAME | Note
-------------|--------------------------|---------
I4-FP1       | postVpre1                | normal
I4-FP2       | postVpre2                | normal
I4-FP3       | postVpre3                | normal
I4-LP1       | recVpre1                 | normal
I4-LP2       | recVpre2                 | normal
I4-LP4       | recVpre3                 | normal
I4-LP5       | recVpre4                 | normal
I4-LP6       | recVpre5                 | normal
I4-RP1       | postVrec1                | inverted
I4-RP2       | postVrec2                | inverted
I4-RP4       | postVrec3                | inverted
I4-RP5       | postVrec4                | inverted
I4-RP6       | postVrec5                | inverted

Similarly, interal names for the timepoints were different:

Timepoint | Internal name
----------|--------------
L-92      | 2021-06
L-44      | 2021-08
L-3       | 2021-09-1pre
R+1       | 2021-09-2post
R+45      | 2021-11
R+82      | 2021-12
R+194     | 2022-03

You may want to modify `05-m6anet-and-methylKit/r/samples.tsv` to point to  the
locations where m6anet outputs can be found.



### 06-plotting/06_01-DEG.ipynb

This is a Jupyter notebook that generates Figure 3b (DEG analysis).  
Note that throughout the notebook, internal comparison and timepoint names are
used, and are converted to the format seen in the paper just before generating
the plot.



### 06-plotting/06_02-m6a.ipynb

This is a Jupyter notebook that generates Figure 3c (m6a analysis).  
Note that throughout the notebook, internal comparison and timepoint names are
used, and are converted to the format seen in the paper just before generating
the plot.

Uses file `06-plotting/python/transcript_offsets.tsv.xz`, which can be generated
from the Gencode v41 GTF file using the script
`06-plotting/python/generate-transcript_offsets.py`.



### 06-plotting/latex/Figure-3.tex

This is a LaTeX file that combines all parts of Figure 3 into one PDF.

Uses file `06-plotting/latex/flowchart-yEd.png`, which was generated using
[yEd](https://www.yworks.com/products/yed) from file
`06-plotting/latex/flowchart-yEd.graphml`.
