# Protein existence and strain Phenotype correlation calculator

Own implementation of bacterial GWAS (genome-wide association study).

Simple correlation calculator, given
1. A phenotype table
   - Test strain with target phenotype
   - Other strains with or without target phenotype
   - (Multiple phenotypes should be supported)
2. Reference proteome (focus strain)
3. Proteome database with all tested strain
   - (Shall this include the focus strain?)

Calculates: Correlation score of all proteins in the reference proteome, for
**how much they correlate with the given phenotype (combinations)**.

With simple Pearson correlation, the result is quite noisy. You can see it
contains information, but how to interpret the result still is a challenge.

Now that [evo2](https://github.com/ArcInstitute/evo2) has been released, I need
to re-evaluate the method used here.

## Step 1 Gather proteome

Make a proteome database, in `fasta` format, for `jackhmmer` to search.

## Step 2 Run jackhmmer

Given the threshold settings in `project_settings.py`, run `jackhmmer` for all
proteins in reference proteome, produce a table:

| protein ID | strain A | strain B | strain C|
| - | - | - | - |
| LC001 | 1 | 1 | 1 |
| LC002 | 1 | 0 | 1 |
| LC003 | 1 | 0 | 0 |

## Step 3 Calculate correlation

Currently Pearson correlation is used.
