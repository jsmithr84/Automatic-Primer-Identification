# Automatic Primer Identification (Primer3)

This repository contains a Python script that **automatically designs PCR primers** to genotype a **SNP at a user-specified chromosome and position** in the *C. elegans* genome.

The script:
1. Extracts a window around the SNP from a genome FASTA.
2. Constrains primer placement so the SNP is **at least 100 bp away** from both ends of the PCR product.
3. Calls **primer3_core** to design a primer pair that yields a **600–800 bp amplicon**.
4. Prints the **best left and right primer sequences** to stdout.

---

## Primer Placement Requirements

The goal is to amplify a region around a SNP for genotyping:

- **Amplicon length:** 600–800 bp  
- **SNP safety margin:** SNP must be **≥ 100 bp** from both primer binding sites (i.e., not near either end of the PCR product)
- **Allowed primer search zones (relative to the extracted 1000 bp window):**
  - **Left primer:** positions **0–400**
  - **Right primer:** positions **600–1000**
- The SNP should sit in the middle region (roughly **400–600**) so primers do not overlap the SNP region.

The script enforces this by building a `SEQUENCE_PRIMER_PAIR_OK_REGION_LIST` constraint for primer3.

---

## Inputs

The script takes **three positional arguments**:

1. `fasta_file` — path to the *C. elegans* genome FASTA file  
2. `chromosome` — chromosome identifier (must match a FASTA header, e.g. `I`, `II`, `III`, `IV`, `V`, `X`, or `chrI` depending on your file)
3. `position` — 1-based genomic position of the SNP (integer)

---

## Output

The script prints two lines:

1. **left primer sequence** (lowercase)
2. **right primer sequence** (lowercase)

Example:
```text
atgctgac...
cctagaaa...
