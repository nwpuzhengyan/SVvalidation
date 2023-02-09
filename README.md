# SVvaliation
## Introduction

Although there are many methods have been developed to detect the Structural variations (SVs) from the genomic sequences, there are few methods used to validate their results. We observed that the results of current several SV callers all have a large number of false positive SVs, and current validation methods are not accurate enough. Hence, we propose a new method-SVvalidation, a more accurate and highly efficient method for validating SVs by using long-read sequencing data. Compared with existing methods, SVvalidation not only has better performance for the SVs in repeat regions and validates whether a SV is homozygous or heterozygous. In addition, SVvalidation has the highest precision and F1-score (improved by over 10\%) in all datasets.

---
## Installation
```
git clone https://github.com/nwpuzhengyan/SVvaliation.git
```
---
## Dependence
    1. python3
	2. pysam
	3. cigar
	4. pyfaidx
	5. argparse
---
## Running
The input files are sorted bam file and a SV bed file.
```
cd dist
SVvalidation <input SV bed file>	<input sorted bam>
```
---
## Input format
The input format is as follows. The first column is chromosome name, the second and third columns are SV start and end position. The fourth column is the SV length and fifth column are SV type.
```
chr1	1925143	1925143	147	INS
chr1	2768484	2768692	208	DEL
chr1	2866089	2866236	147	INV
```
---
## Output format
The output format is as follows. The first three columns are same with input format. The fourth column is the validated SV length and fifth column are validated SV type. homo_INS represents the SV is a homozygous INS and heter_DEL represents the SV is a heterozygous DEL. false_SV represents there is no INV in this location.
```
chr1	1925143	1925143	147	homo_INS
chr1	2768484	2768692	208	heter_DEL
chr1	2866089	2866236	147	false_INV
```
---
## Contact
For advising, bug reporting and requiring help, please contact yan.zheng@nwpu-bioinformatics.com.

