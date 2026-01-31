# Provenance of PlinkArray/biocXqtl materials

## Acquiring and splitting 3202 1KG samples

The link used to get the full hg38 genotyping for 3202 1KG samples is
provided by the plink2 project:

```
https://www.cog-genomics.org/plink/2.0/resources#phase3_1kg
```

Get the pgen resource:

```
wget https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst?dl=1
```

Transform to bed/bim/fam; the log is

```
PLINK v2.0.0-a.6.32LM AVX2 Intel (10 Jan 2026)
Options in effect:
  --make-bed
  --max-alleles 2
  --out all1kg
  --pfile all_hg38

Hostname: vincelargegpu
Working directory: /home/exouser/1KG
Start time: Sun Jan 11 20:47:12 2026

Random number seed: 1768164432
60263 MiB RAM detected, ~58094 available; reserving 30131 MiB for main
workspace.
Using up to 16 threads (change this with --threads).
3202 samples (1603 females, 1598 males, 1 ambiguous; 2583 founders) loaded from
all_hg38.psam.
Note: 2585 nonstandard chromosome codes present.
74929081 out of 75193455 variants loaded from all_hg38.pvar.
2 categorical phenotypes loaded.
74929081 variants remaining after main filters.
Writing all1kg.fam ... done.
Writing all1kg.bim ... done.
Writing all1kg.bed ... done.

End time: Sun Jan 11 20:50:36 2026
```

## Make these usable in R

Example DelayedArray ingestion of chr17:

```
> pl17
<3202 x 2073624> DelayedMatrix object of type "double":
          rs1482117838 rs2039636079 ... rs1186006889 rs1212764665
0_HG00096            0            0   .            0            0
0_HG00097            0            0   .            0            0
0_HG00099            0            0   .            0            0
0_HG00100            0            0   .            0            0
0_HG00101            0            0   .            0            0
      ...            .            .   .            .            .
0_NA21137            0            0   .            0            0
0_NA21141            0            0   .            0            0
0_NA21142            0            0   .            0            0
0_NA21143            0            0   .            0            0
0_NA21144            0            0   .            0            0
```
