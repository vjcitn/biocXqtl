# Provenance of PlinkArray/biocXqtl materials

## Gotcha

Before considering how nice all the following material is, note
```
> sum(duplicated(colnames(pl17)))
[1] 104523
```
so there are duplicated rsids in the columns of the ingested bed.
The associated bim has:
```
15 17 rs1410863556  0 113953   A      G
16 17 rs1410863556  0 113953   C      G
```

We have
```
> bim17[bim17$V2=="rs1007347099",]
       V1           V2 V3       V4                V5      V6
365941 17 rs1007347099  0 12697335               GTA       G
365942 17 rs1007347099  0 12697335             GTATA       G
365943 17 rs1007347099  0 12697335 GTATATATATATATATA       G
365944 17 rs1007347099  0 12697335                 G     GTA
365945 17 rs1007347099  0 12697335                 G   GTATA
365946 17 rs1007347099  0 12697335                 G GTATATA
```

This should be kept in mind when using the approach noted here.


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

## Split by chromosome, sometimes (e.g., chr6) to arbitrary subchromosomes

For chr17:
```
PLINK v2.0.0-a.6.32LM AVX2 Intel (10 Jan 2026)
Options in effect:
  --bfile all1kg
  --chr 17
  --make-bed
  --out all1kg_chr17

Hostname: vincelargegpu
Working directory: /home/exouser/1KG
Start time: Sun Jan 11 21:46:57 2026

Random number seed: 1768168017
60263 MiB RAM detected, ~56126 available; reserving 30131 MiB for main
workspace.
Using up to 16 threads (change this with --threads).
3202 samples (1603 females, 1598 males, 1 ambiguous; 2583 founders) loaded from
all1kg.fam.
Note: 2585 nonstandard chromosome codes present.
2073624 out of 74929081 variants loaded from all1kg.bim.
Note: No phenotype data present.
Writing all1kg_chr17.fam ... done.
Writing all1kg_chr17.bim ... done.
Writing all1kg_chr17.bed ... done.

End time: Sun Jan 11 21:47:06 2026
```

For a subchromosome, we cut up the bim file.  For chr6 it had 3 pieces

```
PLINK v2.0.0-a.6.32LM AVX2 Intel (10 Jan 2026)
Options in effect:
  --bfile all1kg_chr6
  --extract c6a.txt
  --make-bed
  --out all1kg_chr6a

Hostname: vincelargegpu
Working directory: /home/exouser/1KG
Start time: Sun Jan 18 19:52:52 2026

Random number seed: 1768765972
60263 MiB RAM detected, ~57712 available; reserving 30131 MiB for main
workspace.
Using up to 16 threads (change this with --threads).
3202 samples (1603 females, 1598 males, 1 ambiguous; 2583 founders) loaded from
all1kg_chr6.fam.
4315217 variants loaded from all1kg_chr6.bim.
Note: No phenotype data present.
--extract: 1500000 variants remaining.
1500000 variants remaining after main filters.
Writing all1kg_chr6a.fam ... done.
Writing all1kg_chr6a.bim ... done.
Writing all1kg_chr6a.bed ... done.

End time: Sun Jan 18 19:52:55 2026
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
