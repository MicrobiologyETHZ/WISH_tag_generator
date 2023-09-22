# WISH tag generator

![Design and principles of WISH tag construction](https://github.com/MicrobiologyETHZ/WISH_tag_generator/blob/main/Figure1.svg?raw=true)

A program to generate Wild-type Isogenic Standardized Hybrid (WISH) tags for QPCR and Illumina sequencing with the following structure:

```
| Fwd primer | Fwd spacer | Tag start | Tag primer | Rev spacer | Rev primer |
|<---24bp--->|<---24bp--->|<---16bp-->|<---24bp--->|<---8bp---->|<---24bp--->|
```

The design principle is to have a fixed sequence consisting of the Fwd primer, Fwd spacer, Rev spacer and Rev primer, with a variable tag consisting of the Tag start and Tag primer.

QPCR would be performed with the Fwd and Tag primers, whilst Illumina sequencing would be performed with the Fwd and Rev primers.

Individual primers are composed as follows:
* 19bp randomly chosen from [ACGT]
* 3bp randomly chosen from [AT]
* 2bp randomly chosen from [CG]

Many primers are generated then filtered to ensure that:
* GC content is between 45% and 55% inclusive
* the primer is not a palindrome, i.e.: identical in the reverse complement
* no base occurs more than 3 times in a row, i.e.: no AAAA
* no pair of bases occurs more than 3 times in a row, i.e.: no ACACACAC
* no triplet of bases occurs anywhere in the reverse complement (to prevent hairpins)
* there is no alignment with up to 2 mismatches in a specified sequence database
* the melting temperature according to Primer3 is between 59.5 and 62.5

The primers are then appended to the Illumina Nextera transposase adapters and again filtered so that no triplet of bases occurs anywhere in the reverse complement (to prevent hairpins).

Primers were then sorted from least to most penalty according to Primer3. The best 10 each of Fwd and Rev primers were paired in all possible combinations, scored as a pair with Primer3 and the pair with the least penalty chosen for the final constructs.

Spacers are composed of however many base pairs randomly chosen from [ACGT] then filtered to ensure that:
* GC content is between 45% and 55% inclusive
* the spacer is not a palindrome, i.e.: identical in the reverse complement

The final constructs are then once again checked for mispriming with two rounds of Primer3 (Fwd-Rev and Fwd-Tag), then sorted by their Primer3 penalties.

The method does not test final constructs against each other for possible hybridisation.

## Software requirements

The program requires the following software, with the version used in brackets (other versions are untested):

* Python (3.8.6)
* BWA (0.7.17)
* Primer3 (2.4.0)

## Sequence database

The sequence database should be a fasta file, `db/db.fasta` containing all sequences of interest that could be accidentally amplified in your system. It should be indexed with the `bwa index` command.

For the published tags we generated a database file from:

| Genome(s) | Reference |
| --------- | --------- |
| Escherichia coli 8178 | internal assembly |
| Escherichia coli Z1331 | internal assembly |
| Escherichia coli Nissle 1917 | NZ_CP007799.1 |
| Escherichia coli HS | NC_009800.1 |
| Salmonella enterica subsp. enterica serovar Typhimurium str. 14028S | CP001363.1 |
| Salmonella enterica subsp. enterica serovar Typhimurium SL1344 | FQ312003.1 |
| Salmonella enterica subsp. enterica serovar Typhimurium str. SL1344 plasmid pSLT_SL1344 | HE654724.1 |
| Salmonella enterica subsp. enterica serovar Typhimurium str. SL1344 plasmid pCol1B9_SL1344 | NC_017718.1 |
| Salmonella enterica subsp. enterica serovar Typhimurium str. SL1344 plasmid pRSF1010_SL1344 | NC_017719.1 |
| Eubacterium rectale ATCC 33656 | NC_012781.1 |
| Bacteroides thetaiotaomicron VPI-5482 | AE015928.1 |
| Whole genome sequencing of Prevotella isolates from human stool | PRJNA559898 |
| At-LSPHERE genome collection | PRJNA297956 |
| Pseudomonas syringae pv. tomato str. DC3000 | GCF_000007805.1 |
| Arabidopsis thaliana | GCF_000001735.4_TAIR10.1 |
| Mus musculus strain C57BL/6J GRCm38.p6 | GCF_000001635.26 |

