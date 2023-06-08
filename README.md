# tnseq_tags

A program to generate randomised tags for QPCR and Illumina sequencing with the following structure:

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
