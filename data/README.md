# Yeast licensing dataset

The analysis expects the unchanged supplementary spreadsheet used by Dukaj and Rhind (2021):

> Dukaj L and Rhind N (2021) The capacity of origins to load MCM establishes replication timing patterns. *PLOS Genetics* 17, e1009467. Manuscript reference 14.

The file is not redistributed in this repository. Download it from the published study and place it here with the following name:

```text
data/OriDep.xlsx
```

The copy used for the manuscript analysis has these identifying properties:

- worksheet: `Hoja1`
- rows: 421 origins, excluding the header
- columns: 10
- SHA-256: `ae85379437d04e089cd93a011b68ef04b2621e2697f057c13812d23bf1b93c4b`

The analysis uses the genomic interval and chromosome fields plus the two untreated 0 micromolar auxin MCM ChIP-seq replicates. No values in the downloaded spreadsheet are curated or altered before import.
