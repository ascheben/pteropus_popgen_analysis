# Metadata

We generated a metadata table as well as two popmap files used as input for population genetic analyses and plotting scripts. The species name abbreviations SFF (*P. conspicillatus*), BFF (*P. alecto*), IBFF (Indonesian *P. alecto alecto*), LRFF (*P. scapulatus*), GHFF (*P. poliocephalus*)

* `pteropus_metadata.txt`: Here we aggregate all relevant sampling metadata with self-explanatory column headers. The collection date format is DD.MM.YYYY. Note that although the `Sample Identifier` column contains species abbreviations, the `Species` column encodes the final species determination. In six cases, there is a conflict between `Species` and `Sample Identifier` due to re-determinations of species status, and the `Species` column is the final determination but the `Sample Identifier` is retained for consistency with historical IDs.
* `pteropus.popmap`: Population mapping file for all species in the format `<Sample Identifier>\t<Species>\t<Roost>\t<Region>`
* `palecto_pconspicillatus.popmap`: Population mapping file for *P. conspicillatus* and *P. alecto* in the format `<Sample Identifier>\t<Population>`. The populations defined in this file are used for all population genetic analyses of these two species.
