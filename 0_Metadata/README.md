# Metadata

We generated a sample metadata table and two popmap files. They are used as input for the population genetic analyses and plotting scripts in all other modules.

Species abbreviations: SFF (*P. conspicillatus*), BFF (*P. alecto*), IBFF (Indonesian *P. alecto alecto*), LRFF (*P. scapulatus*), GHFF (*P. poliocephalus*).

* `pteropus_metadata.txt`: sampling metadata for all 242 samples that passed quality control, with self-explanatory column headers (Table S1). The collection date format is DD.MM.YYYY.
  * **Species:** the `Species` column is the final species determination. Although `Sample Identifier` contains species abbreviations, six samples were re-determined, so the two columns conflict (e.g. `1510019_LRFF` is *P. poliocephalus*). `Sample Identifier` is kept for consistency with historical IDs and matches the sample names in all VCFs.
  * **Source:** `CSIRO`, `Rocky bat carer` and `Broome bat carer` are extant samples. `WA museum` samples are historical (1989–2003).
  * **Population:** `BFF_SFF_Population_ID` gives the regional population of *P. alecto* and *P. conspicillatus* samples (`Outgroup` = *P. alecto alecto*).
* `pteropus.popmap`: population map for all species, in the format `<Sample Identifier>\t<Species>\t<Roost>\t<Region>`.
* `palecto_pconspicillatus.popmap`: population map for *P. alecto* and *P. conspicillatus* (150 samples), in the format `<Sample Identifier>\t<Population>`. These populations are used for all population genetic analyses of the two species.

### Population codes

Regional populations (`palecto_pconspicillatus.popmap`) and their manuscript names:

| Code | Manuscript region | Species |
|---|---|---|
| `INDO` | Indonesia | *P. alecto* |
| `INDO_outgroup` | Indonesia | *P. alecto alecto* |
| `N_AUS` | North West Australia | *P. alecto* |
| `NQ` | North East Australia | *P. alecto* |
| `E_COAST` | East Australia | *P. alecto* |
| `Wet_tropics` | Wet Tropics Australia | *P. conspicillatus* |
| `PNG` | New Guinea | *P. conspicillatus* |

The TreeMix and fastsimcoal2 analyses use genetic populations defined from the PCA (Figure S7). These codes are used in `4_Treemix` and `5_Fastsimcoal`:

| Code | Manuscript code | Samples |
|---|---|---|
| `BFFINDOplusWAmuseum` (TreeMix: `INDO`) | INDO+3NWA | `INDO` plus the three North West Australian museum samples M53754, M54659 and M54660 |
| `BFFINDO` | INDO | `INDO` only |
| `BFFNAUS` (TreeMix: `N_AUS`) | NWA+BURK | remaining `N_AUS` samples plus all Burketown samples |
| `BFFEAST` | EA+NEA | remaining `NQ` and all `E_COAST` samples |
| `BFFAUS` | NWA+EA+NEA | all Australian *P. alecto* except M53754, M54659 and M54660 |
| `SFF` | WTA+NG | all *P. conspicillatus* (`Wet_tropics` and `PNG`) |
| `IBFF` | – | *P. alecto alecto* (`INDO_outgroup`) |
