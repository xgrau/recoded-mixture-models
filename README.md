# Recoded protein mixture models for IQ-TREE

This repository contains protein mixture models for some common amino-acid recoding schemes (Dayhoff6, Dayhoff9, SR4, SR6, KGB6) in `nexus` format for use in IQ-TREE.

An R script to convert an input alignment to various recoding schemes is also included.

## How to use

The model files are found in `recoded_models/`. It includes all combinations of the following:

* Recoding schemes: Dayhoff6, Dayhoff9, SR4, SR6, and KGB6. All models are encoded as `0-9` morphological characters, except for SR4 (encoded as `ACGT`).
* Protein mixture models: `C10` through `C60`.

Note that the models are unweighted.

To try them, indicate the model name with IQ-TREE `-m` flag, and model definition file with `-mdef`:

```bash
# recode your aminoacid alignment into various recoding schemes:
Rscript recoding_alignments.R test_data/test.fasta

# run iqtree pointing to your desired input fasta (recoded) and the corresponding model file:
iqtree -s test_data/test.recDayhoff6.fasta -m xmC10Dayhoff6 -mdef recoded_models/xmC10Dayhoff6.nex
```

## Included recoding schemes

Here:

```R
recoding_dictionary = list(
  "SR4" = list(
    "A" = "AGNPST",
    "C" = "CHWY",
    "G" = "DEKQR",
    "T" = "FILMV"
  ),
  "Dayhoff6" = list(
    "0" = "AGPST",
    "1" = "DENQ", 
    "2" = "HKR", 
    "3" = "ILMV",
    "4" = "FWY",
    "5" = "C"
  ),
  "Dayhoff9" = list(
    "0" = "DEHNQ",
    "1" = "ILMV",
    "2" = "FY",
    "3" = "AST",
    "4" = "KR",
    "5" = "G",
    "6" = "P",
    "7" = "C",
    "8" = "W"
  ),
  "SR6" = list(
    "0" = "APST",
    "1" = "DENG",
    "2" = "QKR",
    "3" = "MIVL",
    "4" = "WC",
    "5" = "FYH"
  ),
  "KGB6" = list(
    "0" = "AGPS",
    "1" = "DENQHKRT",
    "2" = "MIL",
    "3" = "W",
    "4" = "FY",
    "5" = "CV"
  )
)
```

## References

These models have been prepared following instructions in the IQ-TREE discussion group, [here](https://groups.google.com/g/iqtree/c/j884eSJiugY) and [here](https://groups.google.com/u/1/g/iqtree/c/GXNzbembxlA/m/NGze4jKWAgAJ).

References for the recoding models:

```bash
# references:
# Dayhoff6: M.O. Dayhoff, R.M. Schwartz, B.C. Orcutt, A model of evolutionary change in proteins
# Dayhoff9: Alexandra M Hernandez, Joseph F Ryan,  Six-State Amino Acid Recoding is not an Effective Strategy to Offset Compositional Heterogeneity and Saturation in Phylogenetic Analyses
# SR6: On reduced amino acid alphabets for phylogenetic inference, Mol. Biol. Evol., 24 (2007), pp. 2139-2150
# KGB6: A new criterion and method for amino acid classification J. Theor. Biol., 228 (2004), pp. 97-106
```
