# Manhattan Plotter
Functions to Annotate, Visualize and Interact with GWAS Outputs



#### Runner for Variant Annotator

```bash
#! /bin/bash
#SBATCH  -J VariantAnnotator
#SBATCH --mem=64G
#SBATCH -o VariantAnnotator.out
#############################################################################################
#
#
#       Name: VariantAnnotator_RUNNER.sh
#
#       Description: Runs VariantAnnotator.R
#
#       Date : 10/28/2021
#
#       By: wagu
#
#
##############################################################################################

# Example:
# Rscript -e 'source("/data/nrnb03/REF/SOFTWARE/RCode/VariantAnnotator.R"); VariantAnnotator("[SAIGE SNPLIST.TXT]")'
```



#### Runner for Manhattan Plotter

```bash
#! /bin/bash
#SBATCH  -J ManhattanPlotter
#SBATCH --mem=64G
#SBATCH -o ManhattanPlot.out
#############################################################################################
#
#
#       Name: ManhattanPlotter_RUNNER.sh
#
#       Description: Runs ManhattanPlotter.R
#
#       Date : 10/28/2021
#
#       By: wagu
#
#
##############################################################################################

# Example:
# Rscript -e 'source("/data/nrnb03/REF/SOFTWARE/RCode/ManhattanPlotter.R"); ManhattanPlotter("[SAIGE MERGED OUTPUT.TXT]")'
```



#### Runner for Interactive Manhattan Plotter

```bash
#! /bin/bash
#SBATCH  -J InteractiveManhattanPlotter
#SBATCH --mem=64G
#SBATCH -o InteractiveManhattanPlotter.out
#############################################################################################
#
#
#       Name: InteractiveManhattanPlotter_RUNNER.sh
#
#       Description: Runs InteractiveManhattanPlotter.R
#
#       Date : 10/28/2021
#
#       By: wagu
#
#
##############################################################################################

# Example:
# Rscript -e 'source("/data/nrnb03/REF/SOFTWARE/RCode/InteractiveManhattanPlotter.R"); InteractiveManhattanPlotter("[SAIGE MERGED OUTPUT.TXT]")'
```

