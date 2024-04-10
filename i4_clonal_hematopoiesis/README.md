This folder contains:
- AstronautAnalysis.R
- data/
- utils/

**AstronautAnalysis.R** - contains code required to produce the SpiderPlot figure showing the change in VAF for the Astronaut samples and a similarly age matched cohort without any treatment exposure. 
This script was developed using R4.3.2 and the following libraries:
- tidyr_1.3.1
- data.table_1.15.2
- vcfR_1.15.0

**data/** - contains the required raw data required by the AstronautAnalysis.R code to produce the figure. All data produced within this directory was either manually created or was created through the [ArCH pipeline](https://doi.org/10.1093/bioinformatics/btae121)

**utils/** - contains additional R functions to process and create the timepoints dataframe required in the analysis.
