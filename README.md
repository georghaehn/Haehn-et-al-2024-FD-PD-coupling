Summary of the workflow for Haehn et al. 2024 “Global decoupling of
functional and phylogenetic diversity in plant communities”.
================
03/12/2024

Hähn, Georg J. A., Gabriella Damasceno, Esteban Alvarez-Davila, Isabelle
Aubin, Marijn Bauters, Erwin Bergmeier, Idoia Biurrun, et al. 2024.
“Global Decoupling of Functional and Phylogenetic Diversity in Plant
Communities”. *Nature Ecology & Evolution*:
<https://doi.org/10.1038/s41559-024-02589-0>

### Executive summary

1.  Creating the phylogenetic tree based on the taxonomic backbone of
    sPlot 3.0

2.  Calculating functional and phylogenetic diversity and standardized
    effect size

3.  Prepare the data for statistical analysis

4.  Statistical analysis - relationship of FD and PD, environmental
    drivers

5.  Plotting of the distribution of functional and phylogenetic
    diversity and the results of the models

6.  Evaluate the method

### About this README

This README is made to guide you trough the analysis to calculate
phylogenetic <br> and functional diversity in sPlot 3.0. and to show how
the additional analysis <br> were performed to test

1.  whether patterns of coupling or decoupling dominate at the global
    level,
2.  have regional patterns,
3.  differ between forest and non-forest ecosystems, and
4.  correlated with current and past climatic gradients.

The used scripts are made to run on three different kinds of machines:
<br>

1.  “normal” computer (~16GB of RAM) <br>
2.  RStudio Server (\>50GB of RAM) <br>
3.  HPC Cluster (\>16GB of RAM per core) <br>

In each script you will find a number referring to the recommended
machine, <br> in addition you will find (if needed) the .sh scripts to
run the scripts on a Cluster <br> (please, check if your Cluster match
the properties). <br> The author highly recommend to do so, the scripts
were already optimized to run with <br> the lowest resources possible.

### Reproducing the code

### Step 0: The idea

Code to produce Figure 1.

``` r
source("00.Concept-figure.R", local = knitr::knit_global())
```

### Step 1: The phylogentic tree

Calculate the phylogenetic tree <br>

The $V.Phylomaker$ package was used to build the phylogenetic tree.
Additional binding was done with $phytools$.

``` r
source("01.phylogenetic_tree_calculation.R", local = knitr::knit_global())
```

### Step 2: Calculate the diversity indices

#### 2.1 Prune the tree

Phylogenetic trees were pruned to contain only the species of one data
subset. <br> This was done for each unique dataset (GIVD) in sPlot and
for sPlotOpen. <br>

``` r
source("02.01.phylogenetic_pruning.R", local = knitr::knit_global())
```

#### 2.2 Calculate phylogenetic diversity (RQEP and MPD)

``` bash
sbatch -a 1-2000000:5000 02.02.phylo-calculations.sh
sbatch -a 1-95000:5000 02.02.phylo-calculations.sPlotOpen.sh
```

#### 2.3 Phylogenetic NULL calculation

To calculate standardized effect size of phylogenetic diversity, <br>
for each species richness 499 times random species communities were
sampled <br> and RQEP and MPD calculated based on global, biome, and
phylogenetic cluster specific species pools. <br>

``` bash
sbatch -a 2-718:2 02.03.phylo-calculations-null.sh
sbatch -a 2-718:2 02.03.01.phylo-calculations-null-biome.sh
sbatch -a 2-718:2 02.03.02.phylo-calculations-null-biome-weighted.sh
sbatch -a 2-718:2 02.03.03.phylo-calculations-null-phyl-clust.sh
```

#### 2.4 Calculate functional diversity

``` bash
sbatch -a 1-2000000:5000 02.04.funct-div-calculations.sh
sbatch -a 1-95000:5000 02.04.funct-div-calculations.sPlotOpen.sh
```

#### 2.5 Functional NULL calculation

Analogue to 2.3.

``` bash
sbatch -a 2-718:2 02.05.funct-div-calculations-null.sh
sbatch -a 2-718:2 02.05.01.funct-div-calculations-null-biome.sh
sbatch -a 2-718:2 02.05.02.funct-div-calculations-null-biome-weighted.sh
sbatch -a 2-718:2 02.05.03.funct-div-calculations-null-phyl-clust.sh
```

### Step 3: Prepare the data and add the response variables

#### 3.1 Indicies handling

Receiving the data from the Cluster calculations. Calculating mean and
standard <br> error for the NULL models.

``` r
source("03.01.indices-handling.R", local = knitr::knit_global())
source("03.01.indices-handling.sPlotOpen.R", local = knitr::knit_global())
```

#### 3.2 Explanatory variables

Downloading the climate variables from “CHELSA” and fitting the PCA.
Extracting the climate <br> condition from the last glacial maximum from
“StableClim” <br> The data was available at a 2.5° spatial scale and
were extracted with the <br> *extract* function from the $raster$
package for each vegetation-plot.

Two additional variables were present in the sPlot database:

1.  which plants were recorded in the plot (e.g., all vascular plants,
    only dominant species, all woody plants, only trees)
2.  vegetation type (forest vs. non-forest)

``` r
source("03.02.explanatory-climate-variables.R", local = knitr::knit_global())
```

#### 3.3 Standardized effect size

Calculating the standardized effect size based on the calculated NULL
data and <br> and merging the additional explanatory data to the
dataset.

``` r
source("03.03.SES-calculation_merging.R", local = knitr::knit_global())

source("03.04.SES-calculation_BIOME.R", local = knitr::knit_global())
```

### Step 4: Statistical analysis

#### 4.1 The relationship of functional and phylogeentic diversity

A generalised additive model (GAM) was used to analyse the relationship
between <br> functional and phylogenetic diversity. A GAM is a
generalised linear model in <br> which the linear response depends on
unknown smooth functions of the explanatory <br> variables. To account
for the spatial structure of the data, the spatial <br> coordinates were
included as smooth spherical splines. All GAMs included a basis <br>
penalty smoother spline on the sphere (bs = ”sos”), applied to the
geographic <br> coordinates of every site, thus taking spatial
autocorrelation into account. <br> The model was performed using the
*gam* function from the $mgcv$ package.

The script also includes the code to produce Figure 2 A, Figure S 1 and
Figure S 2.

``` r
source("04.01.FD-PD-GAM.R", local = knitr::knit_global())
```

#### 4.2 Boosted Regression Trees

The parameters of the BRT were set as follows: a tree complexity of five
and <br> a bag fraction of 0.5. The learning rate was set to 0.01 with a
maximum number <br> of 20,000 trees. The BRTs were calculated using the
*gbm.step* routine from the $dismo$ package.

``` bash
sbatch -a 1-2:1 04.02.BRT.sh
source("04.03.BRT-BW.R", local = knitr::knit_global())
```

#### 4.3 Explaining functional and phylogenetic diversity

The variables that were considered as relevant from the BRTs were used
in a GAM <br> as explanatory variables of functional and phylogenetic
diversity. The same syntax <br> as for the relationship between
functional and phylogenetic diversity was used <br> for the GAMs. This
means for each response variable (SES.RQEP, SES.RQEF) <br> a GAM was
performed with the explanatory variables that turned out to be relevant
<br> in the BRTs. The spatial coordinates were included as smooth
spherical splines in <br> the model as explained above.

``` r
source("04.03.GAM-expl-var.R", local = knitr::knit_global())
```

### Step 5: Plotting

#### 5. The distribution of functional and phylogeentic diversity

Functional and phylogenetic variables where plotted as mean for each
grid cell with <br> a size of 863.8 km^2.

``` r
source("05.00.spatial-distribution.R", local = knitr::knit_global())
```

#### 5.1 Relative influence of the explanatory variables.

Code to produce Figure 3 and Figure S 3 - 6. The fitted smooth of the
BRT function was plotted between the $2.5^{th}$% quantile <br> and the
$97.5^{th}$% quantile to remove the very extreme points.

``` r
source("05.01.BRT-summary.R", local = knitr::knit_global())
```

#### 5.2 The GAM results

Code to produce Figure 4 and 5.

``` r
source("05.02.GAM-expl-plots.R", local = knitr::knit_global())
```

### Step 6: Evaluation

Calculation of phylogenetic signals and the relationship of FD and PD
when FD is based on different traits. Code to produce Figure S 7.

``` r
source("06.phylosignals-trait-evaluation.R", local = knitr::knit_global())
```
