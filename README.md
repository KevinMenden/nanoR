# Installation
## clone from git
```
git clone git@git.eu.boehringer.com:menden/nanoR.git
```
## R install 
```
install.packages("nanoR/",repos=NULL,type="source")
```


# Handling and analyzation of NanoString data
This package contains functions for parsing, handling, quality control and analyzation of data created on the
NanoString nCounter system in the form of RCC files.

## File Parsing
It is best to store all the RCC files of interest in a separate directory, or at least with in a directory with no RCC files
from other experiments. The parseRCC() function will then allow to transfer your files into a nano object:
```
nano <- parseRCC(dir = "path/to/RCC/directory/")
```

The nano object is a list containing the raw counts and the header from each RCC file. After normalization, other 
objects will be attached to it, e.g. normalization factors.

## Background Correction and Positive Control Normalization
The first step should be the calculation and subtraction of the background and positive control normalization, i.e. normalizing using the given spike-ins.
By default, the background is calculated using the mean of all background probes, and two standard deviations are added to it. You can adjust this using the "bm" and "sd.factor"
parameters. Similarly, you can choose which method should be used to calculate the positive control factors using the "pcm" parameter in the nsPositiveControlNormalization() function.
By default, the geometric mean is used.

```
nano <- nsBackgroundCorrect(nano)
nano <- nsPositiveControlNormalization(nano)
```


## Quality Control
The three basic plots that shold always be looked at are visualizing the FOV (Field Of View) counts, binding density and positive scaling factors.

```
plotBindingDensities(nano)
plotFOV(nano)
plotPositiveScalingFactors(nano)
```

## Content Normalization

### Decision over normalization method
Currently, three different methods are implemented for content normalization: top100, housekeeping and total (aka global).
The primary choice of the normalizaton method should be dependent on the experimental design. Generally, the following should be considered:

- Housekeeping gene normalization is a good choice if the trust in the housekeepers is high
- Global normalization can be used if only a relatviely small fraction of genes is expected to be differentially expressed
- Top100 normalization is usefull if only few genes of the panel are expressed above threshold (miRNA data)

Some functions have been included to compare different normalization methods. To asses correlation of the housekeepers, you can use:

```
plotHousekeepingGenes(nano)
housekeepingCorrelation(nano)
```

The first function will generate a line chart, which allows to see whether all houesekeeping genes follow the same pattern over all genes. Additionally,
the houseekepingCorrelation function will show a correlation matrix for the housekeeping genes. Housekeeping genes should show high correlation. To compare the
effect of all three normalization methods, you can use the following function:

```
plotNormDistanceRatio(nano, groups)
```

Where groups is a vector that assigns every sample to a specific group. The function then calculates the ratio of the mean inner-group distance to the mean global-distance
for every normalization and for the un-normalized counts. A smaller ratio after normalization implies better clustering of the groups. If groups are expected to be very different,
then this might be of interest.

### Normalization
After a content normalization method has been chosen, the counts can be normalized using the following function (houseekping normalization):
```
nano <- nsNormalize(nano, method="housekeeping")
```


## Downstream Analysis Functionality
After normalization, there exist several functions to visualize and explore the data. The following functions can be used for data visualization:

```
# Define a group vectors (examplary)
groups <- c("Group1", "Group1", "Group1", "Group1", "Group2", "Group2", "Group2", "Group2")
# Define a mean count that genes should have for consideration (default is 3 or 5 for most functions)
my_cutoff = 3

# Plot group-wise boxplot (log count values)
plotGroupBoxplot(nano, groups, title="my_title", countCutoff=my_cutoff)

# Plot scaled Heatmap of expression values
plotHeatmapExpression(nano, groups, countCutoff=my_cutoff)

# Plot a PCA
plotPCA(nano, groups, countCutoff=my_cutoff)

```

To make differential expression analysis, a limma based analysis and the t-test are both implemented.
To be continued ...
