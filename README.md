# inForm2fcs
Utility for working with data tables produced by inForm software (PerkinElmer) for cell segmentation, phenotyping, quantification, and localization of Vectra scans.
Enable flow-style analysis of Vectra tissue scans for inForm segmentation results.
Create FCS files, grouped by "Tissue Category" and "Phenotype" per each "Slide ID".

# Installation


## Install required R packages

If you haven't already, you will need to install these packages:
* flowCore
* Shiny

#### flowCore
At the R command line enter:
```
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("flowCore")
```

#### Shiny
```
install.packages("shiny")
```


#### inForm2fcs
Download inForm2fcs.R at the top of this page.



# Usage

In an R session, navigate to the directory containing inForm2fcs.R, and start it by typing:

```
source("inForm2fcs.R")
runApp("inForm2fcs.R")
```

In the user interface, select a cell segmentation file produced by inForm. These files typically end in "\_cell\_seg\_data.txt"

This will create one FCS file for each Tissue Category, Phenotype, and Slide ID, named accordingly.

Additionally, one FCS file will be created for each slide containing all cells for that slide: slideName\_All.fcs

Exported fields include Mean and Total (if present) for every marker, as well as:
* Cell X Position
* Cell Y Position
* Nucleus Area (pixels)
* Confidence

FCS output files will be sent to the directory containing inForm2fcs.R



## Command Line Usage
It will often be much faster to bypass the Siny interface and create FCS files directly from the command line:
```
make_fcs(fname = "/path/to/_cell_seg_data.txt", outdir = "/path/to/output/directory/")
```
###### fname:  
Full path to input \_cell\_seg\_data.txt file.
###### outdir:  
Directory to write resulting FCS files. This directory must exist before calling make\_fcs(). Defaults to the current working directory.





