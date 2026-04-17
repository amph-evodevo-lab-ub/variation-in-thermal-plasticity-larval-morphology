# Variation in thermal plasticity of larval morphology among crested newt species and their reciprocal hybrids (Salamandridae: Triturus)

This repository includes R scripts and data for analysing the effect of elevated temperature on larval morphology in Triturus newts 
We explored effect of temperature on two morphological structures:
- Head: analysing size and shape variation
- Tail: analysing size and shape variation
- These two structures were analysed independently within two different R scripts ('tirturus_tail_script.R' and 'tirturus_head_script.R')
- The third script ('metamorphosis_barplot_and_fisher_test.R') is associated with analysis of metamorphic data

---

## Input Files

- `metamorphosis_individual.csv` – File defining affiliation of each individual to one of four metamorphic categories
- `metamorphosis_table.csv` – Table with count data used for creating the barplot
- `sliders_tail.csv` – File defining semilandmarks for tail shape data  
- `sliders_head.txt` – File defining semilandmarks for head shape data  
- `triturus_temperature_tail.tps` – Tail landmark coordinates
- `triturus_temperature_head.TPS` – Head landmark coordinates  
- `links_head.txt` – File defining landmark links used in visualization of head shape
- `lm_pairs_head.txt` – File defining pairs of bilateral landmarks for head data

---

## Dependencies

Install the following R libraries before running analyses:

```r
library(RRPP)
library(geomorph)
library(ggplot2)
library(abind)
library(tibble)
library(tidyr)
library(gridExtra)
```

---

## Workflow Overview ('tirturus_tail_script.R' and 'tirturus_head_script.R')

### Data Import & Preprocessing

- Read TPS coordinates
- Link with metadata
- Create classifier
- Adjust landmark orientation for tail data

---

### Procrustes Superimposition & Extracting Symmetric Component

Steps:
- Correct for size, orientation, and translation
- Slide semilandmarks to minimize bending energy
- Extract the symmetric shape component

_Output:_  
- Shape (symmetrized Procrustes coordinates)  
- Size (centroid size)

---

### Size differences
- Creating three-way crossed ANOVA model
- Pairwise comparison of group means
- Creating a boxplot for visualization of size differences

---

### Shape Differences

- Creating three-way crossed ANOVA model
- Pairwise comparison of group means
- Analysis of phenotypic change vectors
- Plotting genotype-specific vectors
- Visualizing genotype-specific shape changes

---

### Testing for dependence of shape data on rates of metamorphosis

- Loading metamorphic data
- Extracting the data from the end of the experiment and creating a new classifier
- One-way ANOVA testing for overall dependence
  ```r
  procD.lm(shape ~ metamo.chr, ...)
  ```
- Two-way crossed ANOVA testing for metamorphosis x genotype interaction
  ```r
  procD.lm(shape ~ metamo.chr * genotype, ...)
  ```
- Genotype-specific ANOVA models
- Visualizing the dependence of shape data on rates of metamorphosis for each genotype

---

## Workflow Overview ('metamorphosis_barplot_and_fisher_test.R')

### Data Import & Preprocessing

- Read metamorphic data
- Transform values from metamorphic table into percentages

---

### Visualisation

- Creating a data frame for the treatment group
- Creating ggplot object for treatment group's barplot
- Creating a data frame for the control group
- Creating ggplot object for control group's barplot
- Plotting
  ```r
  gridExtra::grid.arrange(y, x, nrow = 2)
  ```
- Performing Fisher's exact test for count data

---

## Citation

If using this workflow in your publication, please cite:

> Milic, M., Ivanovic, A., Nikolic, S., Avdalovic, A., Petrovic, T., Prokic, M., Vucic, T. (2025). Variation in thermal plasticity of larval morphology among crested newt species and their reciprocal hybrids (Salamandridae: Triturus) (2025), Under Review

---

## Contact

**For issues within script:**  
Mihajlo Milic  
University of Belgrade – Faculty of Biology  
Email: mihajlo.milic155@gmail.com

**For data requests:**  
Tijana Vucic  
Leiden University – Institute of Biology Leiden, University of Belgrade – Faculty of Biology  
Email: tijana.vucic@bio.bg.ac.rs 
