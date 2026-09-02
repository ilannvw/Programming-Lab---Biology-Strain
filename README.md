# Programming Lab - Biology Strain

## Project Overview

This project investigates how bacteria moves and proliferates under controlled nutrient and stress conditions. Using simulated data from bacterial experiments, the analysis compares how the wild-type (WT) strains differ from the mutant / drug-resistant strains. 

By analysing the cellular counts and the three-dimensional momentum distribution, this project aims to identify the phenotypic variations, the advantages of selective proliferation, and possible kinematic differences that come from genetic variation or environmental stressors.

## Research Questions

1. **What are the counts of each bacterial strain and their statistical uncertainties?**

2. **Is there any asymmetry between the normal and the mutant strain? ($$\sigma_A = \frac{2 \sqrt{N_{\text{WT}} \cdot N_{\text{Mutant}}}}{(N_{\text{WT}} + N_{\text{Mutant}})^{3/2}}$$)**

3. **Is there any asymmetry as a function of the momentum? ($$p = \sqrt{p_x^2 + p_y^2 + p_z^2}$$)**

## Dataset and Input

The data used for the simulation follows a structured format that represent actual experiment runs. 

* **Header Row:**
    * Event ID': Unique identifier for the simulation run
    * Number of bacetria tracked: Total cell count.
* **Data Rows:**
    * Momentum components: $p_x$, $p_y$, $p_z$ 
    * Bacterial ID: Numeric code denoting the bacterial strain. 

## Cloning the Repository
```bash
git clone [https://github.com/ilannvw/Programming-Lab---Biology-Strain.git](https://github.com/ilannvw/Programming-Lab---Biology-Strain.git) 
cd Programming-Lab---Biology-Strain