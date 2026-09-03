# Programming Lab - Biology Strain

A Python programme that reads experimental simulation files to be able to compare the population sizes and movements of normal bacteria compared to their mutants.

## Project Overview

### Context

When bacteria mutate or become resistant to certain antibiotics, their behaviour changes.

In this project, we will analyse computer simulation runs of bacterial cultures. Every run tracks individual cells, recording which strain they belong to, their proliferation and movement speed across three dimensions (x, y, z) under controlled nutrient and stress conditions. 

When comparing normal strains (WT, positive IDs) to mutant strains (negative IDs) across different bacterial species, we want to see which strains dominate and whether movement speed plays a role. 

## Research Questions

1. **What are the counts of each bacterial strain and their statistical uncertainties?**

2. **Is there any asymmetry between the normal and the mutant strain? ($$\sigma_A = \frac{2 \sqrt{N_{\text{WT}} \cdot N_{\text{Mutant}}}}{(N_{\text{WT}} + N_{\text{Mutant}})^{3/2}}$$)**

3. **Is there any asymmetry as a function of the momentum? ($$p = \sqrt{p_x^2 + p_y^2 + p_z^2}$$)**

## Dataset and Input

The data files used for the simulation follows a structured format that represent actual experiment runs. 

* **Header Row:**
    * Event ID': Unique identifier for the simulation run
    * Number of bacteria tracked: Total cell count.
* **Data Rows:**
    * Momentum components: $p_x$, $p_y$, $p_z$ 
    * Bacterial ID: Numeric code denoting the bacterial strain. 

## Cloning the Repository and Setup

### Prequisites

* Python 3.10+
* A terminal on macOS, Linux, or Windows

### Cloning the Repository

```bash
git clone [https://github.com/ilannvw/Programming-Lab---Biology-Strain.git](https://github.com/ilannvw/Programming-Lab---Biology-Strain.git) 
cd Programming-Lab---Biology-Strain
```

### Install Dependencies

Install numpy and Matplotlib 
```bash
python3 -m pip install numpy matplotlib
```
Place the Data files in a folder call data
```bash
mkdir -p data
```
Run the main analysis script
```bash
python3 main.py
```