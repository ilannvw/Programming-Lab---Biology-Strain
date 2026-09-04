# Programming Lab - Biology Strain

A Python programme that reads experimental simulation files to be able to compare the population sizes and movements of normal bacteria compared to their mutants.

---

## Project Overview

### Context

When bacteria undergo mutation or become resistant to certain antibiotics, their behaviour changes. In this project, we will analyse computer simulation runs of bacterial cultures. Every run tracks individual cells, recording which strain they belong to, their proliferation and movement speed across three dimensions $(x, y, z)$ under controlled nutrient and stress conditions. 

When comparing normal strains (WT, positive IDs) to mutant strains (negative IDs) across different bacterial species, we want to see which strains dominate and whether movement speed plays a role. 

---

## Research Questions

1. **What are the counts of each bacterial strain and their statistical uncertainties?**

2. **Is there any asymmetry between the normal and the mutant strain? ($$\sigma_A = \frac{2 \sqrt{N_{\text{WT}} \cdot N_{\text{Mutant}}}}{(N_{\text{WT}} + N_{\text{Mutant}})^{3/2}}$$)**

3. **Is there any asymmetry as a function of the momentum? ($$p = \sqrt{p_x^2 + p_y^2 + p_z^2}$$)**

  To answer this, the total momentum $p$ is computed for each cell, and then grouped in momentum bins. Within each group, the asymmetry $\sigma_A$ between WT and mutant counts is calculated using the formula above. This will make it possible to see whether the balance between WT and mutant strains shifts at different momentum ranges, rather than overall.  

  ---

## Dataset and Input

The data files used for the simulation follow a structured format that represents actual experiment runs. 

* **Header Row:**
  * Event ID: Unique identifier for the simulation run
  * Number of bacteria tracked: Total cell count.
* **Data Rows:**
  * Momentum components: $p_x$, $p_y$, $p_z$ (in units of $10^{-20} kg·m/s$)
  * Bacterial ID: Numeric code denoting the bacterial strain. 

### Example Input File

```
1 29
0.153626 0.0787489 0.189973 -211
0.706237 0.293027 -0.188836 211
-0.195769 0.13282 -0.0321855 -211
-0.28052 0.0953441 0.235989 211
-0.174182 -0.401561 7.22433 211
0.100557 -0.302209 5.4045 -321
-0.155021 -0.185483 0.855223 321
-0.394574 0.00589716 1.70819 -321
-0.252323 -1.00072 2.8426 2212
...
```

Each row follows the format: `p_x p_y p_z Bacterial_ID`

### Bacterial ID Reference

Each row in the input files ends with the ID that identifies the bacterial strain. The mapping used is:

| Bacterial ID | Bacterial strain |
| --- | --- |
| 211 | E. coli WT |
| -211 | E. coli mutant |
| 321 | Bacillus subtilis WT |
| -321 | Bacillus subtilis mutant |
| 2212 | Pseudomonas aeruginosa WT |
| -2212 | Pseudomonas aeruginosa antibiotic-resistant |
| 3122 | Streptococcus pneumoniae |
| -3122 | Capsule-deficient Streptococcus pneumoniae |
| 3312 | Mycobacterium tuberculosis |
| -3312 | Drug-resistant Mycobacterium tuberculosis |
| 3334 | Salmonella enterica |
| -3334 | Salmonella mutant |

---

## Usage

### Prerequisites

* Python 3.10+
* A terminal on macOS, Linux, or Windows

### Cloning the Repository

```bash
  git clone https://github.com/ilannvw/Programming-Lab---Biology-Strain.git
```

### Install Dependencies

Install numpy and Matplotlib 
```bash
python3 -m pip install numpy matplotlib
```
Place the Data files in a folder called data
```bash
mkdir -p data
```
Run the main analysis script
```bash
python3 main.py
```

---

**Course:** Programming (PRA2003) - Maastricht University
**Author:** Ilan Noè
**Date:** September 2026