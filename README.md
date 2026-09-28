# Programming Lab - Biology Strain

A Python programme that reads simulation run data to be able to compare the population sizes and movements of normal bacteria compared to their mutants.

---

## Project Overview

### Context

In this project, we will analyse  simulation run data of bacterial cultures. Every run tracks individual cells, recording which strain they belong to and their momentum across three dimensions $(x, y, z)$ under controlled conditions.  

When comparing the wild-type (WT, positive IDs) to mutant strains (negative IDs) across different bacterial species, we want to see which strains dominate.

---

## Research Questions

1. **What are the average counts of each bacterial strain and their statistical uncertainties?**

*See the Method and Results sections below.*

2. **Is there any asymmetry between the normal and the mutant strain?**

  We quantify the relative population asymmetry ($A$):
  $$A = \frac{N_{\text{WT}} - N_{\text{Mutant}}}{N_{\text{WT}} + N_{\text{Mutant}}}$$

  Then we determine its statistical significance using the error propagation formula: 
  $$\sigma_A = \frac{2 \sqrt{N_{\text{WT}} \cdot N_{\text{Mutant}}}}{(N_{\text{WT}} + N_{\text{Mutant}})^{3/2}}$$

3. **Is there any asymmetry as a function of their momentum? ($p = \sqrt{p_x^2 + p_y^2 + p_z^2}$)**

  The total momentum $p$ is computed for each cell, and then grouped in momentum bins. Within each group, the asymmetry $A$, and the uncertainty $\sigma_A$ between WT and mutant counts are calculated using the formulas above. This will make it possible to see whether the balance between the WT and the mutant strains shifts at different momentum ranges, rather than overall.  

---

## Dataset and Input

The data files used for the simulation follow a structured format that represents actual experiment runs. 

* **Header Row:**
  * Event ID: Unique identifier for the simulation run
  * Number of bacteria tracked: Total cell count.
* **Data Rows:**
  * Momentum components: $p_x$, $p_y$, $p_z$ (in units of $10^{-20}$ kg·m/s)
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

A single file may contain many events (header, that many data rows, header, and so on), the full dataset consists of ten files, `output-Set1.txt` to `output-Set10.txt`, each containing 500,000 events in this format, meaning there are 5,000,000 events in total. 
The file `output-Set0.txt` is a small test file with a single event.

### Bacterial ID Reference

Every row in the data files ends with the ID that identifies the bacterial strain. The mapping used is:

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

## Method

### Counting

Every file is read one line at a time, so the whole 770 MB file never has to be fully stored. The header of each event gives the number of data rows there are, and the bacterial ID at the end of every row is checked against the 12 existing IDs in the table above (WT and mutant species are counted separately). Any other ID in these files will be ignored. Dividing the total count of an ID by the number of events in the file will give the average count per event of that ID. 

### Sub-sampling and statistical uncertainty

The statistical uncertainties are calculated with the sub-sampling method:

1. Every one of the ten files is analysed separately (one sub-sample is 500,000 events), which gives ten average counts per event for every ID. 
2. The central value is the average count per event over all the 5,000,000 events. Since all the ten sub-samples are the same size, this is just the plain mean of the ten sub-samples
3. The statistical uncertainty is the sample standard deviation (n - 1) of the ten sub-sample averages, meaning it is the spread of the sub-sample results.

The data is simulated, so no systematic uncertainty is needed.

--- 

## Results

Average count per event of each bacterial strain (Research Question 1), from the analysis of all 5,000,000 events. Every value is given as central value ± statistical uncertainty. The units of the results are counts per event. 

| Bacterial ID | Bacterial strain | Average count per event |
| --- | --- | --- |
| 211 | E. coli WT | 18.425338 ± 0.027777 (stat) |
| -211 | E. coli mutant | 18.395508 ± 0.026927 (stat) |
| 321 | Bacillus subtilis WT | 2.317445 ± 0.004011 (stat) |
| -321 | Bacillus subtilis mutant | 2.312189 ± 0.004780 (stat) |
| 2212 | Pseudomonas aeruginosa WT | 1.115739 ± 0.001647 (stat) |
| -2212 | Pseudomonas aeruginosa antibiotic-resistant | 1.093689 ± 0.002087 (stat) |
| 3122 | Streptococcus pneumoniae | 0.255466 ± 0.000963 (stat) |
| -3122 | Capsule-deficient Streptococcus pneumoniae | 0.250938 ± 0.000891 (stat) |
| 3312 | Mycobacterium tuberculosis | 0.036428 ± 0.000254 (stat) |
| -3312 | Drug-resistant Mycobacterium tuberculosis | 0.036021 ± 0.000367 (stat) |
| 3334 | Salmonella enterica | 0.001096 ± 0.000038 (stat) |
| -3334 | Salmonella mutant | 0.001064 ± 0.000047 (stat) |

---



## Usage

### Prerequisites

* Python 3.10+
* A terminal on macOS, Linux, or Windows

### Cloning the Repository

```bash
  git clone https://github.com/ilannvw/Programming-Lab---Biology-Strain.git
```

### Optional / Future Dependencies

Currently, `main.py` only uses Python's built-in `math` and `statistics` packages, no external libraries are needed yet. `numpy` and `matplotlib` will be needed for plotting.

```bash
python3 -m pip install numpy matplotlib
```

The data files are expected in a `data` folder; create it if it doesn't exist yet:

```bash
mkdir -p data
```
Windows equivalents:

```bash
mkdir data
```

Eleven data files are used:
- `data/output-Set0.txt`: single-event test file, used to check the momentum per cell calculation.
- `data/output-Set1.txt` to `data/output-Set10.txt`: the ten sub-samples of 500,000 events each. 

Each file is around 770 MB and is not included in this repository. Must be downloaded separately from the surfdrive provided in the course materials and placed in the `data` folder before running the code, otherwise it will not work.

Run the main analysis script:

```bash
python3 main.py
```
Windows equivalents:
```bash
python main.py
```
---

### Output

Running `main.py`:

1. Reads a single event from `data/output-Set0.txt` and prints:
  - the event ID and the total number of bacteria tracked
  - for each bacterium: its bacterial ID and total momentum magnitude $p = \sqrt{p_x^2 + p_y^2 + p_z^2}$

2. Reads all 500,000 events from `data/output-Set6.txt` and prints:
  - the average count of all 12 IDs, as an example of a single sub-sample
  - the statistical uncertainty, calculated from the Poisson statistics on the total count summed across all the events

3. Analyses all ten files separately and prints a table with the average count per event of each of the 12 bacterial IDs over the 5,000,000 events, with the statistical uncertainty from the sub-sampling method (table in the Results section).

---

**Course:** Programming (PRA2003) | Maastricht University
**Author:** Ilan Noè
**Date:** September 2026