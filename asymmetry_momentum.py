import statistics # for the mean and the standard deviation
import matplotlib.pyplot as plt # used to make the plot
# from main.py: the week 2 momentum function and the week 4 functions that analyse every file separately
from main import calculate_momentum, analyse_subsamples, calculate_mean_and_uncertainty

# week 5: asymmetry of one pair (WT and mutant), overall and as a function of momentum, using the sub-sampling method
# every file is a sub-sample and is analysed separately, the uncertainty is the standard deviation of the 10 files

# the pair that is analysed (the WT has the positive ID, the mutant has the negative ID)
# to analyse the other pair, change these four lines (3334, -3334 and the names) and the bin edges below
WT_ID = 2212
MUTANT_ID = -2212
WT_NAME = "Pseudomonas aeruginosa WT"
MUTANT_NAME = "Pseudomonas aeruginosa antibiotic-resistant"

# edges of the momentum bins (units = 10^-20 kg m/s): the first bin is 0 to 1, the second is 1 to 2, and so on
# 3334 / -3334 has very few bacteria, so for that pair use wider bins, for example [0, 2, 4, 6, 8, 12, 20]
BIN_EDGES = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 15, 20]

NUMBER_OF_FILES = 10 # output-Set1.txt to output-Set10.txt

# returns the number of the bin the momentum falls in (0 is the first bin), or None if it falls outside all the bins
def find_bin(momentum, bin_edges):
    for i in range(len(bin_edges) - 1):
        if bin_edges[i] <= momentum < bin_edges[i + 1]: # a momentum exactly on an edge goes to the upper bin
            return i
    return None

# goes through one file line by line and counts how many WT and how many mutant bacteria there are in every momentum bin
# then, it also returns how many fell outside all the bins, to check whether the bins cover the data
def count_momentum_bins(filepath):
    n_bins = len(BIN_EDGES) - 1
    wt_counts = [0] * n_bins # a counter per bin, and they start at 0
    mutant_counts = [0] * n_bins
    outside = 0

    rows_left_in_event = 0 # 0 means the next line will be a header (same logic as previous weeks)

    with open(filepath, "r") as infile:
        for line in infile:
            if rows_left_in_event == 0:
                # this line is a header line
                event_id, n_bacteria = line.split()
                rows_left_in_event = int(n_bacteria)
            else:
                # this is the data row
                parts = line.split()
                bacterial_id = int(parts[3])

                # only count bacteria that are the WT or the mutant of the target pair,
                # skips all the other IDs
                if bacterial_id == WT_ID or bacterial_id == MUTANT_ID:
                    px = float(parts[0])
                    py = float(parts[1])
                    pz = float(parts[2])
                    momentum = calculate_momentum(px, py, pz)
                    bin_index = find_bin(momentum, BIN_EDGES)

                    if bin_index is None:
                        outside += 1
                    elif bacterial_id == WT_ID:
                        wt_counts[bin_index] += 1
                    else:
                        mutant_counts[bin_index] += 1

                rows_left_in_event -= 1 # one fewer data row left before the next header
    return wt_counts, mutant_counts, outside

# calculates the asymmetry A between the WT and the mutant. A = (N_WT - N_mutant) / (N_WT + N_mutant)
# N_WT is the number of WT bacteria and N_mutant is the number of mutant bacteria
# A is 0 when there are equally many of both, and positive if there are more WT than mutant
# if there aren't any bacteria at all (0+0), we would divide by 0 crashing the programme, so
# we return None and the main function skips that bin
def calculate_asymmetry(n_wt, n_mutant):
    total = n_wt + n_mutant
    if total == 0:
        return None
    return (n_wt - n_mutant) / total

def main():
    n_bins = len(BIN_EDGES) - 1

    filepaths = []
    for i in range(1, NUMBER_OF_FILES + 1):
        filepaths.append(f"data/output-Set{i}.txt")

    # part 1: the week 4 results of the pair, from main.py
    # the week 4 functions need a dictionary with the IDs to analyse, here only the WT and the mutant
    pair_ids = {WT_ID: WT_NAME, MUTANT_ID: MUTANT_NAME}
    # for every ID, a list with its average count per event in each file
    subsample_averages = analyse_subsamples(filepaths, pair_ids)

    print("\n| Bacterial ID | Bacterial strain | Average count per event |")
    print("| --- | --- | --- |")
    for bacterial_id, strain_name in pair_ids.items():
        mean, uncertainty = calculate_mean_and_uncertainty(subsample_averages[bacterial_id])
        print(f"| {bacterial_id} | {strain_name} | {mean:.6f} ± {uncertainty:.6f} (stat) |")

    # part 2: overall asymmetry, with all the bacteria of the pair (no bins), one value per file
    # the average count per event can be used for A, because dividing both counts by the same number of events does not change A
    overall_asymmetries = []
    for i in range(NUMBER_OF_FILES):
        overall_asymmetries.append(calculate_asymmetry(subsample_averages[WT_ID][i], subsample_averages[MUTANT_ID][i]))

    # mean of the 10 files, uncertainty is the standard deviation of the 10 files (sub-sampling)
    overall_mean = statistics.mean(overall_asymmetries)
    overall_uncertainty = statistics.stdev(overall_asymmetries)
    print(f"\n{WT_NAME} vs {MUTANT_NAME}")
    print(f"overall asymmetry A = {overall_mean:.6f} ± {overall_uncertainty:.6f} (stat)")
    if overall_uncertainty > 0:
        print(f"significance = {abs(overall_mean) / overall_uncertainty:.2f} standard deviations")

    # part 3: asymmetry per momentum bin
    # one list per momentum bin, each list holds the asymmetry of that bin from every file
    bin_asymmetries = []
    for b in range(n_bins):
        bin_asymmetries.append([])

    # every file is analysed separately
    print("")
    for filepath in filepaths:
        print(f"analysing the momenta in {filepath}")
        wt_counts, mutant_counts, outside = count_momentum_bins(filepath)
        print(wt_counts, mutant_counts) # the WT and the mutant counts per bin, to check that the bins have enough bacteria

        if outside > 0:
            print(f"  {outside} bacteria are outside the bins")

        # asymmetry of every bin in this file
        for b in range(n_bins):
            asymmetry = calculate_asymmetry(wt_counts[b], mutant_counts[b])
            if asymmetry is not None: # a bin without any bacteria in this file is skipped
                bin_asymmetries[b].append(asymmetry)

    # asymmetry per bin: mean and standard deviation over the files, printed in a table for the results section in the README
    centres = [] # middle of every bin, the x-position of the point
    half_widths = [] # half the width of the bin, used for the horizontal bars
    means = []
    uncertainties = []

    print("\n| Momentum bin | Asymmetry A |")
    print("| --- | --- |")
    for b in range(n_bins):
        low = BIN_EDGES[b]
        high = BIN_EDGES[b + 1]
        values = bin_asymmetries[b] # the asymmetries of this bin, one per file

        if len(values) < 2: # the standard deviation needs at least 2 files
            print(f"| {low} - {high} | not enough data |")
        else:
            mean = statistics.mean(values)
            uncertainty = statistics.stdev(values)
            print(f"| {low} - {high} | {mean:.5f} ± {uncertainty:.5f} (stat) |")
            centres.append((low + high) / 2)
            half_widths.append((high - low) / 2)
            means.append(mean)
            uncertainties.append(uncertainty)

    # plot: asymmetry vs momentum, vertical bars = statistical uncertainty, horizontal bars = width of the bin
    plt.figure(figsize=(8, 5))
    plt.errorbar(centres, means, xerr=half_widths, yerr=uncertainties, fmt="o", color="black", capsize=3)
    plt.axhline(0, linestyle="--", color="grey") # A = 0 means no asymmetry
    plt.xlabel("momentum p ($10^{-20}$ kg m/s)")
    plt.ylabel("asymmetry A")
    plt.title(f"{WT_NAME} vs {MUTANT_NAME}")
    plt.tight_layout()
    plt.savefig(f"asymmetry_{WT_ID}.png", dpi=200) # saved so it can be used in the README and the presentation
    plt.show()

if __name__ == "__main__":
    main()