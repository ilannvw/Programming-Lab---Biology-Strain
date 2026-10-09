import statistics # for the mean and the standard deviation
# retrieve from main.py the bacterial IDs and the function that analyses every 
# file separately
from main import get_bacterial_ids, analyse_subsamples
# retrieve from asymmetry_momentum.py the asymmetry function
from asymmetry_momentum import calculate_asymmetry

# total asymmetry of every WT / mutant pair, by using the sub-sampling method
# for that pair, A is calculated for each file, and the result is the mean of the 
# 10 files and the uncertainty is the standard deviation

NUMBER_OF_FILES = 10 # output-Set1.txt to output-Set10.txt

def main():
    ids_of_interest = get_bacterial_ids()

    filepaths = []
    for i in range(1, NUMBER_OF_FILES + 1):
        filepaths.append(f"data/output-Set{i}.txt")

    # every file is read once, giving the average count per event of all 12 IDs in every file
    subsample_averages = analyse_subsamples(filepaths, ids_of_interest)

    print("\n| Pair | Asymmetry A | Significance |")
    print("| --- | --- | --- |")
    for wt_id, wt_name in ids_of_interest.items():
        if wt_id > 0: # only the WT IDs
            mutant_id = -wt_id

        # asymmetry of this pair in every file
            asymmetries = []
            for i in range(NUMBER_OF_FILES):
                asymmetry = calculate_asymmetry(subsample_averages[wt_id][i], subsample_averages[mutant_id][i])
                if asymmetry is not None: # a file without any bacteria of this pair is skipped
                    asymmetries.append(asymmetry)

            # mean and standard deviation of the files (sub-sampling), and the significance = |A| / uncertainty
            mean = statistics.mean(asymmetries)
            uncertainty = statistics.stdev(asymmetries)
            significance = abs(mean) / uncertainty
            print(f"| {wt_name} ({wt_id} / {mutant_id}) | {mean:.6f} \u00b1 {uncertainty:.6f} | {significance:.1f} \u03c3 |")

if __name__ == "__main__":
    main()