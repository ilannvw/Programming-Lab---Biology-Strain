import math # gives access to math.sqrt() for the momentum calculation

# week 2: read a single event and calculate the momentum of each bacterium

# reads a single run/event from the file and returns the event ID and list of bacteria
def read_event(filepath):
    # open the file in read mode
    with open(filepath, "r") as infile:
        lines = infile.readlines() # reads all the lines and puts into a list of strings
        event_id, n_bacteria = lines[0].split() # splitting header into event id, and the amount of bacteria
        event_id = int(event_id) # converts event_id into an integer
        n_bacteria = int(n_bacteria) # converts n_bacteria into an integer

        bacteria = [] # empty list to input each bacterium's data
        # looping over the data rows, skipping the header
        # lines[1 : 1 + n_bacteria] selects the exact amount of rows the file has
        for line in lines [1:1 + n_bacteria]: 
            px, py, pz, bacterial_id = line.split() # splitting row into 4 string pieces and naming them to their corresponding value
            bacteria.append((float(px),float(py),float(pz),int(bacterial_id))) # convert each string into their corresponding type
    return event_id, bacteria # returns both values when function is called

# calculates the total momentum of each bacterium
def calculate_momentum(px, py, pz):
    total = px**2 + py**2 + pz**2 # sum of squares of the momentum components
    if total >= 0: # check total value to make sure square root isn't negative
        p = math.sqrt(total) # value of total momentum
        return p
    else: # this should never run as total value should be always positive, but it protects programme from crashing
        print("invalid input")
        return None
    
# week 3: read a file with many events, selects the bacteria of interest, and calculates 
# the average count and uncertainty of each strain per event.

# returns a dictionary mapping each bacterial ID to its strain name
def get_bacterial_ids():
    ids_of_interest = {
        211: "E. coli WT",
        -211: "E. coli mutant",
        321: "Bacillus subtilis WT",
        -321: "Bacillus subtilis mutant",
        2212: "Pseudomonas aeruginosa WT",
        -2212: "Pseudomonas aeruginosa antibiotic-resistant",
        3122: "Streptococcus pneumoniae",
        -3122: "Capsule-deficient Streptococcus pneumoniae",
        3312: "Mycobacterium tuberculosis",
        -3312: "Drug-resistant Mycobacterium tuberculosis",
        3334: "Salmonella enterica",
        -3334: "Salmonella mutant"
    }
    return ids_of_interest

# streams through the file one line at a time (instead of storing the whole file) and counts how many times each id
# of interest appears, and the total number of events that have been processed
def count_bacteria_per_event(filepath, ids_of_interest):
    counts = {} # starting with empty dictionary and every ID starts at 0
    for bacterial_id in ids_of_interest:
        counts[bacterial_id] = 0

    n_events = 0
    rows_left_in_event = 0 # 0 means the next line will be a header

    with open(filepath, "r") as infile:
        for line in infile:
            if rows_left_in_event == 0:
                # this line is a header line: "event_id n_bacteria"
                event_id, n_bacteria = line.split()
                n_bacteria = int(n_bacteria)
                rows_left_in_event = n_bacteria # expecting this many rows in the event
                n_events += 1
            else:
                # this is the data row, px py pz bacterial_id
                parts = line.split()
                bacterial_id = int(parts[3]) # since only the ID matters
                if bacterial_id in ids_of_interest: # only count bacteria that matches one of the 12 IDs of interest
                    counts[bacterial_id] += 1

                rows_left_in_event -= 1 # one fewer data row left before the next header
    return counts, n_events

# calculates the average count and uncertainty in each event
# for one bacterial ID, given its total count across all events
def calculate_average_and_uncertainty(total_count, n_events):
    average = total_count / n_events
    uncertainty = math.sqrt(total_count) / n_events # Poisson uncertainty on the total, divided by n_events
    return average, uncertainty

def main():

    # week 2: momentum calculation for single event
    filepath = "data/output-Set0.txt"
    event_id, bacteria = read_event(filepath) # reads file
 
    print(f"event ID: {event_id}") # use of f (formatted string), makes it easier to read
    print(f"number of bacteria tracked: {len(bacteria)}\n") 
    # loops over each bacterium and calculates and prints the momentum
    for px, py, pz, bacterial_id in bacteria:
        momentum = calculate_momentum(px,py,pz)  
        print(bacterial_id, momentum)

    # week 3: average count per event

    print("\nweek 3: average bacteria counts per event\n")
    full_filepath = "data/output-Set6.txt"
    ids_of_interest = get_bacterial_ids()

    counts, n_events = count_bacteria_per_event(full_filepath, ids_of_interest)

    print(f"number of events processed: {n_events}\n")

    for bacterial_id, strain_name in ids_of_interest.items():
        total = counts[bacterial_id]
        average, uncertainty = calculate_average_and_uncertainty(total, n_events)
        print(f"{strain_name} ({bacterial_id}): average = {average:.4f} \u00b1 {uncertainty:.4f} per event") # \u00b1 prints the ± sign

 
if __name__ == "__main__":
    main()