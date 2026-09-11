import math # gives access to math.sqrt() for the momentum calculation

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

def main():
    # reads file and prints out each bacteria id and its corresponding total momentum
    filepath = "data/output-Set0.txt"
    event_id, bacteria = read_event(filepath) # reads file
 
    print(f"event ID: {event_id}")
    print(f"number of bacteria tracked: {len(bacteria)}\n") 
    # loops over each bacterium and calculates and prints the momentum
    for px, py, pz, bacterial_id in bacteria:
        momentum = calculate_momentum(px,py,pz)  
        print(bacterial_id, momentum)
 
main()