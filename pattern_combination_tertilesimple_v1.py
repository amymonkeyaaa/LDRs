#printing all headers so can find and tally patterns
from pathlib import Path
import pandas as pd

input_file = Path("/Users/monkeyaaa/Documents/GitHub/LDRs/merged files/hpep_blockiness_merged_v3.txt")

#store every MERGED_GOTERM_SET header
headers = []

with input_file.open("r") as f:
    for line in f:
        line = line.strip()

        #skip blank lines
        if not line:
            continue

        #keep merged header rows
        if line.startswith("#MERGED_GOTERM_SET"):
            headers.append(line)

print("All MERGED_GOTERM_SET headers:")

#print each header (should be 33) to verify headers list
#for header in headers[:33]:
    #print(header)

#make nested list within headers list to store each header and its associated patterns
header_patterns = []

for header in headers:
    col = header.split("\t") #list
    header_patterns.append(col) #append list to header_patterns to make nested list

#verify header index and value
#for i, value in enumerate(header_patterns[0]):
    #print(i, value)

#index 0 = #MERGED_GOTERM_SET
#index 1 = GO term
#index 2 = bias type
#index 3 = blockiness tertile
#index 4 = blockiness HxIyLz count
#index 5 = hpep tertile
#index 6 = hpep HxIyLz count
#index 7 = number of rows in this block (matched entries)
#index 8 = GO term

#use a dictionary to store the details/counts of each unique pattern combination
#key is combination, value is a list of all headers with that combination
#first letter is blockiness tertile, second letter is hpep tertile
pattern_dict = {
    "HH": [],
    "HI": [],
    "HL": [],
    "IH": [],
    "II": [],
    "IL": [],
    "LH": [],
    "LI": [],
    "LL": []
}

for header in header_patterns:
    #get the two tertiles
    blockiness = header[3]
    hpep = header[5]
    #build the key
    combo = blockiness + hpep
    #store complete header in appropriate key list
    pattern_dict[combo].append(header)

#check HH pattern:
#for header in pattern_dict["HH"]:
    #print(header)

output_file = Path("/Users/monkeyaaa/Documents/GitHub/LDRs/merged files/pattern_combination_tertile_v1.txt")
with output_file.open("w") as f:

    for pattern, headers in pattern_dict.items():
        #pattern heading
        f.write(f"#Pattern: {pattern}\n")
        f.write("GO_term\tBias\tn_matched\n")
        #one line per merged header
        for header in headers:
            go_term = header[8]
            bias = header[2]
            n_matched = header[7]
            f.write(f"{go_term}\t{bias}\t{n_matched}\n")

        f.write("\n") #add blank line between patterns
