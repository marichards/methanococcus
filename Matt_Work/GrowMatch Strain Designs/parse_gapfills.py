# Parse the gapfill files from CobraPy GrowMatch

# Import the parser and system packages
import re
import sys
# Import pickle
import pickle

# First, pull out all the different reaction combinations into a list
my_list = []
in_file = open(sys.argv[1],'r')
# Read the file into one big string
file_string = in_file.read()

# Split the file on "Run "
file_split = file_string.split('Run ')

# Turn it into a dictionary for each number and its reactions, minus any e0 ones
solution_dict = {}
for solution in file_split:
    # Test to make sure it has "Number" first
    if re.match('Number',solution):
        # Check if there's any e0s
        if not re.search('cpd\d{5}_e0',solution):
            # Split the record into lines 
            record_tag = solution.split('\n')[0]
            # Recognize that the first element is "Number: #", so grab the number
            record_num = re.search('\d+',record_tag).group(0)  
            # Make the values the reactions
            solution_dict[int(record_num)] = re.findall('rxn\d{5}.+',solution)
    
    
# Create an outfile of records and things
out_file = open('parsed_gapfills.txt','w')

# Write the dictionary
for k in sorted(solution_dict.keys()):
    out_file.write('Run Number: %d\n' %k)
    for value in solution_dict[k]:
        out_file.write('%s\n' %value)
    # Put a space between them
    out_file.write('\n---------------------------------\n')

# And save the dictionary using pickle
pickle.dump(solution_dict, open("rxn_dict.p","wb"))