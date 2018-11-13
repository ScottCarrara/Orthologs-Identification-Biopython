import Bio
import sys
from Bio.Blast import NCBIXML

# The data structure used is the pYthon dictionary

# checking for correct number of arguments
if len(sys.argv) != 4 and len(sys.argv) != 5:
    print("Incorrect number of arguments supplied")
    sys.exit()

# FOR TESTING
# cd_file = open('celegans_drosophila.xml')
# dc_file = open('drosophila_celegans.xml')
# output = open('orthologs.txt', 'w')

# Opening the XML files
cd_file = open( sys.argv[1], 'r' )
dc_file = open( sys.argv[2], 'r' )

# Opening the output text file
output = open( sys.argv[3], 'w' )

#Parsing the crossed XML files using the NCBIXML module
cd_parse = NCBIXML.parse(cd_file)
dc_parse = NCBIXML.parse(dc_file)

# Setting the E_VALUE threshold
E_VALUE_THRESH = 10e-20

# using data structure Dictionary / Hash tables
# This is a dictionary of lists ( But the lists do not affect the Big Oh runtime )
# The key used is the match amino acid sequence and the value is a list containing the alignment title (protein name)
# , amino acid query seq and the amino acid subject(target) sequence
cd_dict = {}

# celegens_drosophila
# looping through the parsed XML file for celegens_drosophila
# Looping through the hits to generate our dictionary
for blast_record in cd_parse:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                cd_dict[ hsp.match+blast_record.query+alignment.title ] = [ blast_record.query, hsp.query, alignment.title, hsp.sbjct ]

# Using the same dictionary as cd_dict type
dc_dict = {}

# drosophila_celegens
# looping through the parsed XML file for drosophila_celegens
# Looping through the hits to generate our dictionary
for blast_record in dc_parse:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                dc_dict[ hsp.match+alignment.title+blast_record.query ] = [ blast_record.query, hsp.query, alignment.title, hsp.sbjct ]

# Printing the column labels for the output file
print( 'Species 1', 'Species 2', file=output, sep='\t' )

unique_dict = {} # dict to hold all matches without repeats
# Checking if the FALSE optional argument was not entered
# If not then, the definition of a reciprocal match is both a HSP from species A to B and species B to A

if len(sys.argv) == 4:
    # iterating through all the HSPs from the celegens_drosophila generated dict
    for key, value in cd_dict.items():
        # if there is a matched seq that is the same or reciprocal homologs then add it to the output file
        if key in dc_dict.keys():
            # dc_dict[key][2] = celegens protein name
            # cd_dict[key][2] = drosophila protein name
            keyCheck = dc_dict[key][2] + cd_dict[key][2]
            # Only add unique protein pairs
            if keyCheck not in unique_dict:
                unique_dict [ keyCheck ] = [ dc_dict[key][2], cd_dict[key][2]]
    # Writing the unique homologous pairs to the putput file
    for key, value in unique_dict.items():
        print( value[0], '\t', value[1], file=output )

# Else the optional FALSE argument might have been entered
# FALSE opt agument means the definition of a reciprocal match is no longer both a HSP from species A to B and
# species B to A, but is defined as a HSP from species A to B or species B to A.
else:
    # if the fourth argument is FALSE then
    if sys.argv[4] == 'FALSE':
        for key, value in cd_dict.items():
            # the key will be concatenation of the query protein name and hit/subject protein name
            # value[0] = celegens protein name
            # value[2] = drosophila protein name
            unique_dict[ value[0] + value[2] ] = [ value[0], value[2] ]

        for key, value in dc_dict.items():
            # value[2] = celegens protein name
            # value[0] = drosophila protein name
            keyCheck = value[2] + value[0]
            if keyCheck not in unique_dict:
                unique_dict[ keyCheck ] = [ value[2], value[0] ]
        # Writing the unique homologous pairs to the putput file
        for key, value in unique_dict.items():
            print( value[0], '\t', value[1], file=output )

    # if the fourth cmd line argument is not FALSE then return an error and quit
    else:
        print("Invalid argument(s) supplied")

output.close()
cd_file.close()
dc_file.close()
