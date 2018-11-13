import Bio
import sys
from Bio.Blast import NCBIXML

# The data structure used is the pYthon dictionary

# checking for correct number of arguments
if len(sys.argv) != 4 and len(sys.argv) != 5:
    print("Incorrect number of arguments supplied")
    quit()

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

# using array data structure to store the HSP information
# celegens_drosophila cross array
cd_array= []

# celegens_drosophila
# looping through the parsed XML file for celegens_drosophila
# Looping through the hits to populate our cd_array
for blast_record in cd_parse:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                cd_array.append( [hsp.match+blast_record.query+alignment.title, blast_record.query, hsp.query, alignment.title, hsp.sbjct ] )

# drosophila_celegens cross array
dc_array = []

# drosophila_celegens
# looping through the parsed XML file for drosophila_celegens
# Looping through the hits to populate our dc_array
for blast_record in dc_parse:
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < E_VALUE_THRESH:
                dc_array.append( [hsp.match+alignment.title+blast_record.query, blast_record.query, hsp.query, alignment.title, hsp.sbjct ] )

# Printing the column labels for the output file
print( 'Species 1', 'Species 2', file=output, sep='\t' )


# Checking if the FALSE optional argument was not entered or number of arguments is 4
# If not then, the definition of a reciprocal match is both a HSP from species A to B and species B to A
if len(sys.argv) == 4:
    unique_set_chk = set()
    # iterating through all the HSPs from the celegens_drosophila generated list
    for hsp_match_cd in cd_array:
        for hsp_match_dc in dc_array:
        # if there is a matched match_seq + query_seq + subject_seq ( match + celegens + drosophila )
        # add it to our set to ensure unique elements
            if hsp_match_cd[0] == hsp_match_dc[0]:
                unique_set_chk.add( hsp_match_cd[0] )

    # looping through match sequences that were present in both crosses
    checkSet = set() # To ensure only unique protein sets are added to this set
    for match_seq in unique_set_chk:
        for hsp_match_cd in cd_array:
            # if there is a reciprocal homologue match
            if match_seq == hsp_match_cd[0]:
                proteins_name_comb = hsp_match_cd[1]+hsp_match_cd[3]
                # if the protein name set is not already present in our set then
                if (proteins_name_comb) not in checkSet:
                    checkSet.add(proteins_name_comb)
                    print( hsp_match_cd[1], '\t', hsp_match_cd[3], file=output )


# Else the optional FALSE argument might have been entered
# FALSE opt agument means the definition of a reciprocal match is no longer both a HSP from species A to B and
# species B to A, but is defined as a HSP from species A to B or species B to A.
else:
    # two sets to hold the unique protein sets from the cd and dc arrays
    unique_cd_chk = set()
    unique_dc_chk = set()
    # if the fourth argument is FALSE then
    if sys.argv[4] == 'FALSE':
        # add all the unique protein name combination matches to unique sets
        for hsp_match_cd in cd_array:
            unique_cd_chk.add( hsp_match_cd[0] )
        for hsp_match_dc in dc_array:
            unique_dc_chk.add( hsp_match_dc[0] )

        # Checkset to ensure the protein name sets are uniquely added to this set
        checkSet = set()
        # looping through match sequences that were present in the unique_cd_chk set
        for match_seq in unique_cd_chk:
            for hsp_match_cd in cd_array:
                if match_seq == hsp_match_cd[0]:
                    # the protein name combination is celegens + drosophila proteins
                    proteins_name_comb = hsp_match_cd[1]+'\t'+hsp_match_cd[3]
                    # only adding unique protein sets
                    if proteins_name_comb not in checkSet:
                        checkSet.add(proteins_name_comb)

        # looping through match sequences that were present in the unique_dc_chk set
        for match_seq in unique_dc_chk:
            for hsp_match_dc in dc_array:
                if match_seq == hsp_match_dc[0]:
                    # the protein name combination is celegens + drosophila proteins
                    proteins_name_comb = hsp_match_dc[3]+'\t'+hsp_match_dc[1]
                    # only adding unique protein sets
                    if proteins_name_comb not in checkSet:
                        checkSet.add(proteins_name_comb)
        # print all the unique proteins to the file
        for proteins in checkSet:
            print( proteins, file=output )

    # if the fourth cmd line argument is not FALSE then return an error and quit
    else:
        print("Invalid argument(s) supplied")

output.close()
cd_file.close()
dc_file.close()
