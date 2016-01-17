
# --------------------------------------------- PART 2 -------------------------------------------------------

# 1) Given a protein FASTA file (filename) returns an integer corresponding to the total number of protein sequences
#having a relavie frequency higher or equal than a given threshold for a given residue.

fasta_file = open("ex.fasta")

def get_sequence_from_FASTA_file(filename, residue, threshold=0.3):
	integer, match, rel_freq = 0, 0, 0
	protein = ""
	
	for line in filename:
		if not line.startswith(">"):
			protein += line
		else:
			if protein: 
				for aa in protein:
					if aa.upper() == residue.upper():
						match+=1
				rel_freq = match / len(protein)
				if rel_freq >= float(threshold): 
					integer+=1
			protein = ""
			match = 0
			rel_freq = 0

	# for the very last protein, because it is stored without entering
	# in the else, so it wouldn't be taken into account
	for aa in protein:
		if aa.upper() == residue.upper():
			match+=1
	rel_freq = match / len(protein)
	if rel_freq >= float(threshold): 
		integer+=1

	return integer

integer = get_sequence_from_FASTA_file(fasta_file, "A")
print("We have found " + str(integer) + " proteins with those conditions")

fasta_file.close()

# Given a protein FASTA file (filename), print on the screen the protein identifier,
# the first N-aminoacids and the last M-aminoacids. The three fields must be separated
# by a tabulator, and one protein by line. 

fasta_file = open("ex.fasta")

def print_sequence_tails(filename, first_n=10, last_m=10):
	protein = ""
	for line in filename:
		if line.startswith(">"):
			if protein:
				print("{0} \t{1} \t{2}".format(header, protein[:first_n], protein[-last_m:]))
				protein = ""
			header = line.strip("\n")
		else:
			protein += line
			protein = protein.strip("\n")
	# for the last protein
	print("{0} \t{1} \t{2}".format(header, protein[:first_n], protein[-last_m:]))
		

print_sequence_tails(fasta_file)


fasta_file.close()

# ----------------------------------- PART 3 -------------------------------------------------

# Create Create a function that, given a protein FASTA file (multiline FASTA file) 
# and a “sub-sequences”file, calculates the proportion of proteins in the FASTA 
# file containing at least once each of the sub-sequences (exactly equal). 
# Save it in an output file  with the specified format, ordered by the proportion value
# (descending order) 

fasta_file = open("ex.fasta")
sub_seq_file = open ("sub_seq.txt")
output = open("out.txt", "w")

def calculate_aminoacid_frequencies(fasta_filename, subsequences_filename, output_filename):
	import pprint
	sequence = ""
	count = 0
	
	gen_list, prot_list, list_0, list_1, list_2 = [], [], [], [], []
	
	# PROTEIN LIST
	for line in fasta_filename:
		if line.startswith(">"):
			if sequence:
				prot_list.append(sequence)
				sequence = ""
		else:
			sequence += line.strip("\n")
	prot_list.append(sequence)

	# GETTING VARIABLES
	for sub_seq in subsequences_filename:
		for prot_seq in prot_list:
			if sub_seq.strip("\n") in prot_seq:
				count += 1
		gen_list.append([sub_seq.strip("\n"), count, count/len(prot_list)])
		count = 0
	
	from operator import itemgetter
	gen_list = sorted(gen_list, key=itemgetter(2))
	sort_list = list(reversed(gen_list))
	
	for i in sort_list:
		list_0.append(i[0])
		list_1.append(i[1])
		list_2.append(i[2])

	
	# PRINTING
	output_filename.write("#Number of proteins: {0} \n#Number of sequences: {1} \n#Subsequences_proportions:\n".format(len(prot_list), len(list_0)))
	for i in range(len(list_0)):
		output_filename.write(list_0[i].ljust(10) + str(list_1[i]).rjust(10) + "\t%.4f\n" %list_2[i])
		# !!!!!!!!!!!!!!!!!
		# according to the required format seen in the pdf, the print should be:
		# output.write(list_0[i] + "\t" + str(list_1[i]).rjust(10) + "\t%.4f\n" %list_2[i])
		# but it's not nice
		
calculate_aminoacid_frequencies(fasta_file, sub_seq_file, output)

fasta_file.close()
output.close()
sub_seq_file.close()
