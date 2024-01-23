###NutriVi Categorizing the nutrition of Viruses###
import pandas as pd
import numpy as np
import sys,argparse, re, os




base_score= {'A': 5,'T':2,'G':5,'C':3}

#count the number of each base for a seqrecord object
def count_bases(fasta_file):
    # Define a dictionary to store the base counts
    base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    
    with open(fasta_file, 'r') as f:
        # Skip the first line, which contains the header
        next(f)
        
        # Loop through the rest of the file, counting bases
        for line in f:
            line = line.rstrip()
            for base in line:
                if base in base_counts:
                    base_counts[base] += 1
    
    # Print the results
    print("A: ", base_counts['A'])
    print("C: ", base_counts['C'])
    print("G: ", base_counts['G'])
    print("T: ", base_counts['T'])


#Get nitrogen content of a sequence
def get_nitrogen_content(fasta_file):
    # Define a dictionary to store the base counts
    base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    
    with open(fasta_file, 'r') as f:
        # Skip the first line, which contains the header
        next(f)
        
        # Loop through the rest of the file, counting bases
        for line in f:
            line = line.rstrip()
            for base in line:
                if base in base_counts:
                    base_counts[base] += 1
    #get nitrogen content by multiplying the base score by the number of each base
    nitrogen_content=2*(base_counts['A']*base_score['A']+base_counts['C']*base_score['C']+base_counts['G']*base_score['G']+base_counts['T']*base_score['T'])
    return nitrogen_content

#Get phosphate content of a sequence
def get_phosphate_content(fasta_file):
    # Define a dictionary to store the base counts
    base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    base_score= {'A': 2,'T':2,'G':2,'C':2}
    with open(fasta_file, 'r') as f:
        # Skip the first line, which contains the header
        next(f)
        
        # Loop through the rest of the file, counting bases
        for line in f:
            line = line.rstrip()
            for base in line:
                if base in base_counts:
                    base_counts[base] += 1
    #get phosphate content by multiplying the base score by the number of each base
    phosphate_content=2*(base_counts['A']*base_score['A']+base_counts['C']*base_score['C']+base_counts['G']*base_score['G']+base_counts['T']*base_score['T'])
    return phosphate_content

#Get carbon content
def get_carbon_content(fasta_file):
    # Define a dictionary to store the base counts
    base_counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
    base_score= {'A': 8,'T':9,'G':9,'C':8}
    with open(fasta_file, 'r') as f:
       #skip if the first line has a > in it
        if '>' in f.readline():
            next(f)
        # Loop through the rest of the file, counting bases
        for line in f:
            line = line.rstrip()
            for base in line:
                if base in base_counts:
                    base_counts[base] += 1
    #get carbon content by multiplying the base score by the number of each base
    carbon_content=2*(base_counts['A']*base_score['A']+base_counts['C']*base_score['C']+base_counts['G']*base_score['G']+base_counts['T']*base_score['T'])
    return carbon_content

#Ratio of nitrogen to carbon
def get_nitrogen_carbon_ratio(fasta_file):
    nitrogen_content=get_nitrogen_content(fasta_file)
    carbon_content=get_carbon_content(fasta_file)
    nitrogen_carbon_ratio=nitrogen_content/carbon_content
    return nitrogen_carbon_ratio



#Protein nitrogen content
def get_protein_nitrogen_content(protein_file):
    #dictionary for amino acid scores
    amino_acids = {
    "A": 0,
    "C": 0,
    "D": 0,
    "E": 0,
    "F": 0,
    "G": 0,
    "H": 2,
    "I": 0,
    "K": 1,
    "L": 0,
    "M": 0,
    "N": 1,
    "P": 1,
    "Q": 1,
    "R": 3,
    "S": 0,
    "U": 0,
    "T": 0,
    "V": 0,
    "W": 1,
    "Y": 0}
    # Define a dictionary to store the base counts
    amino_acid_count= {
    "A": 0,
    "C": 0,
    "D": 0,
    "E": 0,
    "F": 0,
    "G": 0,
    "H": 0,
    "I": 0,
    "K": 0,
    "L": 0,
    "M": 0,
    "N": 0,
    "P": 0,
    "Q": 0,
    "R": 0,
    "S": 0,
    "U": 0,
    "T": 0,
    "V": 0,
    "W": 0,
    "Y": 0}
    with open(protein_file, 'r') as f:
        # Skip the first line, which contains the header
        if '>' in f.readline():
            next(f) 
        # Loop through the rest of the file, counting bases
        for line in f:
            line = line.rstrip()
            for amino_acid in line:
                if amino_acid in amino_acids:
                    amino_acid_count[amino_acid] += 1
    #get protein nitrogen content by multiplying the amino acid score by the number of each amino acid
    protein_nitrogen_content=(amino_acid_count['A']*amino_acids['A']+amino_acid_count['C']*amino_acids['C']+amino_acid_count['D']*amino_acids['D']+amino_acid_count['E']*amino_acids['E']+amino_acid_count['F']*amino_acids['F']+amino_acid_count['G']*amino_acids['G']+amino_acid_count['H']*amino_acids['H']+amino_acid_count['I']*amino_acids['I']+amino_acid_count['K']*amino_acids['K']+amino_acid_count['L']*amino_acids['L']+amino_acid_count['M']*amino_acids['M']+amino_acid_count['N']*amino_acids['N']+amino_acid_count['P']*amino_acids['P']+amino_acid_count['Q']*amino_acids['Q']+amino_acid_count['R']*amino_acids['R']+amino_acid_count['S']*amino_acids['S']+amino_acid_count['U']*amino_acids['U']+amino_acid_count['T']*amino_acids['T']+amino_acid_count['V']*amino_acids['V']+amino_acid_count['W']*amino_acids['W']+amino_acid_count['Y']*amino_acids['Y'])
    return protein_nitrogen_content

def get_protein_sulfur_content(protein_file):
    #dictionary for amino acid scores
    amino_acids = {
    "A": 0,
    "C": 1,
    "D": 0,
    "E": 0,
    "F": 0,
    "G": 0,
    "H": 0,
    "I": 0,
    "K": 0,
    "L": 0,
    "M": 1,
    "N": 0,
    "P": 0,
    "Q": 0,
    "R": 0,
    "S": 0,
    "U": 0,
    "T": 0,
    "V": 0,
    "W": 0,
    "Y": 0}
    # Define a dictionary to store the base counts
    amino_acid_count= {
    "A": 0,
    "C": 0,
    "D": 0,
    "E": 0,
    "F": 0,
    "G": 0,
    "H": 0,
    "I": 0,
    "K": 0,
    "L": 0,
    "M": 0,
    "N": 0,
    "P": 0,
    "Q": 0,
    "R": 0,
    "S": 0,
    "U": 0,
    "T": 0,
    "V": 0,
    "W": 0,
    "Y": 0}
    with open(protein_file, 'r') as f:
        # Skip the first line, which contains the header
        if '>' in f.readline():
            next(f) 
        # Loop through the rest of the file, counting bases
        for line in f:
            line = line.rstrip()
            for amino_acid in line:
                if amino_acid in amino_acids:
                    amino_acid_count[amino_acid] += 1
    #get protein nitrogen content by multiplying the amino acid score by the number of each amino acid
    protein_nitrogen_content=(amino_acid_count['A']*amino_acids['A']+amino_acid_count['C']*amino_acids['C']+amino_acid_count['D']*amino_acids['D']+amino_acid_count['E']*amino_acids['E']+amino_acid_count['F']*amino_acids['F']+amino_acid_count['G']*amino_acids['G']+amino_acid_count['H']*amino_acids['H']+amino_acid_count['I']*amino_acids['I']+amino_acid_count['K']*amino_acids['K']+amino_acid_count['L']*amino_acids['L']+amino_acid_count['M']*amino_acids['M']+amino_acid_count['N']*amino_acids['N']+amino_acid_count['P']*amino_acids['P']+amino_acid_count['Q']*amino_acids['Q']+amino_acid_count['R']*amino_acids['R']+amino_acid_count['S']*amino_acids['S']+amino_acid_count['U']*amino_acids['U']+amino_acid_count['T']*amino_acids['T']+amino_acid_count['V']*amino_acids['V']+amino_acid_count['W']*amino_acids['W']+amino_acid_count['Y']*amino_acids['Y'])
    return protein_nitrogen_content

#Get total protein
def get_total_protein(protein_file):
    amino_acid_count= {
    "A": 0,
    "C": 0,
    "D": 0,
    "E": 0,
    "F": 0,
    "G": 0,
    "H": 0,
    "I": 0,
    "K": 0,
    "L": 0,
    "M": 0,
    "N": 0,
    "P": 0,
    "Q": 0,
    "R": 0,
    "S": 0,
    "U": 0,
    "T": 0,
    "V": 0,
    "W": 0,
    "Y": 0}
    with open(protein_file, 'r') as f:
        # Skip the first line, which contains the header
        if '>' in f.readline():
            next(f) 
        # Loop through the rest of the file, counting bases
        for line in f:
            line = line.rstrip()
            for amino_acid in line:
                if amino_acid in amino_acid_count:
                    amino_acid_count[amino_acid] += 1
    #get protein nitrogen content by multiplying the amino acid score by the number of each amino acid
    protein_total= sum(amino_acid_count.values())
    return protein_total

def get_protein_carbon_content(protein_file):
    #dictionary for amino acid scores
    amino_acids = {
    "A": 1,
    "C": 1,
    "D": 2,
    "E": 3,
    "F": 7,
    "G": 0,
    "H": 4,
    "I": 4,
    "K": 4,
    "L": 4,
    "M": 3,
    "N": 2,
    "P": 3,
    "Q": 3,
    "R": 4,
    "S": 1,
    "U": 1,
    "T": 2,
    "V": 3,
    "W": 8,
    "Y": 7}
    # Define a dictionary to store the base counts
    amino_acid_count= {
    "A": 0,
    "C": 0,
    "D": 0,
    "E": 0,
    "F": 0,
    "G": 0,
    "H": 0,
    "I": 0,
    "K": 0,
    "L": 0,
    "M": 0,
    "N": 0,
    "P": 0,
    "Q": 0,
    "R": 0,
    "S": 0,
    "U": 0,
    "T": 0,
    "V": 0,
    "W": 0,
    "Y": 0}
    with open(protein_file, 'r') as f:
        # Skip the first line, which contains the header
        if '>' in f.readline():
            next(f) 
        # Loop through the rest of the file, counting bases
        for line in f:
            line = line.rstrip()
            for amino_acid in line:
                if amino_acid in amino_acids:
                    amino_acid_count[amino_acid] += 1
    #get protein nitrogen content by multiplying the amino acid score by the number of each amino acid
    protein_nitrogen_content=(amino_acid_count['A']*amino_acids['A']+amino_acid_count['C']*amino_acids['C']+amino_acid_count['D']*amino_acids['D']+amino_acid_count['E']*amino_acids['E']+amino_acid_count['F']*amino_acids['F']+amino_acid_count['G']*amino_acids['G']+amino_acid_count['H']*amino_acids['H']+amino_acid_count['I']*amino_acids['I']+amino_acid_count['K']*amino_acids['K']+amino_acid_count['L']*amino_acids['L']+amino_acid_count['M']*amino_acids['M']+amino_acid_count['N']*amino_acids['N']+amino_acid_count['P']*amino_acids['P']+amino_acid_count['Q']*amino_acids['Q']+amino_acid_count['R']*amino_acids['R']+amino_acid_count['S']*amino_acids['S']+amino_acid_count['U']*amino_acids['U']+amino_acid_count['T']*amino_acids['T']+amino_acid_count['V']*amino_acids['V']+amino_acid_count['W']*amino_acids['W']+amino_acid_count['Y']*amino_acids['Y'])
    return protein_nitrogen_content



####The main program
parser = argparse.ArgumentParser(description='Calculate the nitrogen content of a genome')
parser.add_argument('-i', '--input', help='Input DNA file folder', required=False)
parser.add_argument('-p', '--protein', help='Input protein file folder', required=False)

args = parser.parse_args()

input_DNA_file = args.input
protein_file = args.protein
#Average N per side chain
def Average_Nitrogen(protein_file): 
    avg=get_protein_nitrogen_content(protein_file)/get_total_protein(protein_file)
    return avg
def n_to_carbon(protein_file):
    avg=get_protein_nitrogen_content(protein_file)/get_protein_carbon_content(protein_file)
    return avg

def Average_sulfur(protein_file):
    avg=get_protein_sulfur_content(protein_file)/get_total_protein(protein_file)
    return avg
#check if -i is used
if input_DNA_file:
    #create a dataframe
    df = pd.DataFrame(columns=['Genome', 'Nitrogen_Carbon_Ratio'])
    for file in os.listdir(input_DNA_file):
        fasta_file = os.path.join(input_DNA_file, file)
        nitrogen_carbon_ratio=get_nitrogen_carbon_ratio(fasta_file)
        #append this to a dataframe
        df = df.append({'Genome': file, 'Nitrogen_Carbon_Ratio': nitrogen_carbon_ratio}, ignore_index=True)
    #save the dataframe as a csv
    df.to_csv('Nitrogen_Carbon_Ratio.csv', index=False)

#check if -p is used
if protein_file:
    #create a dataframe
    df = pd.DataFrame(columns=['Genome', 'AVG_N_sidechain','N:C Ratio','AVG_Sulfur_sidechain'])
    for file in os.listdir(protein_file):
        protein_file_new = os.path.join(protein_file, file)
        averagen=Average_Nitrogen(protein_file_new)
        n_to_c= n_to_carbon(protein_file_new)
        avg_sulfur=Average_sulfur(protein_file_new)
        #append this to a dataframe
        df = df.append({'Genome': file, 'AVG_N_sidechain': averagen,'N:C Ratio':n_to_c,'AVG_Sulfur_sidechain':avg_sulfur}, ignore_index=True)
    #save the dataframe as a csv
    df.to_csv('Average_N_per_sidechain.csv', index=False)
