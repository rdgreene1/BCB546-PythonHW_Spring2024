######################## BCB 546X: Python Assignment Details ########################

# ** Your Mission: Complete Python code in a Jupyter Notebook ** #

#-- Functions --#
## 1. Document Dr. X's function with comments and with markdown text in your Jupyter notebook.
## 2. Write a function that translates a string of nucleotides to amino acids based on Dr. X's pseudo-code suggestion.
## 3. Write an alternative translation function.
## 4. Write a function that calculates the molecular weight of each 3 amino acid sequence.
## 5. Write a function that computes the GC-content of each DNA sequence.

#-- In the MAIN part of the script --#
## 6. Add two new columns to the penguin DataFrame: (1) molecular weight and (2) GC content.
## 7. Call your functions from step 3 (or step 2) and step 4 and fill in the new columns in the DataFrame.
## 8. Plot a bar-chart of adult body mass per species. In your description of the graph, provide text that answers these questions: 
#       a. What is the smallest penguin species? 
#       b. What is the geographical range of this species?
## 9. Plot a graph that shows the molecular weight as a function of GC content. 
## 10. Write the entire DataFrame to a new CSV file that includes your new columns.
## 11. BONUS: What other visualizations, functions or tasks would you do with this dataset? Add something interesting for fun. (0.5 additional points if your total score is < 15).

#-- Additional Instructions (points will be deducted if these instructions are not heeded) --#
## ** Do all of this in a Jupyter notebook and push it to a GitHub repository.
## ** Your repository should not contain any files other than those associated with this assignment. 
## ** Read all comments carefully and answer the questions by including information in your Jupyter notebook.
## ** Document all of your code (and Dr. X's code) very thoroughly so that it is clear what you did.
## ** Be sure to cite (by providing URLs or other appropriate citations) information appropriately in your documented notebook.
## ** Commit and push your completed work in the Jupyter notebook to your repository.
## ** Submit the URL to your git repository via Canvas by the end of the day on May 6, 2022.

#-- Disclaimer --#
## Not all of these tasks have been covered in class and you will have to use online resources to find out how to do some of these tasks.


######################## Python Translate Script ########################

## Here's the start of our Python script. Thanks for completing it for me! - Dr. X
## IMPORTANT: install BioPython so that this will work

from Bio import SeqIO
from Bio.Data import CodonTable
import pandas as pd
from Bio.Seq import Seq
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.SeqUtils import GC
import seaborn as sns
import matplotlib.pyplot as plt

#%%%%%%%%%%%%%%%#
### FUNCTIONS ###
#%%%%%%%%%%%%%%%#

## 1 ##
## Dr. X: this gets sequences 
## Please properly document this function in the Jupyter notebook 
## Your descriptions of all functions should contain information about what the function does,
## as well as information about the return types and arguments.

def get_sequences_from_file(fasta_fn):
    """ 
    Extract species name and corresponding DNA sequence from a .fasta file.
        
    Parameters: 
        fasta_fn (.fasta): .fasta file containing DNA sequences to be analyzed.
            
    Returns: 
        sequence_data_dict: Dictionary with species names and DNA sequences stored as a Seq object.
    """
    
    sequence_data_dict = {} # empty dictionary that will be filled with the output of the loop
    for record in SeqIO.parse(fasta_fn, "fasta"): # Loop through the .fasta file and parse out Seq records
        description = record.description.split() # Split the description into multiple arguments
        species_name = description[1] + " " + description[2] # Generate the species name by adding the different arguments
        sequence_data_dict[species_name] = record.seq 
    return(sequence_data_dict) # Return species name and Seq object to the empty dictionary

seq_file = get_sequences_from_file('penguins_cytb.fasta')
print(seq_file)

## 2 ##
####### YOUR STRING-TRANSLATE FUNCTION ########
## Write a function that translates sequences
## All sequences start at codon position 1
## Complete a function that translates using a loop over the string of nucleotides
## Here is  some pseudo-code and suggestions
## feel free to change the function and variable names
# def translate_function(string_nucleotides): 
#     mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"] # this should work using BioPython (be sure to check what this returns)
#     for-loop through every 3rd position in string_nucleotides to get the codon using range subsets
#         # IMPORTANT: if the sequence has a stop codon at the end, you should leave it off
#         # this is how you can retrieve the amino acid: mito_table.forward_table[codon]
#         add the aa to aa_seq_string
#     return(aa_seq_string)

def translate_function(seq_dict):
    """
    Translate a dictionary of Seq objects into amino acid strings using the Vertebrate Mitochondrial codon table.
    
    Parameters:
        seq_dict: Dictionary containing species name as key and DNA sequences as Seq objects as the items.
        
    Return:
        translated_aa_dict: Dictionary containing the translated amino acid strings for each DNA sequence.
    """
    mito_table = CodonTable.unambiguous_dna_by_name["Vertebrate Mitochondrial"] # Load Vertebrate Mitochondrial codon table
    translated_aa_dict = {} # Initiate empty dictionary
    for species, seq_object in seq_dict.items(): # Loop through every sequence in the seq_dict
        translated_aa_list = [] # Initiate empty list
        for i in range(0, len(seq_object), 3): # Loop through the sequence and extract every three nucleotides
            codon = seq_object[i:i+3] # Define codon as every three nucleotides from sequence
            if codon == "TAA" or codon == "TAG" or codon == "AGA" or codon == "AGG": # Identify the stop codons
                break # Break loop if stop codon is encountered
            else:
                amino_acid = mito_table.forward_table[codon] # All other codons will be translated based on mito_table and defined as amino acids
            translated_aa_list.append(str(amino_acid)) # Store amino acids in empty list
        translated_aa_seq = ''.join(translated_aa_list) # Create a string of amino acids by joining amino acids from the list
        translated_aa_dict[species] = translated_aa_seq # Store amino acid strings in the empty dictionary
        
    return translated_aa_dict

translated_seq = translate_function(seq_file)
print(translated_seq)

## 3 ##
####### YOUR ALTERNATIVE FUNCTION ########
## Is there a better way to write the translation function? (Hint: yes there is.) 
## Perhaps using available BioPython library utilities?
## Please also write this function.

def new_translate_function(seq_dict, table=2):
    """
    A simpler translation function to translate a dictionary of Seq objects into amino acid strings using the Vertebrate Mitochondrial codon table.
    
    Parameters:
        seq_dict: Dictionary containing species name as key and DNA sequences as Seq objects as the items.
        table: Codon table identifier for NCBI; 1 is for standard codon table, 2 is for Vertebrate Mitochondrial codon table
        
    Returns:
        translated_seq_dict: Dictionary containing the translated amino acid strings for each DNA sequence.
    """
    translated_seq_dict = {}
    for species, seq_object in seq_dict.items():
        aa_sequence = seq_object.translate(table=table, to_stop = True) # Translates the Seq object from the file using the indicated codon table and stops translation at stop codon
        translated_seq_dict[species] = str(aa_sequence) # Stores amino acids as a string in a new dictionary

    return translated_seq_dict

new_translated_seq = new_translate_function(seq_file)
print(new_translated_seq)

## 4 ##
####### YOUR COUNT AA ANALYSIS FUNCTION ########
## Write a function that calculates the molecular weight of each amino acid sequence.
## For this, you can use some BioPython functions. I think you can use the ProtParam module.
## For more info, check this out: http://biopython.org/wiki/ProtParam
## So you should import the following before defining your function:
# def compute_molecular_weight(aa_seq):
#     # I think the ProtParam functions may require aa_seq to be a string.
#     # It may not work if the amino acid sequence has stop codons.
#     run the ProteinAnalysis() function on aa_seq
    #  return the molecular weight

def molecular_weight(translated_sequence):
    """
    Calculates the molecular weight of an amino acid string within a dictionary.
    
    Parameters:
        translated_sequence (dict): Dictionary of species names as key and string of amino acid sequences as items.
        
    Returns:
        molecular_weight_dict (dict): Dictionary that contains the species names and molecular weight of the respective amino acid string in Daltons.
    
    """
    molecular_weight_dict = {} # Empty dictionary to store molecular weights
    for species, aa_seq in translated_sequence.items(): # Loop through input dictionary selecting the amino acid strings
        prot_analysis = ProteinAnalysis(aa_seq) # Load the amino acid sequence into ProteinAnalysis module
        MW = "%0.2f" % prot_analysis.molecular_weight() # Calculate molecular weight for each amino acid string up to two decimal points
        molecular_weight_dict[species] = MW # Add molecular weights to the empty dictionary
        
    return molecular_weight_dict

protein_MW = molecular_weight(new_translated_seq)
print(protein_MW)

## 5 ##
####### YOUR GC CONTENT ANALYSIS FUNCTION ########
## Write a function that calculates the GC-content (proportion of "G" and "C") of each DNA sequence and returns this value.

def GC_content(seq_dict):
    """
    Calculate the GC content for DNA sequences.
    
    Parameters:
        seq_dict (dict): Dictionary of species names and DNA sequences as Seq objects.
        
    Return:
        GC_content_dict (dict): Dictionary of species names and GC percentage of DNA sequence.
    
    """
    GC_content_dict = {} # Empty dictionary to store species names and GC content
    for species, dna_seq in seq_dict.items(): # Loop through the input dictionary extracting the DNA sequences
        GC_cont = "%0.2f" % GC(dna_seq) # Use the GC module to calculate GC percent of DNA sequences
        GC_content_dict[species] = GC_cont # Store GC concent and species names in empty dictionary
        
    return GC_content_dict

seq_GC_cont = GC_content(seq_file)
print(seq_GC_cont)

#%%%%%%%%%%%%%%#
###   MAIN   ###
#%%%%%%%%%%%%%%#

cytb_seqs = get_sequences_from_file("penguins_cytb.fasta") 

penguins_df = pd.read_csv("penguins_mass.csv") # Includes only data for body mass 
species_list = list(penguins_df.species)

## 6 ## 
## Add two new columns to the penguin DataFrame: (1) molecular weight and (2) GC content.
## Set the value to 'NaN' to indicate that these cells are currently empty.

penguins_df['MW'] = 'NaN'
penguins_df['GC'] = 'NaN'
# Add new columns 'MW' and 'GC' with 'NaN' as the values
penguins_df.head()

## 7 ##
## Write a for-loop that translates each sequence and also gets molecular weight and computes the GC content
## of each translated sequence and adds those data to DataFrame
# for key, value in cytb_seqs.items():
#     aa_seq = nuc2aa_translate_function(value) # whichever function you prefer of #2 or #3
#     get the molecular weight of aa_seq
#     get the GC content of the DNA sequence
#     fill in empty cells in DF that you created above

penguins_df['MW'] = penguins_df['species'].map(molecular_weight(new_translated_seq))
# Add the molcular weight values for each species by mapping the species name and adding molecular weight to the MW column
penguins_df['GC'] = penguins_df['species'].map(GC_content(cytb_seqs))
# Add the GC content values for each species by mapping the species name and adding GC content to the GC column
penguins_df.head()

## 8 ##
## Plot a bar-chart of the mass with the x-axes labeled with species names.
## *Q1* What is the smallest penguin species? 
## *Q2* What is the geographical range of this species?

plt.bar(penguins_df['species'], penguins_df['mass']) # Create a bar graph using species and mass columns
plt.xlabel('species') # Label the x-axis as species
plt.ylabel('mass') # Label the y-axis as mass
plt.title('Mass by Species') # Create a graph title
plt.xticks(rotation=90) # Rotate the x-axis labels to be able to read them

# The smallest penguin species is the Eudyptula minor species.
# It looks like this species geographical range is solely in New Zealand. 

## 9 ##
## Plot a visualization of the molecular weight (y-axis) as a function of GC-content (x-axis).

scatter = sns.scatterplot(data=penguins_df, x='GC', y='MW', hue='species') 
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xlabel('GC')
plt.ylabel('MW')
plt.title('Relationship Between Molecular Weight and GC-content')
plt.xticks(rotation=90)

## 10 ##
## Save the new DataFrame to a file called "penguins_mass_cytb.csv"

penguins_df.to_csv('penguins_mass_cytb.csv', index=False)

## 11 - BONUS ##
## What else can we do with this dataset in Python? 
## Add functions or anything that might be interesting and fun. (optional)

def secondary_structure(translated_sequence):
    """
    Determine the secondary structure proportion of an amino acid sequence.
    
    Parameters:
        translated_sequence (dict): Dictionary that contains a species name and a corresponding amino acid sequence.
        
    Return:
        secondary_structure_dict (dict): Dictionary that contains the species name and three percentage values corresponding to proportion of the sequence in a helix, loop, and sheet, respectively.
    """
    secondary_structure_dict = {} # Empty dictionary to store species name and secondary structure percentages
    for species, aa_seq in translated_sequence.items(): # Loop through the dictionary pulling out the amino acid sequences
        prot_analysis = ProteinAnalysis(aa_seq) # Apply the ProteinAnalysis module to the amino acid sequence
        sec_struc = prot_analysis.secondary_structure_fraction() # Calculate the secondary structure percentages for the amino acid sequence
        secondary_structure_dict[species] = sec_struc # Store the percentages and species names in the empty dictionary
        
    return secondary_structure_dict

protein_sec_struc = secondary_structure(new_translated_seq)
print(protein_sec_struc)

protein_struc_df = pd.DataFrame.from_dict(protein_sec_struc, orient='index').reset_index()
protein_struc_df.columns = ['species'] + ['Helix'] + ['Loop'] + ['Sheet']
protein_struc_df.head()

penguins_merged_df = pd.merge(penguins_df, protein_struc_df[['species', 'Helix', 'Loop', 'Sheet']], on='species', how='left')
print(penguins_merged_df)

sns.scatterplot(data=penguins_merged_df, x='Helix', y='mass')
plt.xticks(rotation=90)
plt.title('Relationship Between Mass and Helix-percentage')

sns.scatterplot(data=penguins_merged_df, x='Loop', y='mass')
plt.xticks(rotation=90)
plt.title('Relationship Between Mass and Loop-percentage')

sns.scatterplot(data=penguins_merged_df, x='Sheet', y='mass')
plt.xticks(rotation=90)
plt.title('Relationship Between Mass and Sheet-percentage')

def aa_count(translated_sequence):
    """
    Calculate the amount of each amino acid for an amino acid string.
    
    Parameters:
        translated_sequence (dict): Dictionary that contains species name and a corresponding amino acid string.
        
    Return:
        aa_count_dict (dict): Dictionary that will have the species names and each amino acid single letter code with their respective count.
    
    """
    aa_count_dict = {} # Empty dictionary to store the species name and amino acid counts
    for species, aa_seq in translated_sequence.items(): # Loop through the input dictionary and pull out all the amino acid strings
        prot_analysis = ProteinAnalysis(aa_seq) # Apply the ProteinAnalysis module to the amino acid strings
        aa_count = prot_analysis.count_amino_acids() # Count each amino acid in the amino acid string
        aa_count_dict[species] = aa_count # Add the amino acid counts to the empty dictionary with their species name
        
    return aa_count_dict

protein_aa_count = aa_count(new_translated_seq)
print(protein_aa_count)

aa_count_df = pd.DataFrame.from_dict(protein_aa_count, orient='index').reset_index()
aa_count_df.columns = ['species'] + aa_count_df.columns[1:].tolist()
aa_count_df.head()

final_penguins_df = pd.merge(penguins_df, aa_count_df, on='species', how='left')
final_penguins_df.head()
