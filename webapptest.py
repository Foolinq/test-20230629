# This program has two elements: first, a LLM-powered diagnosis tool for doctors.
# second, a tool to separate coding sequences of genes into individual codons and amino acids
# The goal is to eventually have a tool which doctors will be able to use to detect early cancer (namely ovarian cancer) using DNA samples from patients


# TODO:
# Change column names in gene tables to gene names
# Find way for gene ID input to actually be gene name input to make more user-friendly
# Find possible symptom database, then use that database to make a scrolldown selection menu instead of text box
# Refine prompt to aim towards cancer detection


import os
import time
from Bio import Entrez, SeqIO
from collections import Counter
from Bio.Seq import Seq
import pandas as pd
import streamlit as st
from llm_processing import process_llm_request

# dictionary to map single-letter codes to full names of amino acids
amino_acid_names = {
    'A': 'Alanine', 'C': 'Cysteine', 'D': 'Aspartic acid', 'E': 'Glutamic acid',
    'F': 'Phenylalanine', 'G': 'Glycine', 'H': 'Histidine', 'I': 'Isoleucine',
    'K': 'Lysine', 'L': 'Leucine', 'M': 'Methionine', 'N': 'Asparagine',
    'P': 'Proline', 'Q': 'Glutamine', 'R': 'Arginine', 'S': 'Serine',
    'T': 'Threonine', 'V': 'Valine', 'W': 'Tryptophan', 'Y': 'Tyrosine',
    '*': 'Stop codon'
}

import requests

def get_gene_id(gene_symbol):
    server = "https://rest.ensembl.org"
    ext = f"/xrefs/symbol/homo_sapiens/{gene_symbol}?content-type=application/json"
    r = requests.get(server + ext, headers={"Content-Type": "application/json"})

    if not r.ok:
        r.raise_for_status()
        return None

    decoded = r.json()
    try:
        return decoded[0]['id']
    except (IndexError, KeyError):
        return None


def count_codons(gene_symbols, debug_mode):
    server = "https://rest.ensembl.org"
    codon_counts = {}
    amino_acid_counts = {}

    for gene_symbol in gene_symbols:
        try:
            gene_id = get_gene_id(gene_symbol)

            ext = f"/overlap/id/{gene_id}?feature=cds;content-type=application/json"
            r = requests.get(server + ext, headers={"Content-Type": "application/json"})

            if not r.ok:
                r.raise_for_status()
                sys.exit()

            decoded = r.json()

            # initialize counters for this gene
            codon_count = Counter()
            amino_acid_count = Counter()

            # find the CDS feature
            cds_count = 0
            for feature in decoded:
                if feature['feature_type'] == "cds":
                    cds_count += 1
                    # extract the CDS sequence
                    sequence = feature['seq']
                    if debug_mode:
                        st.write(f"CDS {cds_count} for gene ID {gene_id}: {sequence}")  # output the CDS

                    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
                    
                    if debug_mode:
                        st.write(f"Codons for gene ID {gene_id}: {codons}")  # output the list of codons

                    # update codon count
                    codon_count.update(codons)

                    # translate codons to amino acids and update count
                    amino_acids = Seq(sequence).translate()
                    amino_acid_count.update(amino_acids)

            if debug_mode:
                st.write(f"Number of CDS for gene ID {gene_id}: {cds_count}")
                st.write(f"Codon counts for gene ID {gene_id}: {codon_count}")
                st.write(f"Amino acid counts for gene ID {gene_id}: {amino_acid_count}")

            # replace single-letter codes with full names
            amino_acid_count = {amino_acid_names.get(k, k): v for k, v in amino_acid_count.items()}

            # store counts for this gene
            codon_counts[gene_id] = codon_count
            amino_acid_counts[gene_id] = amino_acid_count

        except Exception as e:
            if debug_mode:
                st.write(f"Error processing gene symbol {gene_symbol}: {str(e)}")

        # Wait for 0.1 seconds to respect the rate limit of 10 requests per second
        time.sleep(0.1)

    return codon_counts, amino_acid_counts


def get_patient_info():
    """Get patient information from Streamlit widgets and return as a dictionary."""
    patient_info = {}

    patient_info['name'] = st.text_input('Patient Name')
    patient_info['age'] = st.slider('Patient Age', 0, 100, 50)
    patient_info['weight'] = st.slider('Patient Weight (kg)', 30.0, 150.0, 70.0)
    patient_info['height'] = st.slider('Patient Height (m)', 1.0, 2.5, 1.7)
    patient_info['bmi'] = patient_info['weight'] / patient_info['height']**2
    patient_info['sex'] = st.selectbox('Sex at Birth', ['Male', 'Female'])
    patient_info['activity'] = st.selectbox(
        'Recent Physical Activity Level',
        [
            "No physical activity",
            "Light physical activity (e.g., walks, light housework)",
            "Moderate physical activity (e.g., brisk walking, cycling)",
            "Intense physical activity (e.g., running, heavy lifting)",
            "Very intense physical activity (e.g., athlete in training)"
        ]
    )
    patient_info['symptoms'] = st.text_input('Enter Symptoms')
    patient_info['question'] = st.text_input('Enter your question or concern regarding the patient')
    patient_info['dna_file'] = st.file_uploader('Upload DNA Sequence', type=['txt', 'fasta'])

    return patient_info


def display_patient_info(patient_info):
    """Display the patient information."""
    st.write('Patient Information:')
    st.write('Name: ', patient_info['name'])
    st.write('Age: ', patient_info['age'])
    st.write('Weight (kg): ', patient_info['weight'])
    st.write('Height (m): ', patient_info['height'])
    st.write('BMI: ', patient_info['bmi'])
    st.write('Sex at Birth: ', patient_info['sex'])
    st.write('Recent Physical Activity Level: ', patient_info['activity'])
    st.write('Symptoms: ', patient_info['symptoms'])

    if patient_info['dna_file'] is not None:
        dna_sequence = patient_info['dna_file'].getvalue().decode()
        st.write('DNA Sequence: ', dna_sequence)


def process_patient_info(patient_info):
        """Process the patient information and display it."""
        display_patient_info(patient_info)

        # If the DNA file was uploaded, decode it and add it to the dictionary
        if patient_info['dna_file'] is not None:
            patient_info['dna_sequence'] = patient_info['dna_file'].getvalue().decode()

        # Process the Langchain LLM model request
        result = process_llm_request(
            patient_info['question'],
            patient_info['name'],
            patient_info['age'],
            patient_info['weight'],
            patient_info['height'],
            patient_info['bmi'],
            patient_info['sex'],
            patient_info['activity'],
            patient_info['symptoms']
        )
        st.write('LLM Result: ', result)


def get_gene_symbols():
    gene_symbols = st.text_input('Enter Gene Symbols (comma-separated)')
    debug_mode = st.checkbox('Activate Debug Mode')
    return gene_symbols, debug_mode

def process_gene_symbols(gene_symbols, debug_mode):
    gene_symbols_list = [gene_symbol.strip() for gene_symbol in gene_symbols.split(',')]
    codon_counts, amino_acid_counts = count_codons(gene_symbols_list, "parentrdavid@gmail.com", os.getenv('NCBI_API_KEY'), debug_mode)

    # Convert dictionaries to pandas DataFrames
    codon_df = pd.DataFrame.from_dict(codon_counts, orient='index').transpose()
    amino_acid_df = pd.DataFrame.from_dict(amino_acid_counts, orient='index').transpose()

    # Replace NA values with 0 and convert to integer
    codon_df = codon_df.fillna(0).astype(int)
    amino_acid_df = amino_acid_df.fillna(0).astype(int)

    # Add a total count column
    codon_df['Total'] = codon_df.sum(axis=1)
    amino_acid_df['Total'] = amino_acid_df.sum(axis=1)

    # Create two columns for side-by-side layout
    col1, col2 = st.columns(2)

    # Display tables side by side
    with col1:
        st.table(codon_df)
    with col2:
        st.table(amino_acid_df)


def clean_genes_in_excel(input_filepath, output_filepath, email, api_key, debug_mode):
    # Read the Excel file
    df = pd.read_excel(input_filepath)

    # Get a list of gene symbols from the column names
    gene_symbols = df.columns.tolist()

    # Initialize a list to hold gene symbols that have at least one CDS
    valid_genes = []

    # Iterate over the gene symbols
    for gene_symbol in gene_symbols:
        try:
            gene_id = get_gene_id(gene_symbol, email, api_key)

            # Fetch the gene record
            handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()

            # Check if there's a CDS feature in the record
            cds_found = False
            for feature in record.features:
                if feature.type == "CDS":
                    cds_found = True
                    break
            
            if cds_found:
                valid_genes.append(gene_symbol)
            elif debug_mode:
                st.write(f"No CDS feature found for gene symbol {gene_symbol}")

        except Exception as e:
            if debug_mode:
                st.write(f"Error processing gene symbol {gene_symbol}: {str(e)}")

        # Wait for 0.1 seconds to respect the rate limit of 10 requests per second
        time.sleep(0.1)

    # Keep only the columns for valid genes
    df_clean = df[valid_genes]

    # Write the cleaned DataFrame to the output Excel file
    df_clean.to_excel(output_filepath, index=False)


# Streamlit section functions, contain calls to all other functions


def patient_info_section():
    st.write('Please enter the patient information below:')
    patient_info = get_patient_info()
    if st.button('Submit Patient Information'):
        process_patient_info(patient_info)

def gene_analysis_section():
    st.markdown("---")
    st.markdown("# Gene Analysis")
    gene_symbols, debug_mode = get_gene_symbols()
    if st.button('Submit Gene Symbols'):
        process_gene_symbols(gene_symbols, debug_mode)

def excel_cleaning_section():
    st.markdown("---")
    st.markdown("# Clean Gene Columns in Excel")
    debug_mode = st.checkbox('Activate Debug Mode for Excel Cleaning')
    input_file = st.file_uploader('Upload Excel File', type=['xlsx'])
    if st.button('Clean Excel File') and input_file is not None:
        # Save the uploaded file to a temporary file
        input_filepath = 'input.xlsx'
        with open(input_filepath, 'wb') as f:
            f.write(input_file.getbuffer())
        output_filepath = 'output.xlsx'
        clean_genes_in_excel(input_filepath, output_filepath, "parentrdavid@gmail.com", os.getenv('NCBI_API_KEY'), debug_mode)
        # Offer the cleaned file for download
        st.markdown(f'[Download Cleaned Excel File]({output_filepath})')


# Main function below, only contains calls to section functions


def main():
    st.title('Ovarian Cancer Diagnosis Interface')
    patient_info_section()
    gene_analysis_section()
    excel_cleaning_section()

if __name__ == "__main__":
    main()
