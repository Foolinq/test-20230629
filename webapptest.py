# This program has two elements: first, a LLM-powered diagnosis tool for doctors.
# second, a tool to separate coding sequences of genes into individual codons and amino acids
# The goal is to eventually have a tool which doctors will be able to use to detect early cancer (namely ovarian cancer) using DNA samples from patients


# TODO:
# Change column names in gene tables
# Make codon tables side-by-side
# Find way for gene ID input to actually be gene name input to make more user-friendly
# Make those gene names their column names
# Add a total codon and amino acid count column
# Find possible symptom database, then use that database to make a scrolldown selection menu instead of text box
# Change BMI input to weight and height, and calculate BMI using those inputs
# Change "Enter question for the LLM" part to say something along the lines of "Question/concern regarding patient"
# Add a "sex at birth" input (scrolldown two options: male and female)
# Refine prompt to aim towards cancer detection





import streamlit as st
from llm_processing import process_llm_request
from Bio import Entrez, SeqIO
from collections import Counter
import pandas as pd
from Bio.Seq import Seq

# dictionary to map single-letter codes to full names of amino acids
amino_acid_names = {
    'A': 'Alanine', 'C': 'Cysteine', 'D': 'Aspartic acid', 'E': 'Glutamic acid',
    'F': 'Phenylalanine', 'G': 'Glycine', 'H': 'Histidine', 'I': 'Isoleucine',
    'K': 'Lysine', 'L': 'Leucine', 'M': 'Methionine', 'N': 'Asparagine',
    'P': 'Proline', 'Q': 'Glutamine', 'R': 'Arginine', 'S': 'Serine',
    'T': 'Threonine', 'V': 'Valine', 'W': 'Tryptophan', 'Y': 'Tyrosine',
    '*': 'Stop codon'
}

def count_codons(gene_ids, email):
    Entrez.email = email
    codon_counts = {}
    amino_acid_counts = {}

    for gene_id in gene_ids:
        try:
            handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="gb", retmode="text")
            record = SeqIO.read(handle, "genbank")
            handle.close()

            # initialize counters for this gene
            codon_count = Counter()
            amino_acid_count = Counter()

            # find the CDS feature
            for feature in record.features:
                if feature.type == "CDS":
                    # extract the CDS sequence
                    sequence = str(feature.extract(record.seq))

                    st.write(f"CDS for gene ID {gene_id}: {sequence}")  # output the CDS

                    codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]
            
                    st.write(f"Codons for gene ID {gene_id}: {codons}")  # output the list of codons

                    # update codon count
                    codon_count.update(codons)

                    # translate codons to amino acids and update count
                    amino_acids = Seq(sequence).translate()
                    amino_acid_count.update(amino_acids)

            # replace single-letter codes with full names
            amino_acid_count = {amino_acid_names.get(k, k): v for k, v in amino_acid_count.items()}

            # store counts for this gene
            codon_counts[gene_id] = codon_count
            amino_acid_counts[gene_id] = amino_acid_count

        except Exception as e:
            st.write(f"Error fetching data for gene ID {gene_id}: {e}")

    return codon_counts, amino_acid_counts



def main():
    st.title('Ovarian Cancer Diagnosis Interface')

    st.write('Please enter the patient information below:')

    # Text input for patient's name
    patient_name = st.text_input('Patient Name')

    # Slider for patient's age
    patient_age = st.slider('Patient Age', 0, 100, 50)

    # Slider for patient's weight and height
    patient_weight = st.slider('Patient Weight (kg)', 30.0, 150.0, 70.0)
    patient_height = st.slider('Patient Height (m)', 1.0, 2.5, 1.7)
    
    # Calculate BMI
    patient_bmi = patient_weight / patient_height**2

    # Select box for sex at birth
    patient_sex = st.selectbox('Sex at Birth', ['Male', 'Female'])

    # Select box for recent physical activity level
    patient_activity = st.selectbox(
        'Recent Physical Activity Level',
        [
            "No physical activity",
            "Light physical activity (e.g., walks, light housework)",
            "Moderate physical activity (e.g., brisk walking, cycling)",
            "Intense physical activity (e.g., running, heavy lifting)",
            "Very intense physical activity (e.g., athlete in training)"
        ]
    )

    # File uploader for patient's DNA sequence
    dna_file = st.file_uploader('Upload DNA Sequence', type=['txt', 'fasta'])

    # Text input for symptoms
    symptoms = st.text_input('Enter Symptoms')

    # Text input for doctor's question with updated label
    doctor_question = st.text_input('Enter your question or concern regarding the patient')

    if st.button('Submit Patient Information'):
        st.write('Patient Information:')
        st.write('Name: ', patient_name)
        st.write('Age: ', patient_age)
        st.write('Weight (kg): ', patient_weight)
        st.write('Height (m): ', patient_height)
        st.write('BMI: ', patient_bmi)
        st.write('Sex at Birth: ', patient_sex)
        st.write('Recent Physical Activity Level: ', patient_activity)
        st.write('Symptoms: ', symptoms)

        if dna_file is not None:
            # Read the file and store the DNA sequence
            dna_sequence = dna_file.getvalue().decode()
            st.write('DNA Sequence: ', dna_sequence)

        # Process the Langchain LLM model request
        result = process_llm_request(
            doctor_question, 
            patient_name, 
            patient_age, 
            patient_weight,
            patient_height,
            patient_bmi, 
            patient_sex, 
            patient_activity, 
            symptoms
        )
        st.write('LLM Result: ', result)

    # add a horizontal line and a title
    st.markdown("---")
    st.markdown("# Gene Analysis")

    # gene input and submit boxes
    gene_ids = st.text_input('Enter Gene IDs (comma-separated)')

    if st.button('Submit Gene IDs'):
        gene_ids_list = [gene_id.strip() for gene_id in gene_ids.split(',')]
        codon_counts, amino_acid_counts = count_codons(gene_ids_list, "parentrdavid@gmail.com")

        # Convert dictionaries to pandas DataFrames
        codon_df = pd.DataFrame.from_dict(codon_counts, orient='index').transpose()
        amino_acid_df = pd.DataFrame.from_dict(amino_acid_counts, orient='index').transpose()

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



if __name__ == "__main__":
    main()
