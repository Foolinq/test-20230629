import streamlit as st
import biopython as Bio
from llm_processing import process_llm_request
from Bio import Entrez, SeqIO
from collections import Counter

def count_codons(gene_ids, email):
    Entrez.email = email  # Always tell NCBI who you are
    codon_counts = {}

    for gene_id in gene_ids:
        handle = Entrez.efetch(db="nucleotide", id=gene_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        sequence = str(record.seq)
        codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]  # split sequence into codons

        # count codons
        codon_count = Counter(codons)
        codon_counts[gene_id] = codon_count

    return codon_counts

def main():
    st.title('Ovarian Cancer Diagnosis Interface')

    st.write('Please enter the patient information below:')

    # Text input for patient's name
    patient_name = st.text_input('Patient Name')

    # Slider for patient's age
    patient_age = st.slider('Patient Age', 0, 100, 50)

    # Slider for patient's BMI
    patient_bmi = st.slider('Patient BMI', 10.0, 50.0, 25.0)

    # Text input for self-described state of health
    patient_health = st.text_input('Self-described State of Health')

    # File uploader for patient's DNA sequence
    dna_file = st.file_uploader('Upload DNA Sequence', type=['txt', 'fasta'])

    # Text input for symptoms
    symptoms = st.text_input('Enter Symptoms')

    # Text input for gene IDs
    gene_ids = st.text_input('Enter Gene IDs (comma-separated)')

    # Text input for doctor's question
    doctor_question = st.text_input('Enter your question for the Langchain LLM model')

    if st.button('Submit Patient Information'):
        st.write('Patient Information:')
        st.write('Name: ', patient_name)
        st.write('Age: ', patient_age)
        st.write('BMI: ', patient_bmi)
        st.write('Self-described State of Health: ', patient_health)
        st.write('Symptoms: ', symptoms)

        if dna_file is not None:
            # Read the file and store the DNA sequence
            dna_sequence = dna_file.getvalue().decode()
            st.write('DNA Sequence: ', dna_sequence)

        # Process the Langchain LLM model request
        result = process_llm_request(doctor_question, patient_name, patient_age, patient_bmi, patient_health, symptoms)
        st.write('LLM Result: ', result)

    if st.button('Submit Gene IDs'):
        # Split the gene IDs string into a list
        gene_ids_list = [gene_id.strip() for gene_id in gene_ids.split(',')]
        # Count codons for each gene
        codon_counts = count_codons(gene_ids_list, "your_email@example.com")  # replace with your email
        st.write('Codon Counts: ', codon_counts)

if __name__ == "__main__":
    main()
