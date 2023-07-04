import streamlit as st
from llm_processing import process_llm_request

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

    # Text input for Langchain LLM model request
    llm_request = st.text_input('Enter a specific request:')

    if st.button('Submit'):
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
        result = process_llm_request(llm_request, patient_name, patient_age, patient_bmi, patient_health, symptoms)
        st.write('LLM Result: ', result)

if __name__ == "__main__":
    main()
