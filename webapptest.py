import pandas as pd
import streamlit as st
import requests
import concurrent.futures
import time

from sqlalchemy import create_engine, Column, String, ForeignKey
from sqlalchemy.orm import sessionmaker, relationship
from sqlalchemy.ext.declarative import declarative_base

# Define the SQLAlchemy models
Base = declarative_base()

class Gene(Base):
    __tablename__ = 'genes'

    name = Column(String, primary_key=True)
    sequence = Column(String)

# Set up the database
engine = create_engine('postgresql://jqczrvezvwlzer:5980122424475b4fcdc2e539dbb0667d254452a21cd3cc633f9a64a8796da2df@ec2-3-212-70-5.compute-1.amazonaws.com:5432/d202roftknash2')
Session = sessionmaker(bind=engine)

def setup_tables():
    Base.metadata.create_all(engine)

# Global variable for debugging
debug = False

def add_gene_and_sequence(gene_name, sequence, gene_number, total_genes):
    # Start a new session
    session = Session()

    try:
        # Check if the gene already exists
        gene = session.query(Gene).filter_by(name=gene_name).first()
        if gene is not None:
            # If the gene already exists, skip it
            if debug:
                st.write(f"Gene {gene_number} of {total_genes}: {gene_name} - Not added, already exists in the database.")
            return

        # Create a new gene with its sequence
        gene = Gene(name=gene_name, sequence=sequence)

        # Add the gene to the session
        session.add(gene)

        # Commit the session to save the changes to the database
        session.commit()

        if debug:
            st.write(f"Gene {gene_number} of {total_genes}: {gene_name} - Successfully added to the database.")
    except:
        # If an error occurs, roll back the session
        session.rollback()
        raise
    finally:
        # Always close the session when you're done with it
        session.close()





# Function to convert HGNC symbol to Ensembl ID
def symbol_to_id(gene_symbol):
    server = "https://rest.ensembl.org"
    ext = f"/xrefs/symbol/homo_sapiens/{gene_symbol}?content-type=application/json"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
        r.raise_for_status()
    decoded = r.json()
    return (gene_symbol, decoded[0]['id']) if decoded else (gene_symbol, None)

# Function to fetch transcript IDs for a given gene
def fetch_transcripts(ensembl_ids):
    server = "https://rest.ensembl.org"
    ext = "/lookup/id/?content-type=application/json;expand=1"
    data = {"ids": ensembl_ids}
    r = requests.post(server+ext, headers={ "Content-Type" : "application/json"}, json=data)
    if not r.ok:
        r.raise_for_status()
    decoded = r.json()
    return {ensembl_id: [transcript['id'] for transcript in result['Transcript']] for ensembl_id, result in decoded.items() if 'Transcript' in result}

# Function to fetch CDS for a given transcript
def fetch_cds(transcript_id):
    server = "https://rest.ensembl.org"
    ext = f"/sequence/id/{transcript_id}?content-type=text/x-fasta;type=cds"
    r = requests.get(server+ext, headers={ "Content-Type" : "text/x-fasta"})
    if not r.ok:
        if r.status_code == 429 and 'Retry-After' in r.headers:
            retry_after = int(r.headers['Retry-After'])
            time.sleep(retry_after)
            return fetch_cds(transcript_id)
        return None
    seq = "".join(r.text.split("\n")[1:])  # Remove the fasta header
    return seq

import pandas as pd

# Streamlit app
def main():
    global debug

    # Uncomment the line below when you need to setup the tables
    setup_tables()

    st.title('Fetch CDS Sequences')

    # Add debug checkbox
    debug_checkbox = st.checkbox('Enable Debugging')
    debug = debug_checkbox

    # File uploader for Excel file
    file = st.file_uploader('Upload an Excel file:', type=['xlsx'], key="file_uploader")

    # Button to fetch CDS sequences
    fetch_button = st.button('Fetch CDS Sequences')

    if fetch_button:
        if file is not None:
            # Read the Excel file into a DataFrame
            df = pd.read_excel(file)

            # Get the gene symbols from the column headers
            gene_symbols = df.columns.tolist()

            # Progress bar
            progress_bar = st.progress(0)

            # Convert HGNC symbols to Ensembl IDs
            ensembl_ids = {}
            for i, gene_symbol in enumerate(gene_symbols):
                try:
                    symbol, ensembl_id = symbol_to_id(gene_symbol)
                    ensembl_ids[symbol] = ensembl_id
                except Exception as e:
                    st.error(f'Failed to fetch Ensembl ID for {gene_symbol}: {e}')

                # Fetch transcript IDs for each gene
                transcripts = fetch_transcripts([ensembl_id for ensembl_id in ensembl_ids.values() if ensembl_id])

                if ensembl_id in transcripts:
                    transcript_ids = transcripts[ensembl_id]
                    if len(transcript_ids) == 1:  # Only fetch sequence for genes with one transcription
                        transcript_id = transcript_ids[0]
                        cds = fetch_cds(transcript_id)
                        if cds:
                            add_gene_and_sequence(gene_symbol, cds, i+1, len(ensembl_ids))
                        elif debug:
                            st.write(f"Gene {i+1} of {len(ensembl_ids)}: {gene_symbol} - Not added, no CDS found.")
                    elif debug:
                        st.write(f"Gene {i+1} of {len(ensembl_ids)}: {gene_symbol} - Not added, multiple CDS found.")

                # Update the progress bar
                progress_bar.progress((i+1) / len(gene_symbols))

        else:
            st.error('Please upload an Excel file.')

if __name__ == "__main__":
    main()
