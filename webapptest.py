from sqlalchemy import create_engine, Column, String, ForeignKey
from sqlalchemy.orm import sessionmaker, relationship
from sqlalchemy.ext.declarative import declarative_base
import pandas as pd
import streamlit as st
import requests
import concurrent.futures
import time

# Define the SQLAlchemy models
Base = declarative_base()

class Gene(Base):
    __tablename__ = 'genes'

    name = Column(String, primary_key=True)
    transcriptions = relationship('Transcription', back_populates='gene')

class Transcription(Base):
    __tablename__ = 'transcriptions'

    id = Column(String, primary_key=True)
    sequence = Column(String)
    gene_name = Column(String, ForeignKey('genes.name'))
    gene = relationship('Gene', back_populates='transcriptions')

# Set up the database
engine = create_engine('postgresql://xmybdbirwwaocp:d2a20d22a09e4948cda6b9688769effa0a9f0b0393c0af7ba95d04b07b5d6fff@ec2-52-205-45-222.compute-1.amazonaws.com:5432/d9addgatbiak5p')
Session = sessionmaker(bind=engine)

Base.metadata.create_all(engine)

from sqlalchemy import Table, MetaData

metadata = MetaData()

if not Table('genes', metadata, autoload_with=engine).exists():
    Base.metadata.create_all(engine)





def add_gene_and_transcriptions(gene_name, transcriptions):
    # Start a new session
    session = Session()

    try:
        # Check if the gene already exists
        gene = session.query(Gene).filter_by(name=gene_name).first()
        if gene is not None:
            # If the gene already exists, skip it
            return

        # Create a new gene
        gene = Gene(name=gene_name)

        # Add the gene to the session
        session.add(gene)

        # Create a new transcription for each sequence and add it to the session
        for transcription_id, sequence in transcriptions.items():
            transcription = Transcription(id=transcription_id, sequence=sequence, gene=gene)
            session.add(transcription)

        # Commit the session to save the changes to the database
        session.commit()
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
    st.title('Fetch CDS Sequences')

    # File uploader for Excel file
    file = st.file_uploader('Upload an Excel file:', type=['xlsx'])

    # Button to fetch CDS sequences
    fetch_button = st.button('Fetch CDS Sequences')

    if fetch_button:
        if file is not None:
            # Read the Excel file into a DataFrame
            df = pd.read_excel(file)

            # Get the gene symbols from the column headers
            gene_symbols = df.columns.tolist()

            # Convert HGNC symbols to Ensembl IDs
            with concurrent.futures.ThreadPoolExecutor() as executor:
                future_to_symbol = {executor.submit(symbol_to_id, gene_symbol): gene_symbol for gene_symbol in gene_symbols}
                ensembl_ids = {}
                for future in concurrent.futures.as_completed(future_to_symbol):
                    gene_symbol = future_to_symbol[future]
                    try:
                        symbol, ensembl_id = future.result()
                        ensembl_ids[symbol] = ensembl_id
                    except Exception as e:
                        st.error(f'Failed to fetch Ensembl ID for {gene_symbol}: {e}')

            # Fetch transcript IDs for each gene
            transcripts = fetch_transcripts([ensembl_id for ensembl_id in ensembl_ids.values() if ensembl_id])

            # Fetch CDS for each gene
            with concurrent.futures.ThreadPoolExecutor() as executor:
                for gene_symbol, ensembl_id in ensembl_ids.items():
                    try:
                        if ensembl_id in transcripts:
                            transcript_ids = transcripts[ensembl_id]
                            cds_dict = {}
                            future_to_cds = {executor.submit(fetch_cds, transcript_id): transcript_id for transcript_id in transcript_ids}
                            for future in concurrent.futures.as_completed(future_to_cds):
                                transcript_id = future_to_cds[future]
                                try:
                                    cds = future.result()
                                    if cds:
                                        cds_dict[transcript_id] = cds
                                except Exception as e:
                                    st.error(f'Failed to fetch CDS for {gene_symbol} transcript {transcript_id}: {e}')
                            if cds_dict:
                                add_gene_and_transcriptions(gene_symbol, cds_dict)
                    except Exception as e:
                        st.error(f'Failed to fetch CDS for {gene_symbol}: {e}')
        else:
            st.error('Please upload an Excel file.')

if __name__ == "__main__":
    main()
