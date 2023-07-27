import os
import streamlit as st
import pandas as pd
import base64
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy import Column, String, Integer, Text, Float, ForeignKey, UniqueConstraint
from sqlalchemy.ext.declarative import declarative_base
from collections import defaultdict
import json

# Get your database connection from the environment variable
DATABASE_URL = os.environ.get('DBURL')

engine = create_engine(DATABASE_URL)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)
Base = declarative_base()

# Your Gene class
class Gene(Base):
    __tablename__ = 'genes'

    name = Column(String, primary_key=True)
    sequence = Column(String)

# Your PatientGeneExpression class
class PatientGeneExpression(Base):
    __tablename__ = 'patient_gene_expression'

    id = Column(Integer, primary_key=True)
    gene_expressions = Column(Text)

class Codon(Base):
    __tablename__ = 'codons'

    id = Column(Integer, primary_key=True)
    gene_name = Column(String, ForeignKey('genes.name'))
    codon = Column(String)
    amino_acid = Column(String)
    count = Column(Integer)
    __table_args__ = (UniqueConstraint('gene_name', 'codon', name='uc_gene_codon'), )

class CodonIncorporationRate(Base):
    __tablename__ = 'codon_incorporation_rate'

    id = Column(Integer, primary_key=True)
    rates = Column(Text)

class GlobalCodonFrequency(Base):
    __tablename__ = 'global_codon_frequency'

    id = Column(Integer, primary_key=True)
    frequencies = Column(Text)

Base.metadata.create_all(bind=engine)

def main():
    st.title("CSV File Cleaning and Data Storing")

    st.header("Clean CSV File")
    csv_file = st.file_uploader("Upload CSV file for cleaning", type=['csv'])
    if csv_file is not None:
        df = pd.read_csv(csv_file)
        if st.button("Start Cleaning"):
            cleaned_df = clean_csv(df)
            st.dataframe(cleaned_df)
            csv_download_link(cleaned_df)

    st.header("Store Gene Expressions")
    cleaned_csv_file = st.file_uploader("Upload cleaned CSV file for data storing", type=['csv'])
    if cleaned_csv_file is not None:
        cleaned_df = pd.read_csv(cleaned_csv_file)
        if st.button("Store Gene Expressions"):
            store_gene_expressions(cleaned_df)

    st.header("Perform Analysis")
    if st.button("Analyze"):
        calculate_codon_incorporation_and_frequency()
        st.write("Analysis complete.")

def clean_csv(df):
    # Create a session
    session = SessionLocal()

    # Get gene names from 'genes' table in the database
    gene_names = session.query(Gene.name).all()
    gene_names = [gene[0] for gene in gene_names]

    # Keep only columns from dataframe that exist in the 'genes' table
    df = df.loc[:, df.columns.isin(gene_names)]
    
    return df

def store_gene_expressions(df):
    # Create a session
    session = SessionLocal()

    for i, row in df.iterrows():
        gene_expressions = ','.join(map(str, row.tolist()))
        patient = PatientGeneExpression(id=i+1, gene_expressions=gene_expressions)
        session.add(patient)

    session.commit()

def csv_download_link(df):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode() 
    href = f'<a href="data:file/csv;base64,{b64}" download="cleaned.csv">Download cleaned CSV File</a>'
    st.markdown(href, unsafe_allow_html=True)

def calculate_codon_incorporation_and_frequency():
    # Create a session
    session = SessionLocal()

    # Declaration of codon_incorporation_rates dictionary
    codon_incorporation_rates = dict()
    
    patients = session.query(PatientGeneExpression).all()
    genes = session.query(Gene).all()
    codons = session.query(Codon).all()

    for patient in patients:
        patient_expressions = list(map(float, patient.gene_expressions.split(',')))
        for gene, expression in zip(genes, patient_expressions):
            for codon in codons:
                if codon.gene_name == gene.name:
                    if patient.id not in codon_incorporation_rates:
                        codon_incorporation_rates[patient.id] = defaultdict(int)
                    codon_incorporation_rates[patient.id][codon.codon] += codon.count * expression

    for patient_id, rates in codon_incorporation_rates.items():
        incorporation_rate = CodonIncorporationRate(id=patient_id, rates=json.dumps(rates))
        session.add(incorporation_rate)

    session.commit()

if __name__ == "__main__":
    main()
