import streamlit as st
import requests

# Function to fetch CDS for a given gene
def fetch_cds(gene_symbol):
    server = "https://rest.ensembl.org"
    ext = f"/sequence/id/{gene_symbol}?content-type=text/x-fasta;type=cds"
    r = requests.get(server+ext, headers={ "Content-Type" : "text/x-fasta"})
    if not r.ok:
        r.raise_for_status()
    seq = "".join(r.text.split("\n")[1:])  # Remove the fasta header
    return seq

# Streamlit app
def main():
    st.title('Fetch CDS Sequences')

    # Text input for gene symbols
    gene_symbols_input = st.text_input('Enter gene symbols, separated by commas:', 'BRCA2, TP53, BRAF')

    # Button to fetch CDS sequences
    fetch_button = st.button('Fetch CDS Sequences')

    if fetch_button:
        # Split input into list of gene symbols
        gene_symbols = [symbol.strip() for symbol in gene_symbols_input.split(',')]

        # Fetch CDS for each gene
        cds_dict = {}
        for gene_symbol in gene_symbols:
            try:
                cds = fetch_cds(gene_symbol)
                cds_dict[gene_symbol] = cds
            except Exception as e:
                st.error(f'Failed to fetch CDS for {gene_symbol}: {e}')

        # Display CDS dictionary
        st.write(cds_dict)

if __name__ == "__main__":
    main()
