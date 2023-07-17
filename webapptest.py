import streamlit as st
import requests

# Function to convert HGNC symbol to Ensembl ID
def symbol_to_id(gene_symbol):
    server = "https://rest.ensembl.org"
    ext = f"/xrefs/symbol/homo_sapiens/{gene_symbol}?content-type=application/json"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
        r.raise_for_status()
    decoded = r.json()
    return decoded[0]['id'] if decoded else None

# Function to fetch CDS for a given gene
def fetch_cds(ensembl_id):
    server = "https://rest.ensembl.org"
    ext = f"/sequence/id/{ensembl_id}?content-type=text/x-fasta;type=cds"
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
                ensembl_id = symbol_to_id(gene_symbol)
                if ensembl_id:
                    cds = fetch_cds(ensembl_id)
                    cds_dict[gene_symbol] = cds
                else:
                    st.error(f'Failed to find Ensembl ID for {gene_symbol}')
            except Exception as e:
                st.error(f'Failed to fetch CDS for {gene_symbol}: {e}')

        # Display CDS dictionary
        st.write(cds_dict)

if __name__ == "__main__":
    main()
