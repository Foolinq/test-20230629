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

# Function to fetch transcript IDs for a given gene
def fetch_transcripts(ensembl_id):
    server = "https://rest.ensembl.org"
    ext = f"/lookup/id/{ensembl_id}?content-type=application/json;expand=1"
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
        r.raise_for_status()
    decoded = r.json()
    return [transcript['id'] for transcript in decoded['Transcript']]

# Function to fetch CDS for a given transcript
def fetch_cds(transcript_id):
    server = "https://rest.ensembl.org"
    ext = f"/sequence/id/{transcript_id}?content-type=text/x-fasta;type=cds"
    r = requests.get(server+ext, headers={ "Content-Type" : "text/x-fasta"})
    if not r.ok:
        return None
    seq = "".join(r.text.split("\n")[1:])  # Remove the fasta header
    return seq

# Streamlit app
def main():
    st.title('Fetch CDS Sequences')

    # Text input for gene symbols
    gene_symbols_input = st.text_input('Enter gene symbols, separated by commas or tabs:', 'BRCA2, TP53, BRAF')

    # Button to fetch CDS sequences
    fetch_button = st.button('Fetch CDS Sequences')

    if fetch_button:
        # Split input into list of gene symbols
        gene_symbols = [symbol.strip() for symbol in gene_symbols_input.replace('\t', ',').split(',')]

        # Fetch CDS for each gene
        cds_dict = {}
        single_cds_dict = {}
        for gene_symbol in gene_symbols:
            try:
                ensembl_id = symbol_to_id(gene_symbol)
                if ensembl_id:
                    transcript_ids = fetch_transcripts(ensembl_id)
                    cds_count = 0
                    for transcript_id in transcript_ids:
                        try:
                            cds = fetch_cds(transcript_id)
                            if cds:
                                cds_dict[gene_symbol] = cds
                                cds_count += 1
                                if cds_count > 1:
                                    break  # Stop after finding more than one CDS
                        except Exception as e:
                            st.error(f'Failed to fetch CDS for {gene_symbol} transcript {transcript_id}: {e}')
                    if cds_count == 1:
                        single_cds_dict[gene_symbol] = cds_dict[gene_symbol]
                else:
                    st.error(f'Failed to find Ensembl ID for {gene_symbol}')
            except Exception as e:
                st.error(f'Failed to fetch CDS for {gene_symbol}: {e}')

        # Display CDS dictionaries
        st.write('All genes with at least one CDS:', cds_dict)
        st.write('Genes with exactly one CDS:', single_cds_dict)

if __name__ == "__main__":
    main()
