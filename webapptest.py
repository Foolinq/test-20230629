import streamlit as st
import requests
import concurrent.futures
import time

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
        cds_dict = {}
        single_cds_dict = {}
        with concurrent.futures.ThreadPoolExecutor() as executor:
            for gene_symbol, ensembl_id in ensembl_ids.items():
                try:
                    if ensembl_id in transcripts:
                        transcript_ids = transcripts[ensembl_id]
                        cds_count = 0
                        future_to_cds = {executor.submit(fetch_cds, transcript_id): transcript_id for transcript_id in transcript_ids}
                        for future in concurrent.futures.as_completed(future_to_cds):
                            transcript_id = future_to_cds[future]
                            try:
                                cds = future.result()
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
