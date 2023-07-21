from sqlalchemy import create_engine, Column, Integer, String, ForeignKey
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base

# Define the database
DATABASE_URI = 'postgresql://jqczrvezvwlzer:5980122424475b4fcdc2e539dbb0667d254452a21cd3cc633f9a64a8796da2df@ec2-3-212-70-5.compute-1.amazonaws.com:5432/d202roftknash2'
engine = create_engine(DATABASE_URI)
Session = sessionmaker(bind=engine)
Base = declarative_base()

# Define the Gene class
class Gene(Base):
    __tablename__ = 'genes'

    id = Column(Integer, primary_key=True)
    name = Column(String)
    sequence = Column(String)
    expression = Column(Integer)

class Codon(Base):
    __tablename__ = 'codons'

    id = Column(Integer, primary_key=True)
    gene_name = Column(String, ForeignKey('genes.name'))
    codon = Column(String)
    amino_acid = Column(String)

# Define the genetic code
GENETIC_CODE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T', 
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                  
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P', 
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R', 
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A', 
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G', 
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L', 
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W', 
}

def setup_tables():
    Base.metadata.create_all(engine)

def add_codon(gene_name, codon, amino_acid):
    session = Session()

    try:
        new_codon = Codon(gene_name=gene_name, codon=codon, amino_acid=amino_acid)
        session.add(new_codon)
        session.commit()
    except:
        session.rollback()
        raise
    finally:
        session.close()

def process_gene_codons(gene_name):
    session = Session()

    try:
        gene = session.query(Gene).filter_by(name=gene_name).first()

        sequence = gene.sequence
        codons = [sequence[i:i+3] for i in range(0, len(sequence), 3)]

        for codon in codons:
            if len(codon) != 3:
                continue

            amino_acid = GENETIC_CODE.get(codon, 'X')
            add_codon(gene_name, codon, amino_acid)
    except Exception as e:
        print(f'Failed to process {gene_name}: {e}')
    finally:
        session.close()

if __name__ == "__main__":
    setup_tables()
    
    # Query all the genes from your database
    session = Session()
    gene_names = [gene.name for gene in session.query(Gene).all()]
    session.close()

    # Process all the genes
    for gene_name in gene_names:
        process_gene_codons(gene_name)

    st.write("Done!")
