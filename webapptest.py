import base64
import streamlit as st
import pandas as pd
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker
from sqlalchemy import Column, String
from sqlalchemy.ext.declarative import declarative_base

# Your database connection
DATABASE_URL = "postgresql://jqczrvezvwlzer:5980122424475b4fcdc2e539dbb0667d254452a21cd3cc633f9a64a8796da2df@ec2-3-212-70-5.compute-1.amazonaws.com:5432/d202roftknash2"

engine = create_engine(DATABASE_URL)
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)
Base = declarative_base()

# Your Gene class
class Gene(Base):
    __tablename__ = 'genes'

    name = Column(String, primary_key=True)
    sequence = Column(String)

Base.metadata.create_all(bind=engine)

def main():
    st.title("CSV File Cleaning")

    csv_file = st.file_uploader("Upload CSV file", type=['csv'])
    if csv_file is not None:
        df = pd.read_csv(csv_file)

        if st.button("Start Cleaning"):
            cleaned_df = clean_csv(df)
            st.dataframe(cleaned_df)
            csv_download_link(cleaned_df)

def clean_csv(df):
    # Create a session
    session = SessionLocal()

    # Get column names from 'genes' table in the database
    table_columns = [c.name for c in Gene.__table__.columns]
    
    # Drop columns from dataframe that do not exist in the database table
    df = df.loc[:, df.columns.isin(table_columns)]
    
    return df

def csv_download_link(df):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode() 
    href = f'<a href="data:file/csv;base64,{b64}" download="cleaned.csv">Download cleaned CSV File</a>'
    st.markdown(href, unsafe_allow_html=True)

if __name__ == "__main__":
    main()
