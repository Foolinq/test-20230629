import os
from langchain_experimental.sql import SQLDatabaseChain
from langchain import OpenAI, SQLDatabase
import streamlit as st

# Retrieve Heroku Config Vars
OPENAI_API_KEY = os.environ.get('OPENAI_API_KEY')
DBURL = os.environ.get('DBURL')

# Set Up Database Connection
db = SQLDatabase.from_uri(DBURL)

# Set Up LangChain
llm = OpenAI(temperature=0, openai_api_key=OPENAI_API_KEY, model_name='gpt-3.5-turbo')
db_chain = SQLDatabaseChain.from_llm(llm, db)

def query_database(prompt, return_sql=False):
    # Preprocess the prompt
    processed_prompt = preprocess_input(prompt)
    
    db_chain.return_sql = return_sql
    inputs = {"query": processed_prompt}
    response = db_chain(inputs)
    return response['result']

def preprocess_input(user_input):
    # Add context or researched information to the user input
    context = "The data is housed in a PostgreSQL database and consists of gene-related information. There are tables containing gene annotations, detailing the characteristics and functions of genes. Other tables count the codons within coding sequences, providing insights into genetic coding structure. Additionally, there are tables with specific information related to ovarian cancer patients, including gene expression levels and codon incorporation rates. Together, these tables offer a multifaceted view of genetic data and its relevance to ovarian cancer."
    researched_info = ""
    return f"{context} {researched_info} {user_input}"

# Streamlit Interface
def main():
    st.title('Natural Language Query Interface')
    user_input = st.text_input('Please enter your natural language query:')
    
    send_button = st.button('Send')
    check_button = st.button('Check')

    if send_button and user_input:
        response = query_database(user_input)
        st.write(response)

    if check_button and user_input:
        sql_query = query_database(user_input, return_sql=True)
        st.write(f"SQL Query: {sql_query}")

if __name__ == '__main__':
    main()
