import os
from langchain_experimental.sql import SQLDatabaseChain
from langchain import OpenAI, SQLDatabase
from langchain.schema.prompt_template import BasePromptTemplate
import streamlit as st

# Retrieve Heroku Config Vars
OPENAI_API_KEY = os.environ.get('OPENAI_API_KEY')
DBURL = os.environ.get('DBURL')

# Set Up Database Connection
db = SQLDatabase.from_uri(DBURL)

# Set Up LangChain
llm = OpenAI(temperature=0, openai_api_key=OPENAI_API_KEY, model_name='gpt-3.5-turbo')

# Define your custom prompt class
class CustomPrompt(BasePromptTemplate):
    template: str = "This data represents genetic information. The 'codon_incorporation_rate' table contains the incorporation rate of each codon for the genes contained in the 'genes' for 578 patients. The 'codons' table contains the quantity of each codon for a couple hundred genes. The 'patient_gene_expression' table contains the level of expression of each gene in the 'genes' table for each patient."
    input_variables: list = ["query"]

    def format(self, **kwargs):
        return self.template.format(**kwargs)

# Set up LangChain with the custom prompt
custom_prompt = CustomPrompt()
db_chain = SQLDatabaseChain.from_llm(llm, db, prompt=custom_prompt)

# Query Function
def query_database(prompt):
    inputs = {"query": prompt}
    response = db_chain(inputs)
    return response['result']

# Streamlit Interface
def main():
    st.title('Natural Language Query Interface')
    user_input = st.text_input('Please enter your natural language query:')
    if user_input:
        response = query_database(user_input)
        st.write(response)

if __name__ == '__main__':
    main()
