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

# Query Function
def query_database(prompt):
    inputs = {"question": prompt}
    response = db_chain(inputs)
    return response

# Streamlit Interface
def main():
    st.title('Natural Language Query Interface')
    user_input = st.text_input('Please enter your natural language query:')
    if user_input:
        response = query_database(user_input)
        st.write(response)

if __name__ == '__main__':
    main()
