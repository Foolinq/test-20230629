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
def query_database(prompt, return_sql=False):
    db_chain.return_sql = return_sql
    inputs = {"query": prompt}
    response = db_chain(inputs)
    return response['result']

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
