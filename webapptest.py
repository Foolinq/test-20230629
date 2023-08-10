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
db_chain = SQLDatabaseChain.from_llm(llm, db, return_sql=True)  # Set return_sql to True

# Query Function
def query_database(prompt):
    inputs = {"query": prompt}
    response = db_chain(inputs)  # Return SQL without executing
    return response['result']

# Streamlit Interface
def main():
    st.title('Natural Language Query Interface')
    user_input = st.text_input('Please enter your natural language query:')
    send_button = st.button('Send')
    if send_button and user_input:
        sql_query = query_database(user_input)
        st.write(f"SQL Query: {sql_query}")
        confirm_button = st.button('Confirm')
        if confirm_button:
            db_chain.return_sql = False  # Set return_sql to False to execute the query
            response = query_database(user_input)
            st.write(response)

if __name__ == '__main__':
    main()
