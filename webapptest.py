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
    inputs = {"query": prompt}
    response = db_chain(inputs, return_sql=True)  # Return SQL without executing
    return response['result']

# Streamlit Interface
def main():
    st.title('Natural Language Query Interface')
    user_input = st.text_input('Please enter your natural language query:')
    send_button = st.button('Send')  # Send button

    if send_button and user_input:
        sql_query = query_database(user_input)
        st.write(f"The SQL query to be run is: {sql_query}")
        confirm_button = st.button('Confirm')  # Confirm button

        if confirm_button:
            # Execute the SQL query
            inputs = {"query": user_input}
            response = db_chain(inputs)
            st.write(response['result'])

if __name__ == '__main__':
    main()
