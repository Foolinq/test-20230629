import os
from langchain import OpenAI, SQLDatabase
from langchain_experimental.sql import SQLDatabaseChain
import streamlit as st



OPENAI_API_KEY = os.environ.get('OPENAI_API_KEY')
DBURL = os.environ.get('DBURL')




db = SQLDatabase.from_uri(DBURL)




llm = OpenAI(temperature=0, openai_api_key=OPENAI_API_KEY, model_name='gpt-3.5-turbo')
db_chain = SQLDatabaseChain.from_llm(llm, db)




def query_database(prompt):
    return db_chain.call({"question": prompt})



def main():
    st.title('Natural Language Query Interface')
    user_input = st.text_input('Please enter your natural language query:')
    if user_input:
        response = query_database(user_input)
        st.write(response)

if __name__ == '__main__':
    main()
