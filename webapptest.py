import os
from langchain import SQLDatabase, SQLDatabaseChain
from langchain.chat_models import ChatOpenAI
import streamlit as st


API_KEY = os.environ.get('OPENAI_API_KEY')
DBURL = os.environ.get('DBURL')


db = SQLDatabase.from_uri(DBURL)


llm = ChatOpenAI(temperature=0, openai_api_key=API_KEY, model_name='gpt-3.5-turbo')

QUERY = """
Given an input question, first create a syntactically correct postgresql query to run, then look at the results of the query and return the answer.
Use the following format:

Question: "Question here"
SQLQuery: "SQL Query to run"
SQLResult: "Result of the SQLQuery"
Answer: "Final answer here"

{question}
"""

db_chain = SQLDatabaseChain(llm=llm, database=db, verbose=True)


def query_database(prompt):
    question = QUERY.format(question=prompt)
    return db_chain.run(question)


def main():
    st.title('Natural Language Query Interface')
    user_input = st.text_input('Please enter your natural language query:')
    if user_input:
        response = query_database(user_input)
        st.write(response)

if __name__ == '__main__':
    main()
