import os
from langchain import OpenAI, SQLDatabase
from langchain.chains import SQLDatabaseSequentialChain
import streamlit as st



API_KEY = os.environ.get('OPENAI_API_KEY')
DBURL = os.environ.get('DBURL')




db = SQLDatabase.from_uri(DBURL)




llm = OpenAI(temperature=0, openai_api_key=OPENAI_API_KEY, model_name='gpt-3.5-turbo')

PROMPT = """
Given an input question, first create a syntactically correct postgresql query to run,
then look at the results of the query and return the answer.
The question: {question}
"""

db_chain = SQLDatabaseSequentialChain(llm=llm, database=db, verbose=True, top_k=3)



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
