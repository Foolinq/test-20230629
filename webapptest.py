from langchain.agents import create_sql_agent
from langchain.agents.agent_toolkits import SQLDatabaseToolkit
from langchain.sql_database import SQLDatabase
from langchain.llms.openai import OpenAI
from langchain.agents import AgentExecutor
from langchain.agents.agent_types import AgentType
from langchain.chat_models import ChatOpenAI



db = SQLDatabase.from_uri("sqlite:///path/to/your/database.db")
toolkit = SQLDatabaseToolkit(db=db, llm=OpenAI(temperature=0))



agent_executor = create_sql_agent(
    llm=OpenAI(temperature=0),
    toolkit=toolkit,
    verbose=True,
    agent_type=AgentType.ZERO_SHOT_REACT_DESCRIPTION,
)



def query_database(query):
    result = agent_executor.run(query)
    return result



import streamlit as st

def main():
    st.title('Natural Language Query Interface')
    user_input = st.text_input('Please enter your natural language query:')
    if user_input:
        response = query_database(user_input)
        st.write(response)

if __name__ == '__main__':
    main()
