import os
from langchain import LLMChain, OpenAI, PromptTemplate

def process_llm_request(llm_request):
    api_key = os.getenv('OPENAI_API_KEY')
    llm = LLMChain(llm=OpenAI(api_key=api_key), prompt=PromptTemplate(input_variables=["adjective"], template=llm_request))
    # Here you would call the function to process the request with the LLM model
    # result = llm.process_request()
    # return result
