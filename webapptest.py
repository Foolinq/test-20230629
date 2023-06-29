import streamlit as st

st.write("""
# Les bhy sont on?""")

input = st.text_area("","",height=10)


st.write("""
***
""")

if input == "Oui":
  st.header("Ketttt")
elif input == "Non":
  st.header("Yo t off")
else:
  st.header("Comme si que...")
