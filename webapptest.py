import streamlit as st

st.write("""
# Test(icles)""")

input = st.text_area("Ouin", "Les gars sont on?", height=70)


st.write("""
***
""")

st.header("C on?")
input

if input == "Oui":
  st.header("Ketttt")
elif input == "Non":
  st.header("Yo t off")
else:
  st.header("Comme si que...")
