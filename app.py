import streamlit as st
from sequence_utils import load_data, gc_content, reverse_compliment, sequence_transcribe, translate, find_motif

st.title("DNA Sequence Analyzer")

uploaded_file = st.file_uploader("Upload a FASTA file", type=["fasta", "fa"])

if uploaded_file:
    seq = load_data(uploaded_file)
    st.write(f"Sequence length: {len(seq)}")
    st.write(f"GC Content: {gc_content(seq)}%")
    st.write(f"Reverse Complement: {reverse_compliment(seq)}")
    st.write(f"Transcribed RNA: {sequence_transcribe(seq)}")
    st.write(f"Translated Protein: {translate(sequence_transcribe(seq))}")

