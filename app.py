import streamlit as st
from Bio import SeqIO

st.title("DNA Sequence Analyzer")

uploaded_file = st.file_uploader("Upload a FASTA file", type=["fasta", "fa"])

if uploaded_file:
    record = next(SeqIO.parse(uploaded_file, "fasta"))
    seq = record.seq
    st.write(f"Sequence length: {len(seq)}")
    st.write(f"GC Content: {(seq.count('G') + seq.count('C')) / len(seq) * 100:.2f}%")
    st.write(f"Reverse Complement: {seq.reverse_complement()}")
    st.write(f"Transcribed RNA: {seq.transcribe()}")
    st.write(f"Translated Protein: {seq.transcribe().translate()}")
