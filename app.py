import streamlit as st
from sequence_utils import load_data, gc_content, reverse_complement, sequence_transcribe, translate_rna, find_motif

st.title("🧬 DNA Sequence Analyzer")

uploaded_file = st.file_uploader("Upload a FASTA file", type=["fasta", "fa"])

if uploaded_file:
    seq = load_data(uploaded_file)
    
    st.subheader("📄 Raw Sequence:")
    st.code(str(seq))

    # GC Content
    st.subheader("🧪 GC Content")
    gc = gc_content(str(seq))
    st.write(f"{gc}%")

    # Reverse Complement
    st.subheader("🔁 Reverse Complement")
    st.code(str(reverse_complement(seq)))

    # Transcription
    st.subheader("📤 Transcribed RNA")
    rna = sequence_transcribe(seq)
    st.code(str(rna))

    # Translation
    st.subheader("🧬 Protein Translation")
    protein = translate_rna(rna)
    st.code(str(protein))

    # Motif Search
    st.subheader("🔍 Motif Search")
    motif = st.text_input("Enter motif (e.g., ATG):")
    if motif:
        positions = find_motif(str(seq), motif)
        if positions:
            st.write(f"Motif found at positions: {positions}")
        else:
            st.warning("Motif not found.")