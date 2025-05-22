import streamlit as st
from sequence_utils import (
    load_data, gc_content, reverse_complement,
    sequence_transcribe, translate_rna,
    find_motif, trim_to_codon_length
)
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Page setup
st.set_page_config(page_title="DNA Sequence Analyzer", layout="wide")
st.title("🧬 DNA Sequence Analyzer")

# File uploader
uploaded_file = st.file_uploader("📄 Upload a FASTA file", type=["fasta", "fa"])


@st.cache_data
def process_sequence(_seq: Seq) -> dict:
    """
    Analyze a DNA sequence and return GC content, reverse complement,
    transcribed RNA, and protein translation.
    """
    rna = sequence_transcribe(_seq)
    protein = translate_rna(trim_to_codon_length(rna))
    return {
        "gc": gc_content(str(_seq)),
        "rev_comp": reverse_complement(_seq),
        "rna": rna,
        "protein": protein,
    }


if uploaded_file:
    with st.spinner("🔍 Reading and parsing sequences..."):
        try:
            records = load_data(uploaded_file)
        except ValueError as e:
            st.error(f"Error: {e}")
            st.stop()

    if not records:
        st.error("❌ No valid sequences found in the file.")
        st.stop()

    # Select a sequence
    selected_id = st.selectbox(
        "🔢 Select a sequence ID", [rec.id for rec in records[:10]]
    )
    selected_record: SeqRecord = next(rec for rec in records if rec.id == selected_id)
    seq: Seq = selected_record.seq

    # Show raw sequence
    st.subheader(f"🧾 Raw Sequence for {selected_id} (First 200 bp):")
    st.text_area("📜 Sequence (scrollable)", str(seq), height=100)

    # Process the selected sequence
    with st.spinner("🧪 Running sequence analysis..."):
        result = process_sequence(seq)
        print(result)

        # Show results
        st.subheader("🧪 GC Content")
        st.write(f"{gc_content(seq)}%")

        st.subheader("🔁 Reverse Complement")
        st.code(str(reverse_complement(seq)), language="text")
        #st.code(str(result["rev_comp"]), language="text")

        st.subheader("🔤 Transcribed RNA")
        st.code(str(sequence_transcribe(seq)), language="text")
        #st.code(str(result["rna"]), language="text")

        st.subheader("🧬 Protein Translation")
        st.code(str(translate_rna(seq)), language="text")

        # Motif Search
        st.subheader("🔎 Motif Search")
        motif = st.text_input("Enter motif (e.g., ATG or TATA):").upper().strip()

        if motif:
            if not motif.isalpha():
                st.warning("❗ Please enter a valid DNA motif using only letters (A, T, C, G).")
            elif len(motif) > 20:
                st.warning("❗ Motif too long. Try something shorter.")
            else:
                positions = find_motif(str(seq), motif)
                if positions:
                    st.success(f"✅ Motif found at positions: {positions}")
                else:
                    st.warning("🔍 Motif not found.")
