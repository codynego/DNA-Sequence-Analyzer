import streamlit as st
from Bio import SeqIO
from io import StringIO

def load_data(uploaded_file):
    string_data = uploaded_file.read().decode("utf-8")
    fasta_io = StringIO(string_data)
    seq_record = SeqIO.read(fasta_io, "fasta")
    return seq_record.seq  # Only return the sequence string

def gc_content(seq):
    g = seq.count("G")
    c = seq.count("C")
    return round((g + c) / len(seq) * 100, 2)

def reverse_complement(seq):
    return seq.reverse_complement()

def sequence_transcribe(seq):
    return seq.transcribe()

def translate_rna(rna_seq):
    return rna_seq.translate()

def find_motif(seq, motif):
    return [i for i in range(len(seq) - len(motif) + 1) if seq[i:i+len(motif)] == motif]