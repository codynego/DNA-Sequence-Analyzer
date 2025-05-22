from typing import List
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from io import StringIO


def load_data(uploaded_file) -> List[SeqRecord]:
    """
    Load a FASTA file and return a list of SeqRecord objects.
    Handles multiple records in the file.
    """
    try:
        string_data = uploaded_file.read().decode("utf-8")
        fasta_io = StringIO(string_data)
        return list(SeqIO.parse(fasta_io, "fasta"))
    except Exception as e:
        raise ValueError("Failed to read or parse FASTA file.") from e


def trim_to_codon_length(seq: Seq) -> Seq:
    """
    Trim the sequence length so it's divisible by 3.
    Useful for proper translation into codons.
    """
    trimmed_length = len(seq) - (len(seq) % 3)
    return seq[:trimmed_length]


def gc_content(seq: str) -> float:
    """
    Calculate the GC content percentage of a DNA sequence.
    """
    if not seq:
        return 0.0
    g = seq.upper().count("G")
    c = seq.upper().count("C")
    return round((g + c) / len(seq) * 100, 2)


def reverse_complement(seq: Seq) -> Seq:
    """
    Return the reverse complement of a DNA sequence.
    """
    return seq.reverse_complement()


def sequence_transcribe(seq: Seq) -> Seq:
    """
    Transcribe DNA sequence into RNA.
    """
    return seq.transcribe()


def translate_rna(rna_seq: Seq) -> Seq:
    """
    Translate RNA sequence into amino acids (protein).
    """
    try:
        return rna_seq.translate()
    except Exception as e:
        raise ValueError("Invalid RNA sequence for translation.") from e


def find_motif(seq: str, motif: str) -> List[int]:
    """
    Find all start positions (0-based) of a motif in a given sequence.
    """
    seq = seq.upper()
    motif = motif.upper()
    return [i for i in range(len(seq) - len(motif) + 1) if seq[i:i + len(motif)] == motif]
