from Bio import SeqIO

record = next(SeqIO.parse("Homo_sapiens_INS_sequence.fasta", "fasta"))
sequence = record.seq
print(f"Loaded sequence: {sequence}")

def gc_content(seq):
    g = seq.count("G")
    c = seq.count("C")
    return round((g + c) / len(seq) * 100, 2)

def reverse_compliment(seq):
    return seq.reverse_compliment()  

def sequence_transcribe(seq):
    """Convert Dna to Rna

    Args:
        seq (_type_): The dna sequence to transcribe
    """

    return seq.transcribe()

def translate(rna_seq):
    """Making Protein from RNA

    Args:
        seq (_type_): The rna sequence
    """
    return rna_seq.translate()

def find_motif(seq, motif):
    return [i for i in range(len(seq)) if seq[i:i+len(motif)] == motif]