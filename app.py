from Bio import SeqIO

record = next(SeqIO.parse("Homo_sapiens_INS_sequence.fasta", "fasta"))
sequence = record.seq
print(f"Loaded sequence: {sequence}")

def gc_content(seq):
    g = seq.count("G")
    c = seq.count("C")
    return round((g + c) / len(seq) * 100, 2)