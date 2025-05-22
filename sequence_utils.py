from Bio import SeqIO

def load_data(file):
    """Load a fasta file and return the first sequence record.

    Args:
        file_path (str): Path to the fasta file.

    Returns:
        Bio.SeqRecord.SeqRecord: The first sequence record in the fasta file.
    """
    with open(file, "r") as fasta_file:
        seq_record = SeqIO.read(fasta_file, "fasta")
    return seq_record.seq



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