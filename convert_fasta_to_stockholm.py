from bio import SeqIO

records = SeqIO.parse("PF08854.fasta", "fasta")
count = SeqIO.write(records, "PF08854.stockholm", "stockholm")
print("Converted %i records" % count)