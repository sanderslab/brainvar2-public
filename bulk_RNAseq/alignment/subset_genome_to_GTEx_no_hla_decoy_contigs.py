def read_contig_names(filename):
    """
    Read contig names from a file and return them as a set.
    """
    contigs = set()
    with open(filename, 'r') as file:
        for line in file:
            contig = line.strip()
            if contig:
                contigs.add(contig)
    return contigs

def subset_fasta(input_fasta, output_fasta, contig_names):
    """
    Subset a FASTA file to include only lines related to specified contig names.
    """
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        write_flag = False
        for line in infile:
            if line.startswith('>'): 
                header = line.strip()[1:]
                contig_name = header.split()[0]
                if contig_name in contig_names:
                    write_flag = True
                    outfile.write(line)
                else:
                    write_flag = False
            elif write_flag:
                outfile.write(line)

if __name__ == "__main__":
    input_fasta_file = "/home/ec2-user/ref_genome/GRCh38.p13.genome.fa"  # Replace with your input FASTA file
    output_fasta_file = "/home/ec2-user/GRCh38.p13.GTEx.no-alt-hla-decoy.genome.fa"  # Replace with the desired output file
    contig_names_file = "/home/ec2-user/GRCh38.p13.Broad.no_alt_no_decoy.contigs.txt"  # Replace with the file containing contig names to subset
    
    contig_names = read_contig_names(contig_names_file)
    subset_fasta(input_fasta_file, output_fasta_file, contig_names)
    print("Subsetting complete.")