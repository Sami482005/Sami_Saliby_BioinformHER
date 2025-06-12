from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
from Bio import Entrez
import argparse
import os
from retrieving_sequences import getting_full_fasta_records

def find_similar_sequences(input_fasta, output_file, log_file=None):
    '''
    This function will loop through the homo sapiens HBB sequences.
    It will use NCBI's BLAST to search for HBB nucleotide sequences from other
    organisms and returns the sequences of the 5 sequences with the highest similarity

    params:
    -------
    - input_fasta: str, path to the fasta file 
    - output_file: str, path to the text file where the accession IDs of the top 5 hits will be written
    - log_file: str, path to log file (optional)

    output:
    -------
    - sim_acc_IDS: txt file containing the accession IDs of the sequences that match

    return:
    -------
    - sim_acc_id_list: list, list of ids of the matches
    '''
    
    sep="\n-------------------------\n"
    
    id_list = []

    if log_file:
        f_log = open(log_file, 'a')
        f_log.write(sep)
    with open(input_fasta) as handle:
        record = SeqIO.read(handle, format="fasta")

    if log_file:
        f_log.write(f"Running BLASTn for {record.id}...")
    
    # Accessing NCBI's online BLASTn tool to searh for the similar sequence
    result_handle = NCBIWWW.qblast("blastn", "nt", record.format("fasta"), entrez_query="HBB[Gene] NOT Homo sapiens[Organism] NOT Humans [Title]")
    blast_record = NCBIXML.read(result_handle)
    
    top_hits = 20
    if log_file:
        f_log.write(f"\n>>> BLAST results for {record.id} <<<\n")
        for alignment in blast_record.alignments[:top_hits]:
            
            hsp = alignment.hsps[0]
            f_log.write(f"Hit ID: {alignment.accession}\n")
            f_log.write(f"  Title: {alignment.title}\n")
            f_log.write(f"  Length: {alignment.length}\n")
            f_log.write(f"  Identity: {hsp.identities}/{hsp.align_length} ({round(hsp.identities / hsp.align_length * 100, 2)}%)\n")
            f_log.write(f"  E-value: {hsp.expect}\n\n")
        f_log.write(sep)
    
    with open(output_file, 'w') as out_handle:
        for alignment in blast_record.alignments[:top_hits]:
            if "Homo sapiens" not in alignment.title:
                out_handle.write(f"{alignment.accession}\n")
                id_list.append(alignment.accession)
    if log_file:
        f_log.write(f"Top hits saved to {output_file}\n")
        f_log.close()
    print(f"Top hits saved to {output_file}")
    print(f"BLASTn search completed for {record.id}. Results saved to {output_file}.")

    return id_list

def main():
    Entrez.email = input("Enter your email address to make use of NCBI's E-utilities: ").strip()
    if not Entrez.email or '@' not in Entrez.email:
        raise ValueError("You must provide a valid email address for NCBI Entrez access.")
    
    log_path = 'logs/similar_sequences.log'
    output_path = 'data/metadata/similar_sequences.txt'

    parser = argparse.ArgumentParser(description='Retrieve similar Hemoglobin sequences from NCBI')
    parser.add_argument('-o', '--output', type=str, default=output_path, dest='output_file',
                        help="Path to save the similar sequences' accession IDs")
    parser.add_argument('--log', type=str, default=log_path, dest='log_path',
                        help='Path to save the log file')
    parser.add_argument('-i', '--input', dest='fasta_file', 
                        help='Input fasta file containing Human Hemoglobin Gene sequence', required=False)
    parser.add_argument('-d', '--directory', dest='directory',
                        help='Location of output fasta files', required=False)

    args = parser.parse_args()
    output_path = args.output_file
    log_path = args.log_path
    output_directory = args.directory if args.directory else input("Enter the directory to save the fasta files: ").strip()
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)


    if  args.fasta_file:
        input_fasta = args.fasta_file
        f_log = open(log_path, 'w')
        f_log.write('----------------------------\n')
        f_log.write(f'<> accession file provided: {input_fasta}\n')

        #if file doesn't exist exit 1
        if not os.path.exists(input_fasta):
            f_log.write('<> file does not exist, exiting\n')
            f_log.write('----------------------------\n')
            f_log.close()
            exit(1)
        
    else:
        input_fasta = input("Enter the path to the input fasta file containing Human Hemoglobin Gene sequence: ").strip()
        if not os.path.exists(input_fasta):
            f_log.write('<> file does not exist, exiting\n')
            f_log.write('----------------------------\n')
            f_log.close()
            exit(1)
    
    f_log.write(f'<> input fasta file: {input_fasta}\n')
    f_log.write(f'<> output file: {output_path}\n')
    f_log.write(f'<> log file: {log_path}\n')
    f_log.write('----------------------------\n')

    id_list = find_similar_sequences(input_fasta, output_path, log_file=log_path)
    getting_full_fasta_records(id_list, output_directory, log_file=log_path)

if __name__ == "__main__":
    main()