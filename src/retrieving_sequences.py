#!~/miniconda3/envs/bioinformher/bin/python3

'''
This script aims to:
    1- search the nucleotide database in NCBI for HBB gene of humans,
    2- retrieve the sequences in FASTA format.
'''

import os
from Bio import Entrez, SeqIO
import argparse

def retrieving_accession_IDs(species, output_path=None, log_file=None):
    '''
    Takes a species list, searches the NCBI nucleotide database for Human HBB gene sequence,
    and writes the accession IDs to a file.

    params:
    -------
    - species: str, species name to search for
    - output_path: str, path to the output file where accession IDs will be saved 
    - log_file: str, path to the log file where search results will be logged

    return:
    -------
    - ids: list of accession IDs found for the given sequences
    '''

    sep = '\n-------------------------------------------\n'

    if log_file:
        f_log = open(log_file, 'a')
        f_log.write(sep)
    
    id_list = []

    # Will loop through each species and search for Hemoglobin sequences
    
    search = f'("{species}"[Organism] AND HBB[Gene])'
    handle = Entrez.esearch(db='nucleotide', term=search, retmax=1)
    if log_file:
        f_log.write(f'searching for {search} in nucleotide database\n')
    record = Entrez.read(handle)
    handle.close()
    id_list.extend(record['IdList'])
    print(f'there are {len(id_list)} records')
    if log_file:
        f_log.write(f'there are {len(id_list)} records\n')
    
    if not os.path.exists("data/metadata"):
        os.mkdir("data/metadata")
    if not output_path:
        output_path = 'data/metadata/accession_IDs.txt'
    with open(output_path, 'w') as f:
        for id in id_list:
            f.write(id + '\n')
    print(f'Accession ID saved to {output_path}')
    if log_file:
        f_log.write(f'Accession IDs saved to {output_path}\n')
        f_log.write(sep)
        f_log.write('<> ids extracted successfully!')
        f_log.close()
    
    return id_list

def getting_full_fasta_records(id_list, destination_directory, log_file=None):
    '''
    Takes a list of accession IDs and downloads the fasta file corresponding to each
    from NCBI's nucleotide database.
    These files are saved in data/genomes by their id number.
    
    params:
    -------
    - id_list: list, list of accession IDs to retrieve
    - destination_directory: str, directory where the fasta files will be saved
    - log_file: str, path to the log file where search results will be logged
    
    return:
    -------
    - None
    '''
    sep = '\n-------------------------------------------\n'
    if log_file:
        f_log = open(log_file, 'a')
        f_log.write(sep)
    
    if not os.path.exists(destination_directory):
        os.makedirs(destination_directory)

    # Will loop through each id in the list and retrieve the fasta file
    for id in id_list:
        if os.path.exists(f'{destination_directory}/{id}.fasta'):
            if log_file:
                f_log.write(f'File {id}.fasta already exists, skipping download.\n')
            continue
        handle = Entrez.efetch(db='nucleotide', id=id, rettype='fasta', retmode='text')
        record = SeqIO.read(handle, 'fasta')
        with open(f'{destination_directory}/{id}.fasta', 'w') as fa:
            SeqIO.write(record, fa, 'fasta')
            if log_file:
                f_log.write(f'Saved {id}.fasta to {destination_directory}/\n')
        handle.close()
    
    if log_file:
        f_log.write(sep)
        f_log.write('<> fasta records extracted successfully!')
        f_log.close()

def main():
    '''
    Main function to execute the script.
    Requires the user to input their email address for NCBI's E-utilities.
    Searches for Hemoglobin sequences in specified species and retrieves them in FASTA format.
    Checks if an input file with accession IDs is provided, otherwise retrieves them from NCBI.
    Allows the user to run the script with command line arguments to specify output paths for accession IDs and logs.
    
    Args:
    ----
    - output: str, path to save the accession IDs (default: 'data/metadata/accession_IDs.txt')
    - log: str, path to save the log file (default: 'data/logs/retrieving_sequences.log')
    - input: str, input file name containing accession IDs (optional)
    '''
    
    Entrez.email = input("Enter your email address to make use of NCBI's E-utilities: ")
    species = "homo sapiens"
    if not os.path.exists("logs"):
        os.mkdir("logs")
    if not os.path.exists("data"):
        os.mkdir("data")
    log_path = 'logs/retrieving_sequences.log'


    parser = argparse.ArgumentParser(description='Retrieve Hemoglobin sequences from NCBI')
    parser.add_argument('--output', type=str, default='data/metadata/accession_IDs.txt',
                        help='Path to save the accession IDs')
    parser.add_argument('--log', type=str, default=log_path,
                        help='Path to save the log file')
    parser.add_argument('-i', '--input', dest='accession_file', 
                        help='Input file name containing accession IDs', required=False)
    parser.add_argument('-d', dest='destination_directory',
                        help='Directory to save output file' )

    args = parser.parse_args()
    
    if args.accession_file:
        f_log=open(log_path, 'w')
        f_log.write('----------------------------\n')
        f_log.write(f'<> accession file provided: {args.accession_file}\n')

        #if file doesn't exist exit 1
        if not os.path.exists(args.accession_file):
            f_log.write('<> file does not exist, exiting\n')
            f_log.write('----------------------------\n')
            f_log.close()
            exit(1)

        with open(args.accession_file,'r') as f:
            accession_list=f.readlines()
            accession_list=[acc.strip() for acc in accession_list] #stripping the \n

    else:
        f_log=open(log_path, 'w'); f_log.close()
        accession_list = retrieving_accession_IDs(species=species,log_file=log_path)
    if not args.destination_directory:
        destination_directory = 'data/Humans_gene'
    else:
        destination_directory = args.destination_directory
    getting_full_fasta_records(accession_list, destination_directory, log_file=log_path)


if __name__ == "__main__":
    main()
