import subprocess
import random
import os
import pandas as pd
import numpy as np

def clean_up():
    command = ["ls *aln *fq | grep -v combined | rm"]
    subprocess.run(command,check=True)

def combine_reads(read_files, output_file):
    command = ["cat"] + read_files
    with open(output_file, 'w') as outfile:
        subprocess.run(command, stdout=outfile, check=True)

def shuffle_reads_pe(reads_file1, reads_file2, output_file1, output_file2):
    # Read in the reads
    with open(reads_file1, 'r') as f1, open(reads_file2, 'r') as f2:
        reads1 = f1.readlines()
        reads2 = f2.readlines()
    # Make sure we have the same number of reads in both files
    assert len(reads1) == len(reads2)
    # Get the number of reads (each read is 4 lines in the FASTQ format)
    num_reads = len(reads1) // 4
    # Generate a random permutation of the read indices
    permutation = np.random.permutation(num_reads)
    # Apply the permutation to the reads
    shuffled_reads1 = [reads1[i * 4:(i + 1) * 4] for i in permutation]
    shuffled_reads2 = [reads2[i * 4:(i + 1) * 4] for i in permutation]
    # Write the shuffled reads back to the files
    with open(output_file1, 'w') as f1, open(output_file2, 'w') as f2:
        for read in shuffled_reads1:
            f1.writelines(read)
        for read in shuffled_reads2:
            f2.writelines(read)

def create_simulation_config(sample_id, read1_path, read2_path):
    with open('simulation_config.txt', 'a') as f:
        f.write(f"{sample_id}\t{read1_path}\t{read2_path}\n")

def simulate_metagenome(sample_id, genome_files):
    sampled_files = random.sample(genome_files, 100)
    # Generate random relative abundances
    random_abundances = np.random.dirichlet(np.ones(100),size=1)[0]
    # Mapping for genome to its abundance (reads count) and count (total reads)
    genome_abundance = {}
    genome_count = {}
    output_reads_files1 = []
    output_reads_files2 = []
    for file, abundance in zip(sampled_files, random_abundances):
        output_prefix = os.path.basename(file)
        # Scale fold coverage by relative abundance
        fold_coverage = 100 * abundance
        # Calculate genome length
        genome_length = sum(len(line.strip()) for line in open(file) if not line.startswith('>'))
        # Calculate number of reads
        num_reads = fold_coverage * genome_length / 150
        command = ["art_illumina", "-ss", "HS25", "-i", file, "-p", "-l", "150", "-f", str(fold_coverage), "-m", "200", "-s", "10", "-o", output_prefix]
        subprocess.run(command, check=True)
        output_reads_files1.append(output_prefix + '1.fq')
        output_reads_files2.append(output_prefix + '2.fq')
        genome_abundance[file] = abundance
        genome_count[file] = num_reads
    combined_file1 = f'combined_{sample_id}_1.fq'
    combined_file2 = f'combined_{sample_id}_2.fq'
    shuffle_reads1 = f'combined_shuffled_{sample_id}_1.fq'
    shuffle_reads2 = f'combined_shuffled_{sample_id}_2.fq'
    combine_reads(output_reads_files1, combined_file1)
    combine_reads(output_reads_files2, combined_file2)
    shuffle_reads_pe(combined_file1, combined_file2,shuffle_reads1, shuffle_reads2)
    create_simulation_config(sample_id, shuffle_reads1, shuffle_reads2)
    return genome_abundance, genome_count

def simulate_multiple_metagenomes(num_samples, genome_files):
    all_abundances = []
    all_counts = []
    for i in range(num_samples):
        sample_id = f'sample_{i}'
        abundance, count = simulate_metagenome(sample_id, genome_files)
        all_abundances.append(abundance)
        all_counts.append(count)
    # Convert list of dicts to DataFrame
    df_abundances = pd.DataFrame(all_abundances).transpose()
    df_counts = pd.DataFrame(all_counts).transpose()
    # Transpose DataFrame so that viruses are rows and samples are columns
    df_abundances = df_abundances.transpose()
    df_counts = df_counts.transpose()
    # Use sample IDs as column names
    df_abundances.index = [f'sample_{i}' for i in range(num_samples)]
    df_abundances.columns = [x.split('/')[2].replace('.fna','') for x in df_abundances.columns]
    df_counts.index = [f'sample_{i}' for i in range(num_samples)]
    df_counts.columns = [x.split('/')[2].replace('.fna','') for x in df_counts.columns]
    # Write DataFrame to CSV
    df_abundances.to_csv('expected_abundances_matrix.csv', index=True)
    df_counts.to_csv('expected_counts_matrix.csv', index=True)

def main():
    with open('genome_file_paths.txt', 'r') as f:
        genome_files = [line.strip() for line in f]
    simulate_multiple_metagenomes(10, genome_files)
    #clean_up()

if __name__ == "__main__":
    main()
