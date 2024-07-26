import csv
from pyfaidx import Fasta

def NNN_counts(fasta_file):
    all_results = []  # List to store results for all chromosomes
    fasta = Fasta(fasta_file)
    print(fasta[0])
    output_file = f"{fasta_file}.csv"  # Create the output CSV filename based on the input filename

    for chrom_name in fasta.keys():
        results = []  # List to store results for the current chromosome
        seq = str(fasta[chrom_name])
        nnn_start = None
        ctt_count = 0

        for i in range(len(seq)):
            if seq[i] == 'N':
                if nnn_start is None:
                    nnn_start = i
                if i+6 < len(seq) and seq[i+1:i+7] == 'CTTAAG':
                    ctt_count += 1
            elif nnn_start is not None:
                if ctt_count > 0:
                    nnn_end = i
                    nnn_length = nnn_end - nnn_start + 1
                    results.append([chrom_name, nnn_start, nnn_end, nnn_length, ctt_count])
                    ctt_count = 0
                nnn_start = None

        # Append the results for the current chromosome to the list of all results
        all_results.extend(results)

    # Write the results to a CSV file
    with open(output_file, 'w', newline='') as csv_file:
        csv_writer = csv.writer(csv_file)
        csv_writer.writerow(["Chromosome", "Start", "End", "Length", "CTTAAG Count"])
        csv_writer.writerows(all_results)
    return all_results

# Call the function with your FASTA file
NNN_counts("C:\\Users\\Harini Mudradi\\Documents\\NGS_LAB_STUFF\\lg_mapped_lja_CV4m_contigs_to_CV2_molecules.fasta\\lg_mapped_lja_CV4m_contigs_to_CV2_molecules.fasta")







































# import csv
# from pyfaidx import Fasta

# def NNN_counts(fasta_file):
#     results = []  # List to store results
#     fasta = Fasta(fasta_file)
#     print(fasta_file)
#     output_file = f"{fasta_file}.csv"  # Create the output CSV filename based on the input filename

#     for chrom_name in fasta.keys():
#         seq = str(fasta[chrom_name])
#         nnn_start = None
#         ctt_count = 0
#         print(chrom_name)
#         print(seq)

#         for i in range(len(seq)):
#             if seq[i] == 'N':
#                 if nnn_start is None:
#                     nnn_start = i
#                     print("nnn_start:", nnn_start)
#             elif nnn_start is not None:
#                 if seq[i:i+6] == 'CTTAAG':
#                     ctt_count += 1
#                     #print("ctt_count:",ctt_count)
#                 elif seq[i] == 'N':
#                     if ctt_count > 0:
#                         nnn_end = i - 1
#                         print("nnn_end:",nnn_end)
#                         nnn_length = nnn_end - nnn_start + 1 - (6 * ctt_count)
#                         print("nnn_length:" , nnn_length)
#                         if nnn_length > 0:
#                             results.append([chrom_name, nnn_start, nnn_end, nnn_length])
#                             print("results:",results)

#                         ctt_count = 0
#                         nnn_start = i
#                         print("nnn_start:",nnn_start)

#     # Write the results to a CSV file
#     with open(output_file, 'w', newline='') as csv_file:
#         csv_writer = csv.writer(csv_file)
#         print("csv_file_name:",csv_file)
#         csv_writer.writerow(["Chromosome", "Start", "End", "Length"])
#         csv_writer.writerows(results)
#         print("csv_writet:",csv_writer)

#     return results

# # List of input FASTA files
# input_files = [ 
#         r"C:\Users\Harini Mudradi\Documents\NGS_LAB_STUFF\lg_mapped_lja_CV4m_contigs_to_CV2_molecules.fasta\lg_mapped_lja_CV4m_contigs_to_CV2_molecules.fasta",
#          r"C:\Users\Harini Mudradi\Documents\NGS_LAB_STUFF\lg_mapped_lja_mc3_contigs_to_mc1_molecules.fasta\lg_mapped_lja_mc3_contigs_to_mc1_molecules.fasta",
#          r"C:\Users\Harini Mudradi\Documents\NGS_LAB_STUFF\lg_mapped_lja_MZ4.3_to_UMD2a.fasta\lg_mapped_lja_MZ4.3_to_UMD2a.fasta",
#          r"C:\Users\Harini Mudradi\Documents\NGS_LAB_STUFF\lg_mapped_lja_mz4_contigs_to_sample9_molecules.fasta\lg_mapped_lja_mz4_contigs_to_sample9_molecules.fasta"]

# # Process each input FASTA file and create corresponding output CSV files
# for input_file in input_files:
#     results = NNN_counts(input_file)
#     print(f"Processed {input_file}. Saved results to {input_file}.csv")
































#     # Write results to a CSV file
#     with open(output_csv, 'w', newline='') as csvfile:
#         csv_writer = csv.writer(csvfile)
#         # Write a header row
#         csv_writer.writerow(['Chrom_Name', 'NNN_Start', 'NNN_End', 'NNN_Length'])
#         # Write the data rows
#         csv_writer.writerows(results)
#         print("Results written to CSV file.")

# # Example usage:
# fasta_file =  r'C:\Users\Harini Mudradi\Documents\NGS_LAB_STUFF\lg_mapped_lja_CV4m_contigs_to_CV2_molecules.fasta\lg_mapped_lja_CV4m_contigs_to_CV2_molecules.fasta'
# # output_csv = 'output.csv'
# NNN_counts(fasta_file)


# def write_results_to_csv(results, output_csv):
#     with open(output_csv, mode="w", newline="") as csv_file:
#         fieldnames = ["Chromosome", "Start", "End", "Length"]
#         writer = csv.writer(csv_file)
#         writer.writerow(fieldnames)

#         for result in results:
#             writer.writerow(result)

# if __name__ == "__main__":
#     fasta_files = [
#         r"C:\Users\Harini Mudradi\Documents\NGS_LAB_STUFF\lg_mapped_lja_CV4m_contigs_to_CV2_molecules.fasta\lg_mapped_lja_CV4m_contigs_to_CV2_molecules.fasta",
#         r"C:\Users\Harini Mudradi\Documents\NGS_LAB_STUFF\lg_mapped_lja_mc3_contigs_to_mc1_molecules.fasta\lg_mapped_lja_mc3_contigs_to_mc1_molecules.fasta",
#         r"C:\Users\Harini Mudradi\Documents\NGS_LAB_STUFF\lg_mapped_lja_MZ4.3_to_UMD2a.fasta\lg_mapped_lja_MZ4.3_to_UMD2a.fasta",
#         r"C:\Users\Harini Mudradi\Documents\NGS_LAB_STUFF\lg_mapped_lja_mz4_contigs_to_sample9_molecules.fasta\lg_mapped_lja_mz4_contigs_to_sample9_molecules.fasta"
#     ]

#     for fasta_file in fasta_files:
#         results = NNN_counts(fasta_file)
#         output_csv = f"{fasta_file.split('.')[0]}_output.csv"
#         write_results_to_csv(results, output_csv)

#         print(f"Results for {fasta_file} have been written to {output_csv}")