import pandas as pd
import argparse

def load_vcf(filename):
    """
    Load VCF Files into pandas DataFrame

    ARGS:
    filename(str): The path to the VCF File to be loaded

    Returns:
    pd.Dataframe: A Dataframe containing the VCF file data.

    Raises:
    ValueError: If the VCF file does not contain a #CHROM header.

    """
    with open(filename, "r") as f:
        lines = f.readlines()
    
    """
    Finding the index of the #CHROM header line,
    Extracting data starting from the #CHROM header line, and
    Create a DataFrame from the VCF data
    """
        
   
    chrom_index = [i for i, line in enumerate(lines) if line.strip().startswith("#CHROM")]
    
    if not chrom_index:
        raise ValueError(f"The VCF file {filename} does not contain a #CHROM header.")
    
    data = lines[chrom_index[0]:]  
    header = data[0].strip().split("\t")
    informations = [d.strip().split("\t") for d in data[1:]]
    
    vcf = pd.DataFrame(informations, columns=header)
    
    return vcf

def processVcf(file1, file2):
    """
    Main function to load two VCF files, find their intersection based on
    #CHROM, POS, REF, and ALT columns, and print the first five columns of the intersection.
    
    Args:
        file1 (str): The path to the first VCF file.
        file2 (str): The path to the second VCF file.
    """
    vcf1 = load_vcf(file1)
    vcf2 = load_vcf(file2)
    
    # Finding intersection of both dataframes on #CHROM, POS, REF, ALT columns
    common_columns = ["#CHROM", "POS", "REF", "ALT"]
    intersection = pd.merge(vcf1, vcf2, on=common_columns)
    
    # Selecting the first five columns
    first_five_columns = intersection.columns[:5]
    intersection_vcf = intersection[first_five_columns]
    
    #std output
    print("Intersection of both VCF files")
    print(intersection_vcf)

if __name__ == "__main__":
    """

    Setting up argument parser to accept two VCF file paths from the command line,
    Parsing the command-line arguments, and
    Executing the main function with the provided file paths
    """
    parser = argparse.ArgumentParser(description="Load two VCF files using python3 {CMD Line- python3 intersection_vcf.py test_a.vcf test_b.vcf}")
    parser.add_argument("file1", help="First VCF file")
    parser.add_argument("file2", help="Second VCF file")
    
    args = parser.parse_args()
    
    processVcf(args.file1, args.file2)