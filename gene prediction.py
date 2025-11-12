import pandas as pd
import os
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3
import time
import requests
from time import sleep
from Bio import SeqIO
from io import StringIO
import json
from functools import reduce
import csv
from pathlib import Path
from upsetplot import from_contents, UpSet
import seaborn as sns
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
# Define base folder for the genomes
base_folder = '/Users/hala/Documents/Genomes'
# List of country folders to loop through
countries = ['China', 'Saudi_Arabia', 'Pakistan', 'Nepal', 'India']
def read_files():
    countries = ['Pakistan', 'India', 'China', 'Nepal', 'Saudi_Arabia']
    base_folder = '/Users/hala/Documents/Genomes'
    output_folder = os.path.join(base_folder, 'comparison_results_pathways')
    os.makedirs(output_folder, exist_ok=True)

    for country in countries:
        metabolism_file = os.path.join(base_folder, country, f"{country}_metabolism.tsv")

        if os.path.exists(metabolism_file):
            print(f"Processing {country}...")

            # Read the metabolism file
            df = pd.read_csv(metabolism_file, sep='\t', header=None)

            # Add column names
            df.columns = ['Gene_ID', 'Type', 'Length', 'Gene_Name', 'EC_Number', 'COG', 'Product']

            # Filter for rows with EC numbers
            df_filtered = df[df['EC_Number'].notnull() & (df['EC_Number'] != '')]

            # Keep only selected columns
            final_df = df_filtered[['Gene_Name', 'EC_Number', 'Product']]

            # Save cleaned pathways
            save_path = os.path.join(output_folder, f"{country}_pathways.csv")
            final_df.to_csv(save_path, index=False)

            print(f"✅ Saved {country} pathways to {save_path}")
        else:
            print(f"❗ {country} metabolism file not found. Skipping.")

    print("\n✅ All countries processed!")
def add_headers_to_metabolism_files(base_folder, countries):
    """
    Adds standard headers to each country's metabolism.tsv file.

    Parameters:
    - base_folder: str – Path to the 'Genomes' directory.
    - countries: list – List of country folder names.
    """
    column_names = ["locus_tag", "ftype", "length_bp", "gene", "EC_number", "COG", "Product"]

    for country in countries:
        file_path = os.path.join(base_folder, country, f"{country}_metabolism.tsv")

        try:
            # Read file without header
            df = pd.read_csv(file_path, sep="\t", header=None, names=column_names)

            # Save file with headers
            df.to_csv(file_path, sep="\t", index=False)
            print(f"✅ Header added to: {file_path}")

        except FileNotFoundError:
            print(f"❌ File not found: {file_path}")
        except Exception as e:
            print(f"⚠️ Error in {file_path}: {e}")
def count_and_plot_ec_numbers():
    """
    Count and plot common and unique EC numbers
    across countries' pathway files.
    """
    data_folder = '/Users/hala/Documents/Genomes/comparison_results_pathways'
    countries = ['Pakistan', 'India', 'China', 'Nepal', 'Saudi_Arabia']
    country_pathways = {}
    # Load each country's pathways
    for country in countries:
        path = os.path.join(data_folder, f"{country}_pathways.csv")
        if os.path.exists(path):
            df = pd.read_csv(path)
            country_pathways[country] = set(df['EC_Number'].dropna().unique())
        else:
            print(f"❗ {country} pathways file not found!")

    # ------
    # 1. Find COMMON EC Numbers
    # ------
    common_ec_numbers = set.intersection(*country_pathways.values())
    print(f"\n🌟 Common EC Numbers across all countries ({len(common_ec_numbers)} pathways):")
    print(common_ec_numbers)

    # Save common EC numbers
    common_path = os.path.join(data_folder, 'common_pathways.csv')
    pd.DataFrame(list(common_ec_numbers), columns=['Common_EC_Number']).to_csv(common_path, index=False)

    # ------
    # 2. Find UNIQUE EC Numbers for each country
    # ------
    for country in countries:
        others = set.union(*(v for k, v in country_pathways.items() if k != country))
        unique_ec = country_pathways[country] - others

        print(f"\n🎯 Unique pathways for {country} ({len(unique_ec)} pathways):")
        print(unique_ec)

        unique_path = os.path.join(data_folder, f"{country}_unique_pathways.csv")
        pd.DataFrame(list(unique_ec), columns=['Unique_EC_Number']).to_csv(unique_path, index=False)
def fetch_kegg_pathway_info(ec_number):
    """
    Given an EC number, fetch KEGG pathway information using the KEGG API.
    """
    base_url = "http://rest.kegg.jp/link/pathway/ec:"
    try:
        response = requests.get(base_url + ec_number)
        if response.status_code == 200 and response.text:
            lines = response.text.strip().split("\n")
            pathways = [line.split("\t")[1] for line in lines]
            return pathways
        else:
            return []
    except Exception as e:
        print(f"❗ Error fetching KEGG info for {ec_number}: {e}")
        return []
def map_ec_to_pathways():
    """
    Map common EC numbers to KEGG pathways.
    """
    data_folder = '/Users/hala/Documents/Genomes/comparison_results_pathways'
    common_path = os.path.join(data_folder, 'common_pathways.csv')

    if not os.path.exists(common_path):
        print("❗ Common pathways file not found!")
        return

    common_df = pd.read_csv(common_path)
    all_mappings = []

    for ec in common_df['Common_EC_Number']:
        print(f"Fetching pathways for EC {ec}...")
        pathways = fetch_kegg_pathway_info(ec)
        time.sleep(1)  # Be gentle with KEGG servers!

        if pathways:
            for pw in pathways:
                all_mappings.append({'EC_Number': ec, 'KEGG_Pathway': pw})
        else:
            all_mappings.append({'EC_Number': ec, 'KEGG_Pathway': 'No pathway found'})

    mapped_df = pd.DataFrame(all_mappings)
    mapped_path = os.path.join(data_folder, 'common_ec_kegg_mapping.csv')
    mapped_df.to_csv(mapped_path, index=False)

    print(f"\n✅ Mapped EC numbers to KEGG pathways. Saved to {mapped_path}")
    return mapped_df
def categorize_kegg_pathways(mapped_df):
    """
    Categorize mapped KEGG pathways into broader categories
    like Metabolism, Genetic Information Processing, etc.
    """
    # Simple keyword matching (can be expanded)
    category_mapping = {
        'metabolism': 'Metabolism',
        'genetic': 'Genetic Information Processing',
        'environmental': 'Environmental Information Processing',
        'cellular': 'Cellular Processes',
        'organismal': 'Organismal Systems',
        'human': 'Human Diseases',
        'drug': 'Drug Development'
    }

    def assign_category(pathway_id):
        if pd.isna(pathway_id) or pathway_id == 'No pathway found':
            return 'Unknown'
        try:
            url = f"http://rest.kegg.jp/get/{pathway_id}"
            response = requests.get(url)
            if response.status_code == 200:
                data = response.text.lower()
                for keyword, category in category_mapping.items():
                    if keyword in data:
                        return category
                return 'Other'
            else:
                return 'Unknown'
        except:
            return 'Unknown'

    mapped_df['Category'] = mapped_df['KEGG_Pathway'].apply(assign_category)

    # Save categorized results
    data_folder = '/Users/hala/Documents/Genomes/comparison_results_pathways'
    categorized_path = os.path.join(data_folder, 'categorized_common_pathways.csv')
    mapped_df.to_csv(categorized_path, index=False)

    print(f"\n✅ Categorized KEGG pathways. Saved to {categorized_path}")
def process_country_genes():
    # Important: define combined_data **inside** the function
    combined_data = []

    for country in countries:
        country_folder = os.path.join(base_folder, country)

        # Gene list file
        gene_list_file = os.path.join(country_folder, f'{country}_gene_list.tsv')
        if os.path.exists(gene_list_file):
            gene_df = pd.read_csv(gene_list_file, sep='\t')
            gene_df['Country'] = country
            combined_data.append(gene_df)
            print(f'Processed gene list for {country}')

        # Metabolism file
        metabolism_file = os.path.join(country_folder, f'{country}_metabolism.tsv')
        if os.path.exists(metabolism_file):
            metabolism_df = pd.read_csv(metabolism_file, sep='\t')
            metabolism_df['Country'] = country
            combined_data.append(metabolism_df)
            print(f'Processed metabolism data for {country}')

    # Combine all DataFrames
    final_df = pd.concat(combined_data, ignore_index=True)

    # Save the combined data
    output_file = os.path.join(base_folder, 'combined_data.tsv')
    final_df.to_csv(output_file, sep='\t', index=False)
    print(f'Combined data saved to {output_file}')
def create_gene_list(base_folder, countries):
    for country in countries:
        file_path = os.path.join(base_folder, country, f'{country}_gene_list.tsv')
        output_path = os.path.join(base_folder, country, 'gene_list.tsv')

        if not os.path.exists(file_path):
            print(f"Warning: {file_path} does not exist. Skipping {country}.")
            continue

        df = pd.read_csv(file_path, sep='\t')

        if 'gene' not in df.columns:
            print(f"Error: 'gene' column not found in {file_path}")
            continue

        gene_list = df[['gene']].dropna()  # Only keep rows where gene is not empty
        gene_list = gene_list[gene_list['gene'].str.strip() != '']  # Remove blank spaces

        gene_list.to_csv(output_path, sep='\t', index=False)
        print(f"Created gene_list.tsv for {country}")
def find_consensus_and_unique_genes(base_folder, countries):
    country_genes = {}

    # Step 1: Load all genes for each country
    for country in countries:
        file_path = os.path.join(base_folder, country, f'{country}_gene_list.tsv')
        if os.path.exists(file_path):
            df = pd.read_csv(file_path, sep='\t')
            gene_set = set(df['gene'].dropna())  # Remove NaN and make set
            country_genes[country] = gene_set
        else:
            print(f"Warning: {file_path} does not exist. Skipping {country}.")

    if not country_genes:
        print("No gene data found. Exiting.")
        return

    # Step 2: Find consensus genes (intersection of all sets)
    consensus_genes = set.intersection(*country_genes.values())

    # Step 3: Find unique genes for each country
    unique_genes_per_country = {}
    for country, genes in country_genes.items():
        other_countries = set()
        for other_country, other_genes in country_genes.items():
            if other_country != country:
                other_countries.update(other_genes)
        unique_genes = genes - other_countries
        unique_genes_per_country[country] = unique_genes

    # Step 4: Save consensus genes
    consensus_file = os.path.join(base_folder, 'consensus_genes.tsv')
    pd.DataFrame({'gene': list(consensus_genes)}).to_csv(consensus_file, sep='\t', index=False)
    print(f"Consensus genes saved to {consensus_file}")

    # Step 5: Save unique genes for each country
    for country, unique_genes in unique_genes_per_country.items():
        unique_file = os.path.join(base_folder, f'unique_genes_{country}.tsv')
        pd.DataFrame({'gene': list(unique_genes)}).to_csv(unique_file, sep='\t', index=False)
        print(f"Unique genes for {country} saved to {unique_file}")
def process_gene_file(country):
    gene_file = os.path.join(base_folder, country, f'{country}_gene_list.tsv')
    if os.path.exists(gene_file):
        # Load the gene list data into a pandas DataFrame
        gene_data = pd.read_csv(gene_file, sep='\t')
        return gene_data
    else:
        print(f"Warning: {gene_file} does not exist. Skipping {country}.")
        return None
def process_all_countries():
    all_gene_data = []
    for country in countries:
        country_gene_data = process_gene_file(country)
        if country_gene_data is not None:
            all_gene_data.append(country_gene_data)
            print(f"Processed gene list for {country}")

    if all_gene_data:
        # Concatenate all gene data
        combined_gene_data = pd.concat(all_gene_data, ignore_index=True)
        # Save the combined data to a new file
        combined_gene_data.to_csv(os.path.join(base_folder, 'combined_gene_list.tsv'), sep='\t', index=False)
        print("All gene data combined and saved to combined_gene_list.tsv.")
    else:
        print("No gene data to combine.")
if __name__ == '__main__':
    process_all_countries()
def extract_amr_genes(fasta_path, aro_json_path, output_path):
    """
    Extract AMR genes from CARD files using Biopython and aro.json.
    """
    # Load aro.json (list of dictionaries)
    with open(aro_json_path) as f:
        aro_data = json.load(f)

    # Create mapping from ARO accession to gene name
    accession_to_gene = {
        entry["accession"]: entry["name"]
        for entry in aro_data
        if "accession" in entry and "name" in entry
    }

    # Collect gene names based on ARO accessions found in FASTA headers
    gene_names = set()
    for record in SeqIO.parse(fasta_path, "fasta"):
        header = record.description  # Full FASTA header
        for accession in accession_to_gene:
            if accession in header:
                gene_names.add(accession_to_gene[accession])
                break  # Stop after first match for efficiency

    # Save to file
    with open(output_path, "w") as f:
        for gene in sorted(gene_names):
            f.write(gene + "\n")

    print(f"✅ {len(gene_names)} AMR genes written to '{output_path}'")
fasta_path = "/Users/hala/Documents/CARD/protein_fasta_protein_homolog_model.fasta"
aro_json_path = "/Users/hala/Documents/CARD/aro.json"
output_path = "/Users/hala/Documents/CARD/amr_gene_list.txt"
base_folder = "/Users/hala/Documents/Genomes"
amr_gene_file = "/Users/hala/Documents/CARD/amr_gene_list.txt"
def match_amr_genes_to_country_gene_lists(base_folder, amr_gene_file):
    # Load AMR gene list
    with open(amr_gene_file) as f:
        amr_genes = set(line.strip() for line in f if line.strip())

    print(f"✅ Loaded {len(amr_genes)} AMR genes from CARD.")

    # Loop through each country folder inside Genomes/
    for country_folder in os.listdir(base_folder):
        country_path = os.path.join(base_folder, country_folder)
        if not os.path.isdir(country_path):
            continue  # Skip files, only enter directories

        gene_list_filename = f"{country_folder}_gene_list.tsv"
        gene_list_path = os.path.join(country_path, gene_list_filename)

        if not os.path.exists(gene_list_path):
            print(f"⚠️ Gene list file not found: {gene_list_path}")
            continue

        # Load the gene list
        df = pd.read_csv(gene_list_path, sep="\t")

        if "gene" not in df.columns:
            print(f"⚠️ No 'gene' column found in {gene_list_filename}")
            continue

        # ✅ Filter rows where gene name is in the AMR gene list
        matched = df[df["gene"].isin(amr_genes)]

        # Save matched genes to a new file
        output_file = os.path.join(country_path, f"{country_folder}_amr_genes.tsv")
        matched.to_csv(output_file, sep="\t", index=False)

        print(f"✅ {len(matched)} AMR genes matched for {country_folder} and saved to '{output_file}'")
def extract_consensus_amr_genes(genomes_path, countries, output_file):
    gene_sets = []

    for country in countries:
        amr_file = os.path.join(genomes_path, country, f"{country}_amr_genes.tsv")
        if os.path.exists(amr_file):
            df = pd.read_csv(amr_file, sep="\t")
            if "gene" in df.columns:
                genes = set(df["gene"].dropna().str.strip().str.lower())
                gene_sets.append(genes)
            else:
                print(f"⚠️ 'gene' column not found in {amr_file}")
        else:
            print(f"⚠️ File not found: {amr_file}")

    if not gene_sets:
        print("❌ No AMR gene sets were collected. Please check file paths and column names.")
        return

    # Get intersection (consensus genes)
    consensus_genes = set.intersection(*gene_sets)

    # Save result
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    pd.DataFrame(sorted(consensus_genes), columns=["Consensus_AMR_Genes"]).to_csv(output_file, sep="\t", index=False)
    print(f"✅ Consensus AMR genes saved to: {output_file}")
genomes_path = "/Users/hala/Documents/Genomes"
countries = ["China", "India", "Nepal", "Pakistan", "Saudi_Arabia"]
output_file = os.path.join(genomes_path, "gene_comparison", "consensus_amr_genes.tsv")
def find_unique_amr_genes(genomes_base_path, output_csv_path):
    """
    Reads AMR gene files from each country's folder in the given base path,
    identifies unique AMR genes for each country, and saves the result as a CSV.

    Parameters:
    - genomes_base_path (str): Path to the main Genomes folder.
    - output_csv_path (str): Path to save the resulting CSV file.
    """

    country_genes = {}

    for country in os.listdir(genomes_base_path):
        country_folder = os.path.join(genomes_base_path, country)
        gene_file = os.path.join(country_folder, f"{country}_amr_genes.tsv")

        if os.path.isfile(gene_file):
            try:
                df = pd.read_csv(gene_file, sep="\t")

                if "gene" in df.columns:
                    genes = set(df["gene"].dropna().str.strip().str.lower())
                    country_genes[country] = genes
                else:
                    print(f"⚠️ Column 'gene' not found in {gene_file}")
            except Exception as e:
                print(f"❌ Error reading {gene_file}: {e}")
        else:
            print(f"❌ File not found: {gene_file}")

    # Compute unique genes for each country
    output_rows = []
    for country, genes in country_genes.items():
        other_genes = set().union(*(g for c, g in country_genes.items() if c != country))
        unique_genes = genes - other_genes

        for gene in unique_genes:
            output_rows.append({"Country": country, "Unique_AMR_Gene": gene})

    # Save to CSV
    output_df = pd.DataFrame(output_rows)
    if not output_df.empty:
        output_df.to_csv(output_csv_path, index=False)
        print(f"✅ Unique AMR genes saved to: {output_csv_path}")
    else:
        print("⚠️ No unique AMR genes found.")
genomes_folder = "/Users/hala/Documents/Genomes"
output_file = os.path.join(genomes_folder, "unique_amr_genes_per_country.csv")
def map_amr_genes_to_pathways(base_folder, pathways_folder):
    """
    Map AMR genes to pathways using country-specific files.
    AMR genes are in: /Genomes/[Country]/[Country]_amr_genes.tsv
    Pathways are in: /Genomes/comparison_results_pathways/[Country]_pathways.csv
    """
    for country in os.listdir(base_folder):
        country_path = os.path.join(base_folder, country)
        if not os.path.isdir(country_path):
            continue

        amr_genes_file = os.path.join(country_path, f"{country}_amr_genes.tsv")
        pathways_file = os.path.join(pathways_folder, f"{country}_pathways.csv")

        if not os.path.exists(amr_genes_file):
            print(f"⚠️ AMR genes file not found: {amr_genes_file}")
            continue
        if not os.path.exists(pathways_file):
            print(f"⚠️ Pathways file not found: {pathways_file}")
            continue

        amr_df = pd.read_csv(amr_genes_file, sep="\t")
        if "gene" not in amr_df.columns:
            print(f"⚠️ Column 'gene' not found in {amr_genes_file}")
            continue

        pathways_df = pd.read_csv(pathways_file)
        if "Gene_Name" not in pathways_df.columns:
            print(f"⚠️ Column 'Gene_Name' not found in {pathways_file}")
            continue

        # Merge on gene names
        merged_df = pd.merge(amr_df, pathways_df, left_on="gene", right_on="Gene_Name", how="left")

        output_file = os.path.join(country_path, f"{country}_amr_pathways.tsv")
        merged_df.to_csv(output_file, sep="\t", index=False)

        print(f"✅ Mapped AMR genes to pathways for {country} and saved to '{output_file}'")
base_folder = "/Users/hala/Documents/Genomes"
pathways_folder = "/Users/hala/Documents/Genomes/comparison_results_pathways"
def compare_amr_pathways(base_folder, countries):
    output_folder = os.path.join(base_folder, 'comparison_results_pathways')
    os.makedirs(output_folder, exist_ok=True)
    pathway_sets = {}
    print("🔄 Collecting AMR pathways from countries...\n")
    for country in countries:
        file_path = os.path.join(base_folder, country, f'{country}_amr_pathways.tsv')
        if os.path.exists(file_path):
            df = pd.read_csv(file_path, sep='\t')
            if 'EC_number' in df.columns:
                ec_numbers = df['EC_number'].dropna().unique().tolist()
            else:
                ec_numbers = []

            # If EC numbers missing, fallback to Product
            if not ec_numbers and 'Product' in df.columns:
                ec_numbers = df['Product'].dropna().unique().tolist()

            pathway_sets[country] = set(ec_numbers)
            print(f"✅ Loaded {len(ec_numbers)} pathways for {country}")
        else:
            print(f"⚠️ File not found: {file_path}")

    if not pathway_sets:
        print("❌ No valid country data found.")
        return

    all_sets = list(pathway_sets.values())
    common_pathways = set.intersection(*all_sets)
    all_pathways = set.union(*all_sets)

    # Save common pathways
    common_df = pd.DataFrame({'Common_Pathways': sorted(common_pathways)})
    common_df.to_csv(os.path.join(output_folder, 'common_pathways.csv'), index=False)
    print("\n✅ Saved common pathways to common_pathways.csv")

    # Save unique pathways for each country
    for country, pathways in pathway_sets.items():
        unique = pathways - common_pathways
        unique_df = pd.DataFrame({'Unique_Pathways': sorted(unique)})
        unique_df.to_csv(os.path.join(output_folder, f'{country}_AMR_unique_pathways.csv'), index=False)
        print(f"✅ Saved unique pathways for {country} to {country}_AMR_unique_pathways.csv")

    print("\n🎉 Done! All results saved in:", output_folder)
base_folder = '/Users/hala/Documents/Genomes'
countries = ['Pakistan', 'India', 'Nepal', 'China', 'Saudi_Arabia']
def query_kegg_api_for_ec(ec_number):
    """
    Queries the KEGG API using an EC number and returns pathway mappings.
    """
    url = f"http://rest.kegg.jp/link/pathway/ec:{ec_number}"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            lines = response.text.strip().split("\n")
            results = []
            for line in lines:
                parts = line.split("\t")
                if len(parts) == 2:
                    kegg_id = parts[1].replace("path:", "")
                    pathway_name = get_kegg_pathway_name(kegg_id)
                    results.append((ec_number, kegg_id, pathway_name))
            return results
        else:
            return []
    except Exception as e:
        print(f"Error fetching EC {ec_number}: {e}")
        return []
def get_kegg_pathway_name(kegg_id):
    """
    Gets the pathway name for a given KEGG pathway ID.
    """
    url = f"http://rest.kegg.jp/get/{kegg_id}"
    try:
        response = requests.get(url)
        if response.status_code == 200:
            for line in response.text.split("\n"):
                if line.startswith("NAME"):
                    return line.replace("NAME", "").strip()
        return "Unknown Pathway"
    except:
        return "Unknown Pathway"
def map_ec_numbers_to_kegg(input_folder, output_folder):
    """
    Maps EC numbers from CSV files to KEGG pathways and saves results with proper headers.
    """
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    for file in os.listdir(input_folder):
        if file.endswith(".csv"):
            file_path = os.path.join(input_folder, file)
            output_path = os.path.join(output_folder, file.replace(".csv", "_KEGG.tsv"))

            try:
                df = pd.read_csv(file_path)
                if "EC_number" not in df.columns:
                    print(f"⚠️ Skipped (no EC_number): {file}")
                    continue

                ec_list = df["EC_number"].dropna().astype(str).unique()
                all_results = []

                for ec in ec_list:
                    mappings = query_kegg_api_for_ec(ec)
                    all_results.extend(mappings)
                    time.sleep(0.5)  # Respect KEGG server rate limits

                result_df = pd.DataFrame(all_results, columns=["EC_number", "KEGG_Pathway_ID", "Pathway_Name"])

                # Add indicator column
                if "common" in file.lower():
                    result_df["Pathway_Type"] = "Common Pathway"
                elif "unique" in file.lower():
                    result_df["Pathway_Type"] = "Unique Pathway"
                else:
                    result_df["Pathway_Type"] = "Unknown Type"

                result_df.to_csv(output_path, sep="\t", index=False)
                print(f"✅ Saved: {output_path}")

            except Exception as e:
                print(f"❌ Error with {file}: {e}")
# === Run the function ===
input_folder = "/Users/hala/Documents/Genomes/comparison_results_pathways"
output_folder = "/Users/hala/Documents/Genomes/comparison_results_pathways/kegg_mapped"
output_dir = "/Users/hala/Documents/Genomes/amr_kegg_pathways"
def extract_mutation_linked_genes(aro_json_path, output_file):
    # Check if file exists
    if not os.path.exists(aro_json_path):
        print(f"❌ File not found: {aro_json_path}")
        return

    # Keywords to search for mutation-linked mechanisms
    keywords = ['mutation', 'variant', 'substitution', 'snps', 'snp']

    # Load the JSON data
    with open(aro_json_path, 'r') as file:
        try:
            data = json.load(file)
        except json.JSONDecodeError:
            print("❌ Failed to parse JSON. Check if file is properly formatted.")
            return

    # Filter entries
    mutation_genes = []
    for entry in data:
        description = entry.get('description', '').lower()
        if any(keyword in description for keyword in keywords):
            mutation_genes.append({
                'Gene Name': entry.get('name', 'N/A'),
                'ARO Accession': entry.get('accession', 'N/A'),
                'Description': entry.get('description', 'N/A')
            })

    # Save results
    if mutation_genes:
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=mutation_genes[0].keys())
            writer.writeheader()
            writer.writerows(mutation_genes)
        print(f"✅ {len(mutation_genes)} mutation-linked AMR genes saved to: {output_file}")
    else:
        print("⚠️ No mutation-linked genes found in the aro.json file.")
aro_json_path = '/Users/hala/Documents/CARD/aro.json'
output_file = '/Users/hala/Documents/CARD/mutation_linked_genes.csv'
def map_mutation_genes_to_countries(card_path, genomes_path, output_file):
    mutation_file = os.path.join(card_path, "mutation_linked_genes.csv")

    if not os.path.exists(mutation_file):
        raise FileNotFoundError(f"Missing: {mutation_file}")

    mutation_genes_df = pd.read_csv(mutation_file)
    mutation_genes = set(mutation_genes_df["Gene Name"].dropna().str.strip())

    countries = ["China", "India", "Nepal", "Pakistan", "Saudi_Arabia"]
    results = []

    for country in countries:
        gene_list_path = os.path.join(genomes_path, country, f"{country}_gene_list.tsv")
        if os.path.exists(gene_list_path):
            country_df = pd.read_csv(gene_list_path, sep="\t")
            country_genes = set(country_df["gene"].dropna().str.strip())
            shared_genes = mutation_genes & country_genes

            for gene in shared_genes:
                results.append({"Country": country, "Mutation_Gene": gene})

    results_df = pd.DataFrame(results)
    results_df.to_csv(output_file, index=False)
    print(f"✅ Saved: {output_file}")
    return results_df
card_path = "/Users/hala/Documents/CARD"
genomes_path = "/Users/hala/Documents/Genomes"
output_file = os.path.join(card_path, "country_mutation_genes.csv")
def combine_aro_data(card_path, output_file):
    # Load files
    aro_index = pd.read_csv(os.path.join(card_path, "aro_index.tsv"), sep="\t")
    aro_categories = pd.read_csv(os.path.join(card_path, "aro_categories.tsv"), sep="\t")

    # Strip column names (important!)
    aro_index.columns = aro_index.columns.str.strip()
    aro_categories.columns = aro_categories.columns.str.strip()

    # Merge based on ARO Accession
    merged = pd.merge(aro_index, aro_categories, on="ARO Accession", how="left")

    # Show actual merged column names for debug
    print("🧾 Merged columns:", list(merged.columns))

    # Select useful fields if available
    expected_cols = [
        "ARO Name_x",
        "ARO Accession",
        "CARD Short Name",
        "Drug Class",
        "Resistance Mechanism",
        "ARO Category"
    ]

    # Check which columns exist
    missing = [col for col in expected_cols if col not in merged.columns]
    if missing:
        print("⚠️ Missing columns:", missing)
    else:
        # Filter and save only if all columns exist
        result = merged[expected_cols].drop_duplicates()
        result.to_csv(os.path.join(card_path, output_file), index=False)
        print("✅ File saved to:", os.path.join(card_path, output_file))
combine_aro_data("/Users/hala/Documents/CARD", "amr_gene_antibiotics.csv")
def plot_amr_gene_class_distribution(amr_class_counts, title="AMR Genes per Antibiotic Class", output_path="amr_class_distribution.png"):
    """
    Plots a bar chart of AMR gene counts per antibiotic class.

    Parameters:
    - amr_class_counts (dict): Dictionary with antibiotic classes as keys and gene counts as values.
    - title (str): Title for the chart.
    - output_path (str): File path to save the image.
    """
    # Sort the dictionary by values (descending)
    sorted_items = sorted(amr_class_counts.items(), key=lambda x: x[1], reverse=True)
    classes, counts = zip(*sorted_items)
    # Plotting
    plt.figure(figsize=(12, 6))
    bars = plt.bar(classes, counts, color='skyblue', edgecolor='black')
    plt.xticks(rotation=45, ha='right')
    plt.ylabel("Number of AMR Genes")
    plt.title(title)
    plt.tight_layout()
    # Add value labels on bars
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, yval + 0.2, int(yval), ha='center', fontsize=9)
    # Save and show
    plt.savefig(output_path, dpi=300)
    plt.show()
    print(f"📊 Bar chart saved to: {output_path}")
amr_class_counts = {
    "Macrolide": 6,
    "Tetracycline": 4,
    "Fluoroquinolone": 3,
    "Phenicol": 3,
    "Aminoglycoside": 2,
    "Glycylcycline": 2,
    "Penicillin β-lactam": 2,
    "Aminocoumarin": 1,
    "Antibacterial Free Fatty Acids": 1,
    "Diaminopyrimidine": 1,
    "Disinfectants": 1,
    "Nitroimidazole": 1,
    "Peptide": 1,
    "Rifamycin": 1,
    "Streptogramin": 1,
    "Streptogramin B": 1
}

plot_amr_gene_class_distribution(amr_class_counts)
def match_genes_to_mutation_descriptions(gene_folder, mutation_file, output_file):
    # Load mutation-linked genes
    mutation_df = pd.read_csv(mutation_file, sep=",")  # assuming CSV with commas
    mutation_df["description_lower"] = mutation_df["Description"].astype(str).str.lower()

    results = []

    # Traverse country folders
    for country in os.listdir(gene_folder):
        country_path = os.path.join(gene_folder, country)
        gene_file = os.path.join(country_path, f"{country}_gene_list.tsv")

        if os.path.isdir(country_path) and os.path.exists(gene_file):
            try:
                # Load genes with header "gene"
                gene_df = pd.read_csv(gene_file, sep="\t")
                gene_list = gene_df["gene"].dropna().str.lower().unique()

                # Search each gene in the mutation descriptions
                for gene in gene_list:
                    matches = mutation_df[mutation_df["description_lower"].str.contains(gene, na=False)]
                    for _, row in matches.iterrows():
                        results.append({
                            "Country": country,
                            "Matched Gene": gene,
                            "ARO Accession": row["ARO Accession"],
                            "Gene Name (CARD)": row["Gene Name"],
                            "Description": row["Description"]
                        })

            except Exception as e:
                print(f"❌ Error reading {gene_file}: {e}")

    # Save the output
    if results:
        result_df = pd.DataFrame(results)
        result_df.to_csv(output_file, sep="\t", index=False)
        print(f"✅ Mutation-linked gene matches saved to: {output_file}")
    else:
        print("⚠️ No matches found.")
# === Usage ===
gene_folder = "/Users/hala/Documents/Genomes"
mutation_csv = "/Users/hala/Documents/CARD/mutation_linked_genes.csv"
output_csv = "/Users/hala/Documents/CARD/mutation_gene_matches.tsv"
def match_genes_to_mutation_descriptions(gene_folder, mutation_file, output_file):
    """
    Matches AMR genes from each country's gene list to mutation-linked resistance descriptions in CARD data.

    Parameters:
        gene_folder (str): Path to the folder containing country folders with *_gene_list.tsv files.
        mutation_file (str): CSV file from CARD with columns: Gene Name, ARO Accession, Description.
        output_file (str): Path to save the output matches as a TSV.
    """

    # Load mutation-linked CARD genes
    mutation_df = pd.read_csv(mutation_file)
    mutation_df["description_lower"] = mutation_df["Description"].astype(str).str.lower()

    results = []

    for country in os.listdir(gene_folder):
        country_path = os.path.join(gene_folder, country)
        gene_file = os.path.join(country_path, f"{country}_gene_list.tsv")

        if os.path.isdir(country_path) and os.path.exists(gene_file):
            try:
                gene_df = pd.read_csv(gene_file, sep="\t")
                gene_list = gene_df['gene'].dropna().str.lower().unique()

                for gene in gene_list:
                    # Use word boundary matching for accurate results
                    pattern = fr'\b{re.escape(gene)}\b'
                    matches = mutation_df[mutation_df["description_lower"].str.contains(pattern, regex=True, na=False)]

                    for _, row in matches.iterrows():
                        results.append({
                            "Country": country,
                            "Matched Gene": gene,
                            "ARO Accession": row["ARO Accession"],
                            "Gene Name (CARD)": row["Gene Name"],
                            "Description": row["Description"]
                        })

            except Exception as e:
                print(f"❌ Error processing {gene_file}: {e}")

    # Save results
    if results:
        result_df = pd.DataFrame(results)
        result_df.to_csv(output_file, sep="\t", index=False)
        print(f"✅ Found {len(result_df)} mutation-linked gene matches. Results saved to:\n{output_file}")
        print("📊 Summary:")
        print(result_df.groupby("Country")["Matched Gene"].count())
    else:
        print("⚠️ No matches found.")
# === Paths ===
gene_folder = "/Users/hala/Documents/Genomes"
mutation_csv = "/Users/hala/Documents/CARD/mutation_linked_genes.csv"
output_csv = "/Users/hala/Documents/CARD/mutation_gene_matches_refined.tsv"
def extract_consensus_mutation_genes(input_file, output_file):
    """
    Extracts consensus mutation-linked genes present in all countries
    from a mutation match file and saves them to an output file.
    """
    try:
        # Load the mutation gene match data
        df = pd.read_csv(input_file, sep="\t")

        # Standardize gene and country names
        df["Matched Gene"] = df["Matched Gene"].str.strip().str.lower()
        df["Country"] = df["Country"].str.strip()

        # Group by country and collect genes
        country_gene_map = df.groupby("Country")["Matched Gene"].apply(set)

        # Compute intersection across all countries
        if not country_gene_map.empty:
            consensus_genes = set.intersection(*country_gene_map)
            print(f"✅ Found {len(consensus_genes)} consensus mutation-linked genes across all countries.")

            # Save to output file
            pd.DataFrame({"Consensus_Genes": sorted(consensus_genes)}).to_csv(output_file, sep="\t", index=False)
            print(f"📁 Saved to: {output_file}")
        else:
            print("⚠️ No gene data found in the input file.")

    except Exception as e:
        print(f"❌ Error: {e}")
# === Example Usage ===
input_file = "/Users/hala/Documents/CARD/mutation_gene_matches_refined.tsv"
output_file = "/Users/hala/Documents/CARD/consensus_mutation_genes.tsv"
data = {
    "Category": [
        "Consensus Genes",
        "Unique Genes (India)",
        "Unique Genes (Nepal)",
        "Unique Genes (Saudi Arabia)",
        "Unique Genes (China)",
        "Unique Genes (Pakistan)"
    ],
    "Count": [1697, 239, 47, 35, 24, 18]
}
df = pd.DataFrame(data)
# === Plot ===
plt.figure(figsize=(10, 6))
sns.set_style("whitegrid")
palette = ["#2E86AB" if "Consensus" in cat else "#F39C12" for cat in df["Category"]]
sns.barplot(data=df, x="Category", y="Count", palette=palette)
plt.title("Distribution of Consensus and Unique Genes Across Countries", fontsize=14)
plt.xlabel("")
plt.ylabel("Number of Genes")
plt.xticks(rotation=30, ha='right')
plt.tight_layout()
# === Save or Show ===
plt.savefig("/Users/hala/Documents/Genomes/genes_txt/gene_distribution_barplot.png", dpi=300)
plt.show()









