#!/bin/env python3

from Bio import SeqIO
import sys, subprocess
import csv
import matplotlib.pyplot as plt
from collections import Counter



IUPAC_codes = 'ATCGNRYKMSWBDHVatcgnrykmswbdhv-'

def fasta_check(fasta_file):
    print("[+] Check FASTA file...")
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if len(record.id) == 0:
                sys.exit("[!] Error: Empty ID found in the FASTA file.");
            
            print(f"\t[+] Checking Chromosome {record.id} ({len(record.seq):,}): ", end="", flush=True)
            print(f"{record.seq[:25]}...{record.seq[len(record.seq) - 3:]}: ", end="", flush=True)
            
            for base in str(record.seq):
                if base not in IUPAC_codes:
                    sys.exit(f"\n[!] Error: Invalid base '{base}' found in the FASTA file.")
            if len(record.seq) < 1:
                sys.exit("\n[!] Error: Sequence length < 1 found in the FASTA file.")

            print("Done")

    except Exception as e:
        sys.exit(f"[!] Error: {e}")
    print("[-] FASTA file check completed successfully!")


def predict_genes(fasta_path, output_gff):
    print("[+] Predicting genes with Augustus...")
    cmd = f"augustus --species=arabidopsis {fasta_path} > {output_gff}"
    subprocess.run(cmd, shell=True, check=True)
    print("[-] Gene prediction completed!")


def extract_sequences(fasta_path, gff_path, cds_out, prot_out):
    print("[+] Extracting CDS and protein sequences...")
    cmd = f"gffread {gff_path} -g {fasta_path} -x {cds_out} -y {prot_out}"
    subprocess.run(cmd, shell=True, check=True)
    print("[-] Sequence extraction completed!")


def annotate_with_blast(prot_file, blast_out):
    print("[+] Annotating proteins with BLAST...")
    cmd = f"blastp -query {prot_file} -db uniprot_sprot_db -out {blast_out} -evalue 1e-5 -outfmt 6"
    subprocess.run(cmd, shell=True, check=True)
    print("[-] BLAST annotation completed!")



def filter_stress_proteins(blast_out, keywords, filtered_out):
    print("[+] Filtering stress-related proteins...")
    with open(blast_out, "r") as infile, open(filtered_out, "w") as outfile:
        for line in infile:
            if any(kw.lower() in line.lower() for kw in keywords):
                outfile.write(line)
    print("[-] Filtering completed! Stress-related proteins.")

def generate_csv(filtered_blast, csv_out):
    print("[+] Generating CSV summary of stress-related proteins...")
    with open(filtered_blast, "r") as infile, open(csv_out, "w", newline="") as outcsv:
        writer = csv.writer(outcsv)
        writer.writerow(["Gene ID", "Function", "Stress Type"])
        for line in infile:
            cols = line.strip().split("\t")
            if len(cols) >= 2:
                gene_id = cols[0]
                function = cols[1]
                writer.writerow([gene_id, function])
    print("[-] CSV summary generated successfully!")


def generate_graph(csv_file, genome_name, keywords):
    counts = Counter()

    with open(csv_file, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            for key in keywords:
                if key in row["Function"].upper():
                    counts[key] += 1

    plt.figure(figsize=(10, 6))
    plt.bar(counts.keys(), counts.values(), color="skyblue")
    plt.title("Gene Counts by Stress-Related Keyword")
    plt.ylabel("Number of Genes")
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(genome_name + ".stress_genes.png")
    plt.show()



def get_description(function_field):
    function_field = function_field.upper()

    if "HSP" in function_field:
        return "HSP – Heat shock proteins, help protect cells from high temperatures"
    elif "DREB" in function_field:
        return "DREB – Help plants respond to drought and water deficiency"
    elif "LEA" in function_field:
        return "LEA – Protect plant cells from severe dehydration"
    elif "HSF" in function_field:
        return "HSF – Activate heat shock response genes"
    elif "SOS" in function_field:
        return "SOS – Controls salt tolerance via sodium ion export"
    elif "NHX" in function_field:
        return "NHX – Regulate sodium/potassium balance under salt stress"
    elif "RD29" in function_field:
        return "RD29 – Activated during drought and salinity stress"
    elif "AREB" in function_field:
        return "AREB – Regulate ABA hormone response during stress"
    else:
        return "–"

def generate_html(csv_file, genome_name):
    print("[+] Generating HTML summary of stress-related genes...")

    html = """<!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <title>Gene Stress Annotation</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 20px; }
            table { border-collapse: collapse; width: 100%; }
            th, td { border: 1px solid #ccc; padding: 8px; text-align: left; }
            th { background-color: #f2f2f2; }
            tr:hover { background-color: #f9f9f9; }
        </style>
    </head>
    <body>
        <h2>Gene Stress Annotation – Human Readable</h2>
        <table>
            <tr>
                <th>Gene ID</th>
                <th>Matched Protein</th>
                <th>Human Description</th>
                <th>Stress Type</th>
            </tr>
    """

    with open(csv_file, newline="") as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            desc = get_description(row['Function'])
            html += f"""        <tr>
                <td>{row['Gene ID']}</td>
                <td>{row['Function']}</td>
                <td>{desc}</td>
                <td>{row['Stress Type']}</td>
            </tr>\n"""

    html += """    </table>
    </body>
    </html>
    """

    with open(genome_name + ".html", "w") as f:
        f.write(html)
    print("[-] HTML summary generated successfully!")




if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(f"[!] Usage: {sys.argv[0]} <fasta_file>")
    genome_name = sys.argv[1][:sys.argv[1].rfind(".")]
    fasta = sys.argv[1]
    gff3 = genome_name + ".gff3"
    cds = genome_name + ".cds.fasta"
    protein = genome_name + ".proteins.fasta"
    blast_output = genome_name + "_blast_results.tsv"
    filtered_output = genome_name + "_filtered_stress.tsv"
    summary_csv = genome_name + "_stress_summary.csv"


    stress_keywords = ["HSP", "DREB", "LEA", "HSF", "SOS", "NHX", "RD29", "AREB"]

    fasta_check(fasta)
    predict_genes(fasta, gff3)
    extract_sequences(fasta, gff3, cds, protein)
    annotate_with_blast(protein, blast_output)
    filter_stress_proteins(blast_output, stress_keywords, filtered_output)
    generate_csv(filtered_output, summary_csv)
    generate_graph(summary_csv, genome_name, stress_keywords)
    generate_html(summary_csv, genome_name)

    print("[+] All steps completed successfully!")
    print("[+] Check the generated files for results.")
    print('\033[32m', end='')
    print("[+] Summary CSV:", summary_csv)
    print("[+] Summary HTML:", genome_name + ".html")
    print("[+] Stress Gene Graph:", genome_name + ".stress_genes.png")
