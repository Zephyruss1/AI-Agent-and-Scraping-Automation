import csv
import os

import requests
from Bio import Entrez, Medline

# Set your email for NCBI API access
Entrez.email = "your.email@example.com"

# Directory to save PDFs
PDF_DIR = "pubmed_pdfs"
os.makedirs(PDF_DIR, exist_ok=True)


# Search for articles related to a given term
def search_pubmed(query, max_results=18865):
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    handle.close()
    return record["IdList"]


# Fetch article details using PubMed IDs
def fetch_pubmed_details(id_list):
    handle = Entrez.efetch(db="pubmed", id=id_list, rettype="medline", retmode="text")
    records = Medline.parse(handle)
    records = list(records)
    handle.close()
    return records


def download_pdf(pmc_id, title):
    pdf_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmc_id}/pdf/"
    headers = {
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3"
    }
    response = requests.get(pdf_url, headers=headers, stream=True)

    if response.status_code == 200:
        pdf_filename = os.path.join(PDF_DIR, f"{pmc_id}.pdf")
        with open(pdf_filename, "wb") as pdf_file:
            for chunk in response.iter_content(chunk_size=1024):
                pdf_file.write(chunk)
        print(f"PDF downloaded: {pdf_filename}")
        return pdf_filename
    else:
        print(f"No free PDF found for {pmc_id}")
        response.raise_for_status()
        return "No PDF available"


# Extract information and save to CSV
def save_articles_to_csv(records, filename="pubmed_articles.csv"):
    with open(filename, mode="w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(
            [
                "Title",
                "Authors",
                "Source",
                "Year",
                "Month",
                "PubMed URL",
                "PMC Full Text URL",
                "DOI",
                "PDF Path",
            ]
        )

        for record in records:
            title = record.get("TI", "No title available")
            authors = ", ".join(record.get("AU", ["No authors listed"]))
            source = record.get("SO", "No source information")
            pmid = record.get("PMID", "No PMID available")

            # Extract year and month from "DP" (Date of Publication)
            pub_date = record.get("DP", "Unknown")
            pub_parts = pub_date.split()  # Split to extract year & month
            year = pub_parts[0] if pub_parts else "Unknown"
            month = pub_parts[1] if len(pub_parts) > 1 else "Unknown"

            # PubMed link
            pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

            # Check for PMC ID (PubMed Central)
            pmc_url, pdf_path = "No PMC Full Text", "No PDF available"
            # if "PMC" in record:
            #     pmc_id = record["PMC"]
            #     pmc_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmc_id}/"
            #     print("Downloading PDF for:", pmc_id)
            #     pdf_path = download_pdf(pmc_id, title)

            doi = record.get("LID", "No DOI available").split()[
                0
            ]  # Some DOIs have extra info

            writer.writerow(
                [
                    title,
                    authors,
                    source,
                    year,
                    month,
                    pubmed_url,
                    pmc_url,
                    doi,
                    pdf_path,
                ]
            )


if __name__ == "__main__":
    search_term = "microbial kinetics"
    pubmed_ids = search_pubmed(search_term)

    if pubmed_ids:
        pubmed_records = fetch_pubmed_details(pubmed_ids)
        save_articles_to_csv(pubmed_records)
        print("PubMed articles saved to 'pubmed_articles.csv'")
    else:
        print("No articles found.")
