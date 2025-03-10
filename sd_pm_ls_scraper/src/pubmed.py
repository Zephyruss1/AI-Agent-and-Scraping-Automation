import csv
import os

import pandas as pd
import requests
from Bio import Entrez, Medline

# Set your email for NCBI API access
Entrez.email = "your.email@example.com"

# Directory to save PDFs
PDF_DIR = "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output"
os.makedirs(PDF_DIR, exist_ok=True)

CSV_FILE = (
    "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/pubmed_articles.csv"
)


# Search for articles related to a given term
def search_pubmed(query, max_results=25500):
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


def download_pdf():
    df = pd.read_csv(CSV_FILE)
    for index, row in df.iterrows():
        pmc_url = row["PMC Full Text URL"]
        print(f"Visiting paper link {index + 1}/{len(df)}: {pmc_url}")
        if "No PMC Full Text" not in pmc_url:
            headers = {
                "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3"
            }
            response = requests.get(pmc_url, headers=headers, stream=True)

        if response.status_code == 200:
            pdf_filename = os.path.join(PDF_DIR, f"{pmc_url}.pdf")
            with open(pdf_filename, "wb") as pdf_file:
                for chunk in response.iter_content(chunk_size=1024):
                    pdf_file.write(chunk)
            print(f"PDF downloaded: {pdf_filename}")
            return pdf_filename
        else:
            print(f"No free PDF found for {pmc_url}")
            response.raise_for_status()
            return "No PDF available"


# Extract information and save to CSV
def save_articles_to_csv(records, filename=CSV_FILE):
    with open(filename, mode="w", newline="", encoding="utf-8") as file:
        writer = csv.writer(file)
        writer.writerow(
            [
                "Title",
                "Source",
                "Authors",
                "Affiliations",
                "Date",
                "Paper Link",
                "PMC Full Text URL",
                "DOI",
            ]
        )

        for record in records:
            title = record.get("TI", "No title available")
            authors = record.get("FAU", ["No authors available"])
            affiliations = record.get("AD", ["No affiliations available"])
            source = record.get("SO", "No source information")
            pmid = record.get("PMID", "No PMID available")

            # Extract year and month from "DP" (Date of Publication)
            pub_date = record.get("DP", "Unknown")

            # PubMed link
            pubmed_url = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"

            # Check for PMC ID (PubMed Central)
            pmc_url = "No PMC Full Text"
            if "PMC" in record:
                pmc_id = record["PMC"]
                pmc_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmc_id}/"

            doi = record.get("LID", "No DOI available").split()[
                0
            ]  # Some DOIs have extra info

            writer.writerow(
                [
                    title,
                    source,
                    " ".join(authors),
                    " ".join(affiliations),
                    pub_date,
                    pubmed_url,
                    pmc_url,
                    doi,
                ]
            )


if __name__ == "__main__":
    search_term = "microbial kinetics"
    pubmed_ids = search_pubmed(search_term)

    if pubmed_ids:
        pubmed_records = fetch_pubmed_details(pubmed_ids)
        save_articles_to_csv(pubmed_records)
        print("PubMed articles saved to 'pubmed_results.csv'")
    else:
        print("No articles found.")
