import csv
import os
import time

from Bio import Entrez, Medline

# Set your email for NCBI API access
Entrez.email = "your.email@example.com"

# Directory to save PDFs
PDF_DIR = "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output"
os.makedirs(PDF_DIR, exist_ok=True)
CSV_FILE = "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/pubmed_results.csv"


def search_pubmed(query, max_results=25500, date_ranges: list = None):
    if not date_ranges:
        date_ranges = [("1900/01/01", "2050/12/31")]
    """Search PubMed
    Args:
        query (str): Search query
        max_results (int): Maximum number of results to return
        date_ranges (list): List of tuples with start and end dates for the search
    Returns:
        list: List of PubMed IDs
    """
    all_ids = []

    # Define date ranges to split the search

    for start_date, end_date in date_ranges:
        print(f"Searching for records from {start_date} to {end_date}")
        handle = Entrez.esearch(
            db="pubmed",
            term=query,
            retmax=10000,
            mindate=start_date,
            maxdate=end_date,
        )
        record = Entrez.read(handle)
        handle.close()

        chunk_ids = record["IdList"]
        print(f"Found {len(chunk_ids)} records in this date range")
        all_ids.extend(chunk_ids)

    print(f"Total records found across all date ranges: {len(all_ids)}")
    return all_ids[:max_results]  # Limit to max_results if needed


# Fetch article details using PubMed IDs in batches
def fetch_pubmed_details_batch(id_list, batch_size=200):
    """Fetch PubMed details in batches to avoid API limitations
    Args:
        id_list (list): List of PubMed IDs
        batch_size (int): Number of IDs to fetch per batch
    Returns:
        list: List of article records
    """
    total_ids = len(id_list)
    all_records = []

    for start in range(0, total_ids, batch_size):
        end = min(start + batch_size, total_ids)
        batch = id_list[start:end]

        print(f"Fetching details for records {start + 1}-{end} of {total_ids}")

        try:
            handle = Entrez.efetch(
                db="pubmed", id=batch, rettype="medline", retmode="text"
            )
            records = list(Medline.parse(handle))
            handle.close()
            all_records.extend(records)

        except Exception as e:
            print(f"Error fetching batch {start + 1}-{end}: {e}")
            # Wait a bit longer before retrying
            time.sleep(5)
            try:
                handle = Entrez.efetch(
                    db="pubmed", id=batch, rettype="medline", retmode="text"
                )
                records = list(Medline.parse(handle))
                handle.close()
                all_records.extend(records)
            except Exception as e:
                print(f"Failed to fetch batch after retry: {e}")

    return all_records


# Extract information and save to CSV
def save_articles_to_csv(records, filename=CSV_FILE):
    """Save PubMed records to a CSV file
    Args:
        records (list): List of PubMed records
        filename (str): Output CSV filename
    Returns:
        csv
    """
    # Create/open the CSV file and write the header
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

        # Process each record and write to CSV
        record_count = 0
        for record in records:
            try:
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

                # Handle DOI - check different fields where DOI might be stored
                doi = "No DOI available"
                if "LID" in record:
                    doi = record["LID"].split()[0]
                elif "AID" in record:
                    for aid in record["AID"]:
                        if "[doi]" in aid:
                            doi = aid.split()[0]
                            break

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

                record_count += 1
                if record_count % 1000 == 0:
                    print(f"Processed {record_count} records")

            except Exception as e:
                print(f"Error processing record: {e}")
                continue

    print(f"Total records written to CSV: {record_count}")


if __name__ == "__main__":
    # Example usage: search_pubmed_with_date_ranges("bacterial growth AND antibiotics", [("2010/01/01", "2020/12/31")])
    search_term = "microbial kinetics AND CFU"
    date_ranges = [
        ("1900/01/01", "2025/12/31"),
        ("2026/01/01", "2050/12/31"),
    ]
    pubmed_ids = search_pubmed(search_term, date_ranges=date_ranges)
    print(f"Searching with term '{search_term}'")
    if pubmed_ids:
        print(f"Found {len(pubmed_ids)} total PubMed IDs")
        pubmed_records = fetch_pubmed_details_batch(pubmed_ids)
        print(f"Successfully fetched details for {len(pubmed_records)} records")
        save_articles_to_csv(pubmed_records)
        print(f"PubMed articles saved to '{CSV_FILE}'")
    else:
        print("No articles found.")
