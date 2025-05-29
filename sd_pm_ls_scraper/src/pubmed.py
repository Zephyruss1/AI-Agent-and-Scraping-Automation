import csv
import os
import time

from Bio import Entrez, Medline
from dotenv import load_dotenv

load_dotenv()
# Set your email for NCBI API access
Entrez.email = "your.email@example.com"
Entrez.api_key = os.getenv("NCBI_API_KEY")

# Directory to save PDFs and results - consistent with Streamlit app
OUTPUT_DIR = "/tmp/output"
os.makedirs(OUTPUT_DIR, exist_ok=True)
CSV_FILE = os.path.join(OUTPUT_DIR, "pubmed_results.csv")


def search_pubmed(query, max_results=None, date_ranges=None):
    if not date_ranges:
        date_ranges = [("1975/01/01", "2025/12/31")]

    all_ids = set()

    print(f"DEBUG: Original query: '{query}'")
    print(f"DEBUG: Date ranges to search: {date_ranges}")

    for start_date, end_date in date_ranges:
        print(f"Searching for records from {start_date} to {end_date}")

        try:
            # Debug the actual search parameters
            search_params = {
                "db": "pubmed",
                "term": query,
                "retmax": 10000,
                "mindate": start_date,
                "maxdate": end_date,
            }
            print(f"DEBUG: Search parameters: {search_params}")

            handle = Entrez.esearch(**search_params)
            record = Entrez.read(handle)
            handle.close()

            # Debug the response
            print(f"DEBUG: API Response keys: {record.keys()}")
            print(f"DEBUG: Count: {record.get('Count', 'N/A')}")
            print(f"DEBUG: RetMax: {record.get('RetMax', 'N/A')}")
            print(f"DEBUG: RetStart: {record.get('RetStart', 'N/A')}")

            chunk_ids = record["IdList"]
            print(f"Found {len(chunk_ids)} records in this date range")

            # If no results, try a broader search to test API connectivity
            if len(chunk_ids) == 0:
                print("DEBUG: No results found. Testing with a simple query...")
                test_handle = Entrez.esearch(
                    db="pubmed",
                    term="cancer",  # Simple, common term
                    retmax=5,
                )
                test_record = Entrez.read(test_handle)
                test_handle.close()
                print(
                    f"DEBUG: Test query 'cancer' returned {len(test_record['IdList'])} results",
                )

                # Try the original query without date filters
                print("DEBUG: Testing original query without date filters...")
                no_date_handle = Entrez.esearch(db="pubmed", term=query, retmax=10)
                no_date_record = Entrez.read(no_date_handle)
                no_date_handle.close()
                print(
                    f"DEBUG: Query without dates returned {len(no_date_record['IdList'])} results",
                )

            all_ids.update(chunk_ids)

        except Exception as e:
            print(
                f"ERROR: Failed to search for date range {start_date} to {end_date}: {e}",
            )
            continue

    all_ids = list(all_ids)
    print(f"Total unique records found across all date ranges: {len(all_ids)}")
    return all_ids[:max_results] if max_results else all_ids


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
                db="pubmed",
                id=batch,
                rettype="medline",
                retmode="text",
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
                    db="pubmed",
                    id=batch,
                    rettype="medline",
                    retmode="text",
                )
                records = list(Medline.parse(handle))
                handle.close()
                all_records.extend(records)
            except Exception as e:
                print(f"Failed to fetch batch after retry: {e}")

    return all_records


# Extract information and save to CSV
def save_articles_to_csv(records, _keyword: str, filename=None):
    """Save PubMed records to a CSV file
    Args:
        records (list): List of PubMed records
        _keyword (str): Search keyword
        filename (str): Output CSV filename (optional)
    Returns:
        str: Path to saved CSV file or None if no records
    """
    if filename is None:
        filename = CSV_FILE

    # Check if we have any records to save
    if not records or len(records) == 0:
        print("No records to save. Skipping CSV file creation.")
        return None

    # Ensure the directory exists
    os.makedirs(os.path.dirname(filename), exist_ok=True)

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
                "Keyword",
            ],
        )

        # Process each record and write to CSV
        record_count = 0
        for record in records:
            filtered_author_list: list = []
            try:
                title = record.get("TI", "No title available")
                authors = record.get("FAU", ["No authors available"])
                if "No authors available" not in authors:
                    for author in authors:
                        parts = author.split(",")
                        if len(parts) == 2:
                            firstname = parts[1].strip()
                            lastname = parts[0].strip()
                            filtered_author_list.append(f"{firstname} {lastname}")
                        else:
                            filtered_author_list.append(author.strip())

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
                        ", ".join(filtered_author_list),
                        " ".join(affiliations),
                        pub_date,
                        pubmed_url,
                        pmc_url,
                        doi,
                        _keyword,
                    ],
                )

                record_count += 1
                if record_count % 1000 == 0:
                    print(f"Processed {record_count} records")

            except Exception as e:
                print(f"Error processing record: {e}")
                continue

    # Only print success message if we actually wrote records
    if record_count > 0:
        print(f"Total records written to CSV: {record_count}")
        print(f"CSV file saved to: {filename}")
        return filename
    else:
        print("No valid records were processed. Removing empty CSV file.")
        # Remove the file if no records were written (only headers)
        if os.path.exists(filename):
            os.remove(filename)
        return None


def generate_yearly_ranges(start_year: int, end_year: int):
    date_ranges = []
    for year in range(start_year, end_year + 1):
        date_ranges.append((f"{year}/01/01", f"{year}/12/31"))
    return date_ranges


def main_run(link: str, proxies=None):
    """Main function to run PubMed scraping
    Args:
        link (str): PubMed search URL
        proxies (dict): Proxy configuration (optional, for future use)
    Returns:
        str: Path to saved CSV file
    """
    print(f"DEBUG: Input link: {link}")

    # Parse search term from URL
    if "&filter=years" not in link:
        if "term=" in link:
            search_term = link.split("term=")[1].replace("+", " ")
            # Handle other URL parameters that might be present
            if "&" in search_term:
                search_term = search_term.split("&")[0]
        else:
            print("ERROR: Could not extract search term from URL")
            return None
    else:
        search_term = link.split("term=")[1].split("&filter=years")[0].replace("+", " ")

    # URL decode the search term (handle %20, etc.)
    import urllib.parse

    search_term = urllib.parse.unquote(search_term)

    print(f"DEBUG: Extracted search term: '{search_term}'")

    # Parse date range from URL
    if "years." not in link:
        print("No date range specified in the link. Using default date range.")
        start_date = "2005"
        end_date = "2025"
    else:
        try:
            date_part = link.split("years.")[1]
            if "-" in date_part:
                start_date = date_part.split("-")[0]
                end_date = date_part.split("-")[1]
                # Clean up end_date if it has other parameters
                if "&" in end_date:
                    end_date = end_date.split("&")[0]
            else:
                print("ERROR: Invalid date range format in URL")
                start_date = "2005"
                end_date = "2025"
        except Exception as e:
            print(f"ERROR: Could not parse date range: {e}")
            start_date = "2005"
            end_date = "2025"

    print(f"DEBUG: Date range - start: {start_date}, end: {end_date}")

    # Validate date range
    try:
        start_year = int(start_date)
        end_year = int(end_date)
        if start_year > end_year:
            print("ERROR: Start year is greater than end year. Swapping...")
            start_year, end_year = end_year, start_year
        if start_year < 1900 or end_year > 2030:
            print("WARNING: Date range seems unusual")
    except ValueError:
        print("ERROR: Invalid year format. Using defaults.")
        start_year, end_year = 2005, 2025

    print("Final search parameters:")
    print(f"  Search term: '{search_term}'")
    print(f"  Year range: {start_year} to {end_year}")

    date_ranges = generate_yearly_ranges(start_year, end_year)
    print(f"DEBUG: Generated date ranges: {date_ranges}")

    pubmed_ids = search_pubmed(search_term, date_ranges=date_ranges)

    if pubmed_ids:
        print(f"Found {len(pubmed_ids)} total PubMed IDs")
        pubmed_records = fetch_pubmed_details_batch(pubmed_ids)
        print(f"Successfully fetched details for {len(pubmed_records)} records")

        if pubmed_records:
            saved_file = save_articles_to_csv(pubmed_records, _keyword=search_term)
            if saved_file:
                print(f"PubMed articles saved to '{saved_file}'")
                return saved_file
            else:
                print("No valid articles were saved.")
                return None
        else:
            print("No article details could be fetched.")
            return None
    else:
        print("No articles found for the search term.")

        # Suggest alternative search strategies
        print("\nDEBUG: Troubleshooting suggestions:")
        print("1. Try a simpler search term")
        print("2. Check if the search term has results on PubMed website")
        print("3. Try without date restrictions")
        print("4. Check your NCBI API key and email configuration")

        return None


if __name__ == "__main__":
    main_run("https://pubmed.ncbi.nlm.nih.gov/?term=machine+learning+AND+CFU")
