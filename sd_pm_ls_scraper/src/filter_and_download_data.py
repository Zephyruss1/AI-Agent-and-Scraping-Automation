import os
import warnings

import pandas as pd
import pycountry
from dotenv import load_dotenv

load_dotenv()

warnings.filterwarnings("ignore")

# Define base directory path
os.makedirs(
    "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/pdfs/sciencedirect",
    exist_ok=True,
)
os.makedirs(
    "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/pdfs/pubmed",
    exist_ok=True,
)
os.makedirs(
    "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/pdfs/springer",
    exist_ok=True,
)
BASE_DIR = "/root/arxiv-and-scholar-scraping"

PDF_DIR_SCIENCEDIRECT = os.path.join(
    BASE_DIR, "sd_pm_ls_scraper/output/pdfs/sciencedirect"
)
PDF_DIR_PUBMED = os.path.join(BASE_DIR, "sd_pm_ls_scraper/output/pdfs/pubmed")
PDF_DIR_SPRINGER = os.path.join(BASE_DIR, "sd_pm_ls_scraper/output/pdfs/springer")

# Load the API key and institution token
API_KEY = os.getenv("ELSEVIER_API_KEY")
INSTTOKEN = os.getenv("ELSEVIER_INSTTOKEN")

# Try to load each CSV file individually
try:
    sciencedirect_path = os.path.join(
        BASE_DIR, "sd_pm_ls_scraper/output/sciencedirect_results_1_merged.csv"
    )
    sciencedirect_csv = pd.read_csv(sciencedirect_path)
    print("ScienceDirect data loaded successfully.")
except Exception as e:
    print(f"Error loading ScienceDirect data: {e}")

try:
    pubmed_path = os.path.join(BASE_DIR, "sd_pm_ls_scraper/output/pubmed_results_1.csv")
    pubmed_csv = pd.read_csv(pubmed_path)
    print("PubMed data loaded successfully.")
except Exception as e:
    print(f"Error loading PubMed data: {e}")

try:
    springer_path = os.path.join(
        BASE_DIR, "sd_pm_ls_scraper/output/springer_results_1.csv"
    )
    springer_csv = pd.read_csv(springer_path)
    print("SpringerLink data loaded successfully.")
except Exception as e:
    print(f"Error loading SpringerLink data: {e}")


def _load_filtered_csv():
    """
    Load the filtered CSV files for PubMed, ScienceDirect, and SpringerLink.
    """
    try:
        pubmed_df = pd.read_csv(pubmed_path)
        sciencedirect_df = pd.read_csv(sciencedirect_path)
        springer_df = pd.read_csv(springer_path)
        return pubmed_df, sciencedirect_df, springer_df
    except FileNotFoundError as e:
        print(f"File not found: {e}")
        raise
    except Exception as e:
        print(f"Error loading filtered CSV files: {e}")
        raise


def reformat_datetime(*csv_files: str) -> None:
    """
    Reformat the datetime column for each CSV file.
    """
    try:
        for csv_file in csv_files:
            data = pd.read_csv(csv_file)  # Read each CSV file
            data["Date"] = pd.to_datetime(data["Date"], format="mixed", errors="coerce")
            data["Date"] = data["Date"].dt.strftime("%Y-%m-%d")
            data.sort_values(by="Date", ascending=False, inplace=True)

            # Save the reformatted CSV
            filtered_csv_path = csv_file.rsplit(".csv", 1)[0]
            data.to_csv(f"{filtered_csv_path}_tests.csv", index=False)
            print(f"Reformatted CSV saved to {filtered_csv_path}_tests.csv.")
    except Exception as e:
        print(f"Error reformatting datetime: {e}")


def download_pdf() -> str:
    """Download the PDF from the PMC Full Text URL."""
    try:
        import time

        import requests
    except NotImplementedError as e:
        print(f"Error importing requests: {e}")

    pubmed_df, sciencedirect_df, springer_df = _load_filtered_csv()
    # Download PDFs from PubMed
    count = 0
    for _index, row in pubmed_df.iterrows():
        pmc_url = row["PMC Full Text URL"]
        if "articles/" in pmc_url:
            file_name = pmc_url.split("articles/")[1].replace("/", "")
            print(
                f"    [INFO]💭 Visiting paper link {_index + 1}/{len(pubmed_df)}: {pmc_url}"
            )
            if "No PMC Full Text" not in pmc_url:
                headers = {
                    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3"
                }
                response = requests.get(pmc_url, headers=headers, stream=True)
                if response.status_code == 200:
                    # Form the filename correctly
                    pdf_filename = os.path.join(PDF_DIR_PUBMED, f"{file_name}.pdf")
                    with open(pdf_filename, "wb") as pdf_file:
                        for chunk in response.iter_content(chunk_size=1024):
                            pdf_file.write(chunk)
                    print(f"    [INFO]✅ PDF downloaded: {pdf_filename}")
                    count += 1
                elif response.status_code == 403:
                    print(f"    [ERROR] 403 Forbidden: {pmc_url}")
                    time.sleep(5)  # Wait before retrying
                else:
                    print(f"    [INFO]⚠️ No free PDF found for {pmc_url}")
                    response.raise_for_status()
        else:
            print(f"    [INFO]⚠️ Invalid PMC URL format: {pmc_url}")
    print(f"Total downloaded pdf from Pubmed: {count}")

    # Download PDFs from ScienceDirect
    count = 0
    headers = {
        "X-ELS-APIKey": API_KEY,
        "X-ELS-Insttoken": INSTTOKEN,
        "Accept": "application/pdf",
        "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.36",
    }

    BASE_URL: str = "https://api.elsevier.com/content/article/doi/"
    for _index, row in sciencedirect_df.iterrows():
        paper_doi = row["DOI_2"]
        print(f"\nProcessing DOI: {paper_doi}")
        paper_doi = str(paper_doi).strip()
        API_URL = f"{BASE_URL}{paper_doi}"

        try:
            response = requests.get(API_URL, headers=headers, allow_redirects=True)

            if response.status_code == 200:
                if response.headers.get("Content-Type") == "application/pdf":
                    pdf_filename = os.path.join(
                        PDF_DIR_SCIENCEDIRECT, f"{paper_doi.replace('/', '_')}.pdf"
                    )
                    with open(pdf_filename, "wb") as pdf_file:
                        pdf_file.write(response.content)
                    print(f"PDF downloaded successfully: {pdf_filename}")
                    count += 1
                else:
                    print(
                        f"Error: Content-Type is {response.headers.get('Content-Type')}, not PDF."
                    )
            else:
                print(f"Failed to download PDF. Status code: {response.status_code}")
                print(f"Response: {response.text}")
        except requests.exceptions.RequestException as e:
            print(f"An error occurred: {e}")

    print(f"\nTotal PDFs downloaded from ScienceDirect: {count}")

    # Download PDFs from SpringerLink
    count = 0
    for _index, row in springer_df.iterrows():
        pdf_url = row["Link"]
        if "link.springer.com" in pdf_url:
            doi_code = pdf_url.split("article/")[1] + ".pdf"
            pdf_url = f"https://link.springer.com/content/pdf/{doi_code}"
            print(
                f"    [INFO]💭 Visiting paper link {_index + 1}/{len(springer_df)}: {pdf_url}"
            )
            headers = {
                "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3"
            }
            response = requests.get(pdf_url, headers=headers, stream=True)
            if response.status_code == 200:
                # Form the filename correctly
                pdf_filename = os.path.join(
                    PDF_DIR_SPRINGER, f"{doi_code.replace('/', '_')}"
                )
                with open(pdf_filename, "wb") as pdf_file:
                    for chunk in response.iter_content(chunk_size=1024):
                        pdf_file.write(chunk)
                print(f"    [INFO]✅ PDF downloaded: {pdf_filename}")
                count += 1
            else:
                print(f"    [INFO]⚠️ No free PDF found for {pdf_url}")
                response.raise_for_status()
        else:
            print(f"    [INFO]⚠️ Invalid SpringerLink URL format: {pdf_url}")


def amount_of_mentions() -> None:
    """
    Count the number of mentions of each author in the CSV files,
    considering only unique author-DOI combinations.
    """

    # Load all three exploded CSV files
    try:
        df1 = pd.read_csv(
            "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/pubmed_results_1_tests.csv",
            parse_dates=["Date"],
        )
        df2 = pd.read_csv(
            "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/sciencedirect_results_1_merged_tests.csv",
            parse_dates=["Date"],
        )
        df3 = pd.read_csv(
            "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/springer_results_1_tests.csv",
            parse_dates=["Date"],
        )
    except FileNotFoundError as e:
        print(f"Error loading playground CSV files: {e}")
        return

    # Combine all dataframes by stacking them vertically
    combined = pd.concat([df1, df2, df3], ignore_index=True)

    # Verify required columns exist in the combined dataframe
    required_columns = ["Authors", "DOI"]
    if all(col in combined.columns for col in required_columns):
        # Identify unique author-DOI combinations
        unique_author_dois = combined.drop_duplicates(subset=["Authors", "DOI"])

        # Count author mentions based on unique author-DOI combinations
        author_counts = unique_author_dois["Authors"].value_counts().reset_index()
        author_counts.columns = ["Authors", "amount_of_mentions"]

        # Merge the counts back into the original combined dataframe
        combined = pd.merge(combined, author_counts, on="Authors", how="left")
    else:
        print(
            f"Missing required columns: {required_columns}. Available columns are: {combined.columns.tolist()}"
        )

    # Save the cleaned files
    combined.to_csv(
        "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/cleaned_total_results_1_tests.csv",
        index=False,
    )

    print("Processing complete. Cleaned files saved.")


def compare_authors() -> None:
    print("Starting compare_authors function...")
    # Load the combined CSV file
    try:
        combined = pd.read_csv(
            "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/cleaned_total_results_1_tests.csv",
            parse_dates=["Date"],
        )
    except FileNotFoundError as e:
        print(f"Error loading CSV file: {e}")
        return

    # Check if Source column exists to properly identify cross-source duplicates
    if "Source" not in combined.columns:
        print("Warning: 'Source' column not found. Adding default source information.")
        # You may need to add source information if it's not already there
        combined["Source"] = "Unknown"

    # Find authors that appear multiple times
    author_counts = combined["Authors"].value_counts()
    duplicate_authors = author_counts[author_counts >= 2].index.tolist()
    print(f"Found {len(duplicate_authors)} authors that appear multiple times.")

    # Use a more efficient approach: create a list to store filtered rows
    filtered_rows = []

    # Process each author group
    for author, group in combined.groupby("Authors"):
        if author in duplicate_authors:
            # Sort by date (newest first) and take the first (most recent)
            most_recent = group.sort_values("Date", ascending=False).iloc[0:1]
            filtered_rows.append(most_recent)
        else:
            # For non-duplicates, keep all entries
            filtered_rows.append(group)

    # Combine the filtered rows more efficiently
    filtered_data = pd.concat(filtered_rows, ignore_index=True)

    # Save the cleaned files with a more descriptive name
    output_path = "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/cleaned_total_results_1_final.csv"
    filtered_data.to_csv(output_path, index=False)

    print(f"Processing complete. Final cleaned file saved to: {output_path}")
    print(
        f"Records before deduplication: {len(combined)}, after deduplication: {len(filtered_data)}"
    )
    print("compare_authors function completed.")


def extract_university() -> None:
    """
    Extract university names from affiliations and store them in a new column.
    Handles both "University of X" and "X University" patterns.
    """
    print("Starting extract_university function...")
    try:
        combined = pd.read_csv(
            "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/cleaned_total_results_1_tests.csv",
            parse_dates=["Date"],
        )
    except FileNotFoundError as e:
        print(f"Error loading CSV file: {e}")
        return

    # Create new column if it doesn't exist
    if "University_Name" not in combined.columns:
        combined["University_Name"] = None

    for _index, row in enumerate(combined["Affiliations"]):
        if pd.notna(row) and "University" in row:
            # Split the text at "University"
            parts = row.split("University")

            if len(parts) > 1:
                # If there's content after "University" (like "University of X")
                university_part = parts[1].strip()

                # Check if the part after University is empty or just punctuation
                if not university_part or university_part[0] in ",;.":
                    # It's likely "X University" pattern
                    prefix = parts[0].strip()
                    # Get the last part before "University" (likely the university name)
                    prefix_parts = prefix.split(",")
                    university_name = prefix_parts[-1].strip() + " University"
                else:
                    # It's likely "University of X" pattern
                    university_name = "University" + university_part.split(",")[0]

                # Clean up any extra spaces and punctuation
                university_name = university_name.strip()
                combined.at[_index, "University_Name"] = university_name
                print(f"{_index}, Original: {row}")
                print(f"{_index}, Extracted: {university_name}")
                print("---" * 30)

    combined.to_csv(
        "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/all_results.csv"
    )

    print("extract_university function completed.")


def extract_department() -> None:
    print("Starting extract_department function...")
    try:
        combined = pd.read_csv(
            "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/cleaned_total_results_1_tests.csv",
            parse_dates=["Date"],
        )
    except FileNotFoundError as e:
        print(f"Error loading CSV file: {e}")
        return

    for _index, row in enumerate(combined["Affiliations"]):
        if pd.notna(row) and "University" in row:
            department = row.split("University")[0]
            result = department.split(",")[0]
            full_department = result
            combined.at[_index, "Department"] = full_department
            print(_index, full_department)

        combined.to_csv(
            "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/all_results.csv"
        )

    print("Department extraction completed.")


def extract_countries() -> None:
    print("Starting extract_countries function...")
    try:
        combined = pd.read_csv(
            "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/cleaned_total_results_1_tests.csv",
            parse_dates=["Date"],
        )
    except FileNotFoundError as e:
        print(f"Error loading CSV file: {e}")
        return

    for _index, row in enumerate(combined["Affiliations"]):
        if pd.notna(row):
            for country in pycountry.countries:
                if country.name in row:
                    combined.at[_index, "Country"] = country.name
                    print(_index, country.name)
                    break
    else:
        print(_index, "No country found")

    combined.to_csv(
        "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/all_results.csv"
    )
    print("Country extraction completed.")


def filter_authors() -> None:
    """Filtering authors based on the amount of mentions and date."""
    print("Starting filter_authors function...")
    try:
        combined = pd.read_csv(
            "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/cleaned_total_results_1_tests.csv",
            parse_dates=["Date"],
        )
    except FileNotFoundError as e:
        print(f"Error loading CSV file: {e}")
        return

    if "Date" in combined.columns and not pd.api.types.is_datetime64_any_dtype(
        combined["Date"]
    ):
        combined["Date"] = pd.to_datetime(combined["Date"], errors="coerce")

    if "amount_of_mentions" in combined.columns:
        # Sort by Authors and then by amount_of_mentions (descending)
        combined_sorted = combined.sort_values(
            by=["Authors", "amount_of_mentions", "Date"], ascending=[True, False, False]
        )

        # Drop duplicates keeping the one with highest amount_of_mentions (first after sorting)
        combined_deduped = combined_sorted.drop_duplicates(
            subset=["Authors"], keep="first"
        )

        print(f"Original shape: {combined.shape}")
        print(f"Shape after keeping highest mentions: {combined_deduped.shape}")

        # Save the result
        combined_deduped.to_csv(
            "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/all_results.csv",
            index=False,
        )
    print("Filtered authors saved to all_results.csv.")


def explode_csv(*csv_files):
    """
    Explode the CSV files to get the authors in separate rows.
    """
    try:
        for csv_file in csv_files:
            df = pd.read_csv(f"{csv_file}.csv")

            # Split authors on commas and handle spaces
            df["Authors"] = df["Authors"].str.split(r"\s*,\s*")
            df = df.explode("Authors").reset_index(drop=True)
            df["Authors"] = df["Authors"].str.strip()  # Clean whitespace

            df.to_csv(f"{csv_file}.csv", index=False)
            print(f"Exploded CSV saved to {csv_file}.")
    except FileNotFoundError:
        print(f"File not found: {csv_file}")
    except Exception as e:
        print(f"Error exploding CSV files: {e}")


def main():
    print("Starting main execution...")
    # download_pdf()

    csv_file_paths = [
        sciencedirect_path,
        pubmed_path,
        springer_path,
    ]

    print("Reformatting datetime...")
    reformat_datetime(*csv_file_paths)

    # Fix the file paths to match what was actually generated
    csv_file_paths = [
        "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/pubmed_results_1_tests",
        "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/sciencedirect_results_1_merged_tests",
        "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/springer_results_1_tests",
    ]

    # Check if files exist before proceeding
    for path in csv_file_paths:
        if not os.path.exists(f"{path}.csv"):
            print(f"File not found: {path}")
            return  # Stop execution if any file is missing

    print("Exploding CSV files...")
    explode_csv(*csv_file_paths)

    print("Calculating author mentions...")
    amount_of_mentions()

    print("Comparing authors across sources...")
    compare_authors()
    print("Extracting university names...")
    extract_university()
    print("Extracting departments...")
    extract_department()
    print("Extracting countries...")
    extract_countries()
    print("Filtering authors based on mentions and date...")
    filter_authors()
    print("Script execution completed successfully.")


if __name__ == "__main__":
    main()
    print("Program terminated.")
