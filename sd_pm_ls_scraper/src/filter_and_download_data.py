import os
import warnings

import pandas as pd
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
        BASE_DIR, "sd_pm_ls_scraper/output/sciencedirect_results.csv"
    )
    sciencedirect_csv = pd.read_csv(sciencedirect_path)
    print("ScienceDirect data loaded successfully.")
except Exception as e:
    print(f"Error loading ScienceDirect data: {e}")

try:
    pubmed_path = os.path.join(BASE_DIR, "sd_pm_ls_scraper/output/pubmed_results.csv")
    pubmed_csv = pd.read_csv(pubmed_path)
    print("PubMed data loaded successfully.")
except Exception as e:
    print(f"Error loading PubMed data: {e}")

try:
    springer_path = os.path.join(
        BASE_DIR, "sd_pm_ls_scraper/output/springer_results.csv"
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


def compare_authors():
    # Load all three CSV files
    df1 = pd.read_csv(
        "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/pubmed_results_tests.csv",
        parse_dates=["Date"],
    )
    df2 = pd.read_csv(
        "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/sciencedirect_results_tests.csv",
        parse_dates=["Date"],
    )
    df3 = pd.read_csv(
        "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/springer_results_tests.csv",
        parse_dates=["Date"],
    )

    # Add source identifier to each dataframe (for tracking)
    df1["source"] = "pubmed_results"
    df2["source"] = "sciencedirect_results"
    df3["source"] = "springer_results"

    # Combine all dataframes
    combined = pd.concat([df1, df2, df3])

    # Find authors that appear in at least two sources
    author_counts = combined["Authors"].value_counts()
    duplicate_authors = author_counts[author_counts >= 2].index.tolist()

    # Process duplicates - keep only the most recent entry
    filtered_data = []

    # Group by author and process each group
    for author, group in combined.groupby("Authors"):
        if author in duplicate_authors:
            # Sort by date (newest first) and take the first (most recent)
            most_recent = group.sort_values("Date", ascending=False).iloc[0]
            filtered_data.append(most_recent)
        else:
            # For non-duplicates, keep all entries
            filtered_data.extend(group.to_dict("records"))

    # Create new dataframe with filtered data
    result_df = pd.DataFrame(filtered_data)

    # Split back into original files if needed
    result_df1 = result_df[result_df["source"] == "pubmed_results"].drop(
        "source", axis=1
    )
    result_df2 = result_df[result_df["source"] == "sciencedirect_results"].drop(
        "source", axis=1
    )
    result_df3 = result_df[result_df["source"] == "springer_results"].drop(
        "source", axis=1
    )

    # Save the cleaned files
    result_df1.to_csv(
        "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/cleaned_pubmed_results.csv",
        index=False,
    )
    result_df2.to_csv(
        "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/cleaned_sciencedirect_results.csv",
        index=False,
    )
    result_df3.to_csv(
        "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/cleaned_springer_results.csv",
        index=False,
    )

    print("Processing complete. Cleaned files saved.")


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
    # download_pdf()

    csv_file_paths = [
        sciencedirect_path,
        pubmed_path,
        springer_path,
    ]

    reformat_datetime(*csv_file_paths)

    csv_file_paths = [
        "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/pubmed_results_tests",
        "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/sciencedirect_results_tests",
        "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/springer_results_tests",
    ]

    # Check if files exist before proceeding
    for path in csv_file_paths:
        if not os.path.exists(f"{path}.csv"):
            print(f"File not found: {path}")
            return  # Stop execution if any file is missing

    explode_csv(*csv_file_paths)
    compare_authors()


if __name__ == "__main__":
    main()
