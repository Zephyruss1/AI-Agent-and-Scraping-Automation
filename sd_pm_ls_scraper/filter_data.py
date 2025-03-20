import os
import warnings

import pandas as pd

warnings.filterwarnings("ignore")

# Define base directory path
os.makedirs("/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/pdfs", exist_ok=True)
BASE_DIR = "/root/arxiv-and-scholar-scraping"
PDF_DIR = os.path.join(BASE_DIR, "sd_pm_ls_scraper/pdfs")

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

# try:
#     springer_path = os.path.join(BASE_DIR, "sd_pm_ls_scraper/output/springer_results.csv")
#     print("SpringerLink data loaded successfully.")
# except Exception as e:
#     print(f"Error loading SpringerLink data: {e}")


def reformat_datetime(data: pd.DataFrame) -> pd.DataFrame:
    """
    Reformat the datetime column.
    """
    try:
        data["Date"] = pd.to_datetime(data["Date"], format="mixed", errors="coerce")
        data["Date"] = data["Date"].dt.strftime("%Y-%m-%d")
    except Exception as e:
        print(f"Error reformatting datetime: {e}")
    return data.sort_values(by="Date", ascending=False)


pubmed = reformat_datetime(pubmed_csv)
sciencedirect = reformat_datetime(sciencedirect_csv)


def compare_authors(
    data1: pd.DataFrame, data2: pd.DataFrame, output1: str, output2: str
):
    """
    Compare authors between two datasets and keep only the latest paper for each author.
    """

    try:
        # Ensure the Date column is in datetime format
        data1["Date"] = pd.to_datetime(data1["Date"], errors="coerce")
        data2["Date"] = pd.to_datetime(data2["Date"], errors="coerce")

        # Merge the two datasets on the Authors column
        merged = pd.merge(
            data1, data2, on="Authors", suffixes=("_1", "_2"), how="outer"
        )

        # Determine the latest paper for each author
        merged["Latest_Date"] = merged[["Date_1", "Date_2"]].max(axis=1)

        # Filter rows to keep only the latest papers
        data1_filtered = merged[merged["Date_1"] == merged["Latest_Date"]].drop(
            columns=["Date_2", "Latest_Date"]
        )
        data2_filtered = merged[merged["Date_2"] == merged["Latest_Date"]].drop(
            columns=["Date_1", "Latest_Date"]
        )

        # Save the filtered datasets to CSV
        data1_filtered.to_csv(output1, index=False)
        data2_filtered.to_csv(output2, index=False)

        print(f"Filtered data saved to {output1}.")
        print(f"Filtered data saved to {output2}.")
    except Exception as e:
        print(f"Error in compare_authors: {e}")


# Call the function with the loaded DataFrames and output file paths
compare_authors(
    pubmed,
    sciencedirect,
    os.path.join(BASE_DIR, "sd_pm_ls_scraper/output/pubmed_results_filtered.csv"),
    os.path.join(
        BASE_DIR, "sd_pm_ls_scraper/output/sciencedirect_results_filtered.csv"
    ),
)


def _load_filtered_csv():
    """
    Load the filtered CSV files.
    """
    os.makedirs(PDF_DIR, exist_ok=True)

    try:
        pubmed_filtered = pd.read_csv(
            "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/pubmed_results_filtered.csv"
        )
        sciencedirect_filtered = pd.read_csv(
            "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/sciencedirect_results_filtered.csv"
        )
        return pubmed_filtered, sciencedirect_filtered
    except FileNotFoundError as e:
        print(f"Error loading filtered CSV files: {e}")


def download_pdf() -> str:
    """ "Download the PDF from the PMC Full Text URL."""
    try:
        import time

        import requests
    except NotImplementedError as e:
        print(f"Error importing requests: {e}")

    pubmed_df, sciencedirect_df = _load_filtered_csv()

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
                    pdf_filename = os.path.join(PDF_DIR, f"{file_name}.pdf")
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
            # Make the API request with redirect following
            response = requests.get(API_URL, headers=headers, allow_redirects=True)

            # Check if the request was successful
            if response.status_code == 200:
                # Verify that the content type is PDF
                if response.headers.get("Content-Type") == "application/pdf":
                    # Save the PDF to a file
                    filename = paper_doi.replace("/", "_") + ".pdf"
                    pdf_filename = os.path.join(
                        PDF_DIR, f"{paper_doi.replace('/', '_')} + .pdf"
                    )
                    with open(pdf_filename, "wb") as pdf_file:
                        pdf_file.write(response.content)
                    print(f"PDF downloaded successfully: {filename}")
                else:
                    print(
                        f"Error: Content-Type is {response.headers.get('Content-Type')}, not PDF."
                    )
            else:
                print(f"Failed to download PDF. Status code: {response.status_code}")
                print(f"Response: {response.text}")
        except requests.exceptions.RequestException as e:
            print(f"An error occurred: {e}")
        break
    print(f"\nTotal PDFs downloaded from ScienceDirect: {count}")


download_pdf()
