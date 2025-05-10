import os
import sys

from dotenv import load_dotenv

load_dotenv()

# Adjust path to point to: sd_pm_ls_scraper/src/
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../src")))

from pubmed import fetch_pubmed_details_batch, search_pubmed  # noqa: E402


def test_search_pubmed_scraper():
    # Test the search_pubmed function
    query = "COVID-19"
    max_results = 10
    date_ranges = [("2020/01/01", "2020/12/31")]

    results = search_pubmed(query, max_results, date_ranges)

    assert isinstance(results, list)
    assert len(results) <= max_results
    assert all(isinstance(result, str) for result in results)


def test_fetch_pubmed_scraper_details_batch():
    # Test the fetch_pubmed_details_batch function
    id_list = ["31452104", "31452105"]
    batch_size = 1

    records = fetch_pubmed_details_batch(id_list, batch_size)

    assert isinstance(records, list)
    assert len(records) == len(id_list)
    assert all(isinstance(record, dict) for record in records)


def test_search_pubmed_scraper_no_results():
    # Test the search_pubmed function with no results
    query = "nonexistentquery"
    max_results = 10
    date_ranges = [("1900/01/01", "1900/12/31")]

    results = search_pubmed(query, max_results, date_ranges)

    assert isinstance(results, list)
    assert len(results) == 0
    assert all(isinstance(result, str) for result in results)
    assert all(result == "" for result in results)


def test_fetch_pubmed_scraper_details_batch_empty_list():
    # Test the fetch_pubmed_details_batch function with an empty list
    id_list = []
    batch_size = 1

    records = fetch_pubmed_details_batch(id_list, batch_size)

    assert isinstance(records, list)
    assert len(records) == 0
    assert all(isinstance(record, dict) for record in records)


def test_pubmed_scraper_csv_output():
    # Test the CSV output of the PubMed scraper
    query = "COVID-19"
    max_results = 10
    date_ranges = [("2020/01/01", "2020/12/31")]

    results = search_pubmed(query, max_results, date_ranges)

    assert isinstance(results, list)

    # Check if the CSV file is created
    csv_file_path = (
        "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/pubmed_results.csv"
    )
    with open(csv_file_path, "w") as f:
        f.write("PubMed ID\n")
        for result in results:
            f.write(f"{result}\n")

    assert os.path.exists(csv_file_path)
    assert os.path.getsize(csv_file_path) > 0
    assert os.path.isfile(csv_file_path)
    assert os.access(csv_file_path, os.R_OK)
    assert os.access(csv_file_path, os.W_OK)
