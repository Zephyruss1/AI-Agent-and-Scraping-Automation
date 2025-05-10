import os
import sys

from dotenv import load_dotenv

load_dotenv()

# Adjust path to point to: sd_pm_ls_scraper/src/
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../src")))

import pytest  # noqa: E402
from sciencedirect import make_request, save_csv  # noqa: E402


@pytest.fixture
def credentials():
    # Load the API key from the environment variable
    return os.getenv("ELSEVIER_API_KEY"), os.getenv("ELSEVIER_INSTTOKEN")


def test_sciencedirect_scraper_initialization_and_return_dtype(credentials):
    API_KEY, INSTTOKEN = credentials

    keyword = "machine learning"
    start_date = "2022-01-01"
    end_date = "2023-01-01"
    scraper = make_request(
        keyword,
        start_date,
        end_date,
        max_papers=5,
        API_KEY=API_KEY,
        INSTTOKEN=INSTTOKEN,
    )

    assert scraper is not None
    assert isinstance(scraper, list)


def test_sciencedirect_scraper_results_exists(credentials):
    API_KEY, INSTTOKEN = credentials

    keyword = "machine learning"
    start_date = "2022-01-01"
    end_date = "2023-01-01"
    scraper = make_request(
        keyword,
        start_date,
        end_date,
        max_papers=5,
        API_KEY=API_KEY,
        INSTTOKEN=INSTTOKEN,
    )

    assert all(isinstance(entry, dict) for entry in scraper)
    assert all("dc:title" in entry for entry in scraper)
    assert all("prism:doi" in entry for entry in scraper)
    assert all("authors" in entry for entry in scraper)
    assert all("$" in entry for entry in scraper)
    assert all("prism:publicationName" in entry for entry in scraper)
    assert all("prism:coverDate" in entry for entry in scraper)
    assert all("openaccess" in entry for entry in scraper)
    assert all("link" in entry for entry in scraper)


def test_sciencedirect_scraper_csv_output(credentials):
    API_KEY, INSTTOKEN = credentials
    keyword = "machine learning"
    start_date = "2022-01-01"
    end_date = "2023-01-01"
    scraper = make_request(
        keyword,
        start_date,
        end_date,
        max_papers=5,
        API_KEY=API_KEY,
        INSTTOKEN=INSTTOKEN,
    )

    assert isinstance(scraper, list)

    # Check if the CSV file is created
    save_csv(scraper, keyword)

    csv_file_path = "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/sciencedirect_results.csv"
    assert os.path.exists(csv_file_path)
    assert os.path.getsize(csv_file_path) > 0
