import csv
import os
import time

import requests
from dotenv import load_dotenv

load_dotenv()

API_KEY = os.getenv("ELSEVIER_API_KEY")
INSTTOKEN = os.getenv("ELSEVIER_INSTTOKEN")


def make_request(
    search_query: str,
    start_year: str = None,
    end_year: str = None,
    max_papers: int = 18000,
) -> list:
    """
    Fetch ScienceDirect articles with a date range filter and specific article types.

    Args:
        search_query (str): The search term (e.g., "microbial kinetics AND CFU").
        start_year (str): Start date in 'YYYY/MM/DD' format (e.g., "2010/01/01").
        end_year (str): End date in 'YYYY/MM/DD' format (e.g., "2020/12/31").
        max_papers (int): Maximum number of papers to retrieve (default: 18000).

    Returns:
        list: List of search result entries.
    """
    retries: int
    start: int
    count: int = 25

    retries, start = 0, 0
    results = []

    year_range = f"{start_year}-{end_year}"

    params = {
        "query": search_query,
        "start": start,
        "count": count,
        "view": "COMPLETE",
        # "subtype": "article,review",  # Article type filter
    }
    if start_year and end_year:
        params["date"] = year_range

    headers = {
        "Accept": "application/json",
        "X-ELS-APIKey": API_KEY,
        "X-ELS-Insttoken": INSTTOKEN,
    }

    print(f"Using year range filter: date={year_range}")

    while len(results) < max_papers and start < 20000:
        try:
            response = requests.get(
                "https://api.elsevier.com/content/search/sciencedirect",
                params=params,
                headers=headers,
                timeout=10,
            )
            print(f"Request URL: {response.url}")
            print(f"Status code: {response.status_code}")
            if response.status_code == 200:
                data = response.json()
                if "search-results" in data:
                    total_results = int(
                        data["search-results"].get("opensearch:totalResults", 0)
                    )
                    print(f"Total available results: {total_results}")
                    if "entry" in data["search-results"]:
                        entries = data["search-results"]["entry"]
                        results.extend(entries)
                        print(
                            f"Retrieved {len(entries)} results in this batch. Total: {len(results)}"
                        )
                        if len(entries) < count and len(results) < total_results:
                            print(
                                "Received fewer results than requested, but more may be available."
                            )
                        elif len(entries) == 0 or len(results) >= total_results:
                            print("No more results available.")
                            break
                        start += count
                        params["start"] = start  # Update start for next batch
                    else:
                        print("No results found in the response.")
                        break
                else:
                    print("Invalid response format.")
                    break
            elif response.status_code == 429:
                wait_time = 2**retries
                print(f"Received 429, waiting {wait_time} seconds")
                time.sleep(wait_time)
                retries += 1
            else:
                print(f"Unexpected status code: {response.status_code}")
                break
            time.sleep(1)
        except Exception as e:
            print(f"Error: {e}")
            break

    return results[:max_papers]


def save_csv(results: list, _keyword: str) -> None:
    if results:
        os.makedirs(
            "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output", exist_ok=True
        )

        headers = [
            "Title",
            "DOI",
            "Authors",
            "Source",
            "Date",
            "Paper Link",
            "Openaccess",
            "Keyword",
        ]
        with open(
            "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/sciencedirect_results.csv",
            "w",
            newline="",
            encoding="utf-8",
        ) as f:
            writer = csv.writer(f)
            writer.writerow(headers)

            for doc in results:
                title = doc.get("dc:title", "N/A")
                doi = doc.get("prism:doi", "N/A")

                authors_data = doc.get("authors")
                if authors_data and isinstance(authors_data, dict):
                    authors = authors_data.get("author", [])
                    if isinstance(authors, list):
                        author_names = [
                            author.get("$", "Unknown") for author in authors
                        ]
                    elif isinstance(authors, dict):
                        author_names = [authors.get("$", "Unknown")]
                    else:
                        author_names = ["N/A"]
                else:
                    author_names = ["N/A"]
                formatted_authors = ", ".join(author_names)

                journal = doc.get("prism:publicationName", "N/A")
                date = doc.get("prism:coverDate", "N/A")

                open_access = doc.get("openaccess", "N/A")

                links = doc.get("link", [])
                paper_link = next(
                    (link.get("@href", "N/A") for link in links if "@href" in link),
                    "N/A",
                )

                writer.writerow(
                    [
                        title,
                        doi,
                        formatted_authors,
                        journal,
                        date,
                        paper_link,
                        open_access,
                        _keyword,
                    ]
                )

        print(
            f"✅ Successfully wrote {len(results)} results to output/sciencedirect_results.csv"
        )
    else:
        print("⚠ No results to write to CSV.")


if __name__ == "__main__":
    search_term = "microbial kinetics AND CFU"
    start_date = "2010"
    end_date = "2020"
    max_papers = 100

    results = make_request(search_term, start_date, end_date, max_papers)
    print(f"Total articles retrieved: {len(results)}")
    save_csv(results, _keyword=search_term)
