import csv
import os
import time

import requests
from dotenv import load_dotenv

load_dotenv()

api_key = os.getenv("ELSEVIER_API_KEY")


def make_request(search_query: str, max_papers: int = 10) -> list:
    retries, start, count = 0, 0, 25
    results = []
    while len(results) < max_papers and start < 20000:
        params = {
            "query": search_query,
            "start": start,
            "count": count,
            "view": "COMPLETE",
        }

        headers = {"Accept": "application/json", "X-ELS-APIKey": api_key}

        response = requests.get(
            "https://api.elsevier.com/content/search/sciencedirect",
            params=params,
            headers=headers,
            timeout=10,
        )
        print(f"Status code: {response.status_code}")
        try:
            if response.status_code == 200:
                data = response.json()
                if "search-results" in data and "entry" in data["search-results"]:
                    results.extend(data["search-results"]["entry"])
                    print(
                        f"Retrieved {len(data['search-results']['entry'])} results. Total: {len(results)}"
                    )
                else:
                    print("No results found in the response.")
                    break
                start += count
            elif response.status_code == 429:
                wait_time = 2**retries
                print(f"Received 429, waiting {wait_time}")
                time.sleep(wait_time)
                retries += 1
        except Exception as e:
            print(f"Error: {e}")
            break
    return results


def save_csv(results: list) -> None:
    if results:
        os.makedirs("output", exist_ok=True)

        headers = ["Title", "DOI", "Authors", "Journal", "Date", "Paper Link", "Free"]
        with open(
            "output/sciencedirect_results.csv", "w", newline="", encoding="utf-8"
        ) as f:
            writer = csv.writer(f)
            writer.writerow(headers)

            for doc in results:
                title = doc.get("dc:title", "N/A")
                doi = doc.get("prism:doi", "N/A")

                authors = doc.get("authors", {}).get("author", [])
                if isinstance(authors, list):
                    author_names = [author.get("$", "Unknown") for author in authors]
                elif isinstance(authors, dict):
                    author_names = [authors.get("$", "Unknown")]
                else:
                    author_names = ["N/A"]
                formatted_authors = ", ".join(author_names)

                journal = doc.get("prism:publicationName", "N/A")
                date = doc.get("prism:coverDate", "N/A")

                links = doc.get("link", [])
                paper_link = next(
                    (link.get("@href", "N/A") for link in links if "@href" in link),
                    "N/A",
                )

                free = doc.get("openaccess", "N/A")

                writer.writerow(
                    [title, doi, formatted_authors, journal, date, paper_link, free]
                )

        print(
            f"✅ Successfully wrote {len(results)} results to output/sciencedirect_results.csv"
        )
    else:
        print("⚠ No results to write to CSV.")


results = make_request("microbiology", max_papers=100)
save_csv(results)
