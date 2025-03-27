import asyncio
import re
from typing import Dict, List

from playwright.async_api import async_playwright


async def _save_to_csv(results: List[Dict], _keyword: str) -> None:
    """Save scraped Springer search results to a CSV file.

    Writes a list of article dictionaries to a CSV file located at a fixed path, including
    headers and the provided keyword as an additional column. Handles empty results and
    missing CSV module gracefully.

    Args:
        results (List[Dict]): A list of dictionaries containing article details. Each dictionary
            should have keys like 'Title', 'Link', 'Description', 'Authors', and 'Date'.
            Missing keys default to 'null'.
        _keyword (str): The search keyword used to generate the results, added to each row.

    Returns:
        None: The function writes to a file and does not return a value.

    Raises:
        ImportError: If the `csv` module cannot be imported, with a descriptive message.
        Exception: If an error occurs during file writing, with details about the failure.

    Notes:
        - The output file is hardcoded to '/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/springer_results.csv'.
        - Uses UTF-8 encoding to support special characters in article data.
        - Prints success or info messages to the console for user feedback.
        - Asynchronous function, though it performs synchronous I/O operations.
    """
    try:
        import csv
    except ImportError as err:
        raise ImportError("csv module not found!") from err

    if not results:
        print("     [INFO] No results to save!")
        return

    headers = [
        "Title",
        "Link",
        "Description",
        "Authors",
        "Date",
        "Keyword",
    ]

    # Open the file and write the headers and rows
    try:
        with open(
            "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/springer_results.csv",
            "w",
            newline="",
            encoding="utf-8",
        ) as f:
            writer = csv.writer(f)
            writer.writerow(headers)  # Write the headers

            for result in results:
                writer.writerow(
                    [
                        result.get("title", "null"),
                        result.get("link", "null"),
                        result.get("description", "null"),
                        result.get("authors", "null"),
                        result.get("published_date", "null"),
                        _keyword,
                    ]
                )
        print(
            f"✅ Successfully wrote {len(results)} results to output/springer_results.csv"
        )
    except Exception as err:
        raise Exception(
            f"    ❌ [ERROR] An error occurred while saving Excel file: {err}"
        ) from err


class SpringerScraper:
    def __init__(self, keyword: str, max_results: int = None):
        """
        Initializes a SpringerScraper instance with a keyword and constructs the search URL.

        Args:
            keyword (str): The search term to query Springer for articles. Spaces are replaced with '+'
                           to form a valid URL query parameter.
            max_results (int, optional): Maximum number of results to scrape.
        """
        self.keyword = keyword.replace(" ", "+") if " " in keyword else keyword
        self.base_url = "https://link.springer.com/search"
        self.url = f"{self.base_url}?new-search=true&query={self.keyword}&content-type=research"
        self.articles = []
        self.page = None
        self.max_results = max_results
        self.unique_links = set()  # Track unique article links

    async def _check_cookies(self, page):
        """
        Handles the cookie consent dialog on the Springer page, if it appears.

        Args:
            page (Page): The Playwright page object representing the current browser tab.
        """
        try:
            await page.wait_for_selector("xpath=/html/body/dialog", timeout=5000)
            close_button = await page.query_selector(
                "xpath=/html/body/dialog/div/div/div[3]/button"
            )
            if close_button:
                await close_button.click()
        except Exception as e:
            print(f"Cookie dialog not found or already closed: {e}")

    async def searching_and_gathering_papers(self, page):
        """
        Performs a search on Springer and gathers research articles from the search results page.
        """
        # Check if we've reached max results
        if self.max_results and len(self.articles) >= self.max_results:
            print("     [INFO] Maximum results reached. Skipping search.")
            return False

        print("\n📍 Step: Gathering papers!")
        self.page = page  # Store page reference if needed later

        await page.wait_for_load_state("networkidle")

        # Wait for article list to load
        try:
            await page.wait_for_selector(
                'xpath=//*[@id="main"]/div/div/div/div[2]/div[4]/ol', timeout=10000
            )
        except Exception as e:
            print(f"Failed to load article list: {e}")
            return False

        articles = await page.query_selector_all(
            'xpath=//*[@id="main"]/div/div/div/div[2]/div[4]/ol/li'
        )

        if not articles:
            print(
                "No articles found. Check if the selector is correct or if content is loaded."
            )
            return False

        print(f"Found {len(articles)} articles on this page")
        new_articles_found = False

        for i, article in enumerate(articles):
            # Stop if max results is reached
            if self.max_results and len(self.articles) >= self.max_results:
                break

            try:
                # Extract article details
                title_link_element = await article.query_selector("xpath=.//div/h3/a")
                if not title_link_element:
                    continue

                title_link = await title_link_element.get_attribute("href")
                full_link = f"https://link.springer.com{title_link}"

                # Skip if link is already scraped
                if full_link in self.unique_links:
                    continue

                # Extract other details
                title_element = await article.query_selector("xpath=.//div/h3")
                title = (
                    await title_element.text_content()
                    if title_element
                    else "No title found"
                )

                description_element = await article.query_selector("xpath=.//div/p")
                description = (
                    await description_element.text_content()
                    if description_element
                    else "No description found"
                )

                authors_element = await article.query_selector(
                    "xpath=.//div/div[3]/div"
                )
                authors = (
                    await authors_element.text_content()
                    if authors_element
                    else "No authors found"
                )
                authors = re.sub(r"\s+in\s.*$", "", authors)
                authors = authors.replace("...", "").strip()

                published_date_element = await article.query_selector(
                    "xpath=.//div/div[4]/div/span[1]"
                )
                published_date = (
                    await published_date_element.text_content()
                    if published_date_element
                    else "No published date found"
                )

                # Add article and track unique link
                article_info = {
                    "title": title.strip(),
                    "link": full_link,
                    "description": description.strip(),
                    "authors": authors.strip(),
                    "published_date": published_date.strip(),
                }

                self.articles.append(article_info)
                self.unique_links.add(full_link)
                new_articles_found = True

            except Exception as e:
                print(f"Error processing article {i}: {e}")

        return new_articles_found

    async def pagination(self, page, max_pages: int = None) -> None:
        """
        Handle pagination for Springer search results with improved tracking.

        Args:
            page: The Playwright page object used for browser interaction.
            max_pages (int): Maximum number of pages to navigate.
        """
        current_page = 1
        while True:
            # Check stopping conditions
            if (max_pages and current_page > max_pages) or (
                self.max_results and len(self.articles) >= self.max_results
            ):
                print(
                    f"     [INFO] Stopping pagination. "
                    f"Current page: {current_page}, "
                    f"Total articles: {len(self.articles)}"
                )
                break

            try:
                # Construct paginated URL
                paginated_url = f"{self.base_url}?new-search=true&query={self.keyword}&content-type=research&page={current_page}"
                print(f"     [INFO] Navigating to page {current_page}: {paginated_url}")

                # Navigate to the paginated URL
                await page.goto(paginated_url)
                await page.wait_for_load_state("networkidle")

                # Small delay to ensure page is loaded
                await asyncio.sleep(2)

                # Try to gather articles
                new_articles = await self.searching_and_gathering_papers(page)

                # If no new articles found, break pagination
                if not new_articles:
                    print("     [INFO] No new articles found. Stopping pagination.")
                    break

                # Increment page counter
                current_page += 1

            except Exception as err:
                print(f"Error during pagination: {err}")
                break

        print(
            f"     [INFO] Pagination complete. Total pages: {current_page - 1}, Total unique articles: {len(self.articles)}"
        )


async def async_springer():
    async with async_playwright() as pw:
        args = [
            "--no-sandbox",
            "--disable-blink-features=AutomationControlled",
            "--disable-infobars",
        ]

        browser = await pw.chromium.launch(headless=False, args=args)
        context = await browser.new_context(viewport={"width": 1280, "height": 720})
        page = await context.new_page()

        search_term = "bacteria AND Plant pathology AND Insecticides"
        max_papers = 1000

        SpringerScraper_obj = SpringerScraper(search_term, max_results=max_papers)

        try:
            # Navigate to initial search page
            await page.goto(SpringerScraper_obj.url)
            await SpringerScraper_obj._check_cookies(page)

            # Perform pagination
            await SpringerScraper_obj.pagination(page, max_pages=20)

            # Save results
            if SpringerScraper_obj.articles:
                print(
                    f"\nSaving {len(SpringerScraper_obj.articles)} total unique articles to CSV file..."
                )
                await _save_to_csv(SpringerScraper_obj.articles, _keyword=search_term)
            else:
                print("No articles found.")

        except Exception as err:
            await page.screenshot(path="unexpected_error.png")
            print(f"An unexpected error occurred: {err}")

        finally:
            print("\nClosing browser...")
            await browser.close()


if __name__ == "__main__":
    asyncio.run(async_springer())
