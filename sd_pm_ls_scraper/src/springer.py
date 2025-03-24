import asyncio
from typing import Dict, List

from playwright.async_api import async_playwright


async def _save_to_csv(results: List[Dict], _keyword: str) -> None:
    """Save scraped Springer search results to a CSV file.

    Writes a list of article dictionaries to a CSV file located at a fixed path, including
    headers and the provided keyword as an additional column. Handles empty results and
    missing CSV module gracefully.

    Args:
        results (List[Dict]): A list of dictionaries containing article details. Each dictionary
            should have keys like 'title', 'link', 'description', 'authors', and 'published_date'.
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
        "title",
        "link",
        "description",
        "authors",
        "published_date",
        "keyword",
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
    """
    SpringerScraper is an asynchronous web scraper designed to gather research articles from the Springer website
    based on a specified keyword. It handles cookie consent dialogs, connects to the Springer search page, and
    extracts article details such as title, link, description, authors, and published date. This class uses
    Playwright for navigating the webpage and scraping content asynchronously.

    Attributes:
        keyword (str): The keyword used to perform the search on Springer. Spaces in the keyword are replaced
                      with '+' for URL compatibility.
        url (str): The URL constructed using the keyword to query Springer for research articles.
        articles (list): A list of dictionaries containing information about the gathered articles.
        page (Page): A reference to the current Playwright page object.

    Methods:
        __init__(self, keyword: str):
            Initializes the SpringerScraper instance with the given keyword and constructs the search URL.

        async _check_cookies(self, page):
            Handles the cookie consent dialog on the Springer webpage if it appears.

        async connecting_to_springer(self, page):
            Navigates to the Springer search page and handles the cookie consent dialog, ensuring the page
            is fully loaded before proceeding.

        async searching_and_gathering_papers(self, page):
            Performs the search on Springer, waits for the list of articles to load, and extracts details
            from each article found. Adds the article information to the `self.articles` list.
    """

    def __init__(self, keyword: str):
        """
        Initializes a SpringerScraper instance with a keyword and constructs the search URL.

        Args:
            keyword (str): The search term to query Springer for articles. Spaces are replaced with '+'
                           to form a valid URL query parameter.
        """
        self.keyword = keyword.replace(" ", "+") if " " in keyword else keyword
        self.url = f"https://link.springer.com/search?new-search=true&query={self.keyword}&content-type=research"
        self.articles = []
        self.page = None

    async def _check_cookies(self, page):
        """
        Handles the cookie consent dialog on the Springer page, if it appears.

        Args:
            page (Page): The Playwright page object representing the current browser tab.

        Raises:
            Exception: Prints an error message if the cookie dialog is not found or has already been closed.
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

    async def connecting_to_springer(self, page):
        """
        Connects to the Springer search page and handles the cookie consent dialog, if present.

        Args:
            page (Page): The Playwright page object representing the current browser tab.

        Raises:
            Exception: If there is an error navigating to the URL or loading the page.
        """
        print("\n📍 Step 1: Connecting to Springer!")
        try:
            await page.goto(self.url)
            await self._check_cookies(page)
            await page.wait_for_load_state("networkidle")
            print("     [INFO] Connected to Springer successfully!")
        except Exception as err:
            raise Exception(f"Failed to connect to Springer: {err}") from err

    async def searching_and_gathering_papers(self, page):
        """
        Performs a search on Springer and gathers research articles from the search results page.

        The method extracts details from each article, including title, link, description, authors, and
        published date, and appends the gathered information to the `self.articles` list.

        Args:
            page (Page): The Playwright page object representing the current browser tab.

        Raises:
            Exception: Prints an error if the list of articles fails to load.
        """
        print("\n📍 Step 2: Gathering papers!")
        self.page = page  # Store page reference if needed later

        await page.goto(self.url)
        await self._check_cookies(page)
        await page.wait_for_load_state("networkidle")

        # Wait for article list to load
        try:
            await page.wait_for_selector(
                'xpath=//*[@id="main"]/div/div/div/div[2]/div[4]/ol', timeout=10000
            )
        except Exception as e:
            print(f"Failed to load article list: {e}")
            return

        articles = await page.query_selector_all(
            'xpath=//*[@id="main"]/div/div/div/div[2]/div[4]/ol/li'
        )

        if not articles:
            print(
                "No articles found. Check if the selector is correct or if content is loaded."
            )
            return

        print(f"Found {len(articles)} articles")
        print("\n📍 Step 2: Gathering articles!")

        total_found = 0
        for i, article in enumerate(articles):
            try:
                # Extract article details with improved error handling
                title_element = await article.query_selector("xpath=.//div/h3")
                title = (
                    await title_element.text_content()
                    if title_element
                    else "No title found"
                )

                title_link_element = await article.query_selector("xpath=.//div/h3/a")
                title_link = (
                    await title_link_element.get_attribute("href")
                    if title_link_element
                    else None
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

                published_date_element = await article.query_selector(
                    "xpath=.//div/div[4]/div/span[1]"
                )
                published_date = (
                    await published_date_element.text_content()
                    if published_date_element
                    else "No published date found"
                )

                # Print and append article
                print(f"\nFound Article {i + 1}")
                total_found += 1
                self.articles.append(
                    {
                        "title": title.strip(),
                        "link": f"https://link.springer.com{title_link}"
                        if title_link
                        else "No link found",
                        "description": description.strip(),
                        "authors": authors.strip(),
                        "published_date": published_date.strip(),
                    }
                )

            except Exception as e:
                print(f"Error processing article {i}: {e}")
                continue

    async def pagination(self, page, max_pages: int) -> None:
        """Handle pagination for Springer search results.

        Navigates through search result pages by clicking the "Next" button, gathering articles
        from each page until there are no more pages, an error occurs, or the maximum number
        of pages is reached. Updates the `self.articles` list with new articles found on each page.

        Args:
            page: The Playwright page object used for browser interaction.
            max_pages (int): Maximum number of pages to navigate. If 0 or None, no limit is applied.

        Returns:
            None: This method modifies the `self.articles` list in place and does not return a value.

        Raises:
            Exception: Propagates any errors encountered during navigation or article gathering,
                with an error message printed to the console.

        Notes:
            - Uses an XPath selector to locate the "Next" button on the Springer search page.
            - Checks for the `aria-disabled` attribute to determine if there are more pages.
            - Calls `searching_and_gathering_papers` to collect articles from each new page.
            - Prints informational messages about pagination progress and completion.
        """
        paginated = 0
        while True:
            if max_pages and paginated >= max_pages:
                print("     [INFO] Reached maximum pages to navigate.")
                break
            try:
                # Check for next page button
                next_page_button = await page.query_selector(
                    'xpath=//*[@id="main"]/div/div/div/div[2]/nav/ul/li[8]/a'
                )

                if (
                    not next_page_button
                    or await next_page_button.get_attribute("aria-disabled") == "true"
                ):
                    print("     [INFO] No more pages to navigate.")
                    break

                print("     [INFO] Navigating to the next page...")
                await next_page_button.click()
                await page.wait_for_load_state("networkidle")

                # Gather articles from the new page
                await self.searching_and_gathering_papers(page)
                paginated += 1
            except Exception as err:
                print(f"Error during pagination: {err}")
                break
        print("     [INFO] Finished pagination. Total pages navigated:", paginated)


async def async_springer():
    async with async_playwright() as pw:
        args = [
            "--no-sandbox",
            "--disable-blink-features=AutomationControlled",
            "--disable-infobars",
            "--disable-background-timer-throttling",
            "--disable-popup-blocking",
            "--disable-backgrounding-occluded-windows",
            "--disable-renderer-backgrounding",
            "--disable-window-activation",
            "--disable-focus-on-load",
            "--no-first-run",
            "--no-default-browser-check",
            "--no-startup-window",
            "--window-position=960,700",
        ]

        browser = await pw.chromium.launch(headless=False, args=args)
        context = await browser.new_context(viewport={"width": 1280, "height": 720})
        page = await context.new_page()

        all_articles = []

        search_term = "machine learning"
        SpringerScraper_obj = SpringerScraper(search_term)
        try:
            await SpringerScraper_obj.connecting_to_springer(page)
            await SpringerScraper_obj.searching_and_gathering_papers(page)
            await SpringerScraper_obj.pagination(page, max_pages=2)
            all_articles.extend(SpringerScraper_obj.articles)
        except Exception as err:
            await page.screenshot(path="unexpected_error.png")
            raise Exception(f"An unexpected error occurred: {err}") from err

        finally:
            print("\nClosing browser...")
            await browser.close()

            if all_articles:
                print(
                    f"\nSaving {len(all_articles)} total articles saving to CSV file..."
                )
                await _save_to_csv(all_articles, _keyword=search_term)
                print("No articles found for any keywords.")


if __name__ == "__main__":
    asyncio.run(async_springer())
