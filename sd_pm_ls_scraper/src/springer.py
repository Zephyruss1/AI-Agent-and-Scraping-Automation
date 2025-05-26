import asyncio
import re
from typing import Dict, List

import pandas as pd
from playwright.async_api import async_playwright


async def _handle_duplicate_csv_filenames(csv_filepath: str) -> None:
    """Handle the case where the CSV file already exists.

    If the file already exists, append a counter to the filename to create a unique name.
    """
    import os

    if not os.path.exists(csv_filepath):
        return csv_filepath

    base, ext = os.path.splitext(csv_filepath)
    counter = 1
    new_csv_filepath = f"{base}_{counter}{ext}"
    while os.path.exists(new_csv_filepath):
        counter += 1
        new_csv_filepath = f"{base}_{counter}{ext}"
    os.rename(csv_filepath, new_csv_filepath)
    print(f"     [INFO] Renamed existing file to {new_csv_filepath}")
    return new_csv_filepath


def _load_csv(csv_filepath: str):
    """Load a CSV file and return a DataFrame."""
    try:
        return pd.read_csv(csv_filepath)
    except Exception as err:
        raise Exception(f"Error loading CSV file: {err}") from err


async def _save_to_csv(results: List[Dict], _keyword: str) -> str:
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
        - The output file is hardcoded to '/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/springer_results_tests.csv'.
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
        "DOI",
        "Description",
        "Authors",
        "Date",
        "Keyword",
    ]

    # Check if the file already exists and handle duplicates
    file_name = "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/springer_results_tests.csv"
    file_name = await _handle_duplicate_csv_filenames(file_name)

    # Open the file and write the headers and rows
    try:
        with open(
            file_name,
            "a",
            newline="",
            encoding="utf-8",
        ) as f:
            writer = csv.writer(f)
            writer.writerow(headers)  # Write the headers

            for result in results:
                # Extract DOI safely
                link = result.get("link", "null")
                doi = "null"
                if "article/" in link:
                    try:
                        doi = link.split("article/")[1]
                    except IndexError:
                        doi = "null"
                writer.writerow(
                    [
                        result.get("title", "null"),
                        link,
                        doi,
                        result.get("description", "null"),
                        result.get("authors", "null"),
                        result.get("published_date", "null"),
                        _keyword,
                    ],
                )
        print(
            f"✅ Successfully wrote {len(results)} results to output/springer_results_tests.csv",
        )
    except Exception as err:
        raise Exception(
            f"    ❌ [ERROR] An error occurred while saving Excel file: {err}",
        ) from err
    return file_name


class SpringerScraper:
    """
    A class to scrape research articles from Springer using Playwright.
    """

    print("\n📍 Step 1: Connecting to Springer!")
    print("     🌐 [INFO] Trying to Connecting to Springer...")

    def __init__(
        self,
        keyword: str,
        max_results: int = None,
        start_date: str = None,
        end_date: str = None,
    ):
        """
        Initializes a SpringerScraper instance with a keyword and constructs the search URL.

        Args:
            keyword (str): The search term to query Springer for articles. Spaces are replaced with '+'
                           to form a valid URL query parameter.
            max_results (int, optional): Maximum number of results to scrape.
        """
        self.start_date = start_date if start_date else ""
        self.end_date = end_date if end_date else ""
        self.keyword = keyword.replace(" ", "+") if " " in keyword else keyword
        self.base_url = "https://link.springer.com/search"
        self.url = f"{self.base_url}?new-search=true&query={self.keyword}&content-type=research&content-type=article&dateFrom={self.start_date}&dateTo={self.end_date}&sortBy=relevance"
        self.articles = []
        self.page = None
        self.max_results = max_results
        self.unique_links = set()  # Track unique article links
        self.cookies_handled: bool = False

    print("     🌐 ✅ [INFO] Connected to Springer!")

    async def _check_cookies(self, page) -> None:
        """
        Handles the cookie consent dialog on the Springer page, if it appears.

        Args:
            page (Page): The Playwright page object representing the current browser tab.
        """
        if self.cookies_handled is True:
            print("     🍪 [INFO] Cookies dialog already handled.")
            return

        print("\n📍 Step 2: Checking for cookies dialog!")
        try:
            print("     🍪 [INFO] Checking for cookies dialog...")

            await page.wait_for_selector("xpath=/html/body/dialog", timeout=15000)
            close_button = await page.query_selector(
                "xpath=/html/body/dialog/div/div/div[3]/button",
            )
            if close_button:
                print("     🍪 [INFO] Cookies dialog found.")
                await close_button.click()
                print("     🍪 ✅ [INFO] Cookies dialog closed.")
        except Exception as e:
            print(f"    🍪 ❌ [INFO] Cookies dialog not found or already closed: {e}")
        self.cookies_handled = True

    async def custom_popup_handler(self, page, popup) -> None:
        """
        Custom handler for a popup dialog on Springer.

        Args:
            page (Page): The Playwright page object representing the current browser tab.
            popup (Popup): The Playwright popup object representing the dialog.

        """
        try:
            # Custom popup Xpath. -> https://prnt.sc/EzOR7WLejg3O
            await page.wait_for_selector(
                "xpath=///html/body/div[8]/div[2]/div",
                timeout=10000,
            )
            close_button = await popup.query_selector(
                "xpath=///html/body/div[8]/div[2]/div/div[3]/button[2]",
            )
            if close_button:
                await close_button.click()
                print("     [INFO] Custom popup closed.")
        except Exception as e:
            print(f"Error handling custom popup: {e}")

    async def discipline_selector_length(self, page) -> int:
        """
        Retrieves the number of disciplines available in the discipline selector dialog.

        Args:
            page (Page): The Playwright page object representing the current browser tab.

        Returns:
            int: The number of disciplines available.
        """
        try:
            await page.wait_for_selector(
                'xpath=//*[@id="popup-filters"]/div[2]/div/div[7]/div/details/div/ol/li',
                timeout=30000,
            )
            disciplines_list = await page.query_selector(
                'xpath=//*[@id="popup-filters"]/div[2]/div/div[7]/div/details/div/ol/li',
            )
            if disciplines_list:
                print("     ✅ [INFO] Discipline selector list found.")
                disciplines_list_available = await page.query_selector_all(
                    '//*[@id="list-discipline-filter"]/li/div',
                )
                return len(disciplines_list_available)
            else:
                print("     ❌ [INFO] Discipline selector list not found.")
                return 0
        except Exception as e:
            print(
                f"    ❌ [INFO] Discipline selector not found or already closed [disc_length]: {e}",
            )
            return 0

    async def discipline_selector(self, page, start_index: int) -> None:
        """
        Handles the discipline selector dialog on the Springer page.

        Args:
            page (Page): The Playwright page object representing the current browser tab.
            start_index (int): Index of the discipline to click
        """
        try:
            await page.wait_for_selector(
                'xpath=//*[@id="popup-filters"]/div[2]/div/div[7]/div/details/div/ol/li',
                timeout=30000,
            )
            disciplines_list = await page.query_selector(
                'xpath=//*[@id="popup-filters"]/div[2]/div/div[7]/div/details/div/ol/li',
            )
            if disciplines_list:
                print("     ✅ [INFO] Discipline selector list found.")

                disciplines_list_available = await page.query_selector_all(
                    '//*[@id="list-discipline-filter"]/li/div',
                )

                # Only click on the discipline at the specified index
                if start_index < len(disciplines_list_available):
                    try:
                        await disciplines_list_available[start_index].click()
                        print(f"     ✅ [INFO] Discipline {start_index + 1} clicked.")
                        discipline_name = await disciplines_list_available[
                            start_index
                        ].text_content()
                        print(
                            f"     [INFO] Discipline name {discipline_name} selected.",
                        )
                        # ---------------START OF DEBUG--------------------
                        await page.evaluate(
                            "scrollTo(0, 1200);",
                        )  # Scroll to the top of the page
                        await page.screenshot(path="discipline_selector_debug.png")
                        # ---------------END OF DEBUG--------------------
                    except Exception as e:
                        print(
                            f"     ❌ [INFO] Error clicking discipline {start_index + 1}: {e}",
                        )

                # Click the update button
                update_button = await page.query_selector(
                    '//*[@id="popup-filters"]/div[3]/button[2]',
                )
                if update_button:
                    try:
                        await update_button.click()
                        print("     ✅ [INFO] Update button clicked.")
                    except Exception as e:
                        print(f"     ❌ [INFO] Error clicking update button: {e}")
            else:
                print("     ❌ [INFO] Discipline selector list not found.")
                return

        except Exception as e:
            print(f"    ❌ [INFO] Discipline selector not found or already closed: {e}")

    async def get_total_results(self, page) -> int:
        try:
            result_count_locator = page.locator(
                "//span[@data-test='results-data-total']/ancestor::div[@class='app-search-filter__result-count u-hide-at-lg']",
            ).first
            result_count_text = await result_count_locator.inner_text()
            total_result = int(
                result_count_text.split("of")[1]
                .split("results")[0]
                .replace(",", "")
                .strip(),
            )
            print(f"total_result: {total_result}")
            # Check if the locator exists
            if total_result > 0:
                return total_result
            else:
                print("Result count element not found.")
                return 0
        except Exception as e:
            print(f"Error getting result count: {e}")
            return 0  # Explicit return 0 to avoid returning None

    async def searching_and_gathering_papers(self, page) -> bool:
        """
        Performs a search on Springer and gathers research articles from the search results page.
        """
        # Check if we've reached max results
        if self.max_results and len(self.articles) >= self.max_results:
            print("     [INFO] Maximum results reached. Skipping search.")
            return False

        print("\n📍 Step 3: Gathering papers!")
        self.page = page  # Store page reference if needed later

        await page.wait_for_load_state("networkidle")

        # Wait for article list to load
        try:
            found = False
            for i in range(1, 10):
                article_list_xpath = (
                    f'xpath=//*[@id="main"]/div/div/div/div[2]/div[{i}]/ol'
                )
                try:
                    await page.wait_for_selector(
                        article_list_xpath,
                        timeout=10000,
                    )
                    print(f"✅ Found matching element at: {article_list_xpath}")
                    found = True
                    break
                except Exception:
                    continue
            if not found:
                print("❌ No matching element found.")

        except Exception as e:
            print(f"Failed to load article list: {e}")
            return False

        articles = await page.query_selector_all(
            article_list_xpath,
        )

        if not articles:
            print(
                "No articles found. Check if the selector is correct or if content is loaded.",
            )
            return False

        print(f"Found {len(articles)} articles on this page")
        new_articles_found = False
        for i, article in enumerate(articles):
            # Stop if max results is reached
            if self.max_results and len(self.articles) >= self.max_results:
                break

            try:
                # Extract article details inside the for loop in searching_and_gathering_papers function
                title_link_element = await article.query_selector("xpath=.//div/h3/a")
                if not title_link_element:
                    continue

                # Extract the link and validate it
                title_link = await title_link_element.get_attribute("href")
                if title_link:
                    # Ensure the link is properly formed by checking if it starts with '/'
                    if not title_link.startswith("http"):
                        title_link = f"https://link.springer.com{title_link}"

                    # Skip if the link is already scraped (to avoid duplicates)
                    if title_link in self.unique_links:
                        continue

                    # Proceed with the valid, unique link
                    full_link = title_link
                else:
                    print("Skipping article with missing href attribute.")
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
                    "xpath=.//div/div[3]/div",
                )
                authors = (
                    await authors_element.text_content()
                    if authors_element
                    else "No authors found"
                )
                authors = re.sub(r"\s+in\s.*$", "", authors)
                authors = authors.replace("...", "").strip()
                published_date_element = await article.query_selector(
                    "xpath=.//div/div[4]/div/span[1]",
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

    async def pagination(self, page, max_pages: int = None) -> bool:
        """
        Handle pagination for Springer search results with improved tracking.

        Args:
            page: The Playwright page object used for browser interaction.
            max_pages (int): Maximum number of pages to navigate.
        """
        current_page = 1
        found_any_articles = False

        # Get the current URL after filtering
        base_url = page.url

        while True:
            # Check stopping conditions
            if (max_pages and current_page > max_pages) or (
                self.max_results and len(self.articles) >= self.max_results
            ):
                print(
                    f"     [INFO] Stopping pagination. "
                    f"Current page: {current_page}, "
                    f"Total articles: {len(self.articles)}",
                )
                break

            try:
                # Construct paginated URL using current page URL as base
                page_param = "page=" + str(current_page)
                if "page=" in base_url:
                    # Replace existing page parameter
                    paginated_url = re.sub(r"page=\d+", page_param, base_url)
                else:
                    # Add page parameter
                    paginated_url = (
                        base_url + ("&" if "?" in base_url else "?") + page_param
                    )

                print(f"     [INFO] Navigating to page {current_page}: {paginated_url}")

                # Navigate to the paginated URL
                await page.goto(paginated_url)
                await page.wait_for_load_state("networkidle")

                # Try to gather articles
                new_articles = await self.searching_and_gathering_papers(page)
                if new_articles:
                    found_any_articles = True

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
            f"     [INFO] Pagination complete. Total pages: {current_page - 1}, Total unique articles: {len(self.articles)}",
        )
        return found_any_articles

    async def extract_affiliations_and_full_names(
        self,
        page,
        csv_file_path: str,
    ) -> None:
        """Extract affiliations and full names from the Springer article page."""
        # Check cookies only once

        _csv_file = _load_csv(csv_file_path)
        _csv_file.drop_duplicates(inplace=True)
        print("\n📍 Step 3: Extracting affiliations and full names!")

        for _index, row in _csv_file.iterrows():
            paper_link = row["Link"]
            await page.goto(paper_link)
            print(f"    [INFO] Navigating to paper: {paper_link}")
            # await page.wait_for_load_state("networkidle") # FIXME: Having issue with this line about timeout...
            try:
                try:
                    await page.wait_for_selector(
                        'xpath=//*[@id="affiliations"]',
                        timeout=10000,
                    )
                except Exception as e:
                    print(f"Failed to find affiliations section in page: {e}")
                    continue
                # retrieve affiliations full list
                affiliations = await page.query_selector_all(
                    'xpath=//*[@id="Aff1"]/p[1]',
                )

                # reitrieve authors full list
                authors = await page.query_selector_all(
                    '//ul[contains(@class, "c-article-author-list")]//li//a',
                )

                try:
                    for affiliation in affiliations:
                        print("affiliations:", await affiliation.text_content())
                        _csv_file.at[_index, "Affiliations"] = (
                            f"{await affiliation.text_content()}."
                        )

                    authors_list = []
                    for author in authors:
                        author_text = await author.text_content()
                        if (
                            author_text.startswith(("nAff", "na", "ORCID"))
                            or len(author_text) <= 2
                        ):
                            continue
                        print("authors:", author_text.replace("&", ",").strip())
                        authors_list.append(author_text.replace("&", ",").strip())

                    _csv_file.at[_index, "Authors"] = ", ".join(authors_list)
                    print("---" * 30)
                except Exception as e:
                    print(f"Error extracting affiliations or authors: {e}")
                    continue
            except Exception as e:
                print(f"Error extracting affiliations and authors: {e}")
                continue
            finally:
                _csv_file.to_csv(csv_file_path, index=False)


async def async_springer():
    args = [
        "--no-sandbox",
        "--disable-blink-features=AutomationControlled",
        "--disable-infobars",
    ]
    search_term = '(CFU OR "colony forming unit" OR "colony counting") AND ("automated imaging" OR "image analysis" OR "colony counter" OR "colony detection" OR "machine learning" OR robotics)'
    start_date = "2012"
    end_date = "2025"

    async with async_playwright() as pw:
        # In async_springer function
        SpringerScraper_obj = SpringerScraper(
            search_term,
            max_results=None,
            start_date=start_date,
            end_date=end_date,
        )
        browser = await pw.chromium.launch(headless=False, args=args)
        context = await browser.new_context(viewport={"width": 1280, "height": 720})
        page = await context.new_page()
        try:
            if SpringerScraper_obj:
                await page.goto(SpringerScraper_obj.url)
                await SpringerScraper_obj._check_cookies(page)

                total_results = await SpringerScraper_obj.get_total_results(page)
                # First check for total results
                print(f"     [INFO] [Step-1] Checking total results: {total_results}")
                if total_results > 1000:
                    print(
                        f"     ✅ [INFO] [Step-1] After checking total results ({total_results}) are greater than 1000. "
                        "Proceeding with complex filtering.",
                    )
                    for current_start_date in range(
                        int(start_date),
                        int(end_date),
                        step=2,
                    ):
                        current_end_date = current_start_date + 1
                        print(f"{current_start_date} - {current_end_date}")

                        SpringerScraper_obj = SpringerScraper(
                            search_term,
                            max_results=None,
                            start_date=current_start_date,
                            end_date=current_end_date,
                        )
                        await page.goto(SpringerScraper_obj.url)
                        total_results = await SpringerScraper_obj.get_total_results(
                            page,
                        )
                        # Second check for total results
                        print(
                            f"     [INFO] [Step-2] Checking total results: {total_results}",
                        )
                        if total_results > 1000:
                            print(
                                f"     ✅ [INFO] [Step-2] After filtering range years total results ({total_results}) are greater than 1000. "
                                "Proceeding with complex filtering.",
                            )
                            len_of_disciplines = (
                                await SpringerScraper_obj.discipline_selector_length(
                                    page,
                                )
                            )
                            print(
                                f"     [INFO] Number of disciplines: {len_of_disciplines}",
                            )

                            j = 0
                            while j < len_of_disciplines:
                                try:
                                    # Reset articles list and unique links for each discipline
                                    SpringerScraper_obj.articles = []
                                    SpringerScraper_obj.unique_links = set()

                                    print(
                                        f"     [INFO] Processing discipline {j + 1}/{len_of_disciplines}",
                                    )

                                    await SpringerScraper_obj.discipline_selector(
                                        page,
                                        start_index=j,
                                    )

                                    await SpringerScraper_obj.pagination(page)

                                    if SpringerScraper_obj.articles:
                                        print(
                                            f"\nSaving {len(SpringerScraper_obj.articles)} total unique articles to CSV file...",
                                        )
                                    else:
                                        print("No articles found for this discipline.")

                                except Exception as discipline_err:
                                    print(
                                        f"Error processing discipline {j + 1}: {discipline_err}",
                                    )
                                finally:
                                    j += 1
                                    if j < len_of_disciplines:
                                        try:
                                            await page.goto(SpringerScraper_obj.url)
                                        except Exception as nav_err:
                                            print(
                                                f"Error navigating back to search page: {nav_err}",
                                            )
                        else:
                            print(
                                f"     ✅ [INFO] [Step-2] After filtering range years total results ({total_results}) are greater than 1000. "
                                "Proceeding without complex filtering.",
                            )
                            await page.goto(SpringerScraper_obj.url)
                            await SpringerScraper_obj.pagination(page)
                else:
                    print(
                        f"     ✅ [INFO] [Step-1] After checking total results ({total_results}) are greater than 1000. "
                        "Proceeding without complex filtering.",
                    )
                    await page.goto(SpringerScraper_obj.url)
                    await SpringerScraper_obj.pagination(page)

        except Exception as err:
            try:
                if not page.is_closed():
                    await page.screenshot(path="unexpected_error.png")
            except Exception as screenshot_err:
                print(f"Could not take screenshot: {screenshot_err}")
            print(f"An unexpected error occurred: {err}")
        finally:
            try:
                springer_file_path = await _save_to_csv(
                    SpringerScraper_obj.articles,
                    _keyword=search_term,
                )
                # Process affiliations only if articles were found
                await SpringerScraper_obj.extract_affiliations_and_full_names(
                    page,
                    springer_file_path,
                )
            except Exception as affil_err:
                print(f"Error extracting affiliations: {affil_err}")
                # Continue despite affiliation extraction errors
            print("\nClosing browser...")
            await browser.close()


if __name__ == "__main__":
    asyncio.run(async_springer())
