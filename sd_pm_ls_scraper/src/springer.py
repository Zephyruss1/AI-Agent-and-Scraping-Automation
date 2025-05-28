import asyncio
import re
from typing import Dict, List, Optional

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
    """Save scraped Springer search results to a CSV file."""
    try:
        import csv
        import os
    except ImportError as err:
        raise ImportError(
            "    ❌ [ERROR] Could not import the csv module. "
            "Please ensure it is installed in your Python environment.",
        ) from err

    if not results:
        print("     [INFO] No results to save!")
        return None

    # Get current script directory and set paths relative to it
    current_dir = os.path.dirname(os.path.abspath(__file__))
    output_dir = os.path.join(current_dir, "..", "output")
    file_name = os.path.join(output_dir, "springer_results.csv")

    # Ensure output directory exists
    os.makedirs(output_dir, exist_ok=True)

    try:
        with open(file_name, "a", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)

            # Write headers only if file is empty
            if os.path.getsize(file_name) == 0:
                writer.writerow(
                    [
                        "Title",
                        "Link",
                        "DOI",
                        "Description",
                        "Authors",
                        "Date",
                        "Keyword",
                    ],
                )

            # Write article data
            for result in results:
                link = result.get("link", "null")
                doi = link.split("article/")[1] if "article/" in link else "null"

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

        print(f"✅ Successfully wrote {len(results)} results to {file_name}")

        # Verify file was created successfully
        if os.path.exists(file_name) and os.path.getsize(file_name) > 0:
            print(f"📄 CSV file saved at: {file_name}")
            print(f"📊 File size: {os.path.getsize(file_name)} bytes")
            return file_name
        else:
            print("❌ Warning: CSV file may not have saved properly")
            return None

    except Exception as err:
        print(f"❌ Error saving CSV: {err}")
        return None


class SpringerScraper:
    """
    A class to scrape research articles from Springer using Playwright.
    """

    print("\n📍 Step 1: Connecting to Springer!")
    print("     🌐 [INFO] Trying to Connecting to Springer...")

    def __init__(
        self,
        link: str,
        max_results: int = None,
        start_date: str = None,
        end_date: str = None,
    ):
        """
        Initializes a SpringerScraper instance with a link and constructs the search URL.

        Args:
            link (str): The link to the Springer search results page.
            max_results (int, optional): Maximum number of results to scrape.
        """
        self.start_date = start_date if start_date else ""
        self.end_date = end_date if end_date else ""
        self.url = link
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

        try:
            # Wait for the main article container - using the correct data-test attribute
            await page.wait_for_selector('[data-test="darwin-search"]', timeout=10000)

            # Get all article elements directly using the list selector
            articles = await page.query_selector_all(
                'ol[data-test="darwin-search"] > li',
            )

            if not articles:
                print(
                    "No articles found. Check if the selector is correct or if content is loaded.",
                )
                return False

            print(f"Found {len(articles)} articles on this page")
            new_articles_found = False

            for article in articles:
                # Stop if max results is reached
                if self.max_results and len(self.articles) >= self.max_results:
                    break

                try:
                    # Extract title and link
                    title_link_element = await article.query_selector("h3 a")
                    if not title_link_element:
                        continue

                    # Extract the link and validate it
                    title_link = await title_link_element.get_attribute("href")
                    if title_link:
                        if not title_link.startswith("http"):
                            title_link = f"https://link.springer.com{title_link}"

                        # Skip if the link is already scraped
                        if title_link in self.unique_links:
                            continue

                        full_link = title_link
                    else:
                        continue

                    # Extract title
                    title = await title_link_element.text_content()

                    # Extract description
                    description_element = await article.query_selector(
                        "p.c-card__description",
                    )
                    description = (
                        await description_element.text_content()
                        if description_element
                        else "No description found"
                    )

                    # Extract authors
                    authors_element = await article.query_selector(
                        ".c-card__section-header",
                    )
                    authors = (
                        await authors_element.text_content()
                        if authors_element
                        else "No authors found"
                    )
                    authors = re.sub(r"\s+in\s.*$", "", authors)
                    authors = authors.replace("...", "").strip()

                    # Extract published date
                    date_element = await article.query_selector(
                        'span.c-meta__item[data-test="published"]',
                    )
                    published_date = (
                        await date_element.text_content()
                        if date_element
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
                    print(f"Error processing article: {e}")
                    continue

            return new_articles_found

        except Exception as e:
            print(f"Failed to load article list: {e}")
            return False

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

        try:
            _csv_file = _load_csv(csv_file_path)
            _csv_file.drop_duplicates(inplace=True)
            print("\n📍 Step 3: Extracting affiliations and full names!")

            for _index, row in _csv_file.iterrows():
                paper_link = row["Link"]
                try:
                    # Add navigation timeout and options
                    await page.goto(
                        paper_link,
                        timeout=10000,
                        wait_until="domcontentloaded",
                    )
                    print(f"    [INFO] Navigating to paper: {paper_link}")

                    # Add small delay for content loading
                    await page.wait_for_timeout(2000)

                    affiliations = []
                    authors_list = []

                    # Try to get affiliations (if fails, continue with empty list)
                    try:
                        await page.wait_for_selector(
                            'xpath=//*[@id="affiliations"]',
                            timeout=5000,
                        )
                        affiliations = await page.query_selector_all(
                            'xpath=//*[@id="Aff1"]/p[1]',
                        )
                    except Exception:
                        print(f"    ⚠️ No affiliations found for: {paper_link}")

                    # Try to get authors (if fails, continue with empty list)
                    try:
                        authors = await page.query_selector_all(
                            '//ul[contains(@class, "c-article-author-list")]//li//a',
                        )

                        for author in authors:
                            author_text = await author.text_content()
                            if (
                                not author_text.startswith(("nAff", "na", "ORCID"))
                                and len(author_text) > 2
                            ):
                                authors_list.append(
                                    author_text.replace("&", ",").strip(),
                                )
                    except Exception:
                        print(f"    ⚠️ No authors found for: {paper_link}")

                    # Update CSV with whatever data we found
                    if affiliations:
                        affil_texts = []
                        for affiliation in affiliations:
                            try:
                                affil_text = await affiliation.text_content()
                                affil_texts.append(affil_text.strip())
                            except Exception:
                                continue
                        if affil_texts:
                            _csv_file.at[_index, "Affiliations"] = "; ".join(
                                affil_texts,
                            )
                            print(f"    ✅ Added {len(affil_texts)} affiliations")

                    if authors_list:
                        _csv_file.at[_index, "Authors"] = ", ".join(authors_list)
                        print(f"    ✅ Added {len(authors_list)} authors")

                    print("---" * 30)

                except Exception as e:
                    print(f"    ⚠️ Error processing {paper_link}: {str(e)}")
                    continue

                # Save progress after each article
                if _index % 10 == 0:  # Save every 10 articles
                    try:
                        _csv_file.to_csv(csv_file_path, index=False)
                        print(f"    💾 Progress saved at article {_index + 1}")
                    except Exception as save_err:
                        print(f"    ⚠️ Error saving progress: {str(save_err)}")

            # Final save
            try:
                _csv_file.to_csv(csv_file_path, index=False)
                print("    💾 Final save completed")
            except Exception as final_save_err:
                print(f"    ⚠️ Error during final save: {str(final_save_err)}")

        except Exception as e:
            print(f"    ❌ Critical error in extraction process: {str(e)}")


async def async_springer(link: str, proxies: Optional[Dict[str, str]] = None) -> None:
    args = [
        "--no-sandbox",
        "--disable-blink-features=AutomationControlled",
        "--disable-infobars",
    ]
    if "dateFrom=" in link and "dateTo=" in link:
        start_date = link.split("dateFrom=")[1].split("&")[0]
        end_date = link.split("dateTo=")[1].split("&")[0]
        # Set defaults if empty
        if not start_date:
            start_date = "1975"
        if not end_date:
            end_date = "2025"
        print(start_date)
        print(end_date)
    else:
        start_date = "1975"
        end_date = "2025"

    async with async_playwright() as pw:
        # In async_springer function
        SpringerScraper_obj = SpringerScraper(
            link=link,
            max_results=None,
            start_date=start_date,
            end_date=end_date,
        )
        browser = await pw.chromium.launch(headless=False, args=args, proxy=proxies)
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
                        2,
                    ):
                        current_end_date = current_start_date + 1
                        print(f"{current_start_date} - {current_end_date}")

                        SpringerScraper_obj = SpringerScraper(
                            link=link,
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
                        if 1000 < total_results <= 2000:
                            print(
                                f"     ✅ [INFO] [Step-2] After filtering range years total results ({total_results}) are greater than 1000. "
                                "Proceeding with scrape old --> new, new --> old.",
                            )
                            # Scrape both newest first and oldest first sorting
                            for sort_order in ["newestFirst", "oldestFirst"]:
                                sorted_url = page.url.replace(
                                    "sortBy=relevance",
                                    f"sortBy={sort_order}",
                                )
                                await page.goto(sorted_url)
                                await SpringerScraper_obj.pagination(page)
                            page = page.url.replace(
                                f"sortBy={sort_order}",
                                "sortBy=relevance",
                            )

                        elif total_results > 2000:
                            print(
                                f"     ✅ [INFO] [Step-3] After filtering range years total results ({total_results}) are greater than 1000. "
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
                    _keyword=link,
                )
                if springer_file_path:
                    # Process affiliations only if articles were found and CSV was saved
                    await SpringerScraper_obj.extract_affiliations_and_full_names(
                        page,
                        springer_file_path,
                    )
                else:
                    print(
                        "No articles found. CSV not saved, skipping affiliation extraction.",
                    )
            except Exception as affil_err:
                print(f"Error extracting affiliations: {affil_err}")
                # Continue despite affiliation extraction errors
            print("\nClosing browser...")
            await browser.close()


if __name__ == "__main__":
    asyncio.run(
        async_springer(
            link="https://link.springer.com/search?new-search=true&query=machine+learning+AND+micro+biolog&dateFrom=&dateTo=&sortBy=relevance",
        ),
    )
