from typing import List, Dict, Optional
from openpyxl import Workbook
import asyncio
from playwright.async_api import async_playwright
import aiohttp

def _save_excel(author_list: List[Dict], keyword: str) -> None:
    if not author_list:
        print("     [INFO] No author details found!")
        return

    wb = Workbook()
    ws = wb.active
    ws.append(["Authors", "AuthorIDs", "PaperID", "LastPaperLink", "AmountOfMentions", "Keyword"])

    for author in author_list:
        # Extract the list of authors
        authors = author.get("Authors", [])
        if not isinstance(authors, list):
            authors = [authors]  # Ensure it's a list even if it's a single author

        # Extract the list of AuthorIDs
        author_ids = author.get("AuthorIDs", [])
        if not isinstance(author_ids, list):
            author_ids = [author_ids]  # Ensure it's a list even if it's a single ID

        # Ensure the number of authors and AuthorIDs match
        if len(authors) != len(author_ids):
            print(f"     ⚠️ [WARNING] Mismatch between authors and AuthorIDs for paper {author.get('PaperID', 'unknown')}")
            continue

        # Save each author separately
        for i in range(len(authors)):
            ws.append([
                authors[i],  # Single author
                author_ids[i],  # Corresponding AuthorID
                author.get("PaperID", "null"),  # PaperID (repeated for each author)
                ", ".join(author["LastPaperLink"]) if isinstance(author.get("LastPaperLink"), list) else author.get("LastPaperLink", "null"),  # LastPaperLink
                author.get("AmountOfMentions", "null"),  # AmountOfMentions
                keyword  # Keyword
            ])

    try:
        wb.save("arxiv_scraped_data.xlsx")
        print("    ✅ [INFO] Details successfully saved to Excel!")
    except Exception as e:
        print(f"    ❌ [ERROR] An error occurred while saving Excel file: {e}")

class ArxivScraper:
    def __init__(self, keyword: str = "photonic circuits", date_from: Optional[str] = None, date_to: Optional[str] = None):
        self.page = None
        self.keyword = keyword
        self.date_from = date_from
        self.date_to = date_to
        self.papers: List[Dict] = []

    async def pagination(self) -> bool:
        print("     [INFO] Checking for next page...")
        next_page_locator = self.page.locator('xpath=//*[@id="main-container"]/div[2]/nav[1]/a[2]')
        try:
            await next_page_locator.wait_for(state='visible', timeout=10000)
            next_class = await next_page_locator.get_attribute("class")
            if next_class and "disabled" in next_class:
                print("     [INFO] No more pages.")
                return False
            await next_page_locator.click()
            await asyncio.sleep(1)
            return True
        except Exception as e:
            print(f"     ❌ [ERROR] Pagination failed: {e}")
            return False

    async def connect_to_arxiv(self):
        print("\n📍 Step 1: Connecting to Arxiv!")
        if not self.date_from or not self.date_to:
            print("     ⚠️ [INFO] Date range not specified. Searching all papers")
            url = "https://arxiv.org/"
        else:
            url = "https://arxiv.org/search/advanced"

        try:
            await self.page.goto(url, wait_until='load')
            await asyncio.sleep(2)
            print("     ✅ [INFO] Arxiv connection successful!")
        except Exception as e:
            print(f"     ❌ [ERROR] An error occurred during navigation: {e}")

        try:
            if self.date_from or self.date_to:
                search_box = self.page.locator('xpath=//*[@id="terms-0-term"]')
                date_range = self.page.locator('xpath=//*[@id="date-filter_by-3"]')
                await date_range.click()
                date_from_box = self.page.locator('xpath=//*[@id="date-from_date"]')
                date_to_box = self.page.locator('xpath=//*[@id="date-to_date"]')
                await date_from_box.click()
                await date_from_box.fill(self.date_from)
                await asyncio.sleep(1)
                await date_to_box.click()
                await date_to_box.fill(self.date_to)
                await asyncio.sleep(1)
            else:
                search_box = self.page.locator('xpath=//*[@id="header"]/div[2]/form/div/div[1]/input')
            await search_box.fill(self.keyword)
            await search_box.press('Enter')
            print("     [INFO] Searching keywords!")
            await asyncio.sleep(2)
        except Exception as e:
            print(f"     ❌ [ERROR] An error occurred while searching: {e}")

    async def get_paper_ids(self, max_pages_DEBUG_MODE: Optional[int] = None) -> List[Dict]:
        print("\n📍 Step 2: Scraping paper links!")
        page_count = 0
        if max_pages_DEBUG_MODE:
            print(f"     ⚠️ [INFO] Debug mode enabled. Scraping only {max_pages_DEBUG_MODE} pages.")

        while True:
            page_count += 1
            if max_pages_DEBUG_MODE and page_count > max_pages_DEBUG_MODE:
                break

            try:
                paper_entries = self.page.locator('xpath=//*[@id="main-container"]/div[2]/ol/li')
            except Exception as e:
                print(f"     ❌ [ERROR] Waiting for paper entries failed: {e}")
                break

            papers = await paper_entries.all()
            if not papers:
                print("     [INFO] No paper entries found on this page.")
                break

            for paper in papers:
                try:
                    paper_link = await (paper.locator('xpath=//p[@class="list-title is-inline-block"]/a').get_attribute("href"))
                    author_locator = paper.locator('xpath=//p[contains(@class, "authors")]')
                    author_text = await author_locator.text_content()

                    paper_id = None
                    if paper_link:
                        paper_id = paper_link.split("/")[-1]

                    if not paper_id:
                        print(f"     ⚠️ [WARNING] Paper ID missing for one entry, skipping...")
                        continue

                    authors = []
                    if author_text:
                        author_names = author_text.split(":")[1].split(",")
                        authors = [author.strip() for author in author_names]
                        print(authors)
                    self.papers.append({
                        "PaperID": paper_id,
                        "Authors": authors,
                        "AuthorIDs": authors,
                        "LastPaperLink": [],
                        "AmountOfMentions": 0,
                        "Keyword": self.keyword
                    })

                except Exception as e:
                    print(f"     ❌ [ERROR] Failed to process a paper entry: {e}")

            print(f"     [INFO] Scraped page {page_count}")
            if not await self.pagination():
                break

        print("     ✅ [INFO] Paper links scraped successfully!")
        return self.papers

    async def get_author_details(self) -> Optional[List[Dict]]:
        print("\n📍 Step 3: Scraping author details!")
        if not self.papers:
            print("     [INFO] No paper links found!")
            return None

        async with aiohttp.ClientSession() as session:
            for paper in self.papers:
                await asyncio.sleep(2)
                api = f"https://api.semanticscholar.org/v1/paper/arXiv:{paper['PaperID']}?include_unknown_references=true"
                try:
                    async with session.get(api) as response:
                        response.raise_for_status()
                        data = await response.json()

                        if 'authors' not in data:
                            print(f"    [WARNING] No authors found for paper {paper['PaperID']}")
                            continue

                        author_ids = [author.get("authorId", "null") for author in data["authors"]]
                        paper["AuthorIDs"] = author_ids

                        print(f"     [INFO] Processed authors for paper {paper['PaperID']}")
                except aiohttp.ClientResponseError as e:
                    print(f"    ❌ [ERROR] HTTP error for paper {paper['PaperID']}: {e}")
                except KeyError as e:
                    print(f"    ❌ [ERROR] Missing key in response for paper {paper['PaperID']}: {e}")
                except Exception as e:
                    print(f"    ❌ [ERROR] An error occurred for paper {paper['PaperID']}: {e}")
                await asyncio.sleep(1)

        if not any(paper["AuthorIDs"] for paper in self.papers):
            print("     [INFO] No author details found for any papers.")
            return None

        print("    ✅ [INFO] Author details scraped successfully!")
        return self.papers

async def main():
    async with async_playwright() as pw:
        browser = await pw.chromium.launch(headless=False)
        context = await browser.new_context(viewport={'width': 1280, 'height': 720})
        page = await context.new_page()

        # Create the ArxivScraper instance and assign the open page
        arxiv_scraper = ArxivScraper(keyword="photonic circuits", date_from="2022-01-01", date_to="2022-01-31")
        arxiv_scraper.page = page

        # Step 1: Connect and search
        await arxiv_scraper.connect_to_arxiv()

        # Step 2: Get paper IDs (use max_pages_DEBUG_MODE=1 for testing)
        await arxiv_scraper.get_paper_ids(max_pages_DEBUG_MODE=1)

        # Step 3: Get author details
        await arxiv_scraper.get_author_details()

        # Save the scraped details to an Excel file
        _save_excel(arxiv_scraper.papers, arxiv_scraper.keyword)

if __name__ == "__main__":
    asyncio.run(main())