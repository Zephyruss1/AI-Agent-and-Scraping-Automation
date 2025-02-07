from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By
from selenium.webdriver.common.action_chains import ActionChains
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.chrome.options import Options
from webdriver_manager.chrome import ChromeDriverManager
import requests
from typing import List, Optional, Dict
import time

class ArxivScraper:
    def __init__(self, keyword: str = "photonic circuits",
                 date_from: Optional[str] = None, date_to: Optional[str] = None):
        options = Options()
        options.add_argument("--no-sandbox")
        options.add_argument("--disable-dev-shm-usage")
        options.add_argument("--ignore-certificate-errors")
        # options.add_argument("--headless")  # Uncomment to run in headless mode

        self.keyword = keyword
        self.date_from = date_from
        self.date_to = date_to
        self.paper_id_list : List[str] = []
        self.author_list : List[Dict] = []

        self.driver = webdriver.Chrome(service=Service(ChromeDriverManager().install()), options=options)

    def _save_excel(self):
        if not self.author_list:
            print("     [INFO] No author details found!")
            return None

        from openpyxl import Workbook
        wb = Workbook()
        ws = wb.active
        ws.append(["Author",
                   "ID",
                   "PaperID",
                   "LastPaperLink",
                   "LastPaperDate",
                   "AmountOfMentions",
                   "Keyword"])
        for author in self.author_list:
            ws.append([author["Author"],
                       author["ID"],
                       author["PaperID"],
                       author.get("LastPaperLink", "null"),
                       author.get("LastPaperDate", "null"),
                       author.get("AmountOfMentions", "null"),
                       self.keyword
                       ])
        try:
            wb.save(f"arxiv_scraped_data.xlsx")
            print("    ✅ [INFO] Details successfully saved to Excel!")
        except Exception as e:
            print(f"    ❌ [ERROR] An error occurred while saving Excel file: {e}")

    def pagination(self) -> bool:
        print("     [INFO] Checking for next page...")
        next_page_locator = (By.XPATH, '//*[@id="main-container"]/div[2]/nav[2]/a[2]')
        try:
            next_page = WebDriverWait(self.driver, 10).until(EC.element_to_be_clickable(next_page_locator))
            if "disabled" in next_page.get_attribute("class"):
                print("     [INFO] No more pages.")
                return False
            next_page.click()
            time.sleep(1)
            return True
        except Exception as e:
            print(f"     ❌ [ERROR] Pagination failed: {e}")
            return False

    def connect_to_arxiv(self):
        print("\n📍 Step 1: Connecting to Arxiv!")
        if not self.date_from or not self.date_to:
            print("     ⚠️ [INFO] Date range not specified. Searching all papers")
            url = "https://arxiv.org/"
        else:
            url = "https://arxiv.org/search/advanced"

        try:
            self.driver.get(url)
            time.sleep(2)
            print("     ✅ [INFO] Arxiv connection successful!")
        except Exception as e:
            print(f"     ❌ [ERROR] An error occurred: {e}")

        try:
            search_box = None
            if self.date_from or self.date_to:
                search_box = self.driver.find_element(By.XPATH, '//*[@id="terms-0-term"]')
                date_range = self.driver.find_element(By.XPATH, '//*[@id="date-filter_by-3"]')
                ActionChains(self.driver).move_to_element(date_range).click().perform()
                date_from_box = self.driver.find_element(By.XPATH, '//*[@id="date-from_date"]')
                date_to_box = self.driver.find_element(By.XPATH, '//*[@id="date-to_date"]')
                ActionChains(self.driver).move_to_element(date_from_box).click().perform()
                date_from_box.send_keys(self.date_from)
                time.sleep(1)
                ActionChains(self.driver).move_to_element(date_to_box).click().perform()
                date_to_box.send_keys(self.date_to)
                time.sleep(1)
            else:
                search_box = self.driver.find_element(By.XPATH, '//*[@id="header"]/div[2]/form/div/div[1]/input')
            search_box.send_keys(self.keyword)
            search_box.send_keys(Keys.RETURN)
            print("     [INFO] Searching keywords!")
            time.sleep(2)
        except Exception as e:
            print(f"     ❌ [ERROR] An error occurred: {e}")

    def get_paper_ids(self, max_pages_DEBUG_MODE: Optional[int] = None) -> List[str]:
        print("\n📍 Step 2: Scraping paper links!")
        page_count = 0
        while True:
            page_count += 1
            if max_pages_DEBUG_MODE and page_count > max_pages_DEBUG_MODE:
                break

            # Scrape current page
            paper_entries = self.driver.find_elements(By.XPATH, '//*[@id="main-container"]/div[2]/ol/li')
            for paper in paper_entries:
                paper_link = paper.find_element(By.XPATH, './/p[@class="list-title is-inline-block"]/a').get_attribute("href")
                paper_id = paper_link.split("/")[-1]
                self.paper_id_list.append(paper_id)
            print(f"     [INFO] Scraped page {page_count}")

            # Attempt pagination
            if not self.pagination():
                break

        print("     ✅ [INFO] Paper links scraped successfully!")
        return self.paper_id_list

    def get_author_details(self) -> Optional[List[Dict]]:
        print("\n📍 Step 3: Scraping author details!")

        if not self.paper_id_list:
            print("     [INFO] No paper links found!")
            return None

        for paper_id in self.paper_id_list:
            api = f"https://api.semanticscholar.org/v1/paper/arXiv:{paper_id}?include_unknown_references=true"
            try:
                response = requests.get(api)
                response.raise_for_status()  # Raises an HTTPError for bad responses

                data = response.json()
                if 'authors' not in data:
                    print(f"    [WARNING] No authors found for paper {paper_id}")
                    continue

                authors = data["authors"]
                for author in authors:
                    self.author_list.append({
                        "Author": author["name"],
                        "ID": author.get("authorId", "null"),
                        "PaperID": paper_id,
                    })
                print(f"     [INFO] Processed authors for paper {paper_id}")

            except requests.exceptions.HTTPError as e:
                print(f"    ❌ [ERROR] HTTP error for paper {paper_id}: {e}")
            except KeyError as e:
                print(f"    ❌ [ERROR] Missing key in response for paper {paper_id}: {e}")
            except Exception as e:
                print(f"    ❌ [ERROR] An error occurred for paper {paper_id}: {e}")
            time.sleep(1)
        if not self.author_list:
            print("     [INFO] No author details found for any papers.")
            return None

        print("    ✅[INFO] Author details scraped successfully!")
        return self.author_list

    def get_related_papers(self) -> Optional[List[Dict]]:
        print("\n📍 Step 4: Connecting to Scholar and scraping related papers!")
        if not self.author_list:
            print("     [INFO] No author details found!")
            return None

        for author in self.author_list:
            author_id = author.get("ID")
            if author_id == "null":
                print(f"    [INFO] No valid author ID for {author['Author']}. Skipping.")
                continue

            try:
                url = f"https://api.semanticscholar.org/graph/v1/author/{author_id}/papers?fields=url,publicationDate,authors&limit=1000"
                response = requests.get(url)
                response.raise_for_status()
                data = response.json()

                if 'data' not in data or not data['data']:
                    print(f"    ⚠️ [WARNING] No papers found for author {author['Author']}")
                    continue

                paper_info = data['data'][0]
                amount_of_mentions = data['data']
                paper_link = paper_info.get('url', 'null')
                paper_year = paper_info.get('year', 'null')

                print(
                    f"📄 {author['Author']} -> Last Paper Link: {paper_link},"
                    f" Last Paper Date: {paper_year}, Mentions: {len(amount_of_mentions)}")

                for details in self.author_list:
                    if details["Author"] == author["Author"]:  # Ensure we update the correct author
                        details["LastPaperLink"] = paper_link
                        details["LastPaperDate"] = paper_year
                        details["AmountOfMentions"] = len(amount_of_mentions)
                time.sleep(1)
            except requests.exceptions.HTTPError as e:
                print(f"    ❌ [ERROR] HTTP error for author {author['Author']}: {e}")
            except KeyError as e:
                print(f"    ❌ [ERROR] Missing key in response for author {author['Author']}: {e}")
            except Exception as e:
                print(f"    ❌ [ERROR] An error occurred for author {author['Author']}: {e}")

        print("    ✅ [INFO] Related papers scraped successfully!")
        return self.author_list

if __name__ == "__main__":
    scraper = ArxivScraper()
    scraper.connect_to_arxiv()
    scraper.get_paper_ids()
    scraper.get_author_details()
    scraper.driver.quit()
    scraper.get_related_papers()
    import json
    print(json.dumps(scraper.author_list, indent=4))
    scraper._save_excel()
