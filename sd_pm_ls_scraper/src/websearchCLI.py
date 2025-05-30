import asyncio
import os
import re
from dataclasses import dataclass, field
from typing import List, Literal, Optional

import google.auth.transport.requests
import gspread
import pandas as pd
import requests
from browser_use import Agent, Browser, BrowserConfig
from dotenv import load_dotenv
from google.oauth2 import service_account
from langchain_openai import ChatOpenAI
from oauth2client.service_account import ServiceAccountCredentials
from pydantic import BaseModel
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity

load_dotenv()


def load_spreadsheet(spreadsheet_link: str = None) -> pd.DataFrame:
    """Load the spreadsheet from the given link and download as .xlsx."""
    print("     💭 [INFO] Checking spreadsheet link...")
    if spreadsheet_link is None:
        raise ValueError("Spreadsheet link is required.")
    if not isinstance(spreadsheet_link, str):
        raise ValueError("Spreadsheet link must be a string.")

    # Extract spreadsheet ID
    match = re.search(r"/d/([a-zA-Z0-9-_]+)", spreadsheet_link)
    if not match:
        raise ValueError("Invalid Google Sheets URL.")
    spreadsheet_id = match.group(1)

    # Setup credentials
    json_keyfile = (
        "/root/arxiv-and-scholar-scraping/aiagenttest-455613-a4d701f3b9ce.json"
    )
    scope = [
        "https://spreadsheets.google.com/feeds",
        "https://www.googleapis.com/auth/drive",
        "https://www.googleapis.com/auth/drive.readonly",
    ]

    creds = ServiceAccountCredentials.from_json_keyfile_name(json_keyfile, scope)
    client = gspread.authorize(creds)

    # Load sheet data into DataFrame
    sheet = client.open_by_url(spreadsheet_link).sheet1
    data = sheet.get_all_records()
    if data:
        print("       ✅ [INFO] Spreadsheet loaded successfully!")
    else:
        print("       ❌ [INFO] Spreadsheet loading failed!")

    # Download as .xlsx
    # Refresh token using google-auth
    creds2 = service_account.Credentials.from_service_account_file(
        json_keyfile,
        scopes=["https://www.googleapis.com/auth/drive.readonly"],
    )
    auth_req = google.auth.transport.requests.Request()
    creds2.refresh(auth_req)
    token = creds2.token

    # Make download request
    export_url = (
        f"https://docs.google.com/spreadsheets/d/{spreadsheet_id}/export?format=xlsx"
    )
    headers = {"Authorization": f"Bearer {token}"}
    response = requests.get(export_url, headers=headers)

    if response.status_code == 200:
        with open(
            "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/downloaded_sheet.xlsx",
            "wb",
        ) as f:
            f.write(response.content)
        print("     ✅ [INFO] Spreadsheet downloaded as downloaded_sheet.xlsx")
    else:
        print(f"⚠️ Failed to download sheet: {response.status_code}")

    return (data, True) if data else (False, False)


def convert_excel_to_csv() -> pd.DataFrame:
    """Convert the excel file to a CSV file."""
    file_path = (
        "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/downloaded_sheet.xlsx"
    )
    if not file_path.endswith(".xlsx"):
        raise ValueError("File name must end with .xlsx")

    xlsx_file = pd.read_excel(file_path)
    try:
        csv_file = xlsx_file.to_csv(
            file_path.replace(".xlsx", ".csv"),
            index=False,
            header=False,
        )
        print(
            f"    ✅ [INFO] Excel file converted to CSV: {file_path.replace('.xlsx', '.csv')}",
        )
    except Exception as err:
        print("    ❌ [INFO] Failed to convert Excel file to CSV.")
        raise Exception(f"Error: {err}") from err
    return csv_file


def _load_csv(file_name: str) -> object:
    """Load the excel file."""
    csv = pd.read_csv(file_name)
    return csv


@dataclass(frozen=True)
class PerplexityConfig:
    """Perplexity API Config Class."""

    url: str = "https://api.perplexity.ai/chat/completions"
    model: Literal["sonar-pro", "sonar-deep-research"] = "sonar-pro"
    search_context_size: str = "Low"  # Low/Medium/High
    temperature: float = 0.2
    api_key: str = field(default_factory=lambda: os.getenv("PERPLEXITY_API_KEY", ""))

    def get_headers(self) -> dict:
        """Dynamically returning headers."""
        if not self.api_key:
            raise ValueError("PERPLEXITY_API_KEY is not set.")
        return {"Authorization": f"Bearer {self.api_key}"}


# @dataclass(frozen=True)
# class BrowserConfig:
#     """Browser Config Class."""

#     headless: bool = True
#     disable_security: bool = True
#     user_agent: str = "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/58.0.3029.110 Safari/537.3"
#     model: Literal["gpt-4o", "gpt-4o-mini"] = "gpt-4o"
#     temprature: float = 0.2

#     def get_headers(self) -> dict:
#         """Dynamically returning headers."""
#         return {"User-Agent": self.user_agent}


class EmailFormat(BaseModel):
    """Answer Format for Perplexity API."""

    email_adress: str


class JobTitleFormat(BaseModel):
    """ "Job Title Format for Perplexity API."""

    job_title: str


class WebSearch:
    def __init__(
        self,
        name: str,
        csv_file: pd.DataFrame,
        _keyword: str,
        config: Optional[PerplexityConfig] = None,
    ) -> None:
        if config is None:
            config = PerplexityConfig()
        """Implement PerplexityConfig to WebSearch."""
        self.config = config
        self.author_name = name
        self.csv_file = csv_file
        self.keyword_index = _keyword
        self.browser = Browser(
            config=BrowserConfig(
                headless=BrowserConfig.headless,
                disable_security=BrowserConfig.disable_security,
            ),
        )
        import re

        # This will search for text within quotes
        match = re.search(r'"([^"]+)"', self.keyword_index)
        if match:
            self.keyword_index = match.group(1)

    def perplexity_search_for_email(self) -> List[str]:
        """Search for email addresses for the provided author name."""
        print("\n📍 Step 8: [Perplexity] Extracting Email Addresses!")
        payload = {
            "model": self.config.model,
            "search_context_size": self.config.search_context_size,
            "messages": [
                {
                    "role": "system",
                    "content": (
                        "You are a web searcher assistant. The user will provide an author name."
                        f"Your task is: Search email addresses in the {self.keyword_index} field for provided author name."
                        "If you find no email addresses, return 'None'. "
                        "Output only the emails or 'None'—no additional explanations."
                    ),
                },
                {
                    "role": "user",
                    "content": f"{self.author_name} email adress",
                },
            ],
            "temperature": self.config.temperature,
            "response_format": {
                "type": "json_schema",
                "json_schema": {"schema": EmailFormat.model_json_schema()},
            },
        }

        response = requests.post(
            self.config.url,
            headers=self.config.get_headers(),
            json=payload,
        )
        print(f"Prompt: {payload.get('messages')[1].get('content')}")
        if response.status_code == 200:
            try:
                response_json = response.json()
                print(response_json["choices"][0]["message"]["content"])
                list_of_emails = response_json["choices"][0]["message"][
                    "content"
                ].split("\n")
                if list_of_emails and list_of_emails[0] != "None":
                    print(f"    ✅ [INFO] Email addresses listed: {list_of_emails}")
                    return list_of_emails
                else:
                    print("    ❌ [INFO] No email addresses found.")
                    return ["None"]
            except requests.JSONDecodeError:
                print("Error: Response content is not valid JSON")
                return ["None"]
        else:
            print(f"Error: Received status code {response.status_code}")
            return ["None"]

    def perplexity_search_for_job_title(self) -> List[str]:
        """Search for job title for the provided author name."""
        import json

        print("\n📍 Step 8: [Perplexity] Extracting job titles!")
        payload = {
            "model": self.config.model,
            "search_context_size": self.config.search_context_size,
            "messages": [
                {
                    "role": "system",
                    "content": (
                        "You are a web search assistant. The user will provide an author name. "
                        f"Search the {self.keyword_index} field for that author's job title. "
                        "If no job title is found, return 'None'. "
                        "Output only the job title or 'None' with no additional commentary."
                    ),
                },
                {
                    "role": "user",
                    "content": f"{self.author_name} job title",
                },
            ],
            "temperature": self.config.temperature,
            "response_format": {
                "type": "json_schema",
                "json_schema": {"schema": JobTitleFormat.model_json_schema()},
            },
        }

        response = requests.post(
            self.config.url,
            headers=self.config.get_headers(),
            json=payload,
        )
        print(f"Prompt: {payload.get('messages')[1].get('content')}")
        if response.status_code == 200:
            try:
                response_json = response.json()
                raw_content = response_json["choices"][0]["message"]["content"]
                print(f"Raw Response: {raw_content}")

                # Parse the JSON content to extract the job title
                try:
                    parsed_content = json.loads(raw_content)
                    job_title = parsed_content.get("job_title", "None")
                    if job_title != "None":
                        print(f"    ✅ [INFO] Job title listed: [{job_title}]")
                        return [job_title]
                    else:
                        print("    ❌ [INFO] No Job title found.")
                        return ["None"]
                except json.JSONDecodeError:
                    print("Error: Response content is not valid JSON")
                    return ["None"]

            except KeyError:
                print("Error: Unexpected response format")
                return ["None"]
        else:
            print(f"Error: Received status code {response.status_code}")
            return ["None"]

    def browser_use_find_email(self) -> List[str]:
        """Search for email addresses using the browser-use."""
        print("\n📍 Step 8: [Browser-Use] Extracting Email Addresses!")

        agent = Agent(
            task=f"""
            1. Search '{self.author_name} email address in the {self.keyword_index} field' and enter.
            2. Output only the emails or 'None'—no additional explanations.
            """,
            llm=ChatOpenAI(model=BrowserConfig.model),
            browser=self.browser,
            max_failures=1,
            max_actions_per_step=20,
        )
        loop = asyncio.get_event_loop()
        result = loop.run_until_complete(agent.run())

        if hasattr(result, "final_result"):
            res = result.final_result().split("\n")
        else:
            res = ["None"]
        return res

    def browser_use_find_job(self) -> List[str]:
        """Search for job title using the browser-use."""
        print("\n📍 Step 8: [Browser-Use] Extracting job titles!")

        agent = Agent(
            task=f"""
            1. Search '{self.author_name} job title in the {self.keyword_index} field' and enter.
            2. Output only the job title or 'None'—no additional explanations.
            """,
            llm=ChatOpenAI(model=BrowserConfig.model),
            browser=self.browser,
            max_failures=1,
            max_actions_per_step=20,
        )
        loop = asyncio.get_event_loop()
        result = loop.run_until_complete(agent.run())

        if hasattr(result, "final_result"):
            res = result.final_result().split("\n")
        else:
            res = ["None"]
        return res


class FindSimilarity:
    """Find similarity between the author names and emails."""

    def __init__(self, csv_file: pd.DataFrame, csv_file_path: str) -> None:
        """
        Initialize with a DataFrame and its file path.

        :param csv_file: The DataFrame object.
        :param csv_file_path: The path of the CSV file (including the file name).
        """
        self.csv_file = csv_file
        self.csv_file_path = csv_file_path  # Store the file path
        self.headers = [header for header in self.csv_file.columns]

    def preprocess_emails(self, email_list: List) -> Optional[str]:
        """Preprocess the email address to extract email author name."""
        return " ".join([email.split("@")[0].lower() for email in email_list])

    def find_email_author_and_save(self, list_of_emails: List[str]) -> None:
        """Find the email address and author name using Cosine Similarity."""
        print("\n📍 Step 9: Finding Similarity and Saving to Csv File!")

        if "email" not in self.headers:
            self.csv_file["email"] = ""

        for email in list_of_emails:
            if "None" in email or "protected" in email or "*" in email:
                list_of_emails = "None"
                continue

            print(f"🔄🕒💭 Processing email: {email}")
            doc1 = self.preprocess_emails([email])
            match_found = False

            for _index, row in self.csv_file.iterrows():
                author_name = row["Authors"]
                if not author_name:
                    continue

                doc2 = str(author_name)

                vectorizer = TfidfVectorizer(
                    analyzer="char",
                    ngram_range=(1, 2),
                    lowercase=True,
                )
                tfidf_matrix = vectorizer.fit_transform([doc1, doc2])

                cosine_sim = cosine_similarity(tfidf_matrix[0:1], tfidf_matrix[1:2])

                try:
                    if cosine_sim[0][0] > 0.6:
                        print(
                            f"    ✅ [INFO] Match Found: {author_name} | {cosine_sim[0][0]}",
                        )
                        current_email = row.get("email", "")

                        if pd.notna(current_email) and current_email:
                            if email in current_email:
                                print(
                                    "     ⚠️ [INFO] Email already exists in the CSV file. Continuing to next email.",
                                )
                                break
                            else:
                                try:
                                    new_email = f"{current_email}, {email}"
                                    self.csv_file.at[_index, "email"] = new_email
                                    print(
                                        f"     🔄 [INFO] Updating email: {current_email} -> {email}",
                                    )
                                except Exception as err:
                                    raise Exception(f"Error: {err}") from err
                        else:
                            self.csv_file.at[_index, "email"] = email
                            [
                                print(
                                    f"     🔄 [INFO] Adding email: {email} to {author_name}",
                                ),
                            ]
                        match_found = True
                except Exception as err:
                    raise Exception(f"Error: {err}") from err

            if not match_found:
                print(f"    ❌ [INFO] No match found for email: {email}")
                new_row = {"Authors": "", "email": email}
                self.csv_file = pd.concat(
                    [self.csv_file, pd.DataFrame([new_row])],
                    ignore_index=True,
                )
            print("---" * 30)
        # Save the updated DataFrame to the same path
        self.csv_file.to_csv(self.csv_file_path, index=False)
        print(f"✅ [INFO] Results saved to: '{self.csv_file_path}'")

    def find_job_title_and_save(
        self,
        list_of_jobs: List[str],
        current_index: int,
    ) -> None:
        """Find the job title and author name using Cosine Similarity, and update only the current author."""
        print("\n📍 Step 9: Finding Similarity and Saving to Csv File!")

        # Ensure 'job_title' column exists
        if "job_title" not in self.headers:
            self.csv_file["job_title"] = ""

        for job_title in list_of_jobs:
            # Skip if job_title is None or has the string "None"
            if not job_title or "None" in job_title:
                continue

            print(f"🔄🕒💭 Processing job title: {job_title}")

            # Retrieve the current row for the specified author index
            author_name = self.csv_file.at[current_index, "Authors"]
            current_job_title = self.csv_file.at[current_index, "job_title"]

            if pd.notna(current_job_title) and current_job_title:
                if job_title in current_job_title:
                    print(
                        "     ⚠️ [INFO] Job title already exists for the current author.",
                    )
                    continue
                else:
                    # Append the new job title to the existing string
                    new_job_title = f"{current_job_title}, {job_title}"
                    self.csv_file.at[current_index, "job_title"] = new_job_title
                    print(
                        f"     🔄 [INFO] Updating job title: {current_job_title} -> {new_job_title}",
                    )
            else:
                # If no job title exists, simply assign it
                self.csv_file.at[current_index, "job_title"] = job_title
                print(f"     🔄 [INFO] Adding job title: {job_title} to {author_name}")

            print("---" * 10)

        # Save the updated CSV file
        self.csv_file.to_csv(self.csv_file_path, index=False)
        print(f"✅ [INFO] Results saved to: '{self.csv_file_path}'")


def fill_empty_emails_with_search():
    csv_paths = [
        # "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/cleaned_sciencedirect_results.csv",
    ]

    for csv_path in csv_paths:
        _csv = _load_csv(file_name=csv_path)

        headers = [header for header in _csv.columns]
        if "email" not in headers:  # Ensure the column name matches
            _csv["email"] = ""  # Add the 'email' column if missing

        # Step 5: Search for email addresses using the AI Search
        for _index, row in _csv.iterrows():
            keyword = row["Keyword"]

            author_name = row["Authors"]
            email = row["email"]
            if pd.notna(email) and email:
                print(f"Email already exists for author: {author_name}")
                continue
            print(f"Processing author: {author_name}")

            web_search = WebSearch(
                name=str(author_name),
                csv_file=_csv,
                _keyword=keyword,
            )
            list_of_emails = web_search.ai_search()

            # Step 6: Find similarity between the author names and emails
            df_csv = _load_csv(file_name=csv_path)
            similarity_finder = FindSimilarity(
                csv_file=df_csv,
                csv_file_path=f"{csv_path}",
            )
            similarity_finder.find_email_author_and_save(list_of_emails)
            print("-----" * 15)


def fill_empty_jobs_with_search():
    csv_paths = [
        # "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/cleaned_sciencedirect_results.csv",
    ]

    for csv_path in csv_paths:
        _csv = _load_csv(file_name=csv_path)

        headers = [header for header in _csv.columns]
        if "job_title" not in headers:  # Ensure the column name matches
            _csv["job_title"] = ""  # Add the 'job_title' column if missing

    # Step 7: Search for job title addresses using the AI Search
    for _index, row in _csv.iterrows():
        keyword = row["Keyword"]

        author_name = row["Authors"]
        job = row["job_title"]
        if pd.notna(job) and job:
            print(f"Job title already exists for author: {author_name}")
            continue
        print(f"Processing author: {author_name}")

        web_search = WebSearch(name=str(author_name), csv_file=_csv, _keyword=keyword)
        list_of_jobs = web_search.ai_search()

        # Step 8: Find similarity between the author names and emails
        df_csv = _load_csv(file_name=csv_path)
        similarity_finder = FindSimilarity(csv_file=df_csv, csv_file_path=f"{csv_path}")
        similarity_finder.find_job_title_and_save(list_of_jobs, _index)
        print("-----" * 15)


if __name__ == "__main__":
    print("=========================🚀 AI Agent TEST 🚀=========================")
    print("---------------------------------------------------------------")
    print("📍 Step 1: Load your spreadsheet!")
    spreadsheet_link = input(str("Enter the link to your spreadsheet: "))
    print("---------------------------------------------------------------")
    _, condition = load_spreadsheet(spreadsheet_link)
    print("---------------------------------------------------------------")
    print("📍 Step 2: Convert the spreadsheet to CSV!")
    convert_excel_to_csv()
    print("---------------------------------------------------------------")
    print("📍 Step 3: AI Agent 🤖 :")
    model_choice = input(
        str(
            "Enter the model you want to use (1 or 2):\n"
            "1. ChatGPT [Browser-Use]\n"
            "2. Perplexity\n",
        ),
    )
    print("📍 Step 4: AI Agent 🔧 :")
    print("---------------------------------------------------------------")
    if model_choice == "1":
        model_config_choice_gpt = input(
            str(
                "Enter the model you want to use (1 or 2):\n"
                "1. gpt-4o\n"
                "2. gpt-4o-mini\n"
                "->",
            ),
        )
        if model_config_choice_gpt == "1":
            BrowserConfig.model = "gpt-4o"
            print(f"     ✅ [INFO] Using {BrowserConfig.model} model")
        elif model_config_choice_gpt == "2":
            BrowserConfig.model = "gpt-4o-mini"
            print(f"     ✅ [INFO] Using {BrowserConfig.model} model")
        else:
            print("Invalid choice. Exiting.")
    elif model_choice == "2":
        model_config_choice_perplexity = input(
            str(
                "Enter the model you want to use (1 or 2):\n"
                "1. sonar-pro\n"
                "2. sonar-deep-research\n",
            ),
        )
        if model_config_choice_perplexity == "1":
            PerplexityConfig.model = "sonar-pro"
            print(f"     ✅ [INFO] Using {PerplexityConfig.model} model")
        elif model_config_choice_perplexity == "2":
            PerplexityConfig.model = "sonar-deep-research"
            print(f"     ✅ [INFO] Using {PerplexityConfig.model} model")
        else:
            print("Invalid choice. Exiting.")
    else:
        print("Invalid choice. Exiting.")
    print("---------------------------------------------------------------")
    purpose_choice = input(
        str(
            "Enter the purpose of your search (1 or 2):\n"
            "1. Find email addresses\n"
            "2. Find job titles\n",
        ),
    )
    print("---------------------------------------------------------------")
    if purpose_choice == "1":
        print("     ✅ [INFO] Searching for email addresses")
        fill_empty_emails_with_search()
    elif purpose_choice == "2":
        print("     ✅ [INFO] Searching for job titles")
        fill_empty_jobs_with_search()
    else:
        print("Invalid choice. Exiting.")
