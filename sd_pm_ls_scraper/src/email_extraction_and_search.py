import os
import re
from dataclasses import dataclass, field
from typing import List, Optional

import pandas as pd
import pdfplumber
import requests
from dotenv import load_dotenv
from openai import OpenAI
from pydantic import BaseModel
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity

load_dotenv()


def _load_csv(file_name: str) -> object:
    """Load the excel file."""
    csv = pd.read_csv(file_name)
    return csv


class LoadPDF:
    """Load PDFs from a specified path and extract text from them."""

    def __init__(self, pdf_path: str):
        self.pdf_path = pdf_path
        print("\n📍 Step 7: Loading PDFs and Extracting Text!")

    def extract_text_from_pdf(self, pdf_file: str) -> str:
        """Extract text from a single PDF file, handling errors gracefully."""
        all_text = ""
        try:
            with pdfplumber.open(os.path.join(self.pdf_path, pdf_file)) as pdf:
                print(f"     [INFO] Extracting text from {pdf_file}")
                for page in pdf.pages:
                    extracted_text = page.extract_text()
                    if extracted_text:  # Check for empty pages
                        all_text += extracted_text + "\n"
        except Exception as e:
            print(
                f"     ❌ [ERROR] Could not process {pdf_file}. Skipping it. Error: {e}"
            )

        return all_text

    def list_pdf_files(self) -> List[str]:
        """List all PDF files in the specified directory."""
        try:
            pdf_files = [
                file for file in os.listdir(self.pdf_path) if file.endswith(".pdf")
            ]
            print(f"    ✅ [INFO] Found {len(pdf_files)} PDF file(s) in the directory.")
            return pdf_files
        except FileNotFoundError as e:
            print(f"     ❌ [ERROR] Directory not found: {self.pdf_path}. Error: {e}")
            return []

    def __len__(self) -> int:
        """Return the number of PDF files in the directory."""
        return len(self.list_pdf_files())


@dataclass(frozen=True)
class ChatGPTConfig:
    """ChatGPT API Config Class."""

    model: str = "gpt-4o"
    api_key: str = field(default_factory=lambda: os.getenv("OPENAI_API_KEY", ""))
    temperature: float = 0.2


class ExtractEmails:
    """Implement ChatGPTConfig to ExtractEmails."""

    def __init__(self, config: Optional[ChatGPTConfig] = None) -> None:
        if config is None:
            config = ChatGPTConfig()
        self.config = config

    def chatgpt_response(self, pdf_texts: str) -> List[str]:
        """Extract email addresses from the text using ChatGPT."""
        print("\n📍 Step 8: [ChatGPT] Extracting Email Addresses!")
        client = OpenAI()
        client.api_key = self.config.api_key
        completion = client.chat.completions.create(
            model=self.config.model,
            store=True,
            messages=[
                {
                    "role": "system",
                    "content": (
                        "You are a data extraction assistant. The user will provide text from a PDF. "
                        "Your task is to extract and return only the email addresses that appear in the text. "
                        "Do not return any other information. If there are no email addresses, return 'None'."
                    ),
                },
                {"role": "user", "content": pdf_texts[:128000]},
            ],
            temperature=self.config.temperature,
        )

        # Use regex to filter valid emails from the response
        raw_response = completion.choices[0].message.content
        list_of_emails = re.findall(
            r"[a-zA-Z0-9._%+-]+@[a-zA-Z0-9.-]+\.[a-zA-Z]{2,}", raw_response
        )

        if list_of_emails:
            print(f"    ✅ [INFO] Email addresses listed: {list_of_emails}")
            return list_of_emails
        else:
            print("    ❌ [INFO] No email addresses found.")
            return ["None"]


@dataclass(frozen=True)
class PerplexityConfig:
    """Perplexity API Config Class."""

    url: str = "https://api.perplexity.ai/chat/completions"
    model: str = "sonar-pro"
    search_context_size: str = "Low"  # Low/Medium/High
    temperature: float = 0.2
    api_key: str = field(default_factory=lambda: os.getenv("PERPLEXITY_API_KEY", ""))

    def get_headers(self) -> dict:
        """Dynamically returning headers."""
        if not self.api_key:
            raise ValueError("PERPLEXITY_API_KEY is not set.")
        return {"Authorization": f"Bearer {self.api_key}"}


class AnswerFormat(BaseModel):
    """Answer Format for Perplexity API."""

    email_adress: str


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

        if "AND" in self.keyword_index:
            self.keyword_index = self.keyword_index.split("AND")[0].strip()

    def ai_search(self) -> str:
        """Search for email addresses using the AI search."""

        def build_search_pipeline(*methods):
            """Build a search pipeline that tries methods in sequence until one succeeds."""

            def pipeline():
                for method in methods:
                    result = method()
                    if result and result[0] != "None":
                        return result
                return ["None"]

            return pipeline

        def perplexity_search() -> List[str]:
            """Search for email addresses for the provided author name."""
            print("\n📍 Step 8: [Perplexity] Extracting Email Addresses!")
            payload = {
                "model": self.config.model,
                "search_context_size": self.config.search_context_size,
                "messages": [
                    {
                        "role": "system",
                        "content": (
                            "You are a web searcher assistant. The user will provide an author name. "
                            f"Your task is: Search email addresses in the {self.keyword_index} field for provided author name."
                            "If you find no email addresses, return 'None'. "
                            "Output only the emails or 'None'—no additional explanations."
                        ),
                    },
                    {"role": "user", "content": f"{self.author_name} email adress"},
                ],
                "temperature": self.config.temperature,
                "response_format": {
                    "type": "json_schema",
                    "json_schema": {"schema": AnswerFormat.model_json_schema()},
                },
            }

            response = requests.post(
                self.config.url, headers=self.config.get_headers(), json=payload
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

        def browser_use() -> List[str]:
            """Search for email addresses using the browser-use."""
            print("\n📍 Step 8: [Browser-Use] Extracting Email Addresses!")

            try:
                import asyncio

                from browser_use import Agent, Browser, BrowserConfig
                from langchain_openai import ChatOpenAI
            except ImportError as err:
                raise ImportError(
                    "Please install the required packages to run this function."
                ) from err

            browser = Browser(
                config=BrowserConfig(headless=True, disable_security=True)
            )

            agent = Agent(
                task=f"""
                1. Go to Google.com.
                2. Search '{self.author_name} email address in the {self.keyword_index} field' and enter.
                3. Output only the emails or 'None'—no additional explanations.
                """,
                llm=ChatOpenAI(model="gpt-4o"),
                browser=browser,
            )
            loop = asyncio.get_event_loop()
            result = loop.run_until_complete(agent.run())
            return result.final_result().split("\n")

        # Build the search pipeline with the desired methods
        search_pipeline = build_search_pipeline(perplexity_search)

        # Execute the pipeline
        return search_pipeline()


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

            print(f"Processing email: {email}")
            doc1 = self.preprocess_emails([email])
            match_found = False

            for _index, row in self.csv_file.iterrows():
                author_name = row["Authors"]
                if not author_name:
                    continue

                doc2 = str(author_name)

                vectorizer = TfidfVectorizer(
                    analyzer="char", ngram_range=(1, 2), lowercase=True
                )
                tfidf_matrix = vectorizer.fit_transform([doc1, doc2])

                cosine_sim = cosine_similarity(tfidf_matrix[0:1], tfidf_matrix[1:2])

                if cosine_sim[0][0] > 0.6:
                    print(f"Match Found: {author_name} | {cosine_sim[0][0]}")
                    current_email = self.csv_file.at[_index, "email"]
                    if pd.notna(current_email) and current_email:
                        if email not in current_email:
                            if email not in row["email"]:
                                new_email = f"{current_email}, {email}"
                                self.csv_file.at[_index, "email"] = new_email
                    else:
                        self.csv_file.at[_index, "email"] = email
                    match_found = True

            if not match_found:
                print(f"No match found for email: {email}")
                new_row = {"Authors": "", "email": email}
                self.csv_file = pd.concat(
                    [self.csv_file, pd.DataFrame([new_row])], ignore_index=True
                )

        # Save the updated DataFrame to the same path
        self.csv_file.to_csv(self.csv_file_path, index=False)
        print(f"    ✅ [INFO] Results saved to '{self.csv_file_path}'")


def extract_emails_from_pdf():
    # Step 1: Load PDFs from each path and process them one by one
    pdf_paths = [
        "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/pdfs/sciencedirect/",
        "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/pdfs/pubmed/",
        "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/pdfs/springer/",
    ]

    for pdf_path in pdf_paths:
        pdf_loader = LoadPDF(pdf_path=pdf_path)
        list_pdf_files = (
            pdf_loader.list_pdf_files()
        )  # Get list of PDF files in the directory

        print(f"{len(list_pdf_files)} PDFs found in the directory.")

        # Step 2: Determine which CSV to load based on the PDF path
        if "sciencedirect" in pdf_path:
            csv_file_path = "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/sciencedirect_results.csv"
        elif "pubmed" in pdf_path:
            csv_file_path = "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/pubmed_results.csv"
        else:
            csv_file_path = "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/springer_results.csv"

        # Load the CSV file and initialize FindSimilarity class with the file path
        df_csv = _load_csv(file_name=csv_file_path)
        similarity_finder = FindSimilarity(csv_file=df_csv, csv_file_path=csv_file_path)

        # Step 3: Loop through each PDF, extract emails, and process them
        for pdf_file in list_pdf_files:
            pdf_texts = pdf_loader.extract_text_from_pdf(
                pdf_file
            )  # Extract text from this PDF

            # Step 4: Extract email addresses from the text of this PDF
            email_extractor = ExtractEmails()
            list_of_emails = email_extractor.chatgpt_response(pdf_texts)

            # Step 5: Find similarity between author names and emails, and save CSV
            similarity_finder.find_email_author_and_save(list_of_emails)
            print("-----" * 15)


def fill_empty_emails_with_search():
    csv_paths = [
        "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/sciencedirect_results.csv",
        "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/pubmed_results.csv",
        "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/springer_results.csv",
    ]

    for csv_path in csv_paths:
        _csv = _load_csv(file_name=csv_path)

        headers = [header for header in _csv.columns]
        if "email" not in headers:  # Ensure the column name matches
            _csv["email"] = ""  # Add the 'email' column if missing

        if "sciencedirect" in csv_path:
            csv_file_path = "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/sciencedirect_results.csv"
        elif "pubmed" in csv_path:
            csv_file_path = "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/pubmed_results.csv"
        else:
            csv_file_path = "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/springer_results.csv"

        # Step 5: Search for email addresses using the AI Search
        for _index, row in _csv.iterrows():
            keyword = ""
            if "sciencedirect" in csv_path:
                keyword = row["Keyword_2"]
            elif "pubmed" in csv_path:
                keyword = row["Keyword_1"]
            elif "springer" in csv_path:
                keyword = row["Keyword"]

            author_name = row["Authors"]

            print(f"Processing author: {author_name}")

            web_search = WebSearch(
                name=str(author_name), csv_file=_csv, _keyword=keyword
            )
            list_of_emails = web_search.ai_search()

            # Step 6: Find similarity between the author names and emails
            df_csv = _load_csv(file_name=csv_file_path)
            similarity_finder = FindSimilarity(
                csv_file=df_csv, csv_file_path=csv_file_path
            )
            similarity_finder.find_email_author_and_save(list_of_emails)
            print("-----" * 15)


if __name__ == "__main__":
    # extract_emails_from_pdf()
    fill_empty_emails_with_search()
