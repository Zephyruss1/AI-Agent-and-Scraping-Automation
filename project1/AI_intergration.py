import pdfplumber
from openai import OpenAI
from dotenv import load_dotenv
import os

load_dotenv()

def extract_text_from_pdf(pdf_path):
    text = ""
    with pdfplumber.open(pdf_path) as pdf:
        for page in pdf.pages:
            text += page.extract_text() + "\n"
    return text

pdf_text = extract_text_from_pdf("pdfs/1908.09519.pdf")

client = OpenAI()
client.api_key = os.getenv("OPENAI_API_KEY")
completion = client.chat.completions.create(
    model="gpt-4o",
    store=True,
    messages = [
    {
        "role": "system",
        "content": (
            "You are a data extraction assistant. The user will provide the text from a PDF. "
            "Identify and list all email addresses that appear in the text (i.e., strings containing '@'). "
            "If you find no email addresses, return 'None'. "
            "Output only the emails or 'None'—no additional explanations."
        )
    },
    {
        "role": "user",
        "content": pdf_text
    }
    ],
    temperature=0.2,
)

# Extract and print the message content
message_content = completion.choices[0].message.content
print(message_content)

from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity

doc1 = 'shifan.xu@yale.edu'
doc2 = "Shifan Xu"

vectorizer = TfidfVectorizer(analyzer='char', ngram_range=(1, 2), lowercase=True)
tfidf_matrix = vectorizer.fit_transform([doc1.split("@")[0], doc2])

# Calculate cosine similarity
cosine_sim = cosine_similarity(tfidf_matrix[0:1], tfidf_matrix[1:2])
print(f"Cosine Similarity (Character n-grams): {cosine_sim[0][0]}")

if cosine_sim > 0.6:
    print("Match")


import requests
from pydantic import BaseModel
from dotenv import load_dotenv
import os

load_dotenv()

class AnswerFormat(BaseModel):
    email_adress: str

url = "https://api.perplexity.ai/chat/completions"
headers = {"Authorization": f"Bearer {os.getenv('PERPLEXITY_API_KEY')}"}
payload = {
    "model": "sonar-pro",
    "messages": [
        {"role": "system",
         "content":(
            "You are a web searcher assistant. The user will provide an author name."
            "Your task is: Search email addresses for provided author name."
            "If you find no email addresses, return 'None'. "
            "Output only the emails or 'None'—no additional explanations."
         ),
         },
        {"role": "user", "content": "Kohr Holger"},
    ],
    "response_format": {
        "type": "json_schema",
        "json_schema": {"schema": AnswerFormat.model_json_schema()},
    },
}

response = requests.post(url, headers=headers, json=payload)

if response.status_code == 200:
    try:
        response_json = response.json()
        print(response_json["choices"][0]["message"]["content"])
    except requests.JSONDecodeError:
        print("Error: Response content is not valid JSON")
else:
    print(f"Error: Received status code {response.status_code}")