import os
import sys

# Adjust path to point to: sd_pm_ls_scraper/src/
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "../src")))

from email_extraction_and_search import ExtractEmails


def test_extract_with_regex_on_sample_text():
    extractor = ExtractEmails()
    text = "Contact us at test@example.com and info@company.org."
    emails = extractor.extract_with_regex(text)

    assert len(emails) == 2
    expected = ["test@example.com", "info@company.org"]
    assert emails == expected, f"Expected {expected}, but got {emails}"


def test_extract_with_regex_on_empty_text():
    extractor = ExtractEmails()
    text = ""
    emails = extractor.extract_with_regex(text)

    assert len(emails) == 1  # "None" is 1
    expected = ["None"]
    assert emails == expected, f"Expected {expected}, but got {emails}"


def test_extract_with_regex_on_no_email_text():
    extractor = ExtractEmails()
    text = "This is a test string without any email addresses."
    emails = extractor.extract_with_regex(text)

    assert len(emails) == 1  # "None" is 1
    expected = ["None"]
    assert emails == expected, f"Expected {expected}, but got {emails}"


def test_extract_with_regex_on_return_dtype():
    extractor = ExtractEmails()
    text = "Contact us at test@example.com and info@company.org."
    emails = extractor.extract_with_regex(text)
    assert len(emails) == 2
    assert isinstance(emails, list), f"Expected a list, but got {type(emails)}"


def test_extract_with_regex_on_invalid_email():
    extractor = ExtractEmails()
    text = "Contact us at test@.com and info@company..org."
    emails = extractor.extract_with_regex(text)

    assert len(emails) == 1  # "None" is 1
    expected = ["info@company.org"]
    assert emails == expected, f"Expected {expected}, but got {emails}"


def test_extract_with_regex_on_complex_email():
    extractor = ExtractEmails()
    text = "Contantact us at test@@example.com and info@company..org"
    emails = extractor.extract_with_regex(text)
    assert len(emails) == 2
    expected = ["test@example.com", "info@company.org"]
    assert emails == expected, f"Expected {expected}, but got {emails}"
