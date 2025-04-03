import pandas as pd
import streamlit as st
from utils.websearch_utils import (
    chatgpt_fill_empty_emails_with_search,
    chatgpt_fill_empty_jobs_with_search,
    convert_excel_to_csv,
    load_spreadsheet,
    perplexity_fill_empty_emails_with_search,
    perplexity_fill_empty_jobs_with_search,
)


def main():
    st.title("Web Search Tool")
    st.sidebar.header("Configuration")

    spreadsheet_link = st.sidebar.text_input("Enter the link to your spreadsheet:")

    if st.sidebar.button("Load Spreadsheet"):
        if spreadsheet_link:
            try:
                _, condition = load_spreadsheet(spreadsheet_link)
                st.success("Spreadsheet loaded successfully!")
            except Exception as e:
                st.error(f"Error loading spreadsheet: {e}")
            try:
                convert_excel_to_csv()
                st.success("Spreadsheet converted to CSV successfully!")
            except Exception as e:
                st.error(f"Error converting spreadsheet to CSV: {e}")
        else:
            st.warning("Please enter a valid spreadsheet link.")

    model_choice = st.sidebar.selectbox(
        "Select AI Model", ["ChatGPT [Browser-Use]", "Perplexity"]
    )

    if model_choice == "ChatGPT [Browser-Use]":
        gpt_model = st.sidebar.selectbox("Select GPT Model", ["gpt-4o", "gpt-4o-mini"])
    else:
        perplexity_model = st.sidebar.selectbox(
            "Select Perplexity Model", ["sonar-pro", "sonar-deep-research"]
        )

    purpose_choice = st.sidebar.selectbox(
        "Select Purpose", ["Find email addresses", "Find job titles"]
    )

    if model_choice == "Perplexity":
        st.sidebar.subheader("Custom System Prompt (Optional)")
        custom_system_prompt = st.sidebar.text_area(
            "Custom System Prompt for Perplexity:",
            placeholder="You are a web search assistant. The user will provide an author name. "
            "Search the {self.keyword_index} field for that author's job title. "
            "If no job title is found, return 'None'. "
            "Output only the job title or 'None' with no additional commentary.",
        )
        if custom_system_prompt:
            st.sidebar.success("Custom System Prompt added successfully!")
    else:
        custom_system_prompt = None

    if model_choice == "ChatGPT [Browser-Use]":
        st.sidebar.subheader("Custom Prompt (Optional)")
        custom_prompt = st.sidebar.text_area(
            "Custom Prompt for ChatGPT:",
            placeholder="1. Search '{self.author_name} email address in the {self.keyword_index} field' and enter."
            "2. Output only the emails or 'None'—no additional explanations.",
        )
        if custom_prompt:
            st.sidebar.success("Custom prompt added successfully!")
    elif model_choice == "Perplexity":
        st.sidebar.subheader("Custom Prompt (Optional)")
        custom_prompt = st.sidebar.text_area(
            "Custom Prompt for ChatGPT:", placeholder="{self.author_name} email adress"
        )
        if custom_prompt:
            st.sidebar.success("Custom prompt added successfully!")
    else:
        custom_prompt = None
    st.sidebar.markdown(
        '<p style="color: gray; font-size: 12px;">This prompt will guide the AI in its search process.</p>'
        '<p style="color: gray; font-size: 12px;">Please add always {self.author_name} and {self.keyword_index} to mention</p>',
        unsafe_allow_html=True,
    )

    if st.sidebar.button("Start Search"):
        output_csv_path = "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/downloaded_sheet.csv"
        output_xlsx_path = "/root/arxiv-and-scholar-scraping/sd_pm_ls_scraper/output/downloaded_sheet.xlsx"

        if model_choice == "ChatGPT [Browser-Use]":
            if purpose_choice == "Find email addresses":
                chatgpt_fill_empty_emails_with_search(gpt_model, prompt=custom_prompt)
                st.success("Searching for email addresses...")
            elif purpose_choice == "Find job titles":
                chatgpt_fill_empty_jobs_with_search(gpt_model, prompt=custom_prompt)
                st.success("Searching for job titles...")
        elif model_choice == "Perplexity":
            if purpose_choice == "Find email addresses":
                perplexity_fill_empty_emails_with_search(
                    perplexity_model,
                    system_prompt=custom_system_prompt,
                    prompt=custom_prompt,
                )
                st.success("Searching for email addresses...")
            elif purpose_choice == "Find job titles":
                perplexity_fill_empty_jobs_with_search(
                    perplexity_model,
                    system_prompt=custom_system_prompt,
                    prompt=custom_prompt,
                )
                st.success("Searching for job titles...")
        else:
            st.warning("Please select a valid purpose.")

        # Convert CSV to XLSX
        try:
            df = pd.read_csv(output_csv_path)
            with pd.ExcelWriter(output_xlsx_path, engine="openpyxl") as writer:
                df.to_excel(writer, index=False, sheet_name="Sheet1")
            st.success("File converted to XLSX format successfully!")
        except Exception as e:
            st.error(f"Error converting file to XLSX: {e}")

        # Add download buttons for both CSV and XLSX files
        try:
            with open(output_csv_path, "rb") as csv_file:
                st.download_button(
                    label="Download Processed File (CSV)",
                    data=csv_file,
                    file_name="processed_file.csv",
                    mime="text/csv",
                )
            with open(output_xlsx_path, "rb") as xlsx_file:
                st.download_button(
                    label="Download Processed File (XLSX)",
                    data=xlsx_file,
                    file_name="processed_file.xlsx",
                    mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                )
        except FileNotFoundError:
            st.error(
                "Processed file not found. Please ensure the process completed successfully."
            )


if __name__ == "__main__":
    main()
