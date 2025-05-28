import asyncio
import io
import json
import os
import random
import subprocess
import sys
import zipfile

import streamlit as st

sys.path.append(
    os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))),
)

from utils.app_utils import extract_proxy_format

from arxiv_project.src import arxiv
from sd_pm_ls_scraper.src import pubmed, sciencedirect, springer

os.system("playwright install")


def main():
    st.title("Research Article Scraper")
    st.subheader("Provide search links for any scraper(s) you want to run:")

    arxiv_link = st.text_input("ArXiv search link:", key="arxiv")
    pubmed_link = st.text_input("PubMed search link:", key="pubmed")
    scidir_link = st.text_input("ScienceDirect search link:", key="scidir")
    springer_link = st.text_input("Springer search link:", key="springer")

    st.sidebar.header("Settings")
    uploaded_file = st.sidebar.file_uploader("Upload proxies file (.txt)", type=["txt"])

    # Initialize proxy variables
    random_proxy_for_springer = None
    random_proxy_for_arxiv = None
    random_proxy_for_pubmed = None
    random_proxy_for_scidir = None

    if uploaded_file is not None:
        try:
            # Decode bytes and split lines
            content = uploaded_file.read().decode("utf-8")
            proxy_list = [line.strip() for line in content.splitlines() if line.strip()]
            if not proxy_list:
                st.sidebar.warning("No valid proxies found in the uploaded file.")
                formatted_proxy_list = []
            elif any(":" not in proxy for proxy in proxy_list):
                st.sidebar.warning(
                    "Proxies should be in the format 'ip:port:username:password'.",
                )
                formatted_proxy_list = []
            else:
                formatted_proxy_list = extract_proxy_format(proxy_list)
                st.sidebar.success(f"Loaded {len(proxy_list)} proxies.")
        except Exception as e:
            st.sidebar.error(f"Failed to load proxies: {e}")
            formatted_proxy_list = []
    else:
        st.sidebar.info("No proxies file uploaded. Continuing without proxies.")
        formatted_proxy_list = []

    if formatted_proxy_list:
        random_proxy_for_springer = random.choice(formatted_proxy_list)
        st.sidebar.write("Springer Proxy")
        st.sidebar.write(
            f"Using proxy: {random_proxy_for_springer['server']} with username: {random_proxy_for_springer['username']}, password: {random_proxy_for_springer['password']}",
        )
        st.sidebar.write("---" * 10)

        random_proxy_for_arxiv = random.choice(formatted_proxy_list)
        st.sidebar.write("ArXiv Proxy")
        st.sidebar.write(
            f"Using proxy: {random_proxy_for_arxiv['server']} with username: {random_proxy_for_arxiv['username']}, password: {random_proxy_for_arxiv['password']}",
        )
        st.sidebar.write("---" * 10)

        random_proxy_for_pubmed = random.choice(formatted_proxy_list)
        st.sidebar.write("PubMed Proxy")
        st.sidebar.write(
            f"Using proxy: {random_proxy_for_pubmed['server']} with username: {random_proxy_for_pubmed['username']}, password: {random_proxy_for_pubmed['password']}",
        )
        st.sidebar.write("---" * 10)

        random_proxy_for_scidir = random.choice(formatted_proxy_list)
        st.sidebar.write("ScienceDirect Proxy")
        st.sidebar.write(
            f"Using proxy: {random_proxy_for_scidir['server']} with username: {random_proxy_for_scidir['username']}, password: {random_proxy_for_scidir['password']}",
        )
        st.sidebar.write("---" * 10)
    else:
        st.sidebar.write("No valid proxies available. Running without proxies.")

    if st.button("Scrape All"):
        # ArXiv
        st.header("ArXiv Results")
        if arxiv_link:
            with st.spinner("Scraping ArXiv..."):
                if (
                    hasattr(arxiv, "async_main_run")
                    and "proxies" in arxiv.async_main_run.__code__.co_varnames
                ):
                    asyncio.run(
                        arxiv.async_main_run(
                            arxiv_link,
                            proxies=random_proxy_for_arxiv,
                        ),
                    )
                else:
                    asyncio.run(arxiv.async_main_run(arxiv_link))
            st.success("ArXiv scraping completed.")
        else:
            st.info("No ArXiv link provided.")

        # PubMed
        st.header("PubMed Results")
        if pubmed_link:
            with st.spinner("Scraping PubMed..."):
                # Pass proxy if available, otherwise None
                if (
                    hasattr(pubmed, "main_run")
                    and "proxies" in pubmed.main_run.__code__.co_varnames
                ):
                    pubmed.main_run(pubmed_link, proxies=random_proxy_for_pubmed)
                else:
                    pubmed.main_run(pubmed_link)
            st.success("PubMed scraping completed.")
        else:
            st.info("No PubMed link provided.")

        # ScienceDirect
        st.header("ScienceDirect Results")
        if scidir_link:
            with st.spinner("Scraping ScienceDirect..."):
                # Pass proxy if available, otherwise None
                if (
                    hasattr(sciencedirect, "scrape_from_link")
                    and "proxies" in sciencedirect.scrape_from_link.__code__.co_varnames
                ):
                    sciencedirect.scrape_from_link(
                        scidir_link,
                        proxies=random_proxy_for_scidir,
                    )
                else:
                    sciencedirect.scrape_from_link(scidir_link)
            st.success("ScienceDirect scraping completed.")
        else:
            st.info("No ScienceDirect link provided.")

        # Springer
        st.header("Springer Results")
        if springer_link:
            with st.spinner("Scraping Springer..."):
                script_path = os.path.abspath(__file__)
                cmd = [
                    "xvfb-run",
                    "-a",
                    sys.executable,
                    script_path,
                    "--springer",
                    springer_link,
                ]

                # Add proxy argument if available
                if random_proxy_for_springer:
                    # You might need to serialize the proxy dict to pass as argument
                    proxy_json = json.dumps(random_proxy_for_springer)
                    cmd.extend(["--proxy-json", proxy_json])

                subprocess.run(cmd, check=True)
            st.success("Springer scraping completed.")
        else:
            st.info("No Springer link provided.")
    try:
        # concataned_file_path = ""
        # with open(concataned_file_path, "rb") as csv_file:
        #     st.download_button(
        #         label="Download Concatened Final Results (XLSX)",
        #         data=csv_file,
        #         file_name="processed_file.csv",
        #         mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
        #     )

        # Download 4 unique files in one button as a zip
        # Define your 4 unique file paths
        file_paths = [
            "arxiv_project/output/arxiv_results.xlsx",
            "sd_pm_ls_scraper/output/pubmed_results.csv",
            "sd_pm_ls_scraper/output/sciencedirect_results.csv",
            "sd_pm_ls_scraper/output/springer_results.csv",
        ]

        # Create an in-memory zip file
        zip_buffer = io.BytesIO()
        with zipfile.ZipFile(zip_buffer, "w") as zip_file:
            for file_path in file_paths:
                try:
                    with open(file_path, "rb") as f:
                        zip_file.writestr(os.path.basename(file_path), f.read())
                except FileNotFoundError as err:
                    print("File not found:", err)
                    continue
        zip_buffer.seek(0)
        st.download_button(
            label="Download All Results (ZIP)",
            data=zip_buffer,
            file_name="all_results.zip",
            mime="application/zip",
        )
    except FileNotFoundError:
        st.error(
            "Processed file not found. Please ensure the process completed successfully.",
        )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--springer", type=str, help="Springer link")
    parser.add_argument(
        "--proxies",
        type=str,
        help="Path to proxies file (.txt)",
        default=None,
    )
    parser.add_argument(
        "--proxy-json",
        type=str,
        help="Proxy configuration as JSON string",
        default=None,
    )
    args = parser.parse_args()

    if args.springer:
        # Load proxies if provided
        proxy_list = []
        random_proxy_for_springer = None

        if args.proxy_json:
            # Use proxy from JSON argument (from Streamlit)
            random_proxy_for_springer = json.loads(args.proxy_json)
        elif args.proxies and os.path.exists(args.proxies):
            # Load from file (direct command line usage)
            with open(args.proxies, "r", encoding="utf-8") as f:
                proxy_list = [line.strip() for line in f if line.strip()]
            formatted_proxy_list = extract_proxy_format(proxy_list)
            random_proxy_for_springer = (
                random.choice(formatted_proxy_list) if formatted_proxy_list else None
            )

        asyncio.run(
            springer.async_springer(args.springer, proxies=random_proxy_for_springer),
        )
    else:
        main()
