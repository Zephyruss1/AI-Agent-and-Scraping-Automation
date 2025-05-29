import asyncio
import io
import json
import os
import random
import shutil
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

# Define consistent output directory
OUTPUT_DIR = "/tmp/output"
os.makedirs(OUTPUT_DIR, exist_ok=True)


def check_file_exists(file_path):
    """Check if file exists and has content"""
    if os.path.exists(file_path) and os.path.getsize(file_path) > 0:
        return True
    return False


def get_available_files():
    """Get list of available result files"""
    file_mappings = {
        "arxiv_results.xlsx": os.path.join(OUTPUT_DIR, "arxiv_results.xlsx"),
        "pubmed_results.csv": os.path.join(OUTPUT_DIR, "pubmed_results.csv"),
        "sciencedirect_results.csv": os.path.join(
            OUTPUT_DIR,
            "sciencedirect_results.csv",
        ),
        "springer_results.csv": os.path.join(OUTPUT_DIR, "springer_results.csv"),
    }

    # Also check current directory for pubmed (since your script saves there)
    pubmed_current = "pubmed_results.csv"
    if check_file_exists(pubmed_current):
        file_mappings["pubmed_results.csv"] = pubmed_current

    available_files = {}
    for name, path in file_mappings.items():
        if check_file_exists(path):
            available_files[name] = path

    return available_files


def create_download_zip(available_files):
    """Create a zip file with all available results"""
    if not available_files:
        return None

    zip_buffer = io.BytesIO()
    with zipfile.ZipFile(zip_buffer, "w", zipfile.ZIP_DEFLATED) as zip_file:
        for filename, filepath in available_files.items():
            try:
                zip_file.write(filepath, filename)
                st.write(f"✅ Added {filename} to zip")
            except Exception as e:
                st.write(f"⚠️ Could not add {filename}: {str(e)}")

    zip_buffer.seek(0)
    return zip_buffer


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
        random_proxy_for_arxiv = random.choice(formatted_proxy_list)
        random_proxy_for_pubmed = random.choice(formatted_proxy_list)
        random_proxy_for_scidir = random.choice(formatted_proxy_list)

        st.sidebar.write("Using proxies for all scrapers")
    else:
        st.sidebar.write("No valid proxies available. Running without proxies.")

    # Scraping section
    if st.button("Scrape All"):
        scraping_results = {}

        # ArXiv
        st.header("ArXiv Results")
        if arxiv_link:
            with st.spinner("Scraping ArXiv..."):
                try:
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
                    scraping_results["arxiv"] = True
                except Exception as e:
                    st.error(f"ArXiv scraping failed: {str(e)}")
                    scraping_results["arxiv"] = False
        else:
            st.info("No ArXiv link provided.")

        # PubMed
        st.header("PubMed Results")
        if pubmed_link:
            with st.spinner("Scraping PubMed..."):
                try:
                    if (
                        hasattr(pubmed, "main_run")
                        and "proxies" in pubmed.main_run.__code__.co_varnames
                    ):
                        result_file = pubmed.main_run(
                            pubmed_link,
                            proxies=random_proxy_for_pubmed,
                        )
                    else:
                        result_file = pubmed.main_run(pubmed_link)

                    if result_file:
                        # Move file to output directory if it exists in current directory
                        current_file = "pubmed_results.csv"
                        target_file = os.path.join(OUTPUT_DIR, "pubmed_results.csv")
                        if os.path.exists(current_file) and current_file != target_file:
                            shutil.move(current_file, target_file)

                        st.success("PubMed scraping completed successfully.")
                        scraping_results["pubmed"] = True
                    else:
                        st.warning(
                            "PubMed scraping completed but no results were found.",
                        )
                        scraping_results["pubmed"] = False
                except Exception as e:
                    st.error(f"PubMed scraping failed: {str(e)}")
                    scraping_results["pubmed"] = False
        else:
            st.info("No PubMed link provided.")

        # ScienceDirect
        st.header("ScienceDirect Results")
        if scidir_link:
            with st.spinner("Scraping ScienceDirect..."):
                try:
                    if (
                        hasattr(sciencedirect, "scrape_from_link")
                        and "proxies"
                        in sciencedirect.scrape_from_link.__code__.co_varnames
                    ):
                        sciencedirect.scrape_from_link(
                            scidir_link,
                            proxies=random_proxy_for_scidir,
                        )
                    else:
                        sciencedirect.scrape_from_link(scidir_link)
                    st.success("ScienceDirect scraping completed.")
                    scraping_results["sciencedirect"] = True
                except Exception as e:
                    st.error(f"ScienceDirect scraping failed: {str(e)}")
                    scraping_results["sciencedirect"] = False
        else:
            st.info("No ScienceDirect link provided.")

        # Springer
        st.header("Springer Results")
        if springer_link:
            with st.spinner("Scraping Springer..."):
                try:
                    script_path = os.path.abspath(__file__)
                    use_xvfb = shutil.which("xvfb-run") is not None
                    print(f"use_xvfb: {use_xvfb}")
                    if use_xvfb:
                        cmd = [
                            "xvfb-run",
                            "-a",
                            sys.executable,
                            script_path,
                            "--springer",
                            springer_link,
                        ]
                    else:
                        cmd = [
                            sys.executable,
                            script_path,
                            "--springer",
                            springer_link,
                        ]

                    if random_proxy_for_springer:
                        proxy_json = json.dumps(random_proxy_for_springer)
                        cmd.extend(["--proxy-json", proxy_json])

                    subprocess.run(cmd, check=True)
                    st.success("Springer scraping completed.")
                    scraping_results["springer"] = True
                except Exception as e:
                    st.error(f"Springer scraping failed: {str(e)}")
                    scraping_results["springer"] = False
        else:
            st.info("No Springer link provided.")

        # Show summary
        st.header("Scraping Summary")
        for scraper, success in scraping_results.items():
            if success:
                status = "✅ Success - Results found"
            elif success is False:
                status = "⚠️ Completed but no results found"
            else:
                status = "❌ Failed"
            st.write(f"{scraper.title()}: {status}")

    # Download section - always show available files
    st.header("Download Results")

    # Check for available files
    available_files = get_available_files()

    if available_files:
        st.write("Available files:")
        for filename, filepath in available_files.items():
            file_size = os.path.getsize(filepath)
            st.write(f"📄 {filename} ({file_size:,} bytes)")

        # Individual file downloads
        st.subheader("Download Individual Files")
        for filename, filepath in available_files.items():
            try:
                with open(filepath, "rb") as f:
                    file_data = f.read()
                    st.download_button(
                        label=f"Download {filename}",
                        data=file_data,
                        file_name=filename,
                        mime="application/octet-stream",
                        key=f"download_{filename}",
                    )
            except Exception as e:
                st.error(f"Could not prepare {filename} for download: {str(e)}")

        # Zip download
        st.subheader("Download All Files")
        zip_buffer = create_download_zip(available_files)
        if zip_buffer:
            st.download_button(
                label=f"Download All Results (ZIP) - {len(available_files)} files",
                data=zip_buffer,
                file_name="all_research_results.zip",
                mime="application/zip",
                key="download_all_zip",
            )
    else:
        st.info(
            "No result files found. Run the scrapers first to generate downloadable files.",
        )

        # Debug info
        st.subheader("Debug Information")
        st.write(f"Looking for files in: {OUTPUT_DIR}")
        if os.path.exists(OUTPUT_DIR):
            files_in_output = os.listdir(OUTPUT_DIR)
            st.write(f"Files in output directory: {files_in_output}")
        else:
            st.write("Output directory does not exist")

        # Check current directory too
        current_files = [f for f in os.listdir(".") if f.endswith((".csv", ".xlsx"))]
        if current_files:
            st.write(f"CSV/Excel files in current directory: {current_files}")


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
        proxy_list = []
        random_proxy_for_springer = None

        if args.proxy_json:
            random_proxy_for_springer = json.loads(args.proxy_json)
        elif args.proxies and os.path.exists(args.proxies):
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
