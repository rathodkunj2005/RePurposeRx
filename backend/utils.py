import logging
import json
import requests
import aiohttp
import asyncio
import xml.etree.ElementTree as ET

# Initialize logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


# Synchronous API Fetcher
def fetch_data_from_api(url):
    """
    Generic function to fetch JSON data from an API URL with improved error handling.
    """
    try:
        response = requests.get(url)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.HTTPError as http_err:
        logging.error(f"HTTP error fetching data from {url}: {http_err}")
    except requests.exceptions.RequestException as req_err:
        logging.error(f"Request error fetching data from {url}: {req_err}")
    except json.JSONDecodeError as json_err:
        logging.error(f"JSON decode error from {url}: {json_err}")
    return None

# Asynchronous API Fetcher
async def fetch_data_from_api_async(session, url):
    """
    Asynchronous function to fetch JSON data from an API URL.
    """
    try:
        async with session.get(url) as response:
            response.raise_for_status()
            return await response.json()
    except aiohttp.ClientError as err:
        logging.error(f"Client error fetching data from {url}: {err}")
    return None

async def fetch_all_data_async(urls):
    """
    Fetch data from multiple URLs asynchronously.
    """
    async with aiohttp.ClientSession() as session:
        tasks = [fetch_data_from_api_async(session, url) for url in urls]
        return await asyncio.gather(*tasks)

# Extract PDF URL from PubMed Article
def get_pdf_url(article):
    """
    Extracts the free full-text PDF link from a PubMed article summary if a PubMed Central (PMC) ID is available.
    """
    articleids = article.get("articleids", [])
    for id_entry in articleids:
        if id_entry.get("idtype", "").lower() == "pmc":
            pmc_id = id_entry.get("value")
            if not pmc_id.startswith("PMC"):
                pmc_id = f"PMC{pmc_id}"
            pdf_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmc_id}/pdf/"
            return pdf_url
    return None

# Download PDF
def download_pdf(pdf_url, save_path):
    """
    Downloads a PDF from the given URL and saves it to the specified path.
    """
    try:
        response = requests.get(pdf_url)
        response.raise_for_status()
        with open(save_path, "wb") as f:
            f.write(response.content)
        logging.info(f"Downloaded PDF from {pdf_url} to {save_path}")
        return save_path
    except Exception as e:
        logging.error(f"Error downloading PDF from {pdf_url}: {e}")
        return None
