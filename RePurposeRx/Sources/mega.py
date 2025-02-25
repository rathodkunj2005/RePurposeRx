import json
import logging
import requests
from ollama import chat
import xml.etree.ElementTree as ET
import os
import concurrent.futures
from config import (
    METADATA_EXTRACTION_MODEL,
    RESPONSE_GENERATION_MODEL,
    PUBCHEM_BASE_URL,
    CHEMBL_BASE_URL,
    UNIPROT_BASE_URL,
    ENSEMBL_BASE_URL,
    RXNAV_BASE_URL,
    KEGG_BASE_URL,
    REACTOME_BASE_URL,
    HPA_BASE_URL,
    PUBMED_EUTILS_BASE_URL,
    OMIM_API_BASE_URL,
    DISGENET_API_BASE_URL,
    CLINICALTRIALS_GOV_BASE_URL,
    EUROPE_PMC_BASE_URL,
    TCGA_GDC_API_BASE_URL,
    FDA_API_BASE_URL,
    DRUG_NAME_BLACKLIST,
    PUBCHEM_PROPERTIES,
    OMIM_API_KEY,
    DISGENET_API_KEY
)


# Initialize logging for better error tracking
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# ------------------------------
# Configuration & Constants
# ------------------------------
PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
CHEMBL_BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/molecule/search"
UNIPROT_BASE_URL = "https://rest.uniprot.org/uniprotkb"
ENSEMBL_BASE_URL = "https://rest.ensembl.org"
RXNAV_BASE_URL = "https://mor.nlm.nih.gov/RxNav"
KEGG_BASE_URL = "https://rest.kegg.jp"
REACTOME_BASE_URL = "https://reactome.org/ContentService"
HPA_BASE_URL = "https://www.proteinatlas.org/api"
PUBMED_EUTILS_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
OMIM_API_BASE_URL = "https://api.omim.org/api"
DISGENET_API_BASE_URL = "https://www.disgenet.org/api"
CLINICALTRIALS_GOV_BASE_URL = "https://clinicaltrials.gov/api/query"
EUROPE_PMC_BASE_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
TCGA_GDC_API_BASE_URL = "https://api.gdc.cancer.gov"
FDA_API_BASE_URL = "https://api.fda.gov/drug"

DRUG_NAME_BLACKLIST = ["primary targets", "targets", "target", "pathway", "information"]
PUBCHEM_PROPERTIES = "MolecularWeight,IUPACName,CanonicalSMILES"
OMIM_API_KEY = "3FZX0uYkSVmNO4Ako4d8yg"  # Existing OMIM API Key
DISGENET_API_KEY = "87dec4c5-14f3-43ec-8432-9b16bc1a0b26"  # New DisGeNET API Key

# ------------------------------
# Helper Function: API Fetcher (Improved Error Handling & Logging)
# ------------------------------
def fetch_data_from_api(url):
    """
    Generic function to fetch JSON data from an API URL with improved error handling.
    """
    try:
        response = requests.get(url)
        response.raise_for_status()  # Raise HTTPError for bad responses (4xx or 5xx)
        return response.json()
    except requests.exceptions.HTTPError as http_err:
        logging.error(f"HTTP error fetching data from {url}: {http_err}")
    except requests.exceptions.RequestException as req_err:
        logging.error(f"Request error fetching data from {url}: {req_err}")
    except json.JSONDecodeError as json_err:
        logging.error(f"JSON decode error from {url}: {json_err}")
    return None

# ------------------------------
# Metadata Extraction Function (Modularized)
# ------------------------------
def extract_metadata(user_query, model_name=METADATA_EXTRACTION_MODEL):
    """
    Extracts relevant metadata from the user query using a specified LLM model.
    """
    metadata_prompt = f"""
    You are an expert at extracting relevant entity names, data, and related keywords from a query.
    Please analyze the following user query and extract the information into the following JSON structure.
    For each category, output a list of names or pieces of data that are relevant. If no relevant data exists for a category, use an empty list.
    Additionally, identify related keywords that might help answer the user query and include them in the JSON.
    Output only the JSON object with no additional commentary.

    {{
      "Drug Information": [ "names or data" ],
      "Gene Information": [ "names or data" ],
      "Target and Pathway Information": [ "names or data" ],
      "Disease Information": [ "names or data" ],
      "Clinical Trial Data": [ "names or data" ],
      "Literature and Evidence": [ "names or data" ],
      "Omics Data": [ "names or data" ],
      "Other Important Resources": {{
          "Patents": [ "names or data" ],
          "FDA": [ "names or data" ],
          "European Medicines Agency": [ "names or data" ]
      }},
      "General Query": true/false, // Indicates if the query is general and can be answered using the METADATA model without scraping or deepseek
      "Related Keywords": [ "keywords that might help answer the user query" ]
    }}

    User query: {user_query}
    """

    try:
        metadata_response = chat(
            model=model_name,
            messages=[{'role': 'user', 'content': metadata_prompt}],
            stream=False,
        )
        metadata = json.loads(metadata_response['message']['content'])
        logging.info("Metadata extracted successfully.")
        return metadata
    except json.JSONDecodeError as json_err:
        logging.error(f"Error parsing metadata JSON: {json_err}")
        return {}
    except Exception as e:
        logging.error(f"Error during metadata extraction: {e}")
        return {}


# ------------------------------
# Data Validation Function
# ------------------------------
def is_valid_drug_name(name, blacklist=DRUG_NAME_BLACKLIST):
    """
    Heuristic to check if a name is likely a drug name.
    """
    for word in blacklist:
        if word in name.lower():
            return False
    return True

# ------------------------------
# PubChem Data Retrieval Functions (Modularized)
# ------------------------------
def fetch_pubchem_info(drug_name, base_url=PUBCHEM_BASE_URL, properties=PUBCHEM_PROPERTIES, output_format="JSON"):
    """
    Query PubChem for drug information.

    Parameters:
      drug_name (str): The name of the drug.
      base_url (str): The base URL for the PubChem API.
      properties (str): A comma-separated list of properties to extract.
      output_format (str): The desired output format ("JSON" or "TXT"). Default is "JSON".

    Returns:
      dict or str: Returns a parsed JSON dictionary if output_format is "JSON",
                   or a plain text string if output_format is "TXT". In case of failures,
                   returns None.
    """
    output_format_upper = output_format.upper()
    url = f"{base_url}/compound/name/{drug_name}/property/{properties}/{output_format_upper}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        if output_format_upper == "JSON":
            return response.json()
        elif output_format_upper == "TXT":
            return response.text
        else:
            logging.error(f"Unsupported output format: {output_format}")
            return None
    except requests.exceptions.HTTPError as http_err:
        logging.error(f"HTTP error fetching data from {url}: {http_err}")
    except requests.exceptions.RequestException as req_err:
        logging.error(f"Request error fetching data from {url}: {req_err}")
    except json.JSONDecodeError as json_err:
        logging.error(f"JSON decode error from {url}: {json_err}")
    return None

# ------------------------------
# ChEMBL Data Retrieval Functions (Modularized)
# ------------------------------
def fetch_chembl_info(drug_name, base_url=CHEMBL_BASE_URL):
    """
    Query ChEMBL for drug information and return text data.
    """
    params = {'q': drug_name, 'format': 'txt'}
    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        return response.text
    except requests.exceptions.HTTPError as http_err:
        logging.error(f"HTTP error fetching ChEMBL data for '{drug_name}': {http_err}")
    except requests.exceptions.RequestException as req_err:
        logging.error(f"Request error fetching ChEMBL data for '{drug_name}': {req_err}")
    return None

# ------------------------------
# UniProt Data Retrieval Functions (Modularized)
# ------------------------------
def fetch_uniprot_id(gene_name, base_url=UNIPROT_BASE_URL):
    """
    Fetch UniProt ID for a gene name.
    """
    url = f"{base_url}/search?query={gene_name}&fields=accession&format=json"
    data = fetch_data_from_api(url)
    if data and 'results' in data and data['results']:
        return data['results'][0]['primaryAccession']
    logging.warning(f"UniProt ID not found for {gene_name}") # Changed to warning
    return None

def fetch_and_extract_uniprot_function(uniprot_id, base_url=UNIPROT_BASE_URL):
    """
    Fetch UniProt data and extract the 'Function' section.
    """
    url = f"{base_url}/{uniprot_id}.txt"
    try:
        response = requests.get(url)
        response.raise_for_status()
        function_section = []
        capture = False
        for line in response.text.split('\n'):
            if line.startswith("CC   -!- FUNCTION:"):
                capture = True
                function_section.append(line.replace("CC   -!- FUNCTION:", "").strip())
            elif capture and line.startswith("CC       "):
                function_section.append(line.replace("CC       ", "").strip())
            elif capture:
                break
        if function_section:
            return ' '.join(function_section)
        logging.warning("'Function' section not found in UniProt data for {uniprot_id}") # Changed to warning
    except requests.exceptions.HTTPError as http_err:
        logging.error(f"HTTP error fetching UniProt data for {uniprot_id}: {http_err}")
    except requests.exceptions.RequestException as req_err:
        logging.error(f"Request error fetching UniProt data for {uniprot_id}: {req_err}")
    return None

def get_uniprot_protein_info(gene_name="FGFR2", base_url=UNIPROT_BASE_URL):
    """Fetch protein information from UniProt using gene name."""
    search_url = f"{base_url}/search"
    params = {
        'query': f'gene_name:{gene_name}',
        'format': 'json'
    }
    try:
        response = requests.get(search_url, params=params)
        response.raise_for_status()
        data = response.json()
        if data and 'results' in data and data['results']:
            protein_info = data['results'][0] # Take the first result
            return protein_info
        else:
            return "No protein information found for this gene name."
    except requests.exceptions.RequestException as e:
        return f"Error fetching UniProt data: {e}"

# ------------------------------
# Ensembl Data Retrieval Functions (Modularized)
# ------------------------------
def fetch_ensembl_gene_id(gene_name, base_url=ENSEMBL_BASE_URL):
    """
    Fetch Ensembl gene ID for a gene name.
    """
    url = f"{base_url}/xrefs/symbol/homo_sapiens/{gene_name}?content-type=application/json"
    data = fetch_data_from_api(url)
    if data and len(data) > 0:
        return data[0]['id']
    logging.warning(f"Ensembl gene ID not found for {gene_name}") # Changed to warning
    return None

def fetch_ensembl_gene_info(gene_id, base_url=ENSEMBL_BASE_URL):
    """
    Fetch detailed gene information from Ensembl using gene ID.
    """
    url = f"{base_url}/lookup/id/{gene_id}?expand=1&content-type=application/json"
    data = fetch_data_from_api(url)
    if data:
        return {
            'Gene name': data.get('display_name'),
            'Description': data.get('description'),
            'Organism': data.get('species'),
            'Biotype': data.get('biotype'),
            'Location': f"{data.get('seq_region_name')}:{data.get('start')}-{data.get('end')}",
            'Strand': data.get('strand'),
            'Ensembl Link': f"https://www.ensembl.org/Homo_sapiens/Gene/Summary?g={gene_id}"
        }
    logging.warning(f"Ensembl gene details not found for {gene_id}") # Changed to warning
    return {}

# ------------------------------
# RxNav Data Retrieval Functions (Modularized) - Function provided but not used in core logic currently
# ------------------------------
def fetch_drug_info(drug_name, base_url=RXNAV_BASE_URL):
    """
    Fetch drug information from NIH RxNav API.
    """
    url = f"{base_url}/search?searchBy=String&searchTerm={drug_name}"
    response = requests.get(url)
    if response.status_code == 200:
        return url # Consider returning more structured data if API allows
    logging.warning(f"Failed to fetch drug information for {drug_name}, status code {response.status_code}") # Changed to warning
    return None

# ------------------------------
# Gene Ontology (GO) Data Retrieval
# ------------------------------
def get_gene_ontology_info(gene_symbol="FGFR2"):
    """Fetch Gene Ontology terms associated with a gene symbol."""
    base_url = "https://www.ebi.ac.uk/QuickGO/services/annotation/search" # EBI QuickGO API
    params = {
        'geneSymbol': gene_symbol,
        'format': 'json'
    }
    try:
        response = requests.get(base_url, params=params)
        response.raise_for_status()
        data = response.json()
        if data and 'results' in data:
            return data['results']
        else:
            return "No GO terms found for this gene symbol."
    except requests.exceptions.RequestException as e:
        return f"Error fetching Gene Ontology data: {e}"

# ------------------------------
# KEGG Pathway Data Retrieval
# ------------------------------
def get_kegg_pathways_for_gene(gene_name="FGFR2", organism="hsa", base_url=KEGG_BASE_URL): # hsa for Homo sapiens
    """Fetch KEGG pathways associated with a gene."""
    link_url = f"{base_url}/link/pathway"
    url = f"{link_url}/{organism}:{gene_name}"

    try:
        response = requests.get(url)
        response.raise_for_status()
        pathway_ids = response.text.strip().split('\n')
        pathway_list = []
        if pathway_ids and pathway_ids != ['']: # Check for empty response
            for pathway_id in pathway_ids:
                pathway_kegg_id = pathway_id.split(":")[1] # Extract pathway ID
                pathway_info = get_kegg_pathway_info(pathway_kegg_id, base_url=base_url)
                if pathway_info:
                    pathway_list.append(pathway_info)
        return pathway_list
    except requests.exceptions.RequestException as e:
        return f"Error fetching KEGG pathway IDs: {e}"

def get_kegg_pathway_info(pathway_id, base_url=KEGG_BASE_URL):
    """Fetch detailed information for a specific KEGG pathway."""
    get_url = f"{base_url}/get"
    url = f"{get_url}/{pathway_id}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        return response.text # Raw text format, can be parsed further if needed
    except requests.exceptions.RequestException as e:
        return f"Error fetching KEGG pathway info for {pathway_id}: {e}"

# ------------------------------
# Reactome Pathway Data Retrieval
# ------------------------------
def get_reactome_pathways_for_gene(gene_name="FGFR2", base_url=REACTOME_BASE_URL):
    """Fetch Reactome pathways associated with a gene."""
    search_url = f"{base_url}/search/query"
    params = {
        'query': gene_name,
        'species': 'Homo sapiens',
        'facet': 'pathway'
    }
    try:
        response = requests.get(search_url, params=params)
        response.raise_for_status()
        data = response.json()
        if data and 'results' in data:
            return data['results']
        else:
            return "No Reactome pathways found for this gene."
    except requests.exceptions.RequestException as e:
        return f"Error fetching Reactome data: {e}"

# ------------------------------
# Human Protein Atlas Expression Data Retrieval
# ------------------------------
def get_human_protein_atlas_expression(gene_name="FGFR2", base_url=HPA_BASE_URL):
    """Fetch gene expression data from Human Protein Atlas."""
    search_url = f"{base_url}/search_download.php"
    params = {
        'query': gene_name,
        'format': 'json',
        'columns': 'g', # Gene info columns
        'compress': 'no'
    }
    try:
        response = requests.get(search_url, params=params)
        response.raise_for_status()
        data = response.json()
        if data and isinstance(data, list) and len(data) > 0: # Check if list and not empty
            return data[0] # Return the first gene entry
        else:
            return "No expression data found for this gene in Human Protein Atlas."
    except requests.exceptions.RequestException as e:
        return f"Error fetching Human Protein Atlas data: {e}"

# ------------------------------
# MeSH Term Retrieval via PubMed E-utilities
# ------------------------------
def search_pubmed_mesh(term="Bladder Cancer", base_url=PUBMED_EUTILS_BASE_URL):
    """Search PubMed for MeSH terms related to a disease."""
    search_url = f"{base_url}/esearch.fcgi"
    params = {
        'db': 'mesh',
        'term': term,
        'retmode': 'xml',
        'retmax': 10 # Limit to first 10 for example
    }
    try:
        response = requests.get(search_url, params=params)
        response.raise_for_status()
        xml_data = response.text
        root = ET.fromstring(xml_data)
        mesh_ids = [id_element.text for id_element in root.findall('.//Id')]
        return mesh_ids
    except requests.exceptions.RequestException as e:
        return f"Error searching PubMed MeSH: {e}"
    except ET.ParseError as e:
        return f"Error parsing XML response from PubMed MeSH: {e}"

# ------------------------------
# OMIM Data Retrieval
# ------------------------------
def get_omim_gene_disorders(gene_symbol="FGFR2", api_key=OMIM_API_KEY, base_url=OMIM_API_BASE_URL):
    """
    Fetches OMIM entries related to a gene symbol using the OMIM API.

    Args:
        gene_symbol (str): The gene symbol to search for (default "FGFR2").
        api_key (str): Your OMIM API Key.
        base_url (str): The base URL for the OMIM API.

    Returns:
        list or str: A list of entry items if found, otherwise a message indicating no entries were found.
    """
    search_url = f"{base_url}/entry/search"
    params = {
         'search': f'gene:{gene_symbol}',
         'format': 'json',
         'apiKey': api_key
    }
    try:
         response = requests.get(search_url, params=params)
         response.raise_for_status()
         data = response.json()
         # OMIM returns results under the 'omim' key. Adjust this logic if the response structure changes.
         if data and 'omim' in data and 'searchResponse' in data['omim'] and 'entryList' in data['omim']['searchResponse']:
             return data['omim']['searchResponse']['entryList']
         else:
             return "No OMIM entries found for this gene symbol."
    except requests.exceptions.RequestException as e:
         logging.error(f"Error fetching OMIM data: {e}")
         return f"Error fetching OMIM data: {e}"

# ------------------------------
# DisGeNET Data Retrieval
# ------------------------------
def get_disgenet_gene_disease_associations(gene_ncbi_id="351", page_number=0,
                                             base_url="https://api.disgenet.com/api/v1",
                                             api_key=DISGENET_API_KEY):
    """
    Fetch gene-disease associations from DisGeNET for a gene using its NCBI gene ID.
    
    Args:
        gene_ncbi_id (str): The NCBI gene identifier (default "351").
        page_number (int): The page number for pagination (default 0 for the first page).
        base_url (str): The base URL for the DisGeNET API (default points to v1).
        api_key (str): Your DisGeNET API key.
        
    Returns:
        dict or str: The JSON response from the DisGeNET API if successful, 
                     or an error message if something goes wrong.
    """
    import time  # Ensure 'time' is imported if not already

    endpoint = f"{base_url}/gda/summary"
    params = {
        "gene_ncbi_id": str(gene_ncbi_id),
        "page_number": str(page_number)
    }
    headers = {
        "Authorization": api_key,
        "accept": "application/json"
    }
    try:
        response = requests.get(endpoint, params=params, headers=headers, verify=False)
        if response.status_code == 429:
            # Handle rate limiting by waiting before retrying.
            retry_after = int(response.headers.get("x-rate-limit-retry-after-seconds", 60))
            logging.warning(f"Rate limit reached. Retrying after {retry_after} seconds.")
            time.sleep(retry_after)
            response = requests.get(endpoint, params=params, headers=headers, verify=False)
        response.raise_for_status()
        data = response.json()
        if data:
            return data
        else:
            return "No gene-disease associations found for the given gene in DisGeNET."
    except requests.exceptions.RequestException as e:
        return f"Error fetching DisGeNET data: {e}"

# ------------------------------
# ClinicalTrials.gov Data Retrieval
# ------------------------------
def search_clinicaltrials_gov(query_term="Bladder Cancer FGFR2 mutation", output_format="TXT", base_url=CLINICALTRIALS_GOV_BASE_URL):
    """
    Search ClinicalTrials.gov for trials matching a query term.

    Parameters:
      query_term (str): The search expression to query ClinicalTrials.gov.
      output_format (str): The desired output format. Use "JSON" to get the raw JSON
                           (as returned by the API) or "TXT" to get a formatted text summary.
                           Default is "TXT".
      base_url (str): The base URL for the ClinicalTrials.gov API.

    Returns:
      If output_format is "JSON", returns the JSON list of trial data.
      If output_format is "TXT", returns a multi-line text summary of the clinical trials.
      In case of errors or no results, returns an appropriate message.
    """
    search_url = f"{base_url}/study_fields"
    params = {
        'expr': query_term,  # Search expression
        'fields': 'NCTId,BriefTitle,Condition,InterventionName,OverallStatus,StudyType',  # Fields to return
        'fmt': 'json',  # Response format from the API
        'min_rnk': 1,
        'max_rnk': 10  # Limit to first 10 for example
    }
    try:
        response = requests.get(search_url, params=params)
        response.raise_for_status()
        data = response.json()
        if data and 'StudyFieldsResponse' in data and 'StudyFieldsList' in data['StudyFieldsResponse']:
            trials = data['StudyFieldsResponse']['StudyFieldsList']
            if output_format.upper() == "TXT":
                # Build a formatted text summary of the trial data
                text_output = ""
                for idx, trial in enumerate(trials, 1):
                    nct_ids = ", ".join(trial.get("NCTId", [])) if isinstance(trial.get("NCTId", []), list) else trial.get("NCTId", "")
                    brief_titles = ", ".join(trial.get("BriefTitle", [])) if isinstance(trial.get("BriefTitle", []), list) else trial.get("BriefTitle", "")
                    conditions = ", ".join(trial.get("Condition", [])) if isinstance(trial.get("Condition", []), list) else trial.get("Condition", "")
                    interventions = ", ".join(trial.get("InterventionName", [])) if isinstance(trial.get("InterventionName", []), list) else trial.get("InterventionName", "")
                    overall_status = ", ".join(trial.get("OverallStatus", [])) if isinstance(trial.get("OverallStatus", []), list) else trial.get("OverallStatus", "")
                    study_type = ", ".join(trial.get("StudyType", [])) if isinstance(trial.get("StudyType", []), list) else trial.get("StudyType", "")

                    text_output += f"Trial {idx}:\n"
                    text_output += f"  NCT ID: {nct_ids}\n"
                    text_output += f"  Brief Title: {brief_titles}\n"
                    text_output += f"  Condition(s): {conditions}\n"
                    text_output += f"  Intervention(s): {interventions}\n"
                    text_output += f"  Overall Status: {overall_status}\n"
                    text_output += f"  Study Type: {study_type}\n\n"
                return text_output.strip()
            else:
                # Return raw JSON list of clinical trial data if requested
                return trials
        else:
            return "No clinical trials found matching the query."
    except requests.exceptions.RequestException as e:
        return f"Error fetching ClinicalTrials.gov data: {e}"

# ------------------------------
# Europe PMC Literature Search
# ------------------------------
def search_europe_pmc(query_term, output_format="TXT", base_url=EUROPE_PMC_BASE_URL, download_dir="europe_pmc_pdfs"):
    """
    Search Europe PMC for articles related to a query.

    Parameters:
      query_term (str): Search keyword(s) for literature retrieval.
      output_format (str): The format of the returned data. 
                           "TXT" will return a formatted text summary of the articles,
                           "PDF" will attempt to download available PDFs of the articles.
                           Default is "TXT".
      base_url (str): The base URL for the Europe PMC API.
      download_dir (str): Directory to save downloaded PDFs (used when output_format is "PDF").

    Returns:
      For "TXT": A multi-line formatted string summarizing the articles.
      For "PDF": A dictionary mapping article identifiers to the local file paths of downloaded PDFs,
                 or a message if no PDFs are available.
      In case of errors or no results, returns an appropriate message.
    """
    # Build the Europe PMC search URL and parameters
    search_url = f"{base_url}"
    params = {
        'query': query_term,
        'format': 'json',
        'pageSize': 10  # Limit to first 10 for example
    }
    try:
        response = requests.get(search_url, params=params)
        response.raise_for_status()
        data = response.json()
        if data and 'resultList' in data and 'result' in data['resultList']:
            articles = data['resultList']['result']
            if output_format.upper() == "TXT":
                # Build a formatted text summary of each article
                text_output = ""
                for idx, article in enumerate(articles, 1):
                    title = article.get("title", "No Title")
                    authors = article.get("authorString", "No Authors")
                    journal = article.get("journalTitle", "No Journal Info")
                    pub_year = article.get("pubYear", "No Publication Year")
                    doi = article.get("doi", "No DOI")
                    abstract = article.get("abstractText", "No Abstract Available")
                    # Collect available full text links, if any
                    full_text_links = ""
                    if "fullTextUrlList" in article:
                        links = article["fullTextUrlList"].get("fullTextUrl", [])
                        if links:
                            full_text_links = " | ".join([link.get("url", "") for link in links])
                    
                    text_output += f"Article {idx}:\n"
                    text_output += f"  Title: {title}\n"
                    text_output += f"  Authors: {authors}\n"
                    text_output += f"  Journal: {journal}\n"
                    text_output += f"  Publication Year: {pub_year}\n"
                    text_output += f"  DOI: {doi}\n"
                    text_output += f"  Abstract: {abstract}\n"
                    if full_text_links:
                        text_output += f"  Full Text Links: {full_text_links}\n"
                    text_output += "\n"
                return text_output.strip()
            
            elif output_format.upper() == "PDF":
                # Create the download directory if it doesn't exist
                if not os.path.exists(download_dir):
                    os.makedirs(download_dir)
                downloaded_files = {}
                for article in articles:
                    # Prefer a unique identifier: PMID or else use the general article ID.
                    article_id = article.get("pmid") or article.get("id") or f"article_{articles.index(article)}"
                    pdf_url = None
                    if "fullTextUrlList" in article:
                        links = article["fullTextUrlList"].get("fullTextUrl", [])
                        # Look for a link with a documentStyle of 'pdf'
                        for link in links:
                            if link.get("documentStyle", "").lower() == "pdf":
                                pdf_url = link.get("url")
                                break
                    if pdf_url:
                        save_path = os.path.join(download_dir, f"{article_id}.pdf")
                        downloaded_path = download_pdf(pdf_url, save_path)
                        if downloaded_path:
                            downloaded_files[article_id] = downloaded_path
                        else:
                            downloaded_files[article_id] = "Failed to download PDF."
                    else:
                        downloaded_files[article_id] = "No PDF available."
                return downloaded_files
            else:
                return "Unsupported output format. Please choose 'TXT' or 'PDF'."
        else:
            return "No Europe PMC articles found for this query."
    except requests.exceptions.RequestException as e:
        return f"Error searching Europe PMC: {e}"

# ------------------------------
# TCGA Project Data Retrieval
# ------------------------------
def get_tcga_projects(disease_type, base_url=TCGA_GDC_API_BASE_URL):
    """Fetch TCGA projects, optionally filtered by disease type."""
    search_url = f"{base_url}/projects"
    params = {
        'filters': json.dumps({
            'op': 'and',
            'content': [
                {
                    'op': 'in',
                    'content': {
                        'field': 'disease_type',
                        'value': [disease_type]
                    }
                }
            ]
        }),
        'fields': 'project_id,name,disease_type,primary_site',
        'format': 'json',
        'size': 10 # Limit to first 10 for example
    }
    try:
        response = requests.get(search_url, params=params)
        response.raise_for_status()
        data = response.json()
        if data and 'data' in data and 'hits' in data['data']:
            return data['data']['hits']
        else:
            return "No TCGA projects found matching the disease type."
    except requests.exceptions.RequestException as e:
        return f"Error fetching TCGA project data: {e}"

# ------------------------------
# FDA Drug Approvals Search
# ------------------------------
def search_fda_drug_approvals(search_term="bladder cancer", base_url=FDA_API_BASE_URL):
    """Search FDA Drug Approvals API for drugs related to a term."""
    search_url = f"{base_url}/drugsfda.json"
    params = {
        'search': f'indications_and_usage:{search_term}',
        'limit': 10 # Limit to first 10 for example
    }
    try:
        response = requests.get(search_url, params=params)
        response.raise_for_status()
        data = response.json()
        if data and 'results' in data:
            return data['results']
        else:
            return "No FDA drug approvals found matching the search term."
    except requests.exceptions.RequestException as e:
        return f"Error fetching FDA drug approval data: {e}"


# ------------------------------
# Data Processing and Response Generation (Modularized)
# ------------------------------
def process_drug_data(metadata):
    """
    Processes drug names from metadata to fetch PubChem and ChEMBL data concurrently.
    """
    try:
        pubchem_data_all = {}
        chembl_data_all = {}

        if "Drug Information" in metadata and isinstance(metadata["Drug Information"], list):
            drugs = [drug.strip() for drug in metadata["Drug Information"] if isinstance(drug, str) and drug and is_valid_drug_name(drug)]
            if drugs:
                with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
                    # Submit both PubChem and ChEMBL tasks concurrently.
                    future_pubchem = {executor.submit(fetch_pubchem_info, drug): drug for drug in drugs}
                    future_chembl = {executor.submit(fetch_chembl_info, drug): drug for drug in drugs}

                    for future in concurrent.futures.as_completed(future_pubchem):
                        drug = future_pubchem[future]
                        result = future.result()
                        pubchem_data_all[drug] = result
                        if result:
                            logging.info(f"PubChem data fetched for '{drug}'.")
                        else:
                            logging.warning(f"No PubChem data found for '{drug}'.")

                    for future in concurrent.futures.as_completed(future_chembl):
                        drug = future_chembl[future]
                        result = future.result()
                        chembl_data_all[drug] = result
                        if result and result.get('molecules'):
                            logging.info(f"ChEMBL data fetched for '{drug}'.")
                        else:
                            logging.warning(f"No ChEMBL data found for '{drug}'.")
            else:
                logging.info("No valid drug names found in metadata.")
        else:
            logging.info("No 'Drug Information' field or invalid format in metadata.")

        return pubchem_data_all, chembl_data_all
    except Exception as e:
        logging.error(f"Error in process_drug_data: {str(e)}")
        return {}, {}


def process_gene_data(metadata):
    """
    Processes gene names from metadata to fetch UniProt and Ensembl data concurrently.
    """
    gene_data_all = {}

    if "Gene Information" in metadata and isinstance(metadata["Gene Information"], list):
        genes = [gene.strip() for gene in metadata["Gene Information"] if isinstance(gene, str) and gene]
        if genes:
            with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
                futures = {}
                for gene in genes:
                    # Submit independent requests concurrently for each gene.
                    futures[gene] = {
                        'uniprot_id': executor.submit(fetch_uniprot_id, gene),
                        'ensembl_id': executor.submit(fetch_ensembl_gene_id, gene),
                        'uniprot_protein_data': executor.submit(get_uniprot_protein_info, gene),
                        'go_data': executor.submit(get_gene_ontology_info, gene),
                        'kegg_pathways': executor.submit(get_kegg_pathways_for_gene, gene),
                        'reactome_pathways': executor.submit(get_reactome_pathways_for_gene, gene),
                        'hpa_expression_data': executor.submit(get_human_protein_atlas_expression, gene)
                    }
                # Collect each gene's results.
                for gene, futs in futures.items():
                    uniprot_id = futs['uniprot_id'].result()
                    # Fetch the UniProt function (depends on uniprot_id)
                    uniprot_function = fetch_and_extract_uniprot_function(uniprot_id) if uniprot_id else None
                    ensembl_id = futs['ensembl_id'].result()
                    ensembl_info = fetch_ensembl_gene_info(ensembl_id) if ensembl_id else {}
                    gene_data_all[gene] = {
                        "UniProt ID": uniprot_id,
                        "UniProt Function": uniprot_function,
                        "Ensembl ID": ensembl_id,
                        "Ensembl Info": ensembl_info,
                        "UniProt Protein Data": futs['uniprot_protein_data'].result(),
                        "Gene Ontology Data": futs['go_data'].result(),
                        "KEGG Pathways": futs['kegg_pathways'].result(),
                        "Reactome Pathways": futs['reactome_pathways'].result(),
                        "Human Protein Atlas Expression": futs['hpa_expression_data'].result(),
                    }
                    logging.info(f"Gene data fetched for '{gene}'.")
        else:
            logging.info("No valid gene names found in metadata.")
    else:
        logging.info("No 'Gene Information' field or invalid format in metadata.")

    return gene_data_all

def process_disease_data(metadata):
    """
    Processes disease names from metadata to fetch MeSH, OMIM, and DisGeNET data concurrently.
    """
    disease_data_all = {}
    if "Disease Information" in metadata and isinstance(metadata["Disease Information"], list):
        diseases = [disease.strip() for disease in metadata["Disease Information"] if isinstance(disease, str) and disease]
        if diseases:
            with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
                futures = {}
                for disease in diseases:
                    futures[disease] = {
                        'mesh_terms': executor.submit(search_pubmed_mesh, disease),
                        'omim_entries': executor.submit(get_omim_gene_disorders, disease),
                        'disgenet_associations': executor.submit(get_disgenet_gene_disease_associations, disease)
                    }
                for disease, futs in futures.items():
                    disease_data_all[disease] = {
                        "MeSH Terms": futs['mesh_terms'].result(),
                        "OMIM Entries": futs['omim_entries'].result(),
                        "DisGeNET Associations": futs['disgenet_associations'].result(),
                    }
                    logging.info(f"Disease data fetched for '{disease}'.")
        else:
            logging.info("No valid disease names found in metadata.")
    else:
        logging.info("No 'Disease Information' field or invalid format in metadata.")
    return disease_data_all

def process_clinical_trial_data(metadata):
    """
    Processes clinical trial queries from metadata and fetches data concurrently.
    """
    clinical_trial_data_all = {}
    if "Clinical Trial Data" in metadata and isinstance(metadata["Clinical Trial Data"], list):
        queries = [query.strip() for query in metadata["Clinical Trial Data"] if isinstance(query, str) and query]
        if queries:
            with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
                future_trials = {executor.submit(search_clinicaltrials_gov, query): query for query in queries}
                for future in concurrent.futures.as_completed(future_trials):
                    query = future_trials[future]
                    clinical_trial_data_all[query] = future.result()
                    logging.info(f"Clinical trial data fetched for query: '{query}'.")
        else:
            logging.info("No valid clinical trial queries found in metadata.")
    else:
        logging.info("No 'Clinical Trial Data' field or invalid format in metadata.")
    return clinical_trial_data_all

def process_literature_data(metadata):
    """
    Processes literature queries from metadata and fetches Europe PMC data concurrently.
    """
    literature_data_all = {}
    if "Literature and Evidence" in metadata and isinstance(metadata["Literature and Evidence"], list):
        queries = [query.strip() for query in metadata["Literature and Evidence"] if isinstance(query, str) and query]
        if queries:
            with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
                future_lit = {executor.submit(search_europe_pmc, query): query for query in queries}
                for future in concurrent.futures.as_completed(future_lit):
                    query = future_lit[future]
                    literature_data_all[query] = future.result()
                    logging.info(f"Literature data fetched for query: '{query}'.")
        else:
            logging.info("No valid literature queries found in metadata.")
    else:
        logging.info("No 'Literature and Evidence' field or invalid format in metadata.")
    return literature_data_all

def process_omics_data(metadata):
    """Processes disease names from metadata to fetch TCGA project data."""
    omics_data_all = {}
    if "Omics Data" in metadata and isinstance(metadata["Omics Data"], list):
        diseases = [disease.strip() for disease in metadata["Omics Data"] if isinstance(disease, str) and disease]
        if diseases:
            for disease in diseases:
                tcga_projects = get_tcga_projects(disease)
                omics_data_all[disease] = tcga_projects
                logging.info(f"TCGA project data fetched for disease: '{disease}'.")
        else:
            logging.info("No valid disease names found in omics metadata.")
    else:
        logging.info("No 'Omics Data' field or invalid format in metadata.")
    return omics_data_all

def process_fda_data(metadata):
    """Processes disease names from metadata to fetch FDA drug approvals data."""
    fda_data_all = {}
    if "Other Important Resources" in metadata and isinstance(metadata["Other Important Resources"], dict) and "FDA" in metadata["Other Important Resources"] and isinstance(metadata["Other Important Resources"]["FDA"], list):
        search_terms = [term.strip() for term in metadata["Other Important Resources"]["FDA"] if isinstance(term, str) and term]
        if search_terms:
            for term in search_terms:
                fda_approvals = search_fda_drug_approvals(term)
                fda_data_all[term] = fda_approvals
                logging.info(f"FDA drug approvals data fetched for term: '{term}'.")
        else:
            logging.info("No valid FDA search terms found in metadata.")
    else:
        logging.info("No 'FDA' terms or invalid format in 'Other Important Resources' metadata.")
    return fda_data_all


def generate_response(user_query, metadata, model_name=RESPONSE_GENERATION_MODEL):
    """
    Generates a response based on the type of query.
    """
    if metadata.get("General Query", False):
        # If it's a general query, use METADATA_MODEL for response
        metadata_response = chat(
            model=METADATA_EXTRACTION_MODEL,
            messages=[{'role': 'user', 'content': user_query}],
            stream=False,
        )
        response = metadata_response['message']['content']
        logging.info(f"Generated response for general query using {METADATA_EXTRACTION_MODEL}.")
    else:
        # If not a general query, proceed with DeepSeek response generation
        pubchem_data_all = metadata.get("PubChem Data", {})
        chembl_data_all = metadata.get("ChEMBL Data", {})
        gene_data_all = metadata.get("Gene Information", {})
        disease_data_all = metadata.get("Disease Information", {})
        clinical_trial_data_all = metadata.get("Clinical Trial Data", {})
        literature_data_all = metadata.get("Literature and Evidence", {})
        omics_data_all = metadata.get("Omics Data", {})
        fda_data_all = metadata.get("FDA", {})

        final_query = (
            user_query +
            "\n\nPubChem Data:\n" + json.dumps(pubchem_data_all, indent=2) +
            "\n\nChEMBL Data:\n" + json.dumps(chembl_data_all, indent=2) +
            "\n\nGene/Protein Data:\n" + json.dumps(gene_data_all, indent=2) +
            "\n\nDisease Data:\n" + json.dumps(disease_data_all, indent=2) +
            "\n\nClinical Trial Data:\n" + json.dumps(clinical_trial_data_all, indent=2) +
            "\n\nLiterature Data:\n" + json.dumps(literature_data_all, indent=2) +
            "\n\nOmics Data (TCGA Projects):\n" + json.dumps(omics_data_all, indent=2) +
            "\n\nFDA Drug Approvals Data:\n" + json.dumps(fda_data_all, indent=2) +
            "\n\nYou are a drug repurposing scientist. Please analyze the data above in detail and provide a comprehensive answer, focusing on potential drug repurposing opportunities for the specific patient case described in the user query."
        )

        print("\nDeepSeek Query:\n")
        print(final_query)

        try:
            stream = chat(
                model=model_name,
                messages=[{'role': 'user', 'content': final_query}],
                stream=True,
            )

            print("\nDeepSeek Response:\n")
            full_response = ""
            for chunk in stream:
                response_chunk = chunk['message']['content']
                print(response_chunk, end='', flush=True)
                full_response += response_chunk
            response = full_response
            logging.info("Generated DeepSeek R1 7b response for specific query.")
        except Exception as e:
            logging.error(f"Error during DeepSeek response generation: {e}")
            response = "Error generating response."

    return response

def get_pdf_url(article):
    """
    Extracts the free full-text PDF link from a PubMed article summary if a PubMed Central (PMC) ID is available.

    Args:
        article (dict): A dictionary representing a PubMed article summary.

    Returns:
        str or None: The URL to the PDF if available, otherwise None.
    """
    articleids = article.get("articleids", [])
    for id_entry in articleids:
        if id_entry.get("idtype", "").lower() == "pmc":
            pmc_id = id_entry.get("value")
            # Ensure the PMC ID starts with 'PMC'
            if not pmc_id.startswith("PMC"):
                pmc_id = f"PMC{pmc_id}"
            pdf_url = f"https://www.ncbi.nlm.nih.gov/pmc/articles/{pmc_id}/pdf/"
            return pdf_url
    return None

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

def search_pubmed(query, retmax=10):
    """
    Search PubMed for articles related to the given query using the NCBI E-utilities.

    This function performs two steps:
      1. Uses the ESearch API to retrieve a list of PubMed IDs matching the query.
      2. Uses the ESummary API to retrieve detailed summaries for each PubMed ID.

    Args:
        query (str): The search query.
        retmax (int): Maximum number of articles to retrieve (default is 10).

    Returns:
        list: A list of dictionaries containing PubMed article summaries.
    """
    # Step 1: ESearch API
    esearch_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    esearch_params = {
        "db": "pubmed",
        "term": query,
        "retmode": "json",
        "retmax": retmax
    }
    try:
        esearch_resp = requests.get(esearch_url, params=esearch_params)
        esearch_resp.raise_for_status()
        esearch_data = esearch_resp.json()
    except Exception as e:
        logging.error(f"Error during PubMed ESearch: {e}")
        return []

    id_list = esearch_data.get("esearchresult", {}).get("idlist", [])
    if not id_list:
        logging.info("No PubMed IDs found for the query.")
        return []

    # Step 2: ESummary API
    esummary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
    esummary_params = {
        "db": "pubmed",
        "id": ",".join(id_list),
        "retmode": "json"
    }
    try:
        esummary_resp = requests.get(esummary_url, params=esummary_params)
        esummary_resp.raise_for_status()
        esummary_data = esummary_resp.json()
    except Exception as e:
        logging.error(f"Error during PubMed ESummary: {e}")
        return []

    result = esummary_data.get("result", {})
    uids = result.get("uids", [])
    articles = [result[uid] for uid in uids if uid in result]
    return articles

def process_literature_data(user_query, retmax=10, mode="TXT"):
    """
    Searches PubMed for articles related to the query and processes the literature data based on the specified mode.

    Args:
        user_query (str): The search query.
        retmax (int): Maximum number of articles to retrieve (default is 10).
        mode (str): The mode to process literature data. Can be "TXT" for text summaries or "PDF" for PDF downloads.

    Returns:
        str or dict: Returns a string containing a nicely formatted text summary in TXT mode, or a dictionary mapping article identifiers to file paths in PDF mode.
    """
    articles = search_pubmed(user_query, retmax)
    if not articles:
        logging.info("No PubMed articles found.")
        return {}

    if mode == "TXT":
        # Process articles for text summaries
        text_summaries = [article.get("summary", "") for article in articles]
        return "\n\n".join(text_summaries)
    elif mode == "PDF":
        # Process articles for PDF downloads
        if not os.path.exists("pubmed_pdfs"):
            os.makedirs("pubmed_pdfs")
        
        # Build a mapping of PubMed ID to PDF URL.
        pdf_links = {}
        for article in articles:
            pdf_url = get_pdf_url(article)
            if pdf_url:
                pmid = article.get("uid", "unknown")
                pdf_links[pmid] = pdf_url

        downloaded_files = {}
        with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
            future_to_pmid = {
                executor.submit(download_pdf, pdf_url, os.path.join("pubmed_pdfs", f"{pmid}.pdf")): pmid
                for pmid, pdf_url in pdf_links.items()
            }
            for future in concurrent.futures.as_completed(future_to_pmid):
                pmid = future_to_pmid[future]
                file_path = future.result()
                if file_path:
                    downloaded_files[pmid] = file_path
        return downloaded_files
    else:
        logging.error(f"Invalid mode specified: {mode}")
        return {}

def download_pubmed_pdfs(query, retmax=10, save_dir="pubmed_pdfs"):
    """
    Searches PubMed for articles using the given query, extracts free full-text PDF links, and downloads the PDFs.

    Parameters:
      query (str): The search query for PubMed.
      retmax (int): The maximum number of PubMed articles to process (default is 10).
      save_dir (str): Directory to save downloaded PDFs.

    Returns:
      dict: A dictionary mapping each PubMed article's unique identifier (uid) to the local file path of the downloaded PDF,
            or a message indicating if a PDF was not available or if the download failed.
    """
    # Search PubMed for articles matching the query
    articles = search_pubmed(query, retmax)
    if not articles:
        logging.info("No PubMed articles found for the query.")
        return {}

    # Ensure the directory to store PDFs exists
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)

    downloaded_files = {}
    for article in articles:
        # Extract PDF URL (if available) using the helper function
        pdf_url = get_pdf_url(article)
        article_id = article.get("uid", f"article_{articles.index(article)}")
        if pdf_url:
            save_path = os.path.join(save_dir, f"{article_id}.pdf")
            downloaded_path = download_pdf(pdf_url, save_path)
            if downloaded_path:
                downloaded_files[article_id] = downloaded_path
            else:
                downloaded_files[article_id] = "Failed to download PDF."
        else:
            downloaded_files[article_id] = "No PDF available."
    return downloaded_files

# ------------------------------
# Main Script Execution
# ------------------------------
if __name__ == "__main__":
    # Get the User Query (for testing, using a hardcoded query)
    user_query = ""

    # Extract Metadata (optional, for more structured queries later)
    metadata = extract_metadata(user_query)
    print("\nExtracted Metadata:")
    print(json.dumps(metadata, indent=2))

    # Process Drug Data (optional, if needed for drug-specific details)
    pubchem_data_all, chembl_data_all = process_drug_data(metadata)
    print("\nPubChem Data Summary:")
    print(json.dumps(pubchem_data_all, indent=2))
    print("\nChEMBL Data Summary:")
    print(json.dumps(chembl_data_all, indent=2))

    # Process Gene Data (optional, if needed for gene-specific details)
    gene_data_all = process_gene_data(metadata)
    print("\nGene Data Summary:")
    print(json.dumps(gene_data_all, indent=2))

    # Process Disease Data (optional, if needed for disease-specific details)
    disease_data_all = process_disease_data(metadata)
    print("\nDisease Data Summary:")
    print(json.dumps(disease_data_all, indent=2))

    # Process Clinical Trial Data (optional, if needed to pre-fetch trials)
    clinical_trial_data_all = process_clinical_trial_data(metadata)
    print("\nClinical Trial Data Summary:")
    print(json.dumps(clinical_trial_data_all, indent=2))

    # Process Literature Data (optional, if needed to pre-fetch literature)
    literature_data_all = process_literature_data(metadata)
    print("\nLiterature Data Summary:")
    print(json.dumps(literature_data_all, indent=2))

    # Process Omics Data (optional, if needed for omics data)
    omics_data_all = process_omics_data(metadata)
    print("\nOmics Data Summary (TCGA Projects):")
    print(json.dumps(omics_data_all, indent=2))

    # Process FDA Data (optional, if needed for FDA approvals)
    fda_data_all = process_fda_data(metadata)
    print("\nFDA Drug Approvals Data Summary:")
    print(json.dumps(fda_data_all, indent=2))

    # Download PDFs from PubMed (existing function)
    downloaded_pdfs = download_pubmed_pdfs(user_query, retmax=10, save_dir="pubmed_pdfs")
    print("\nDownloaded PubMed PDFs:")
    print(json.dumps(downloaded_pdfs, indent=2))

    # Generate DeepSeek Response (optional; remains unchanged)
    deepseek_response = generate_response(user_query, metadata)

    print("\n--- End of Script ---")