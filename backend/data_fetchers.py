import requests
import json
import logging
import xml.etree.ElementTree as ET
from config import (
    PUBCHEM_BASE_URL, CHEMBL_BASE_URL, UNIPROT_BASE_URL, ENSEMBL_BASE_URL,
    RXNAV_BASE_URL, KEGG_BASE_URL, REACTOME_BASE_URL, HPA_BASE_URL,
    PUBMED_EUTILS_BASE_URL, OMIM_API_BASE_URL, DISGENET_API_BASE_URL,
    CLINICALTRIALS_GOV_BASE_URL, EUROPE_PMC_BASE_URL, TCGA_GDC_API_BASE_URL,
    FDA_API_BASE_URL, OMIM_API_KEY, DISGENET_API_KEY, PUBCHEM_PROPERTIES
)
from utils import fetch_data_from_api

# PubChem Data Retrieval
def fetch_pubchem_info(drug_name):
    """
    Query PubChem for drug information.
    """
    url = f"{PUBCHEM_BASE_URL}/compound/name/{drug_name}/property/{PUBCHEM_PROPERTIES}/JSON"
    return fetch_data_from_api(url)

# ChEMBL Data Retrieval
def fetch_chembl_info(drug_name):
    """
    Query ChEMBL for drug information.
    """
    params = {'q': drug_name, 'format': 'json'}
    url = f"{CHEMBL_BASE_URL}"
    try:
        response = requests.get(url, params=params)
        response.raise_for_status()
        return response.json()
    except requests.exceptions.HTTPError as http_err:
        logging.error(f"HTTP error fetching ChEMBL data for '{drug_name}': {http_err}")
    except requests.exceptions.RequestException as req_err:
        logging.error(f"Request error fetching ChEMBL data for '{drug_name}': {req_err}")
    except json.JSONDecodeError as json_err:
        logging.error(f"JSON decode error from ChEMBL for '{drug_name}': {json_err}")
    return None

# UniProt Data Retrieval
def fetch_uniprot_id(gene_name):
    """
    Fetch UniProt ID for a gene name.
    """
    url = f"{UNIPROT_BASE_URL}/search?query={gene_name}&fields=accession&format=json"
    data = fetch_data_from_api(url)
    if data and 'results' in data and data['results']:
        return data['results'][0]['primaryAccession']
    logging.warning(f"UniProt ID not found for {gene_name}")
    return None

def fetch_and_extract_uniprot_function(uniprot_id):
    """
    Fetch UniProt data and extract the 'Function' section.
    """
    url = f"{UNIPROT_BASE_URL}/{uniprot_id}.txt"
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
        logging.warning("'Function' section not found in UniProt data for {uniprot_id}")
    except requests.exceptions.HTTPError as http_err:
        logging.error(f"HTTP error fetching UniProt data for {uniprot_id}: {http_err}")
    except requests.exceptions.RequestException as req_err:
        logging.error(f"Request error fetching UniProt data for {uniprot_id}: {req_err}")
    return None

def get_uniprot_protein_info(gene_name):
    """Fetch protein information from UniProt using gene name."""
    search_url = f"{UNIPROT_BASE_URL}/search"
    params = {
        'query': f'gene_name:{gene_name}',
        'format': 'json'
    }
    try:
        response = requests.get(search_url, params=params)
        response.raise_for_status()
        data = response.json()
        if data and 'results' in data and data['results']:
            protein_info = data['results'][0]
            return protein_info
        else:
            return "No protein information found for this gene name."
    except requests.exceptions.RequestException as e:
        return f"Error fetching UniProt data: {e}"

# Ensembl Data Retrieval
def fetch_ensembl_gene_id(gene_name):
    """
    Fetch Ensembl gene ID for a gene name.
    """
    url = f"{ENSEMBL_BASE_URL}/xrefs/symbol/homo_sapiens/{gene_name}?content-type=application/json"
    data = fetch_data_from_api(url)
    if data and len(data) > 0:
        return data[0]['id']
    logging.warning(f"Ensembl gene ID not found for {gene_name}")
    return None

def fetch_ensembl_gene_info(gene_id):
    """
    Fetch detailed gene information from Ensembl using gene ID.
    """
    url = f"{ENSEMBL_BASE_URL}/lookup/id/{gene_id}?expand=1&content-type=application/json"
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
    logging.warning(f"Ensembl gene details not found for {gene_id}")
    return {}

# RxNav Data Retrieval
def fetch_drug_info(drug_name):
    """
    Fetch drug information from NIH RxNav API.
    """
    url = f"{RXNAV_BASE_URL}/search?searchBy=String&searchTerm={drug_name}"
    response = requests.get(url)
    if response.status_code == 200:
        return url
    logging.warning(f"Failed to fetch drug information for {drug_name}, status code {response.status_code}")
    return None

# Gene Ontology Data Retrieval
def get_gene_ontology_info(gene_symbol):
    """Fetch Gene Ontology terms associated with a gene symbol."""
    base_url = "https://www.ebi.ac.uk/QuickGO/services/annotation/search"
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

# KEGG Pathway Data Retrieval
def get_kegg_pathways_for_gene(gene_name, organism="hsa"):
    """Fetch KEGG pathways associated with a gene."""
    link_url = f"{KEGG_BASE_URL}/link/pathway"
    url = f"{link_url}/{organism}:{gene_name}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        pathway_ids = response.text.strip().split('\n')
        pathway_list = []
        if pathway_ids and pathway_ids != ['']:
            for pathway_id in pathway_ids:
                pathway_kegg_id = pathway_id.split(":")[1]
                pathway_info = get_kegg_pathway_info(pathway_kegg_id)
                if pathway_info:
                    pathway_list.append(pathway_info)
        return pathway_list
    except requests.exceptions.RequestException as e:
        return f"Error fetching KEGG pathway IDs: {e}"

def get_kegg_pathway_info(pathway_id):
    """Fetch detailed information for a specific KEGG pathway."""
    get_url = f"{KEGG_BASE_URL}/get"
    url = f"{get_url}/{pathway_id}"
    try:
        response = requests.get(url)
        response.raise_for_status()
        return response.text
    except requests.exceptions.RequestException as e:
        return f"Error fetching KEGG pathway info for {pathway_id}: {e}"

# Reactome Pathway Data Retrieval
def get_reactome_pathways_for_gene(gene_name):
    """Fetch Reactome pathways associated with a gene."""
    search_url = f"{REACTOME_BASE_URL}/search/query"
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

# Human Protein Atlas Data Retrieval
def get_human_protein_atlas_expression(gene_name):
    """Fetch gene expression data from Human Protein Atlas."""
    search_url = f"{HPA_BASE_URL}/search_download.php"
    params = {
        'query': gene_name,
        'format': 'json',
        'columns': 'g',
        'compress': 'no'
    }
    try:
        response = requests.get(search_url, params=params)
        response.raise_for_status()
        data = response.json()
        if data and isinstance(data, list) and len(data) > 0:
            return data[0]
        else:
            return "No expression data found for this gene in Human Protein Atlas."
    except requests.exceptions.RequestException as e:
        return f"Error fetching Human Protein Atlas data: {e}"

# MeSH Term Retrieval via PubMed E-utilities
def search_pubmed_mesh(term):
    """Search PubMed for MeSH terms related to a disease."""
    search_url = f"{PUBMED_EUTILS_BASE_URL}/esearch.fcgi"
    params = {
        'db': 'mesh',
        'term': term,
        'retmode': 'xml',
        'retmax': 10
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

# OMIM Data Retrieval
def get_omim_gene_disorders(gene_symbol):
    """
    Fetches OMIM entries related to a gene symbol using the OMIM API.
    """
    search_url = f"{OMIM_API_BASE_URL}/entry/search"
    params = {
         'search': f'gene:{gene_symbol}',
         'format': 'json',
         'apiKey': OMIM_API_KEY
    }
    try:
         response = requests.get(search_url, params=params)
         response.raise_for_status()
         data = response.json()
         if data and 'omim' in data and 'searchResponse' in data['omim'] and 'entryList' in data['omim']['searchResponse']:
             return data['omim']['searchResponse']['entryList']
         else:
             return "No OMIM entries found for this gene symbol."
    except requests.exceptions.RequestException as e:
         logging.error(f"Error fetching OMIM data: {e}")
         return f"Error fetching OMIM data: {e}"

# DisGeNET Data Retrieval
def get_disgenet_gene_disease_associations(gene_ncbi_id, page_number=0):
    """
    Fetch gene-disease associations from DisGeNET for a gene using its NCBI gene ID.
    """
    import time
    endpoint = f"{DISGENET_API_BASE_URL}/gda/summary"
    params = {
        "gene_ncbi_id": str(gene_ncbi_id),
        "page_number": str(page_number)
    }
    headers = {
        "Authorization": DISGENET_API_KEY,
        "accept": "application/json"
    }
    try:
        response = requests.get(endpoint, params=params, headers=headers, verify=False)
        if response.status_code == 429:
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

# ClinicalTrials.gov Data Retrieval
def search_clinicaltrials_gov(query_term):
    """Search ClinicalTrials.gov for trials matching a query term."""
    search_url = f"{CLINICALTRIALS_GOV_BASE_URL}/study_fields"
    params = {
        'expr': query_term,
        'fields': 'NCTId,BriefTitle,Condition,InterventionName,OverallStatus,StudyType',
        'fmt': 'json',
        'min_rnk': 1,
        'max_rnk': 10
    }
    try:
        response = requests.get(search_url, params=params)
        response.raise_for_status()
        data = response.json()
        if data and 'StudyFieldsResponse' in data and 'StudyFieldsList' in data['StudyFieldsResponse']:
            return data['StudyFieldsResponse']['StudyFieldsList']
        else:
            return "No clinical trials found matching the query."
    except requests.exceptions.RequestException as e:
        return f"Error fetching ClinicalTrials.gov data: {e}"

# Europe PMC Literature Search
def search_europe_pmc(query_term):
    """Search Europe PMC for articles related to a query."""
    search_url = f"{EUROPE_PMC_BASE_URL}"
    params = {
        'query': query_term,
        'format': 'json',
        'pageSize': 10
    }
    try:
        response = requests.get(search_url, params=params)
        response.raise_for_status()
        data = response.json()
        if data and 'resultList' in data and 'result' in data['resultList']:
            return data['resultList']['result']
        else:
            return "No Europe PMC articles found for this query."
    except requests.exceptions.RequestException as e:
        return f"Error searching Europe PMC: {e}"

# TCGA Project Data Retrieval
def get_tcga_projects(disease_type):
    """Fetch TCGA projects, optionally filtered by disease type."""
    search_url = f"{TCGA_GDC_API_BASE_URL}/projects"
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
        'size': 10
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

# FDA Drug Approvals Search
def search_fda_drug_approvals(search_term):
    """Search FDA Drug Approvals API for drugs related to a term."""
    search_url = f"{FDA_API_BASE_URL}/drugsfda.json"
    params = {
        'search': f'indications_and_usage:{search_term}',
        'limit': 10
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
