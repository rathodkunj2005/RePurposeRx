import json
import logging
import concurrent.futures
from data_fetchers import (
    fetch_pubchem_info, fetch_chembl_info, fetch_uniprot_id, fetch_and_extract_uniprot_function,
    get_uniprot_protein_info, fetch_ensembl_gene_id, fetch_ensembl_gene_info, get_gene_ontology_info,
    get_kegg_pathways_for_gene, get_reactome_pathways_for_gene, get_human_protein_atlas_expression,
    search_pubmed_mesh, get_omim_gene_disorders, get_disgenet_gene_disease_associations,
    search_clinicaltrials_gov, search_europe_pmc, get_tcga_projects, search_fda_drug_approvals
)
from utils import get_pdf_url, download_pdf, search_pubmed, download_pubmed_pdfs


# Process Drug Data
def process_drug_data(metadata):
    """
    Processes drug names from metadata to fetch PubChem and ChEMBL data concurrently.
    """
    pubchem_data_all = {}
    chembl_data_all = {}

    if "Drug Information" in metadata and isinstance(metadata["Drug Information"], list):
        drugs = [drug.strip() for drug in metadata["Drug Information"] if isinstance(drug, str) and drug]
        if drugs:
            with concurrent.futures.ThreadPoolExecutor(max_workers=10) as executor:
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

# Process Gene Data
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
                    futures[gene] = {
                        'uniprot_id': executor.submit(fetch_uniprot_id, gene),
                        'ensembl_id': executor.submit(fetch_ensembl_gene_id, gene),
                        'uniprot_protein_data': executor.submit(get_uniprot_protein_info, gene),
                        'go_data': executor.submit(get_gene_ontology_info, gene),
                        'kegg_pathways': executor.submit(get_kegg_pathways_for_gene, gene),
                        'reactome_pathways': executor.submit(get_reactome_pathways_for_gene, gene),
                        'hpa_expression_data': executor.submit(get_human_protein_atlas_expression, gene)
                    }
                for gene, futs in futures.items():
                    uniprot_id = futs['uniprot_id'].result()
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

# Process Disease Data
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

# Process Clinical Trial Data
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

# Process Literature Data
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

# Process Omics Data
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

# Process FDA Data
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
