import json
import logging
from ollama import chat
from config import METADATA_EXTRACTION_MODEL, RESPONSE_GENERATION_MODEL
from data_processors import (
    process_drug_data, process_gene_data, process_disease_data, process_clinical_trial_data,
    process_literature_data, process_omics_data, process_fda_data
)
from utils import download_pubmed_pdfs


# Initialize logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Metadata Extraction Function
def extract_metadata(user_query, model_name=METADATA_EXTRACTION_MODEL):
    """
    Extracts relevant metadata from the user query using a specified LLM model.
    """
    metadata_prompt = f"""
    You are an expert at extracting relevant entity names and data from a query.
    Please analyze the following user query and extract the information into the following JSON structure.
    For each category, output a list of names or pieces of data that are relevant. If no relevant data exists for a category, use an empty list.
    Output only the JSON object with no additional commentary. Also consider potentially important related entities that might help answer the user query, put that into the json.

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
      }}
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

# Generate DeepSeek Response
def generate_deepseek_response(user_query, pubchem_data_all, chembl_data_all, gene_data_all, disease_data_all, clinical_trial_data_all, literature_data_all, omics_data_all, fda_data_all, model_name=RESPONSE_GENERATION_MODEL):
    """
    Generates a detailed response using DeepSeek based on fetched data.
    """
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
        return full_response

    except Exception as e:
        logging.error(f"Error during DeepSeek response generation: {e}")
        return "Error generating response."

# Main Script Execution
if __name__ == "__main__":
    # Get the User Query (for testing, using a hardcoded query)
    user_query = "FGFR2 mutation bladder cancer drug repurposing"

    # Extract Metadata
    metadata = extract_metadata(user_query)
    print("\nExtracted Metadata:")
    print(json.dumps(metadata, indent=2))

    # Process Drug Data
    pubchem_data_all, chembl_data_all = process_drug_data(metadata)
    print("\nPubChem Data Summary:")
    print(json.dumps(pubchem_data_all, indent=2))
    print("\nChEMBL Data Summary:")
    print(json.dumps(chembl_data_all, indent=2))

    # Process Gene Data
    gene_data_all = process_gene_data(metadata)
    print("\nGene Data Summary:")
    print(json.dumps(gene_data_all, indent=2))

    # Process Disease Data
    disease_data_all = process_disease_data(metadata)
    print("\nDisease Data Summary:")
    print(json.dumps(disease_data_all, indent=2))

    # Process Clinical Trial Data
    clinical_trial_data_all = process_clinical_trial_data(metadata)
    print("\nClinical Trial Data Summary:")
    print(json.dumps(clinical_trial_data_all, indent=2))

    # Process Literature Data
    literature_data_all = process_literature_data(metadata)
    print("\nLiterature Data Summary:")
    print(json.dumps(literature_data_all, indent=2))

    # Process Omics Data
    omics_data_all = process_omics_data(metadata)
    print("\nOmics Data Summary (TCGA Projects):")
    print(json.dumps(omics_data_all, indent=2))

    # Process FDA Data
    fda_data_all = process_fda_data(metadata)
    print("\nFDA Drug Approvals Data Summary:")
    print(json.dumps(fda_data_all, indent=2))

    # Download PDFs from PubMed
    downloaded_pdfs = download_pubmed_pdfs(user_query, retmax=10, save_dir="pubmed_pdfs")
    print("\nDownloaded PubMed PDFs:")
    print(json.dumps(downloaded_pdfs, indent=2))

    # Generate DeepSeek Response
    deepseek_response = generate_deepseek_response(user_query, pubchem_data_all, chembl_data_all, gene_data_all, disease_data_all, clinical_trial_data_all, literature_data_all, omics_data_all, fda_data_all)

    print("\n--- End of Script ---")
