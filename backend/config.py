
import os


METADATA_EXTRACTION_MODEL = os.getenv("METADATA_EXTRACTION_MODEL",'llama3.2:latest')
RESPONSE_GENERATION_MODEL = os.getenv("RESPONSE_GENERATION_MODEL",'deepseek-r1:7b')
# API Base URLs
PUBCHEM_BASE_URL = os.getenv("PUBCHEM_BASE_URL", "https://pubchem.ncbi.nlm.nih.gov/rest/pug")
CHEMBL_BASE_URL = os.getenv("CHEMBL_BASE_URL", "https://www.ebi.ac.uk/chembl/api/data/molecule/search")
UNIPROT_BASE_URL = os.getenv("UNIPROT_BASE_URL", "https://rest.uniprot.org/uniprotkb")
ENSEMBL_BASE_URL = os.getenv("ENSEMBL_BASE_URL", "https://rest.ensembl.org")
RXNAV_BASE_URL = os.getenv("RXNAV_BASE_URL", "https://mor.nlm.nih.gov/RxNav")
KEGG_BASE_URL = os.getenv("KEGG_BASE_URL", "https://rest.kegg.jp")
REACTOME_BASE_URL = os.getenv("REACTOME_BASE_URL", "https://reactome.org/ContentService")
HPA_BASE_URL = os.getenv("HPA_BASE_URL", "https://www.proteinatlas.org/api")
PUBMED_EUTILS_BASE_URL = os.getenv("PUBMED_EUTILS_BASE_URL", "https://eutils.ncbi.nlm.nih.gov/entrez/eutils")
OMIM_API_BASE_URL = os.getenv("OMIM_API_BASE_URL", "https://api.omim.org/api")
DISGENET_API_BASE_URL = os.getenv("DISGENET_API_BASE_URL", "https://www.disgenet.org/api")
CLINICALTRIALS_GOV_BASE_URL = os.getenv("CLINICALTRIALS_GOV_BASE_URL", "https://clinicaltrials.gov/api/query")
EUROPE_PMC_BASE_URL = os.getenv("EUROPE_PMC_BASE_URL", "https://www.ebi.ac.uk/europepmc/webservices/rest/search")
TCGA_GDC_API_BASE_URL = os.getenv("TCGA_GDC_API_BASE_URL", "https://api.gdc.cancer.gov")
FDA_API_BASE_URL = os.getenv("FDA_API_BASE_URL", "https://api.fda.gov/drug")

# API Keys
OMIM_API_KEY = os.getenv("OMIM_API_KEY", "3FZX0uYkSVmNO4Ako4d8yg")
DISGENET_API_KEY = os.getenv("DISGENET_API_KEY", "87dec4c5-14f3-43ec-8432-9b16bc1a0b26")

# Constants
DRUG_NAME_BLACKLIST = ["primary targets", "targets", "target", "pathway", "information"]
PUBCHEM_PROPERTIES = "MolecularWeight,IUPACName,CanonicalSMILES"
