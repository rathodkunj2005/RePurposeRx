from flask import Flask, request, jsonify
from flask_cors import CORS
import main
import logging
from functools import wraps
import time
import concurrent.futures
from werkzeug.serving import WSGIRequestHandler

app = Flask(__name__)
CORS(app)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

# Increase the timeout for the WSGI server
WSGIRequestHandler.timeout = 3000  # 5 minutes

def process_data_concurrently(metadata):
    """Process all data sources concurrently using ThreadPoolExecutor"""
    with concurrent.futures.ThreadPoolExecutor(max_workers=8) as executor:
        # Submit all tasks
        drug_future = executor.submit(main.process_drug_data, metadata)
        gene_future = executor.submit(main.process_gene_data, metadata)
        disease_future = executor.submit(main.process_disease_data, metadata)
        clinical_trial_future = executor.submit(main.process_clinical_trial_data, metadata)
        literature_future = executor.submit(main.process_literature_data, metadata)
        omics_future = executor.submit(main.process_omics_data, metadata)
        fda_future = executor.submit(main.process_fda_data, metadata)

        # Get results
        try:
            pubchem_data_all, chembl_data_all = drug_future.result(timeout=600)
            gene_data_all = gene_future.result(timeout=600)
            disease_data_all = disease_future.result(timeout=600)
            clinical_trial_data_all = clinical_trial_future.result(timeout=600)
            literature_data_all = literature_future.result(timeout=600)
            omics_data_all = omics_future.result(timeout=600)
            fda_data_all = fda_future.result(timeout=600)
            
            return (pubchem_data_all, chembl_data_all, gene_data_all, disease_data_all,
                   clinical_trial_data_all, literature_data_all, omics_data_all, fda_data_all)
        except concurrent.futures.TimeoutError:
            raise TimeoutError("Data processing timed out")

def error_handler(f):
    @wraps(f)
    def wrapper(*args, **kwargs):
        try:
            return f(*args, **kwargs)
        except TimeoutError as e:
            logger.error(f"Timeout in {f.__name__}: {str(e)}")
            return jsonify({
                "error": "Request timed out",
                "details": "The operation took too long to complete. Please try a more specific query."
            }), 504
        except Exception as e:
            logger.error(f"Error in {f.__name__}: {str(e)}", exc_info=True)
            return jsonify({
                "error": "An internal error occurred",
                "details": str(e)
            }), 500
    return wrapper

@app.route("/deepseek", methods=["POST"])
@error_handler
def deepseek_query():
    data = request.get_json()
    if not data or "query" not in data:
        return jsonify({"error": "Missing query in JSON body"}), 400

    query = data["query"]
    logger.info(f"Processing query: {query}")

    try:
        # Extract metadata
        metadata = main.extract_metadata(query)
        logger.info("Metadata extracted successfully")

        # Process all data sources concurrently
        (pubchem_data_all, chembl_data_all, gene_data_all, disease_data_all,
         clinical_trial_data_all, literature_data_all, omics_data_all,
         fda_data_all) = process_data_concurrently(metadata)
        
        logger.info("All data sources processed successfully")

        # Generate response
        deepseek_response = main.generate_deepseek_response(
            query,
            pubchem_data_all,
            chembl_data_all,
            gene_data_all,
            disease_data_all,
            clinical_trial_data_all,
            literature_data_all,
            omics_data_all,
            fda_data_all
        )
        
        return jsonify({
            "success": True,
            "deepseek_response": deepseek_response,
            "timestamp": time.time()
        })
        
    except Exception as e:
        logger.error(f"Error processing query: {str(e)}", exc_info=True)
        return jsonify({
            "error": "Failed to process query",
            "details": str(e)
        }), 500

if __name__ == "__main__":
    # Use threaded=True to handle multiple requests concurrently
    app.run(host="0.0.0.0", port=5000, debug=True, threaded=True)