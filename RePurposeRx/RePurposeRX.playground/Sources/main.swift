import Foundation

// MARK: - Constants & Configuration

let PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
let CHEMBL_BASE_URL = "https://www.ebi.ac.uk/chembl/api/data/molecule/search"
let UNIPROT_BASE_URL = "https://rest.uniprot.org/uniprotkb"
let ENSEMBL_BASE_URL = "https://rest.ensembl.org"
let RXNAV_BASE_URL = "https://mor.nlm.nih.gov/RxNav"
let KEGG_BASE_URL = "https://rest.kegg.jp"
let REACTOME_BASE_URL = "https://reactome.org/ContentService"
let HPA_BASE_URL = "https://www.proteinatlas.org/api"
let PUBMED_EUTILS_BASE_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
let OMIM_API_BASE_URL = "https://api.omim.org/api"
let DISGENET_API_BASE_URL = "https://www.disgenet.org/api"
let CLINICALTRIALS_GOV_BASE_URL = "https://clinicaltrials.gov/api/query"
let EUROPE_PMC_BASE_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
let TCGA_GDC_API_BASE_URL = "https://api.gdc.cancer.gov"
let FDA_API_BASE_URL = "https://api.fda.gov/drug"

let DRUG_NAME_BLACKLIST = ["primary targets", "targets", "target", "pathway", "information"]
let PUBCHEM_PROPERTIES = "MolecularWeight,IUPACName,CanonicalSMILES"
let OMIM_API_KEY = "3FZX0uYkSVmNO4Ako4d8yg"
let DISGENET_API_KEY = "87dec4c5-14f3-43ec-8432-9b16bc1a0b26"

// MARK: - Helper Functions

/// Generic function to fetch data from a URL.
func fetchData(from url: URL) async throws -> Data {
    let (data, response) = try await URLSession.shared.data(from: url)
    guard let httpResponse = response as? HTTPURLResponse, httpResponse.statusCode == 200 else {
        throw URLError(.badServerResponse)
    }
    return data
}

/// A helper to fetch and parse JSON from a URL.
func fetchJSON(from url: URL) async throws -> Any {
    let data = try await fetchData(from: url)
    return try JSONSerialization.jsonObject(with: data, options: [])
}

// MARK: - Metadata Extraction

/// Extracts metadata from the user query using (for example) an LLM call.
/// In this stubbed version, we return a dummy dictionary.
func extractMetadata(userQuery: String) async -> [String: Any] {
    let dummyMetadata: [String: Any] = [
        "Drug Information": [],
        "Gene Information": [],
        "Target and Pathway Information": [],
        "Disease Information": [],
        "Clinical Trial Data": [],
        "Literature and Evidence": [],
        "Omics Data": [],
        "Other Important Resources": [
            "Patents": [],
            "FDA": [],
            "European Medicines Agency": []
        ],
        "General Query": false,
        "Related Keywords": []
    ]
    print("Metadata extracted successfully.")
    return dummyMetadata
}

/// Checks if the given name is likely a drug name, using a blacklist.
func isValidDrugName(_ name: String) -> Bool {
    for word in DRUG_NAME_BLACKLIST {
        if name.lowercased().contains(word) {
            return false
        }
    }
    return true
}

// MARK: - PubChem & ChEMBL Data Retrieval

/// Queries PubChem for drug information.
func fetchPubchemInfo(drugName: String, outputFormat: String = "JSON") async -> Any? {
    let outputFormatUpper = outputFormat.uppercased()
    let urlStr = "\(PUBCHEM_BASE_URL)/compound/name/\(drugName)/property/\(PUBCHEM_PROPERTIES)/\(outputFormatUpper)"
    guard let url = URL(string: urlStr) else {
        print("Invalid URL: \(urlStr)")
        return nil
    }
    do {
        let data = try await fetchData(from: url)
        if outputFormatUpper == "JSON" {
            let jsonObj = try JSONSerialization.jsonObject(with: data, options: [])
            return jsonObj
        } else if outputFormatUpper == "TXT" {
            return String(data: data, encoding: .utf8)
        } else {
            print("Unsupported output format: \(outputFormat)")
            return nil
        }
    } catch {
        print("Error fetching PubChem info for \(drugName): \(error)")
        return nil
    }
}

/// Queries ChEMBL for drug information and returns text data.
func fetchChemblInfo(drugName: String) async -> String? {
    guard var components = URLComponents(string: CHEMBL_BASE_URL) else {
        print("Invalid URL components for ChEMBL")
        return nil
    }
    components.queryItems = [
        URLQueryItem(name: "q", value: drugName),
        URLQueryItem(name: "format", value: "txt")
    ]
    guard let url = components.url else {
        print("Invalid URL for ChEMBL info")
        return nil
    }
    do {
        let data = try await fetchData(from: url)
        return String(data: data, encoding: .utf8)
    } catch {
        print("Error fetching ChEMBL info for \(drugName): \(error)")
        return nil
    }
}

// MARK: - UniProt & Ensembl Data Retrieval

/// Fetches UniProt ID for a given gene name.
func fetchUniprotId(geneName: String) async -> String? {
    guard var components = URLComponents(string: "\(UNIPROT_BASE_URL)/search") else { return nil }
    components.queryItems = [
        URLQueryItem(name: "query", value: geneName),
        URLQueryItem(name: "fields", value: "accession"),
        URLQueryItem(name: "format", value: "json")
    ]
    guard let url = components.url else { return nil }
    do {
        let data = try await fetchData(from: url)
        if let jsonObj = try JSONSerialization.jsonObject(with: data, options: []) as? [String: Any],
           let results = jsonObj["results"] as? [[String: Any]],
           let first = results.first,
           let primaryAccession = first["primaryAccession"] as? String {
            return primaryAccession
        }
        print("UniProt ID not found for \(geneName)")
        return nil
    } catch {
        print("Error fetching UniProt data for \(geneName): \(error)")
        return nil
    }
}

/// Fetch UniProt data and extract the "Function" section.
func fetchAndExtractUniProtFunction(uniprotId: String) async -> String? {
    let urlStr = "\(UNIPROT_BASE_URL)/\(uniprotId).txt"
    guard let url = URL(string: urlStr) else { return nil }
    do {
        let data = try await fetchData(from: url)
        guard let text = String(data: data, encoding: .utf8) else { return nil }
        var functionSection: [String] = []
        var capture = false
        for line in text.components(separatedBy: "\n") {
            if line.hasPrefix("CC   -!- FUNCTION:") {
                capture = true
                let trimmed = line.replacingOccurrences(of: "CC   -!- FUNCTION:", with: "").trimmingCharacters(in: .whitespaces)
                functionSection.append(trimmed)
            } else if capture && line.hasPrefix("CC       ") {
                let trimmed = line.replacingOccurrences(of: "CC       ", with: "").trimmingCharacters(in: .whitespaces)
                functionSection.append(trimmed)
            } else if capture {
                break
            }
        }
        if !functionSection.isEmpty {
            return functionSection.joined(separator: " ")
        }
        print("Function section not found in UniProt data for \(uniprotId)")
        return nil
    } catch {
        print("Error fetching UniProt function for \(uniprotId): \(error)")
        return nil
    }
}

/// Fetches UniProt protein information using a gene name.
func getUniprotProteinInfo(geneName: String = "FGFR2") async -> Any {
    guard var components = URLComponents(string: "\(UNIPROT_BASE_URL)/search") else { return "Invalid URL" }
    components.queryItems = [
        URLQueryItem(name: "query", value: "gene_name:\(geneName)"),
        URLQueryItem(name: "format", value: "json")
    ]
    guard let url = components.url else { return "Invalid URL" }
    do {
        let data = try await fetchData(from: url)
        if let jsonObj = try JSONSerialization.jsonObject(with: data, options: []) as? [String: Any],
           let results = jsonObj["results"] as? [[String: Any]],
           let proteinInfo = results.first {
            return proteinInfo
        }
        return "No protein information found for gene \(geneName)."
    } catch {
        return "Error fetching UniProt data: \(error)"
    }
}

/// Fetches Ensembl gene ID for a given gene name.
func fetchEnsemblGeneId(geneName: String) async -> String? {
    let urlStr = "\(ENSEMBL_BASE_URL)/xrefs/symbol/homo_sapiens/\(geneName)?content-type=application/json"
    guard let url = URL(string: urlStr) else { return nil }
    do {
        let data = try await fetchData(from: url)
        if let jsonArray = try JSONSerialization.jsonObject(with: data, options: []) as? [[String: Any]],
           let first = jsonArray.first,
           let id = first["id"] as? String {
            return id
        }
        print("Ensembl gene ID not found for \(geneName)")
        return nil
    } catch {
        print("Error fetching Ensembl gene ID for \(geneName): \(error)")
        return nil
    }
}

/// Fetches detailed gene information from Ensembl.
func fetchEnsemblGeneInfo(geneId: String) async -> [String: Any] {
    let urlStr = "\(ENSEMBL_BASE_URL)/lookup/id/\(geneId)?expand=1&content-type=application/json"
    guard let url = URL(string: urlStr) else { return [:] }
    do {
        let data = try await fetchData(from: url)
        if let jsonObj = try JSONSerialization.jsonObject(with: data, options: []) as? [String: Any] {
            var info: [String: Any] = [:]
            info["Gene name"] = jsonObj["display_name"]
            info["Description"] = jsonObj["description"]
            info["Organism"] = jsonObj["species"]
            info["Biotype"] = jsonObj["biotype"]
            if let seqRegionName = jsonObj["seq_region_name"], let start = jsonObj["start"], let end = jsonObj["end"] {
                info["Location"] = "\(seqRegionName):\(start)-\(end)"
            }
            info["Strand"] = jsonObj["strand"]
            info["Ensembl Link"] = "https://www.ensembl.org/Homo_sapiens/Gene/Summary?g=\(geneId)"
            return info
        }
        print("No Ensembl gene details found for \(geneId)")
        return [:]
    } catch {
        print("Error fetching Ensembl gene info for \(geneId): \(error)")
        return [:]
    }
}

// MARK: - RxNav Data Retrieval (Stub)

/// Fetches drug information from RxNav.
func fetchDrugInfo(drugName: String) async -> String? {
    let urlStr = "\(RXNAV_BASE_URL)/search?searchBy=String&searchTerm=\(drugName)"
    guard let url = URL(string: urlStr) else { return nil }
    do {
        _ = try await fetchData(from: url)
        // In a real implementation, parse the response.
        return urlStr
    } catch {
        print("Failed to fetch drug information for \(drugName), error: \(error)")
        return nil
    }
}

// MARK: - Gene Ontology Data Retrieval

/// Fetch Gene Ontology terms associated with a gene symbol.
func getGeneOntologyInfo(geneSymbol: String = "FGFR2") async -> Any {
    let baseUrl = "https://www.ebi.ac.uk/QuickGO/services/annotation/search"
    guard var components = URLComponents(string: baseUrl) else { return "Invalid URL" }
    components.queryItems = [
        URLQueryItem(name: "geneSymbol", value: geneSymbol),
        URLQueryItem(name: "format", value: "json")
    ]
    guard let url = components.url else { return "Invalid URL" }
    do {
        let data = try await fetchData(from: url)
        if let jsonObj = try JSONSerialization.jsonObject(with: data, options: []) as? [String: Any],
           let results = jsonObj["results"] {
            return results
        }
        return "No GO terms found for gene \(geneSymbol)."
    } catch {
        return "Error fetching Gene Ontology data: \(error)"
    }
}

// MARK: - KEGG Pathway Data Retrieval

/// Fetch KEGG pathways associated with a gene.
func getKeggPathwaysForGene(geneName: String = "FGFR2", organism: String = "hsa") async -> [String] {
    let linkUrl = "\(KEGG_BASE_URL)/link/pathway"
    let urlStr = "\(linkUrl)/\(organism):\(geneName)"
    guard let url = URL(string: urlStr) else { return [] }
    do {
        let data = try await fetchData(from: url)
        guard let text = String(data: data, encoding: .utf8) else { return [] }
        let lines = text.split(separator: "\n")
        var pathwayList: [String] = []
        for line in lines {
            let parts = line.split(separator: "\t")
            if parts.count >= 2 {
                let pathwayIdFull = parts[1]
                let pathwayId = pathwayIdFull.replacingOccurrences(of: "path:", with: "")
                let pathwayInfo = await getKeggPathwayInfo(pathwayId: pathwayId)
                pathwayList.append(pathwayInfo)
            }
        }
        return pathwayList
    } catch {
        print("Error fetching KEGG pathways for gene \(geneName): \(error)")
        return []
    }
}

/// Fetch detailed information for a specific KEGG pathway.
func getKeggPathwayInfo(pathwayId: String) async -> String {
    let urlStr = "\(KEGG_BASE_URL)/get/\(pathwayId)"
    guard let url = URL(string: urlStr) else { return "Invalid pathway URL" }
    do {
        let data = try await fetchData(from: url)
        return String(data: data, encoding: .utf8) ?? "No data"
    } catch {
        return "Error fetching KEGG pathway info for \(pathwayId): \(error)"
    }
}

// MARK: - Reactome Pathway Data Retrieval

/// Fetch Reactome pathways associated with a gene.
func getReactomePathwaysForGene(geneName: String = "FGFR2") async -> Any {
    let urlStr = "\(REACTOME_BASE_URL)/search/query"
    guard var components = URLComponents(string: urlStr) else { return "Invalid URL" }
    components.queryItems = [
        URLQueryItem(name: "query", value: geneName),
        URLQueryItem(name: "species", value: "Homo sapiens"),
        URLQueryItem(name: "facet", value: "pathway")
    ]
    guard let url = components.url else { return "Invalid URL" }
    do {
        let data = try await fetchData(from: url)
        if let jsonObj = try JSONSerialization.jsonObject(with: data, options: []) as? [String: Any],
           let results = jsonObj["results"] {
            return results
        }
        return "No Reactome pathways found for \(geneName)."
    } catch {
        return "Error fetching Reactome data: \(error)"
    }
}

// MARK: - Human Protein Atlas Data Retrieval

/// Fetch gene expression data from Human Protein Atlas.
func getHumanProteinAtlasExpression(geneName: String = "FGFR2") async -> Any {
    let urlStr = "\(HPA_BASE_URL)/search_download.php"
    guard var components = URLComponents(string: urlStr) else { return "Invalid URL" }
    components.queryItems = [
        URLQueryItem(name: "query", value: geneName),
        URLQueryItem(name: "format", value: "json"),
        URLQueryItem(name: "columns", value: "g"),
        URLQueryItem(name: "compress", value: "no")
    ]
    guard let url = components.url else { return "Invalid URL" }
    do {
        let data = try await fetchData(from: url)
        if let jsonObj = try JSONSerialization.jsonObject(with: data, options: []) as? [[String: Any]],
           !jsonObj.isEmpty {
            return jsonObj[0]
        }
        return "No expression data found for gene \(geneName)."
    } catch {
        return "Error fetching Human Protein Atlas data: \(error)"
    }
}

// MARK: - MeSH Term Retrieval via PubMed E-utilities

/// Searches PubMed for MeSH terms related to the specified term.
func searchPubmedMesh(term: String = "Bladder Cancer") async -> [String] {
    guard var components = URLComponents(string: "\(PUBMED_EUTILS_BASE_URL)/esearch.fcgi") else { return [] }
    components.queryItems = [
        URLQueryItem(name: "db", value: "mesh"),
        URLQueryItem(name: "term", value: term),
        URLQueryItem(name: "retmode", value: "xml"),
        URLQueryItem(name: "retmax", value: "10")
    ]
    guard let url = components.url else { return [] }
    do {
        let data = try await fetchData(from: url)
        guard let xmlText = String(data: data, encoding: .utf8) else { return [] }
        var meshIds: [String] = []
        let pattern = "<Id>(.*?)</Id>"
        if let regex = try? NSRegularExpression(pattern: pattern, options: []) {
            let nsRange = NSRange(xmlText.startIndex..<xmlText.endIndex, in: xmlText)
            let matches = regex.matches(in: xmlText, options: [], range: nsRange)
            for match in matches {
                if let range = Range(match.range(at: 1), in: xmlText) {
                    meshIds.append(String(xmlText[range]))
                }
            }
        }
        return meshIds
    } catch {
        print("Error searching PubMed MeSH: \(error)")
        return []
    }
}

// MARK: - OMIM Data Retrieval

/// Fetches OMIM entries related to a gene symbol.
func getOmimGeneDisorders(geneSymbol: String = "FGFR2") async -> Any {
    guard var components = URLComponents(string: "\(OMIM_API_BASE_URL)/entry/search") else { return "Invalid URL" }
    components.queryItems = [
        URLQueryItem(name: "search", value: "gene:\(geneSymbol)"),
        URLQueryItem(name: "format", value: "json"),
        URLQueryItem(name: "apiKey", value: OMIM_API_KEY)
    ]
    guard let url = components.url else { return "Invalid URL" }
    do {
        let data = try await fetchData(from: url)
        if let jsonObj = try JSONSerialization.jsonObject(with: data, options: []) as? [String: Any],
           let omim = jsonObj["omim"] as? [String: Any],
           let searchResponse = omim["searchResponse"] as? [String: Any],
           let entryList = searchResponse["entryList"] {
            return entryList
        }
        return "No OMIM entries found for gene \(geneSymbol)."
    } catch {
        print("Error fetching OMIM data: \(error)")
        return "Error fetching OMIM data: \(error)"
    }
}

// MARK: - DisGeNET Data Retrieval

/// Fetches gene-disease associations from DisGeNET.
func getDisgenetGeneDiseaseAssociations(geneNcbiId: String = "351", pageNumber: Int = 0) async -> Any {
    guard var components = URLComponents(string: "\(DISGENET_API_BASE_URL)/api/v1") else { return "Invalid URL" }
    components.queryItems = [
        URLQueryItem(name: "gene_ncbi_id", value: geneNcbiId),
        URLQueryItem(name: "page", value: "\(pageNumber)")
    ]
    guard let url = components.url else { return "Invalid URL" }
    var request = URLRequest(url: url)
    request.addValue(DISGENET_API_KEY, forHTTPHeaderField: "api-key")
    do {
        let (data, _) = try await URLSession.shared.data(for: request)
        let jsonObj = try JSONSerialization.jsonObject(with: data, options: [])
        return jsonObj
    } catch {
        print("Error fetching DisGeNET data: \(error)")
        return "Error fetching DisGeNET data: \(error)"
    }
}

// MARK: - Data Processing Functions

/// Processes drug data concurrently.
func processDrugData(metadata: [String: Any]) async -> (pubchemDataAll: [String: Any], chemblDataAll: [String: Any]) {
    var pubchemDataAll: [String: Any] = [:]
    var chemblDataAll: [String: Any] = [:]
    if let drugInfo = metadata["Drug Information"] as? [String], !drugInfo.isEmpty {
        let drugs = drugInfo.filter { isValidDrugName($0) }
        await withTaskGroup(of: (String, Any?).self) { group in
            for drug in drugs {
                group.addTask {
                    let result = await fetchPubchemInfo(drugName: drug)
                    return (drug, result)
                }
            }
            for await (drug, result) in group {
                pubchemDataAll[drug] = result ?? "No data"
            }
        }
        await withTaskGroup(of: (String, Any?).self) { group in
            for drug in drugs {
                group.addTask {
                    let result = await fetchChemblInfo(drugName: drug)
                    return (drug, result)
                }
            }
            for await (drug, result) in group {
                chemblDataAll[drug] = result ?? "No data"
            }
        }
    } else {
        print("No valid drug names found in metadata.")
    }
    return (pubchemDataAll, chemblDataAll)
}

/// Processes gene data concurrently.
func processGeneData(metadata: [String: Any]) async -> [String: Any] {
    var geneDataAll: [String: Any] = [:]
    if let genes = metadata["Gene Information"] as? [String], !genes.isEmpty {
        await withTaskGroup(of: (String, [String: Any]).self) { group in
            for gene in genes {
                group.addTask {
                    var geneData: [String: Any] = [:]
                    if let uniprotId = await fetchUniprotId(geneName: gene) {
                        geneData["UniProt ID"] = uniprotId
                        geneData["UniProt Function"] = await fetchAndExtractUniProtFunction(uniprotId: uniprotId)
                    }
                    if let ensemblId = await fetchEnsemblGeneId(geneName: gene) {
                        geneData["Ensembl ID"] = ensemblId
                        geneData["Ensembl Info"] = await fetchEnsemblGeneInfo(geneId: ensemblId)
                    }
                    geneData["UniProt Protein Data"] = await getUniprotProteinInfo(geneName: gene)
                    geneData["Gene Ontology Data"] = await getGeneOntologyInfo(geneSymbol: gene)
                    geneData["KEGG Pathways"] = await getKeggPathwaysForGene(geneName: gene)
                    geneData["Reactome Pathways"] = await getReactomePathwaysForGene(geneName: gene)
                    geneData["Human Protein Atlas Expression"] = await getHumanProteinAtlasExpression(geneName: gene)
                    return (gene, geneData)
                }
            }
            for await (gene, geneData) in group {
                geneDataAll[gene] = geneData
            }
        }
    } else {
        print("No valid gene names found in metadata.")
    }
    return geneDataAll
}

/// Processes disease data concurrently.
func processDiseaseData(metadata: [String: Any]) async -> [String: Any] {
    var diseaseDataAll: [String: Any] = [:]
    if let diseases = metadata["Disease Information"] as? [String], !diseases.isEmpty {
        await withTaskGroup(of: (String, [String: Any]).self) { group in
            for disease in diseases {
                group.addTask {
                    var data: [String: Any] = [:]
                    data["MeSH Terms"] = await searchPubmedMesh(term: disease)
                    data["OMIM Entries"] = await getOmimGeneDisorders(geneSymbol: disease)
                    data["DisGeNET Associations"] = await getDisgenetGeneDiseaseAssociations(geneNcbiId: disease)
                    return (disease, data)
                }
            }
            for await (disease, data) in group {
                diseaseDataAll[disease] = data
            }
        }
    } else {
        print("No valid disease names found in metadata.")
    }
    return diseaseDataAll
}

/// Processes clinical trial data concurrently.
func processClinicalTrialData(metadata: [String: Any]) async -> [String: Any] {
    var clinicalTrialDataAll: [String: Any] = [:]
    if let queries = metadata["Clinical Trial Data"] as? [String], !queries.isEmpty {
        await withTaskGroup(of: (String, Any).self) { group in
            for query in queries {
                group.addTask {
                    // Simulate clinical trials retrieval.
                    let dummy = "Clinical trial data for \(query)"
                    return (query, dummy)
                }
            }
            for await (query, data) in group {
                clinicalTrialDataAll[query] = data
            }
        }
    } else {
        print("No valid clinical trial queries found in metadata.")
    }
    return clinicalTrialDataAll
}

/// Processes literature data concurrently.
func processLiteratureData(metadata: [String: Any]) async -> [String: Any] {
    var literatureDataAll: [String: Any] = [:]
    if let queries = metadata["Literature and Evidence"] as? [String], !queries.isEmpty {
        await withTaskGroup(of: (String, Any).self) { group in
            for query in queries {
                group.addTask {
                    // Simulate literature search.
                    let dummy = "Literature data for \(query)"
                    return (query, dummy)
                }
            }
            for await (query, data) in group {
                literatureDataAll[query] = data
            }
        }
    } else {
        print("No valid literature queries found in metadata.")
    }
    return literatureDataAll
}

/// Processes omics data.
func processOmicsData(metadata: [String: Any]) async -> [String: Any] {
    var omicsDataAll: [String: Any] = [:]
    if let diseases = metadata["Omics Data"] as? [String], !diseases.isEmpty {
        for disease in diseases {
            // Simulate TCGA project retrieval.
            let dummy = "TCGA project data for \(disease)"
            omicsDataAll[disease] = dummy
        }
    } else {
        print("No valid omics data found in metadata.")
    }
    return omicsDataAll
}

/// Processes FDA approvals data.
func processFdaData(metadata: [String: Any]) async -> [String: Any] {
    var fdaDataAll: [String: Any] = [:]
    if let resources = metadata["Other Important Resources"] as? [String: Any],
       let fdaTerms = resources["FDA"] as? [String], !fdaTerms.isEmpty {
        for term in fdaTerms {
            // Simulate FDA data retrieval.
            let dummy = "FDA approvals data for \(term)"
            fdaDataAll[term] = dummy
        }
    } else {
        print("No valid FDA terms found in metadata.")
    }
    return fdaDataAll
}

// MARK: - PDF Downloading & PubMed Literature Processing

/// Extracts the free full-text PDF URL from a PubMed article summary.
func getPdfUrl(article: [String: Any]) -> String? {
    if let articleIds = article["articleids"] as? [[String: Any]] {
        for idEntry in articleIds {
            if let idtype = idEntry["idtype"] as? String, idtype.lowercased() == "pmc",
               let value = idEntry["value"] as? String {
                var pmc_id = value
                if !pmc_id.hasPrefix("PMC") {
                    pmc_id = "PMC" + pmc_id
                }
                return "https://www.ncbi.nlm.nih.gov/pmc/articles/\(pmc_id)/pdf/"
            }
        }
    }
    return nil
}

/// Downloads a PDF from the given URL and saves it locally.
func downloadPdf(pdfUrl: String, savePath: String) async -> String? {
    guard let url = URL(string: pdfUrl) else { return nil }
    do {
        let data = try await fetchData(from: url)
        let fileURL = URL(fileURLWithPath: savePath)
        try data.write(to: fileURL)
        print("Downloaded PDF from \(pdfUrl) to \(savePath)")
        return savePath
    } catch {
        print("Error downloading PDF from \(pdfUrl): \(error)")
        return nil
    }
}

/// Searches PubMed for articles related to the query using E-utilities.
func searchPubmed(query: String, retmax: Int = 10) async -> [[String: Any]] {
    // Step 1: ESearch
    guard var esearchComponents = URLComponents(string: "\(PUBMED_EUTILS_BASE_URL)/esearch.fcgi") else { return [] }
    esearchComponents.queryItems = [
        URLQueryItem(name: "db", value: "pubmed"),
        URLQueryItem(name: "term", value: query),
        URLQueryItem(name: "retmode", value: "json"),
        URLQueryItem(name: "retmax", value: "\(retmax)")
    ]
    guard let esearchUrl = esearchComponents.url else { return [] }
    do {
        let data = try await fetchData(from: esearchUrl)
        if let jsonObj = try JSONSerialization.jsonObject(with: data, options: []) as? [String: Any],
           let esearchResult = jsonObj["esearchresult"] as? [String: Any],
           let idList = esearchResult["idlist"] as? [String],
           !idList.isEmpty {
            
            // Step 2: ESummary
            guard var esummaryComponents = URLComponents(string: "\(PUBMED_EUTILS_BASE_URL)/esummary.fcgi") else { return [] }
            esummaryComponents.queryItems = [
                URLQueryItem(name: "db", value: "pubmed"),
                URLQueryItem(name: "id", value: idList.joined(separator: ",")),
                URLQueryItem(name: "retmode", value: "json")
            ]
            guard let esummaryUrl = esummaryComponents.url else { return [] }
            let summaryData = try await fetchData(from: esummaryUrl)
            if let summaryJson = try JSONSerialization.jsonObject(with: summaryData, options: []) as? [String: Any],
               let result = summaryJson["result"] as? [String: Any],
               let uids = result["uids"] as? [String] {
                var articles: [[String: Any]] = []
                for uid in uids {
                    if let article = result[uid] as? [String: Any] {
                        articles.append(article)
                    }
                }
                return articles
            }
        }
        return []
    } catch {
        print("Error during PubMed search: \(error)")
        return []
    }
}

/// Processes PubMed literature data and returns either a text summary or downloaded PDFs.
func processLiteratureData(userQuery: String, retmax: Int = 10, mode: String = "TXT") async -> Any {
    let articles = await searchPubmed(query: userQuery, retmax: retmax)
    if articles.isEmpty {
        print("No PubMed articles found.")
        return [:]
    }
    if mode.uppercased() == "TXT" {
        let summaries = articles.compactMap { $0["summary"] as? String }
        return summaries.joined(separator: "\n\n")
    } else if mode.uppercased() == "PDF" {
        let fileManager = FileManager.default
        let saveDir = "pubmed_pdfs"
        if !fileManager.fileExists(atPath: saveDir) {
            try? fileManager.createDirectory(atPath: saveDir, withIntermediateDirectories: true, attributes: nil)
        }
        var pdfLinks: [String: String] = [:]
        for article in articles {
            if let pdfUrl = getPdfUrl(article: article) {
                let pmid = article["uid"] as? String ?? "unknown"
                pdfLinks[pmid] = pdfUrl
            }
        }
        var downloadedFiles: [String: String] = [:]
        for (pmid, pdfUrl) in pdfLinks {
            let savePath = "\(saveDir)/\(pmid).pdf"
            if let downloadedPath = await downloadPdf(pdfUrl: pdfUrl, savePath: savePath) {
                downloadedFiles[pmid] = downloadedPath
            } else {
                downloadedFiles[pmid] = "Failed to download PDF."
            }
        }
        return downloadedFiles
    } else {
        print("Invalid mode specified: \(mode)")
        return [:]
    }
}

/// Searches PubMed and downloads available PDFs.
func downloadPubmedPdfs(query: String, retmax: Int = 10, saveDir: String = "pubmed_pdfs") async -> [String: String] {
    let articles = await searchPubmed(query: query, retmax: retmax)
    if articles.isEmpty {
        print("No PubMed articles found for query.")
        return [:]
    }
    let fileManager = FileManager.default
    if !fileManager.fileExists(atPath: saveDir) {
        try? fileManager.createDirectory(atPath: saveDir, withIntermediateDirectories: true, attributes: nil)
    }
    var downloadedFiles: [String: String] = [:]
    for article in articles {
        let articleId = article["uid"] as? String ?? "article_\(articles.firstIndex { ($0 as AnyObject) === (article as AnyObject) } ?? 0)"
        if let pdfUrl = getPdfUrl(article: article) {
            let savePath = "\(saveDir)/\(articleId).pdf"
            if let downloadedPath = await downloadPdf(pdfUrl: pdfUrl, savePath: savePath) {
                downloadedFiles[articleId] = downloadedPath
            } else {
                downloadedFiles[articleId] = "Failed to download PDF."
            }
        } else {
            downloadedFiles[articleId] = "No PDF available."
        }
    }
    return downloadedFiles
}

// MARK: - Response Generation

/// Generates a response based on the type of query and appends all processed data.
func generateResponse(userQuery: String, metadata: [String: Any]) async -> String {
    if let isGeneral = metadata["General Query"] as? Bool, isGeneral {
        return "General response for query: \(userQuery)"
    } else {
        // In a complete implementation, you would merge data from multiple sources.
        // Here we simulate by appending available metadata.
        let pubchemData = metadata["PubChem Data"] as? [String: Any] ?? [:]
        let chemblData = metadata["ChEMBL Data"] as? [String: Any] ?? [:]
        let geneData = metadata["Gene Information"] as? [String: Any] ?? [:]
        let diseaseData = metadata["Disease Information"] as? [String: Any] ?? [:]
        let clinicalTrialData = metadata["Clinical Trial Data"] as? [String: Any] ?? [:]
        let literatureData = metadata["Literature and Evidence"] as? [String: Any] ?? [:]
        let omicsData = metadata["Omics Data"] as? [String: Any] ?? [:]
        let fdaData = metadata["FDA"] as? [String: Any] ?? [:]
        
        var response = userQuery
        response += "\n\nPubChem Data:\n\(pubchemData)"
        response += "\n\nChEMBL Data:\n\(chemblData)"
        response += "\n\nGene/Protein Data:\n\(geneData)"
        response += "\n\nDisease Data:\n\(diseaseData)"
        response += "\n\nClinical Trial Data:\n\(clinicalTrialData)"
        response += "\n\nLiterature Data:\n\(literatureData)"
        response += "\n\nOmics Data:\n\(omicsData)"
        response += "\n\nFDA Data:\n\(fdaData)"
        response += "\n\nYou are a drug repurposing scientist. Please analyze the above data and provide a comprehensive answer."
        return response
    }
}

// MARK: - Main Script Execution

@main
struct MainScript {
    static func main() async {
        // For testing purposes, the user query is set to an empty string—you can provide a default value.
        let userQuery = ""
        
        // Extract Metadata
        let metadata = await extractMetadata(userQuery: userQuery)
        print("\nExtracted Metadata:")
        if let metaDataJSON = try? JSONSerialization.data(withJSONObject: metadata, options: .prettyPrinted),
           let metaDataString = String(data: metaDataJSON, encoding: .utf8) {
            print(metaDataString)
        }
        
        // Process Drug Data
        let (pubchemData, chemblData) = await processDrugData(metadata: metadata)
        print("\nPubChem Data Summary:")
        print(pubchemData)
        print("\nChEMBL Data Summary:")
        print(chemblData)
        
        // Process Gene Data
        let geneData = await processGeneData(metadata: metadata)
        print("\nGene Data Summary:")
        print(geneData)
        
        // Process Disease Data
        let diseaseData = await processDiseaseData(metadata: metadata)
        print("\nDisease Data Summary:")
        print(diseaseData)
        
        // Process Clinical Trial Data
        let clinicalTrialData = await processClinicalTrialData(metadata: metadata)
        print("\nClinical Trial Data Summary:")
        print(clinicalTrialData)
        
        // Process Literature Data (from metadata)
        let literatureData = await processLiteratureData(metadata: metadata)
        print("\nLiterature Data Summary:")
        print(literatureData)
        
        // Process Omics Data (e.g., TCGA projects)
        let omicsData = await processOmicsData(metadata: metadata)
        print("\nOmics Data Summary (TCGA Projects):")
        print(omicsData)
        
        // Process FDA Data
        let fdaData = await processFdaData(metadata: metadata)
        print("\nFDA Drug Approvals Data Summary:")
        print(fdaData)
        
        // Download PubMed PDFs (using the PubMed literature search)
        let downloadedPdfs = await downloadPubmedPdfs(query: userQuery, retmax: 10, saveDir: "pubmed_pdfs")
        print("\nDownloaded PubMed PDFs:")
        if let pdfsJSON = try? JSONSerialization.data(withJSONObject: downloadedPdfs, options: .prettyPrinted),
           let pdfsString = String(data: pdfsJSON, encoding: .utf8) {
            print(pdfsString)
        }
        
        // Generate DeepSeek Response
        let deepseekResponse = await generateResponse(userQuery: userQuery, metadata: metadata)
        print("\nDeepSeek Response:")
        print(deepseekResponse)
        print("\n--- End of Script ---")
    }
} 