import Vapor
import Foundation

// MARK: - Request and Response Models

/// Structure representing the incoming JSON payload.
struct DeepseekQueryRequest: Content {
    let query: String
}

/// Structure representing the outgoing JSON response.
struct DeepseekResponse: Content {
    let success: Bool
    let deepseekResponse: String
    let timestamp: Double
}

// MARK: - Concurrent Data Processing Stub

/**
 Processes all data sources concurrently based on the extracted metadata.
 
 This function simulates calling several helper methods (like processDrugData,
 processGeneData, etc.) concurrently. In this stubbed version the individual
 functions return empty dictionaries (or tuples) so that you can later expand
 them with real API functionality.
 */
func processDataConcurrently(metadata: [String: Any]) async throws -> (
    pubchemDataAll: Any,
    chemblDataAll: Any,
    geneDataAll: Any,
    diseaseDataAll: Any,
    clinicalTrialDataAll: Any,
    literatureDataAll: Any,
    omicsDataAll: Any,
    fdaDataAll: Any
) {
    // Call the stubbed methods defined in your Mega helper.
    async let drugData = Mega.processDrugData(metadata: metadata) // Returns (pubchem, chembl)
    async let geneData = Mega.processGeneData(metadata: metadata)
    async let diseaseData = Mega.processDiseaseData(metadata: metadata)
    async let clinicalTrialData = Mega.processClinicalTrialData(metadata: metadata)
    async let literatureData = Mega.processLiteratureData(metadata: metadata)
    async let omicsData = Mega.processOmicsData(metadata: metadata)
    async let fdaData = Mega.processFDAData(metadata: metadata)
    
    // Retrieve drug data tuple (pubchem, chembl)
    let (pubchemDataAll, chemblDataAll) = await drugData
    let geneDataAll = await geneData
    let diseaseDataAll = await diseaseData
    let clinicalTrialDataAll = await clinicalTrialData
    let literatureDataAll = await literatureData
    let omicsDataAll = await omicsData
    let fdaDataAll = await fdaData
    
    return (pubchemDataAll, chemblDataAll, geneDataAll, diseaseDataAll, clinicalTrialDataAll, literatureDataAll, omicsDataAll, fdaDataAll)
}

// MARK: - Routes

/**
 Registers application routes. The `/deepseek` endpoint accepts a POST
 request containing a JSON body with a "query" field. It then:
 
 1. Extracts metadata from the query.
 2. Processes all related data concurrently.
 3. Generates a DeepSeek response.
 
 If an error occurs (or a timeout happens in any task), an HTTP 500 error will be returned.
 */
func routes(_ app: Application) throws {
    app.post("deepseek") { req async throws -> DeepseekResponse in
        do {
            // Decode the incoming JSON payload.
            let data = try req.content.decode(DeepseekQueryRequest.self)
            let query = data.query
            req.logger.info("Processing query: \(query)")
            
            // Extract metadata using the Mega helper.
            let metadata = await Mega.extractMetadata(userQuery: query)
            req.logger.info("Metadata extracted successfully.")
            
            // Process all data sources concurrently.
            // This simulates the behavior of Python's process_data_concurrently.
            _ = try await processDataConcurrently(metadata: metadata)
            req.logger.info("Concurrent data processing completed.")
            
            // Generate DeepSeek response.  
            let deepseekResponse = await Mega.generateResponse(userQuery: query, metadata: metadata)
            req.logger.info("Generated DeepSeek response.")
            
            return DeepseekResponse(
                success: true,
                deepseekResponse: deepseekResponse,
                timestamp: Date().timeIntervalSince1970
            )
        } catch {
            req.logger.error("Error processing deepseek query: \(error)")
            throw Abort(.internalServerError, reason: "Failed to process query: \(error.localizedDescription)")
        }
    }
}

// MARK: - Application Configuration

public func configure(_ app: Application) throws {
    try routes(app)
}

// MARK: - Main Entry Point

@main
struct ServerApp {
    static func main() async throws {
        var env = try Environment.detect()
        try LoggingSystem.bootstrap(from: &env)
        let app = Application(env)
        defer { app.shutdown() }
        try configure(app)
        try await app.run()
    }
} 