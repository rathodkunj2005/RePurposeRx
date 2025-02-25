import Foundation

class DeepSeekClient: ObservableObject {
    @Published var deepSeekResponse: String?
    let baseURL: URL

    init(baseURL: URL? = nil) {
        self.baseURL = baseURL ?? URL(string: "http://10.0.0.186:5000")! // Replace with your server URL
    }

    func fetchDeepSeekResponse(
        query: String,
        timeoutInterval: TimeInterval = 3000,
        completion: @escaping (Result<String, Error>) -> Void
    ) {
        let endpoint = baseURL.appendingPathComponent("deepseek")


        var request = URLRequest(url: endpoint)
        request.httpMethod = "POST"
        request.setValue("application/json", forHTTPHeaderField: "Content-Type")
        request.timeoutInterval = timeoutInterval

        let jsonBody = ["query": query]
        do {
            let bodyData = try JSONSerialization.data(withJSONObject: jsonBody)
            request.httpBody = bodyData

            let task = URLSession.shared.dataTask(with: request) { data, response, error in
                if let error = error {
                    let deepSeekError: DeepSeekError
                    if (error as NSError).code == NSURLErrorTimedOut {
                        deepSeekError = .timeout
                    } else if (error as NSError).code == NSURLErrorCannotConnectToHost {
                        deepSeekError = .serverNotReachable
                    } else {
                        deepSeekError = .networkError(error)
                    }
                    completion(.failure(deepSeekError))
                    return
                }

                guard let httpResponse = response as? HTTPURLResponse,
                      (200...299).contains(httpResponse.statusCode) else {
                    completion(.failure(DeepSeekError.invalidResponse))
                    return
                }

                guard let data = data else {
                    completion(.failure(DeepSeekError.invalidResponse))
                    return
                }

                do {
                    if let jsonResponse = try JSONSerialization.jsonObject(with: data) as? [String: Any],
                       let deepseekResponse = jsonResponse["deepseek_response"] as? String {
                        DispatchQueue.main.async {
                            self.deepSeekResponse = deepseekResponse
                            completion(.success(deepseekResponse))
                        }
                    } else {
                        completion(.failure(DeepSeekError.decodingError))
                    }
                } catch {
                    completion(.failure(DeepSeekError.decodingError))
                }
            }
            task.resume()
        } catch {
            completion(.failure(DeepSeekError.networkError(error)))
        }
    }
}

enum DeepSeekError: LocalizedError {
    case serverNotReachable
    case timeout
    case invalidResponse
    case decodingError
    case networkError(Error)

    var errorDescription: String? {
        switch self {
        case .serverNotReachable:
            return "Unable to reach the DeepSeek server. Please ensure the server is running at http://10.0.0.186:5000" // Update if needed
        case .timeout:
            return "Request timed out. Please try a more specific query."
        case .invalidResponse:
            return "Received invalid response from server"
        case .decodingError:
            return "Error processing server response"
        case .networkError(let error):
            return "Network error: \(error.localizedDescription)"
        }
    }
}
