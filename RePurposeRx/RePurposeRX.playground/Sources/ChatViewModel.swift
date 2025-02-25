import SwiftUI
import Foundation

@MainActor
class ChatViewModel: ObservableObject {
    @Published var messages: [ChatMessage] = []
    @Published var inputText = ""
    @Published var isLoading = false
    @Published var errorMessage: String?
    @Published var response: String = ""
    
    private let deepSeekClient = DeepSeekClient()
    
    func sendQuery() {
        guard !inputText.isEmpty else { return }
        
        let userMessage = ChatMessage(content: inputText, isUser: true)
        messages.append(userMessage)
        
        
        let processingMessage = ChatMessage(
            content: "Processing your query...",
            isUser: false,
            status: .processing  // Uses the MessageStatus from ChatMessage.swift
        )
        messages.append(processingMessage)
        
        isLoading = true
        
        deepSeekClient.fetchDeepSeekResponse(query: inputText) { [weak self] result in
            DispatchQueue.main.async {
                self?.messages.removeLast() // Remove processing message
                
                switch result {
                case .success(let response):
                    let agentMessage = ChatMessage(content: response, isUser: false)
                    self?.messages.append(agentMessage)
                case .failure(let error):
                    let errorContent = "Error: \(error.localizedDescription)"
                    self?.errorMessage = errorContent
                    let errorMessage = ChatMessage(
                        content: errorContent,
                        isUser: false,
                        status: .error  // Again, using the unique MessageStatus definition
                    )
                    self?.messages.append(errorMessage)
                }
                self?.inputText = ""
                self?.isLoading = false
            }
        }
    }

    /// Begins streaming a response from the API endpoint.
    func streamResponse() {
        self.response = ""
        self.isLoading = true

        guard let url = URL(string: "http://10.0.0.186:5000/deepseek") else {
            print("Invalid URL")
            self.isLoading = false
            return
        }

        var request = URLRequest(url: url)
        request.httpMethod = "POST"
        // Configure additional headers or HTTP body as needed.
        
        Task {
            do {
                // Start streaming the response using URLSession's async bytes API.
                let (bytes, _) = try await URLSession.shared.bytes(for: request)
                for try await byte in bytes {
                    // Here we're assuming that the API sends the response one byte at a time.
                    // In practice, you might receive larger chunks (Data) that you can convert to a String at once.
                    let chunk = String(decoding: [byte], as: UTF8.self)
                    self.response.append(chunk)
                }
            } catch {
                self.response.append("\n[Error: \(error.localizedDescription)]")
            }
            self.isLoading = false
        }
    }
}
