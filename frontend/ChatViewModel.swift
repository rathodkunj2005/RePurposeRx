import SwiftUI

class ChatViewModel: ObservableObject {
    @Published var messages: [ChatMessage] = []
    @Published var inputText = ""
    @Published var isLoading = false
    @Published var errorMessage: String?
    
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
}
