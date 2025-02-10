import Foundation

enum MessageStatus {
    case normal
    case processing
    case error
}

struct ChatMessage: Identifiable {
    let id = UUID()
    let content: String
    let isUser: Bool
    var status: MessageStatus = .normal

    
    init(content: String, isUser: Bool, status: MessageStatus = .normal) {
        self.content = content
        self.isUser = isUser
        self.status = status
    }
}
