import SwiftUI

struct ChatView: View {
    @StateObject private var viewModel = ChatViewModel()

    var body: some View {
        VStack(alignment: .leading, spacing: 16) {
            // Display the streaming response using your ParsedTextView for markdown rendering.
            ParsedTextView(text: viewModel.response)
                .animation(.default, value: viewModel.response)
                .padding()
            
            // Show a progress view while streaming is in progress.
            if viewModel.isLoading {
                ProgressView("Streaming response...")
                    .padding()
            } else {
                Button("Stream Response") {
                    viewModel.streamResponse()
                }
                .padding()
            }
        }
        .navigationTitle("Chat Streaming")
    }
}

struct ChatView_Previews: PreviewProvider {
    static var previews: some View {
        ChatView()
    }
} 