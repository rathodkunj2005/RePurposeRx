import Foundation
import SwiftUI

struct ContentView: View {
    @StateObject private var viewModel = ChatViewModel()
    @State private var showSidebar = false

    var body: some View {
        NavigationView {
            ZStack {
                // Sidebar
                if showSidebar {
                    SidebarView()
                        .frame(width: 280)
                        // Removed transition as a potential fix for type-checking issue
                        //.transition(.move(edge: .leading))
                }
                
                // Main Chat View
                VStack(spacing: 0) {
                    // **Scrollable Chat Messages**
                    ScrollViewReader { proxy in
                        ScrollView {
                            LazyVStack(spacing: 12) {
                                ForEach(viewModel.messages) { message in
                                    // Assign ChatBubble to a local variable to simplify expression
                                    let bubble = ChatBubble(message: message, aiName: "Drug Repurposing AI") // Added aiName here
                                    bubble
                                        .id(message.id) // Each message needs a unique ID
                                }
                            }
                            .padding()
                        }
                        .onAppear {
                            if let lastMessage = viewModel.messages.last {
                                proxy.scrollTo(lastMessage.id, anchor: .bottom)
                            }
                        }
                        .onChange(of: viewModel.messages.count) { _ in
                            if let lastMessage = viewModel.messages.last {
                                withAnimation {
                                    proxy.scrollTo(lastMessage.id, anchor: .bottom)
                                }
                            }
                        }
                    }
                    
                    // **Error Message**
                    if let errorMessage = viewModel.errorMessage {
                        Text("Error: \(errorMessage)")
                            .foregroundColor(.red)
                            .padding()
                    }
                    
                    // **Input Bar**
                    VStack(spacing: 0) {
                        Divider()
                        HStack(spacing: 12) {
                            TextField("Message Drug Repurposing...", text: $viewModel.inputText)
                                .textFieldStyle(PlainTextFieldStyle())
                                .padding(.horizontal, 12)
                                .padding(.vertical, 8)
                                .background(Color(.systemGray6))
                                .cornerRadius(8)
                                .disabled(viewModel.isLoading)
                            
                            Button(action: viewModel.sendQuery) {
                                Image(systemName: "paperplane.fill")
                                    .foregroundColor(.white)
                                    .frame(width: 32, height: 32)
                                    .background(
                                        viewModel.inputText.isEmpty || viewModel.isLoading
                                            ? Color.blue.opacity(0.5)
                                            : Color.blue
                                    )
                                    .cornerRadius(16)
                            }
                            .disabled(viewModel.inputText.isEmpty || viewModel.isLoading)
                        }
                        .padding(.horizontal)
                        .padding(.vertical, 12)
                        .background(Color(.systemBackground))
                    }
                }
                .navigationTitle("Drug Repurposing Assistant")
                .navigationBarTitleDisplayMode(.inline)
                .toolbar {
                    ToolbarItem(placement: .navigationBarLeading) {
                        Button(action: { withAnimation { showSidebar.toggle() } }) {
                            Image(systemName: "line.3.horizontal")
                        }
                    }
                }
            }
        }
    }
}

enum SidebarOptions: String, CaseIterable {
    case home = "Home"
    case history = "History"
    case settings = "Settings"

    var imageName: String {
        switch self {
        case .home:
            return "house.fill"
        case .history:
            return "clock.fill"
        case .settings:
            return "gearshape.fill"
        }
    }
}

struct SidebarButton: View {
    var option: SidebarOptions
    @Binding var selectedOption: SidebarOptions
    
    var body: some View {
        Button(action: {
            selectedOption = option
            print("Sidebar option selected: \(option)")
            // Action for each sidebar option (e.g., navigate to different views)
        }) {
            HStack {
                Image(systemName: option.imageName)
                    .font(.title2)
                Text(option.rawValue)
                    .font(.body)
                    .foregroundColor(.primary)
            }
            .padding(.vertical, 12)
            .frame(maxWidth: .infinity, alignment: .leading)
        }
        .buttonStyle(PlainButtonStyle())
        .background(selectedOption == option ? Color.blue.opacity(0.2) : Color.clear) // Highlight selected option
        .cornerRadius(8)
    }
}
