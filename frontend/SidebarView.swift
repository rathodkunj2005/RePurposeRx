import SwiftUI

struct ChatItem: Identifiable {
    let id = UUID()
    let title: String
    let date: Date
    var previewText: String
}
struct SidebarView: View {
    @Environment(\.colorScheme) var colorScheme
    @State private var selectedOption: SidebarOptions = .home
    @State private var isSettingsActive = false
    @State private var chatHistory = ["Chat 1", "Chat 2", "Chat 3"] // Example data
    @State private var newChatName = ""
    @State private var isEditingChatName = false
    @State private var renameChatIndex: Int? = nil
    @State private var isAddingNewChat = false
    @State private var newChatTextFieldText = ""

    var body: some View {
        VStack(alignment: .leading) {
            // App Logo or Title
            HStack {
                Image(systemName: "house.fill")
                    .font(.headline)
                Text("DeepSeek App")
                    .font(.headline)
            }
            .padding(.bottom, 20)

            // Navigation Links
            ForEach(SidebarOptions.allCases, id: \.self) { option in
                SidebarButton(option: option, selectedOption: $selectedOption)
            }

            Divider()
                .padding(.vertical)

            // Chat History Section
            HStack {
                Text("Chat History")
                    .font(.subheadline)
                    .foregroundColor(.gray)
                Spacer()
                Button(action: { isAddingNewChat = true }) {
                    Image(systemName: "plus")
                        .font(.body)
                }
            }
            .padding(.bottom, 8)

            // New Chat Creation TextField - shown when isAddingNewChat is true
            if isAddingNewChat {
                HStack {
                    TextField("New Chat Name", text: $newChatTextFieldText)
                        .textFieldStyle(.roundedBorder)
                    Button("Create") {
                        if !newChatTextFieldText.isEmpty {
                            chatHistory.append(newChatTextFieldText)
                            newChatTextFieldText = ""
                            isAddingNewChat = false
                        }
                    }
                    Button("Cancel") {
                        isAddingNewChat = false
                    }
                }
                .padding(.bottom, 8)
            }

            // Chat History List
            ForEach(chatHistory.indices, id: \.self) { index in
                HStack {
                    if renameChatIndex == index {
                        TextField("Rename Chat", text: $newChatName, onCommit: {
                            if !newChatName.isEmpty {
                                chatHistory[index] = newChatName
                            }
                            renameChatIndex = nil
                            isEditingChatName = false
                        })
                        .textFieldStyle(.roundedBorder)
                    } else {
                        Button(action: {
                            // Action for selecting chat
                            print("Selected chat: \(chatHistory[index])")
                        }) {
                            HStack {
                                Image(systemName: "message.circle.fill")
                                    .font(.body)
                                Text(chatHistory[index])
                                    .font(.body)
                                    .foregroundColor(.primary)
                            }
                        }
                        .padding(.vertical, 4)
                        .frame(maxWidth: .infinity, alignment: .leading)
                    }
                    Spacer()
                    // Rename Button - Always visible but action changes based on editing state
                    Button(action: {
                        if renameChatIndex == index {
                            renameChatIndex = nil // Cancel renaming
                            isEditingChatName = false
                        } else {
                            renameChatIndex = index // Start renaming this index
                            newChatName = chatHistory[index] // Pre-fill with current name
                            isEditingChatName = true
                        }
                    }) {
                        Image(systemName: renameChatIndex == index ? "xmark.circle.fill" : "pencil")
                            .font(.body)
                    }
                    .buttonStyle(PlainButtonStyle()) // Ensure button is plain

                    Button(action: {
                        chatHistory.remove(at: index)
                    }) {
                        Image(systemName: "trash.fill")
                            .font(.body)
                    }
                    .buttonStyle(PlainButtonStyle()) // Ensure button is plain
                }
                .padding(.vertical, 2)
            }

            Spacer()

            Button("Settings") {
                isSettingsActive = true
            }
            .onTapGesture {
                selectedOption = .settings // Update selected option in sidebar
            }
            .padding(.bottom)
            .sheet(isPresented: $isSettingsActive) {
                SettingView(isSettingsActive: $isSettingsActive)
            }
        }
        .padding(20)
        .frame(maxWidth: .infinity, alignment: .leading)
        .background(Color(.systemBackground))
        // Deprecation fix for onChange
        .onChange(of: selectedOption) { newValue in
            print("Selected option changed to: \(newValue)")
            // Handle option changes here, e.g., navigate to different views
        }
    }
}
