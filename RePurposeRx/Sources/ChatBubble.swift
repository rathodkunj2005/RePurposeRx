import SwiftUI

struct ChatBubble: View {
    let message: ChatMessage
    @Environment(\.colorScheme) var colorScheme
    let aiName: String // Add AI name as a parameter


    // MARK: - Constants for Styling
    private let avatarSize: CGFloat = 36
    private let avatarIconSize: CGFloat = 18
    private let avatarShadowRadius: CGFloat = 4
    private let avatarShadowYOffset: CGFloat = 2
    private let bubbleCornerRadius: CGFloat = 18
    private let bubbleHorizontalPadding: CGFloat = 14
    private let bubbleVerticalPadding: CGFloat = 10
    private let bubbleMaxWidthRatio: CGFloat = 0.75
    private let messageFontSize: CGFloat = 16
    private let aiNameFontSize: CGFloat = 14
    private let interItemSpacing: CGFloat = 12
    private let bubbleVerticalSpacing: CGFloat = 4
    private let horizontalPadding: CGFloat = 16
    private let verticalPadding: CGFloat = 4

    // Avatar System Names (for clarity and potential slight efficiency)
    private let userAvatarSystemName = "person.circle.fill"
    private let aiAvatarSystemName = "brain.head.profile"


    private var backgroundColor: Color {
        if message.isUser {
            return Color.blue.opacity(0.2)
        } else {
            return colorScheme == .dark ? Color(.systemGray5) : Color(.systemGray6)
        }
    }

    private var avatarBackgroundColor: Color {
        if message.isUser {
            return .blue
        } else {
            return .green
        }
    }

    var body: some View {
        HStack(alignment: .top, spacing: interItemSpacing) { // Single HStack for layout
            if !message.isUser {
                AvatarView(backgroundColor: avatarBackgroundColor, systemImageName: aiAvatarSystemName)
            } else {
                Spacer() // Push user bubble to the right
            }

            VStack(alignment: message.isUser ? .trailing : .leading, spacing: bubbleVerticalSpacing) {
                if !message.isUser {
                    Text(aiName) // Use configurable AI name
                        .font(.system(size: aiNameFontSize, weight: .medium))
                        .foregroundColor(.secondary)
                }
                ParsedTextView(text: message.content) // Assume ParsedTextView exists
                    .font(.system(size: messageFontSize))
                    .foregroundColor(.primary)
                    .fixedSize(horizontal: false, vertical: true)
                    .padding(.vertical, bubbleVerticalPadding)
                    .padding(.horizontal, bubbleHorizontalPadding)
                    .background(backgroundColor)
                    .clipShape(RoundedRectangle(cornerRadius: bubbleCornerRadius, style: .continuous))
            }
            .frame(maxWidth: UIScreen.main.bounds.width * bubbleMaxWidthRatio, alignment: message.isUser ? .trailing : .leading)

            if message.isUser {
                AvatarView(backgroundColor: avatarBackgroundColor, systemImageName: userAvatarSystemName)
            } else {
                Spacer() // Push AI bubble to the left if needed, though VStack already aligns left
            }
        }
        .padding(.horizontal, horizontalPadding)
        .padding(.vertical, verticalPadding)
    }
}

struct AvatarView: View {
    let backgroundColor: Color
    let systemImageName: String

    // MARK: - Constants for AvatarView (if needed to customize AvatarView further)
    private let avatarSize: CGFloat = 36 // Re-define if AvatarView styling needs to be independent
    private let avatarIconSize: CGFloat = 18
    private let avatarShadowRadius: CGFloat = 4
    private let avatarShadowYOffset: CGFloat = 2


    var body: some View {
        Circle()
            .fill(backgroundColor)
            .frame(width: avatarSize, height: avatarSize)
            .overlay(
                Image(systemName: systemImageName)
                    .font(.system(size: avatarIconSize))
                    .foregroundColor(.white)
            )
            .shadow(color: backgroundColor.opacity(0.3), radius: avatarShadowRadius, x: 0, y: avatarShadowYOffset)
    }
}


// Preview Provider for testing
struct ChatBubble_Previews: PreviewProvider {
    static var previews: some View {
        VStack(spacing: 20) {
            ChatBubble(message: ChatMessage(content: """
                    Hello! How can I help you with drug repurposing today?
                    <think>
                        This is some internal **debugging** information.
                    </think>
                    Here's more explanation.
                    """, isUser: false), aiName: "Drug Repurposing AI") // Pass AI Name
            ChatBubble(message: ChatMessage(content: "I'm looking for information about FGFR2 inhibitors for cancer treatment.", isUser: true), aiName: "Drug Repurposing AI") // Pass AI Name
        }
        .padding()
        .previewLayout(.sizeThatFits)
    }
}
