import SwiftUI

/// This enum describes the parts of a message.
enum MessageSegment {
    case normal(String)
    case think(String)
}

struct ParsedTextView: View {
    let text: String
    
    /// Parses the string into segments based on <think>...</think> tags.
    private var segments: [MessageSegment] {
        parseSegments(from: text)
    }
    
    var body: some View {
        VStack(alignment: .leading, spacing: 4) {
            ForEach(Array(segments.enumerated()), id: \.offset) { _, segment in
                switch segment {
                case .normal(let content):
                    // Use Swift's markdown conversion so that text between ** ** is rendered as bold.
                    if let attrText = try? AttributedString(markdown: content) {
                        Text(attrText)
                    } else {
                        Text(content)
                    }
                case .think(let content):
                    // Wrap the <think> content in a collapsible DisclosureGroup.
                    DisclosureGroup("Thinking (click to expand)") {
                        if let attrText = try? AttributedString(markdown: content) {
                            Text(attrText)
                        } else {
                            Text(content)
                        }
                    }
                }
            }
        }
    }
    
    /// Uses a regular expression to break up the text into normal and think segments.
    private func parseSegments(from text: String) -> [MessageSegment] {
        var segments: [MessageSegment] = []
        let pattern = "<think>([\\s\\S]*?)</think>"
        
        guard let regex = try? NSRegularExpression(pattern: pattern, options: []) else {
            return [.normal(text)]
        }
        
        let nsText = text as NSString
        var currentIndex = 0
        let matches = regex.matches(in: text, options: [], range: NSRange(location: 0, length: nsText.length))
        
        for match in matches {
            let matchRange = match.range
            
            // Add any text before the <think> block as normal.
            if matchRange.location > currentIndex {
                let normalPart = nsText.substring(with: NSRange(location: currentIndex, length: matchRange.location - currentIndex))
                segments.append(.normal(normalPart))
            }
            
            // Extract the content within <think>...</think>.
            let thinkRange = match.range(at: 1)
            let thinkPart = nsText.substring(with: thinkRange)
            segments.append(.think(thinkPart))
            
            currentIndex = matchRange.location + matchRange.length
        }
        
        // Append any remaining text as normal.
        if currentIndex < nsText.length {
            let remaining = nsText.substring(from: currentIndex)
            segments.append(.normal(remaining))
        }
        
        return segments
    }
}
