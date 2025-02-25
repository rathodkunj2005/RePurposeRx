// swift-tools-version: 5.9
import PackageDescription

let package = Package(
    name: "DeepSeekApp",
    platforms: [
        .iOS(.v16),
        .macOS(.v13)
    ],
    products: [
        .iOSApplication(
            name: "DeepSeekApp",
            targets: ["DeepSeekApp"],
            bundleIdentifier: "com.yourname.DeepSeekApp",
            teamIdentifier: "YOUR_TEAM_ID",
            displayVersion: "1.0",
            bundleVersion: "1",
            appIcon: .placeholder(icon: .default),
            accentColor: .presetColor(.blue),
            supportedDeviceFamilies: [
                .pad,
                .phone
            ],
            supportedInterfaceOrientations: [
                .portrait,
                .landscapeRight,
                .landscapeLeft,
                .portraitUpsideDown(.when(deviceFamilies: [.pad]))
            ]
        )
    ],
    targets: [
        .executableTarget(
            name: "DeepSeekApp",
            path: "Sources"
        )
    ]
) 