//
//  SettingView.swift
//  DeepSeekApp
//
//  Created by Kunj Rathod on 2/9/25.
//


import SwiftUI
#if os(iOS) || os(tvOS) || targetEnvironment(macCatalyst)
import UIKit
#endif


struct SettingView: View {
    @Binding var isSettingsActive: Bool

    var body: some View {
        NavigationView {
            Form {
                Section(header: Text("General")) {
                    Text("App Version 1.0")
                }
                Section(header: Text("About")) {
                    Text("RePurposeRX App is designed to...")
                }
            }
            .navigationTitle("Settings")
            .toolbar {
                ToolbarItem(placement: .navigationBarTrailing) {
                    Button("Done") {
                        isSettingsActive = false
                    }
                }
            }
        }
    }
}
