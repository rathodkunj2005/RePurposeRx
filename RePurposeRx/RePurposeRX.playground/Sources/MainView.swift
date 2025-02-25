import SwiftUI

struct MainView: View {
    var body: some View {
        VStack {
            // Display the logo at the top of the screen
            LogoView()
            
            // Other UI components
            Text("Welcome to My App!")
                .font(.title)
                .padding()
        }
        .padding()
    }
}

struct MainView_Previews: PreviewProvider {
    static var previews: some View {
        MainView()
    }
} 