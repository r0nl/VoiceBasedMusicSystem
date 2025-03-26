# Voice-Controlled Music Player Using HMM in C++

This project is a C++ console application that simulates a voice-controlled music player. The user can navigate and play music tracks organized in playlists using a Hidden Markov Model (HMM)-based voice recognition system. The application uses the Windows Multimedia API for audio input, performs basic HMM-based speech processing, and allows users to select playlists and songs based on spoken commands.

## Features
- Voice-controlled selection of playlists and songs
- Directory-based organization of playlists and tracks
- Basic Hidden Markov Model (HMM) implementation for speech digit recognition
- Integration with Windows Media Player for playback

## Project Structure
- `src\`: Directory for Main source file.
- `src\src.cpp`: Main source file containing the console application code and for handling HMM model generation and loading.
- `music\`: Directory structure for playlists and songs. Ensure `music` is in the same directory as `src` folder.
- `data\`: Directory structure for saved recordings. Ensure `data` is in the same directory as `src` folder.

## Requirements - Refer `requirement.txt` file.

## Setup Instructions

1. **Extract the Zip Folder**
   
2. **Open terminal in `SPEECH` folder**

3. **Compile the Program**
   - Using **Visual Studio**:
     - Open the project in Visual Studio and build the solution.
   - Using **GCC with MinGW**:
     ```
     g++ src.cpp -o music_player.exe -lwinmm
     ```

4. **Run the Program**
   ```
   ./music_player.exe
   ```

## Usage Instructions

1. Start the program. You will be prompted to either train the model or test it.
   - Training will initialize the HMM model with existing `.txt` data files.
   - Testing will prompt you to navigate playlists and songs using voice commands directly.

2. **Playlist and Song Selection**
   - The application will display available playlists. Press `Enter` to start recording.
   - The system records a short audio sample to interpret your voice command using the HMM model.
   - After a playlist is selected, the program will list songs within that playlist. Use voice commands again to select a song.

3. **Playback Control**
   - Songs will be played through Windows Media Player.
   - After playback, you can choose to continue or exit the program.

## Troubleshooting
- **Audio Recording Issues**: Ensure that your microphone is set up correctly in Windows settings.
- **File Not Found Errors**: Verify the directory structure and that all paths are correctly set relative to `src.cpp`.
- **Compilation Errors**: Make sure all required libraries are installed, and that `winmm.lib` is linked correctly.

## Additional Information
This project is designed as a single-file C++ program, but due to its complexity, separating HMM model handling into another file (`hmm_model.cpp`) may improve readability and modularity.

## Contributors
- Project Author: Rohan Dayal
- Roll No - 244101039
```
