# Noise-suppression-for-speech-signals

## What's in this folder?

This folder contains all the codes implementing the Wiener filter method and the Spectral subtraction method
for speech enhancement. There are separate folders containing the codes for each method. The dataset for 
evaluation is in Test-Data folder. Make sure the test data is in the same directory as the code when running the algorithms.

The results on the Test-Data have been uploaded in the Results folder for each algorithm. 

## How to run the algorithm on custom data?

Navigate to the python notebook for each algorithm. There is a denoise function that takes in 3 arguments:
the clean file, the noisy file and the noisy file SNR respectively. 

Use:
denoise(cleanfile,noisefile,10)

Output:
A wavfile of the reconstructed output is generated. The spectrogram is plotted for the clean, noisy and reconstructed
files. This image is also saved in the local directory.
