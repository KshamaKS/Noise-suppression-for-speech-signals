import numpy as np
import soundfile as sf
import librosa as lr

def nearest_harmonic_band(k, f0, N):
    # given a STFT band k and the estimated fundamental frequency,
    # this function returns the nearest STFT harmonic band k_ and harmonic frequency fh
    w = 2*np.pi*k/N   # frequency of current band
    h = round(w/f0)   # hth harmonic is the closest harmonic
    k_ = int(h*N/(2*np.pi))   # STFT band that hth harmonic is closest to
    return [k_, f0*h]

def window_phase(w, N):
    # this function calculates the phase of a length-N hamming filter at frequency w
    if w == 0:
        return 0.54*2*np.pi
    if w == 2*np.pi/(N-1) or w == 2*np.pi - 2*np.pi/(N-1):
        return -0.46*np.pi
    return 0

def improve_phase(filename):
    data, fs = sf.read(filename)
    print(filename)
    N = int(0.032*fs)   # each frame has a length of 32 ms
    L = int(0.016*fs)   # 16 ms time lag between adjacent frames
    num_frames = int((len(data) - N + L) / L)
    print(N, L, len(data), num_frames)

    hamming = np.hamming(N)
    print(data)

    Y = np.zeros((N, num_frames), dtype=complex)   # STFT
    YB = np.zeros((N, num_frames), dtype=complex)   # baseband STFT
    Phi_SB = np.zeros((N, num_frames))   # baseband phase after improvement
    Y_rec = np.zeros((N, num_frames), dtype=complex)   # STFT after improvement
    data_rec = np.zeros(len(data))    # reconstructed speech

    f0 = lr.yin(data, 65, 800, fs, N, None, L)

    for l in np.arange(num_frames):

        windowed_data = hamming * data[l*L : l*L+N]   # apply hamming window to the frame

        if f0[l] <= 100 or f0[l] > 400:   # this frame is unvoiced
            data_rec[l*L : l*L+N] = data_rec[l*L : l*L+N] + windowed_data

        else:   # this frame is voiced
            
            Y[:,l] = np.fft.fft(windowed_data)   # STFT

            demodulator = np.zeros(N, dtype=complex)
            for k in np.arange(N):
                demodulator[k] = np.exp(-1j*(2*np.pi*k/N)*l*L)

            YB[:,l] = Y[:,l] * demodulator   # shift the STFT coefficients to baseband

            for k in np.arange(N):
                [k_h, fh] = nearest_harmonic_band(k, f0[l], N)
                Phi_SB[k,l] = np.angle(YB[k_h,l]) - (k-k_h)*2*np.pi/N*l*L
                + window_phase(k-fh*N/(2*np.pi),N) - window_phase(k_h-fh*N/(2*np.pi),N)   # equation 12
                
                # combine amplitude and improved phase, and do upconversion
                Y_rec[k,l] = np.abs(YB[k,l])*np.exp(1j*Phi_SB[k,l])*np.exp(1j*(2*np.pi*k/N)*l*L)
            
            data_rec[l*L : l*L+N] = data_rec[l*L : l*L+N] + np.fft.ifft(Y_rec[:,l])

        print(l)
    
    sf.write('SA1_result.WAV', data_rec, fs)

improve_phase('SA1.WAV')


