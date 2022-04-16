function him = plot_phase_diff(xx,fsamp,Lsect)
%plot_phase_diff   plot phase difference between adjacent frames as an image
%
%  usage:   him = plotspec(xx,fsamp,Lsect,DBrange)
%      him = handle to the image object
%       xx = input signal
%    fsamp = sampling rate
%    Lsect = section length (integer, power of 2 is a good choice)
%              amount of data to Fourier analyze at one time

% 16-Feb-2013 J McClellan, created from plotspec.m
% 11-Apr-2022 Y Teng, created from plotspecDB.m

if( nargin<3 )
	Lsect = 256;  %- default section length is 256
end
if( nargin<2 )
	disp('PLOTPHASEDIFF: Sampling Frequency defaulting to 8000 Hz')
	fsamp = 8000;
end
if( length(xx)<1000 )
	warning('PLOTPHASEDIFF: Signal length must be greater than 1000 to get a reasonable spectrogram')
end
Lfft = Lsect;
Noverlap = round(Lsect/2);  %-- overlap defaults to 50%
[B,F,T] = spectgr(xx,Lfft,fsamp,Lsect,Noverlap);
B_baseband = zeros(length(F), length(T));
for l = 1:length(T)
    B_baseband(:,l) = B(:,l).*(exp(-1j*(2*pi*(1:length(F))/length(F))*l*Lsect).');
end
B_ph = angle(B_baseband);
B_ph_diff = zeros(length(F), length(T));
for i = 2:length(T)
    for j = 1: length(F)
        if (abs(B(j,i)) < 0.02) 
            B_ph_diff(j,i) = 0;
        elseif B_ph(j,i) > B_ph(j,i-1)
            B_ph_diff(j,i) = B_ph(j,i) - B_ph(j,i-1);
        else
            B_ph_diff(j,i) = B_ph(j,i) + 2*pi - B_ph(j,i-1);
        end
    end
end
% B_ph_diff(:,2:length(T)) = B_ph(:,2:length(T))-B_ph(:,1:length(T)-1);
him = imagesc(T,F,B_ph_diff);
axis xy
colormap(1-gray)   %-- use colormap(jet) if you like bright colors !