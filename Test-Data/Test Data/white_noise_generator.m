snr_ = zeros(1,10);
for i = 1:10
    [x, fs] = audioread(['SA1_', int2str(i), '.WAV']);
    SNR = 10;
    y = randn(size(x))*std(x)/db2mag(SNR);

    s = x + y;
    audiowrite(['SA1_', int2str(i), '_10DB.WAV'], s, fs);

    snr_(1,i) = snr(x,y)
end
%% 
[x_result, fs] = audioread('SA1_result.WAV');
soundsc(x_result,fs)
y_result = x_result-x;
snr(x_result,y_result)
%% 
for i = 1:10
    [x, fs] = audioread(['SA1_', int2str(i), '.WAV']);
    audiowrite(['SA1_', int2str(i), '.WAV'], x, fs);
end