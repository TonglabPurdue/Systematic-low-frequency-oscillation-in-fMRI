function y=filtf(input_signal,lpass, hpass,freq)
L=length(input_signal);
N=floor(L/10);
signal_padding = padarray(input_signal',[0 N],'symmetric','both');
% input_signal=input_signal.*hanning(L)
[b,a] = butter(3, [lpass*2/freq, hpass*2/freq], 'bandpass');
y_padding = filtfilt(b, a, signal_padding');
y=y_padding(N+1:N+L);