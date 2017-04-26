f = 1.5*1000; % freq in Hz
tL = 12000; %tone length (time) in ? 12000 is about a minute
Fs = 20500; % sampling frequency
t = 0:0.001:tL;
y = sin(2*pi*f*t);

filename = 'srcsignal1500.wav';
audiowrite(filename,y,Fs)

% amp=10; 
% fs=20500;  % sampling frequency
% duration=10;
% freq=1500;
% values=0:1/fs:duration;
% a=amp*sin(2*pi* freq*values);
% %sound(a)
% filename = 'srcsignal1500.wav';
% audiowrite(filename,y,Fs)
