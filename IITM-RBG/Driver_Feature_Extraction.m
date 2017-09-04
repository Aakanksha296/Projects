cd /Users/Dave/Documents/MATLAB/EMG/mcode/;
load('Driver1_2nd5mins.mat');
data = val;

corrected_data = detrend(data);
sample_size = size(corrected_data);
sample_length = length(corrected_data);
fs = 1926;
f_nyq = fs/2;
T_period = 5.192308e-4;
T_total = T_period * sample_length;


plot( 1:sample_length, corrected_data)
xlabel('Sample Number')
ylabel('Detrended Data')
rectified_data=abs(corrected_data);
figure(2)
plot(1:sample_length, rectified_data)
xlabel('Sample Number')
ylabel('Rectified Data')
mean_data = mean(rectified_data);
%d = designfilt('bandpassiir', ...
%               'CutoffFrequency1', 20, 'CutoffFrequency2', 5000, ...
%               'SampleRate', 50000, 'DesignMethod', 'equiripple');
% 
% filter_data=filter(emg_healthym_butter,rectified_data);
% figure(3)
% plot(filter_data, 'r-')
% 
% MVE = max(rectified_data); %NOT USEFUL AS ATUAL MAX IS KNOWN - 0.0029
% 
%  display(MVE);

MVE=0.0055;

Gap = zeros(1, sample_length);
Gap_Count = 0;
for j=1:sample_length
   if (rectified_data(j)<(MVE*0.01))
       Gap(j) = rectified_data(j);
       Gap_Count=Gap_Count+1;
   end
    
end

y = abs(rectified_data-mean_data);    

movsum = zeros(1, sample_length);
z_ma = zeros(1, sample_length);
for n= 151:sample_length-150
     sum = 0;
    for i= (n-150):(n+150)
        sum = sum + y(i);
  
    end    
    movsum(n) = sum;
    z_ma(n) = movsum(n)/301;
   
end
figure(3)
plot (1:sample_length, y)
xlabel('Sample Number')
ylabel('Data Value')
hold on
plot(1:sample_length, z_ma, 'r')
legend('Subtracted Mean Data', 'Average Rectified Value (Moving Window)')

figure(4)
plot(1:sample_length, rectified_data)
xlabel('Sample Number')
ylabel('Data Value')
hold on
plot(1:sample_length, z_ma, 'r')
legend('Rectified Data', 'Average Rectified Value (Moving Window)')
% figure(5)
% plot(envelope(rectified_data, 300, 'rms'));

sq_y = y.^2;
rmssum = zeros(1, sample_length);
z_rms = zeros(1,sample_length);
for n= 151:sample_length-150
     sum = 0;
    for i= (n-150):(n+150)
        sum = sum + sq_y(i);
  
    end    
    rmssum(n) = sum;
    z_rms(n) = sqrt(rmssum(n)/301);
end

figure(5)
plot(1:sample_length, rectified_data)
xlabel('Sample Number')
ylabel('Data Value')
hold on
plot(1:sample_length, z_rms, 'r')
legend('Rectified Data', 'Root Mean Squared Value (Moving Window)')

%% Butterworth Filter

f_cutoff = 5;
[b,a]=butter(4,f_cutoff*1.116/f_nyq);

filt_data = filtfilt(b,a, rectified_data);

figure(6)
plot(1:sample_length, rectified_data)
xlabel('Sample Number')
ylabel('Data Value')
hold on
plot(1:sample_length, filt_data, 'r')
legend('Rectified Data', 'Butterworth Filter')

figure(7)
plot(1:sample_length, rectified_data)
xlabel('Sample Number')
ylabel('Data Value')
hold on
plot(1:sample_length, filt_data, 'r')
hold on
plot(1:sample_length, z_rms, 'g')
legend('Rectified Data', 'RMS')
hold on
plot(1:sample_length, z_ma, 'y')
legend('Rectified Data', 'Butterworth Filter', 'RMS', 'Moving Avg')

%% Band Pass Filter

[b1,a1]=butter(2,[5/f_nyq, 500/f_nyq], 'bandpass');

filt_data1 = filtfilt(b1,a1, rectified_data);

figure(8)
plot(1:sample_length, rectified_data)
xlabel('Sample Number')
ylabel('Data Value')
hold on
plot(1:sample_length, filt_data1, 'r')
legend('Rectified Data', 'Butterworth Filter')

%% FFT
Sig_FFT = fft(filt_data1);
P2 = abs(Sig_FFT/sample_length); %2 sided spectrum
P1 = P2(1:sample_length/2+1); % Single sided spectrum
P1(2:end-1) = 2*P1(2:end-1);

f = fs*(0:(sample_length/2))/sample_length;
figure(9)
plot(f,P1)
title('Single-Sided Amplitude Spectrum of Filtered Signal')
xlabel('f (Hz)')
ylabel('|P1(f)|')


Y= abs(fft(rectified_data));
Y(1) = [];
power = abs(Y(1:sample_length/2)).^2;
freq = (1:sample_length/2)/(sample_length/2)*f_nyq;
figure(10)
plot(freq,power), grid on
xlabel('Sample number (in Frequency)')
ylabel('Power spectrumen');
title({'Single-sided Power spectrum' ...
    ' (Frequency in shown on a log scale)'});
axis tight

xfft = fft(rectified_data-mean_data);
freqs=0:fs/sample_length:f_nyq;
figure(11)
plot(freqs,abs(xfft(1:sample_length/2+1)))
title('Single-Sided Amplitude Spectrum of Raw Signal w/ Subtracted Mean')
xlabel('f (Hz)')
ylabel('|P1(f)|')


Pxx = xfft.*conj(xfft);
figure(12)
plot(freqs,abs(Pxx(1:sample_length/2+1)))
title('Single-Sided Power Spectrum of Raw Signal w/ Subtracted Mean')
xlabel('f (Hz)')
ylabel('Power')


mean_freq = meanfreq(filt_data1, fs); %Signal is filtered, no fft needed

display(mean_freq)

median_freq = medfreq(filt_data1, fs); %Signal is filtered, no fft needed
display(median_freq)

figure(13)
meanfreq(filt_data1, fs)

figure(14)
medfreq(filt_data1, fs)

figure(15)
plot(Gap)

display(Gap_Count)

Gap_Freq = Gap_Count/T_total
