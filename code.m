clc;
clear;
% reading signal 1
[audio1,sample_freq1]=audioread('signal1.m4a');
audio1=audio1(:,1)+audio1(:,2);

% plotting signal 1 in time domain
timeAxis1=linspace(0,length(audio1)/sample_freq1,length(audio1));
figure(1);
subplot(2,1,1);
plot(timeAxis1,audio1);
title('Before Modulation signal 1');
xlabel('time');
ylabel('amplitude');

% plotting signal 1 in frequency domain
N1=size(audio1,1);
freq1=linspace(0,sample_freq1,N1);
audiomag1=abs(fft(audio1,N1));
subplot(2,1,2);
plot(freq1(1:N1),audiomag1(1:N1));
title('Before Modulation signal 1');
xlabel('frequency');
ylabel('magnitude');

%reading signal 2
[audio2,sample_freq2]=audioread('signal2.m4a');
audio2=audio2(:,1)+audio2(:,2);

% plotting signal 2 in time domain
timeAxis2=linspace(0,length(audio2)/sample_freq2,length(audio2));
figure(2);
subplot(2,1,1);
plot(timeAxis2,audio2);
title('Before Modulation signal 2');
xlabel('time');
ylabel('amplitude');

% plotting signal 2 in frequency domain
N2=size(audio2,1);
freq2=linspace(0,sample_freq2,N2);
audiomag2=abs(fft(audio2,N2));
subplot(2,1,2);
plot(freq2(1:N2),audiomag2(1:N2));
title('Before Modulation signal 2');
xlabel('frequency');
ylabel('magnitude');

%reading signal 3
[audio3,sample_freq3]=audioread('signal3.m4a');
audio3=audio3(:,1)+audio3(:,2);

% plotting signal 3 in time domain
timeAxis3=linspace(0,length(audio3)/sample_freq3,length(audio3));
figure(3);
subplot(2,1,1);
plot(timeAxis3,audio3);
title('Before Modulation signal 3');
xlabel('time');
ylabel('amplitude');

% plotting signal 2 in frequency domain
N3=size(audio3,1);
freq3=linspace(0,sample_freq3,N3);
audiomag3=abs(fft(audio3,N3));
subplot(2,1,2);
plot(freq3(1:N3),audiomag3(1:N3));
title('Before Modulation signal 3');
xlabel('frequency');
ylabel('magnitude');


%resampling signal 1 ( forcing three signals to have the same sampling frequency to be able to modulate % them )
resample_freq = 260000;
[p, q] = rat(resample_freq/sample_freq1);
resampled_signal1 = resample(audio1, p, q);

%resampling signal 2
[p, q] = rat(resample_freq/sample_freq2);
resampled_signal2 = resample(audio2, p, q);

%resampling signal 3
[p, q] = rat(resample_freq/sample_freq3);
resampled_signal3 = resample(audio3, p, q);

% modulation
% ( reasonable frequencies were chosen for carriers , same frequency were chosen for signal 2
% and signal 3 as we will use quam modulation to modulate them )
f1 = 50000; 
f2 = 100000;

timescale1 = (0:length(resampled_signal1) - 1) * (1/resample_freq);
carrier1 = cos(2*pi*(f1) * timescale1);

timescale2 = (0:length(resampled_signal2) - 1) * (1/resample_freq);
carrier2 = cos(2*pi*(f2) * timescale2);

timescale3 = (0:length(resampled_signal3) - 1) * (1/resample_freq);
carrier3 = sin(2*pi*(f2) * timescale3);

% modulating carriers
out1 = resampled_signal1 .* carrier1';
out2 = resampled_signal2 .* carrier2';
out3 = resampled_signal3 .* carrier3';

% getting modulated signals length
mslen1 = length(out1);
mslen2 = length(out2);
mslen3 = length(out3);

% plotting modulated signal 1 in frequency domain
fmsignal1 = (resample_freq/mslen1)  * (-mslen1/2: mslen1/2 - 1);
fmmsignal1 = abs(fft(out1));
figure(4);
plot(fmsignal1, fftshift(fmmsignal1));
xlabel("frequency");
ylabel("magnitude");
title("modulated signal 1 ");

% plotting modulated signal 2 in frequency domain
fmsignal2 = (resample_freq/mslen2)  * (-mslen2/2: mslen2/2 - 1);
fmmsignal2 = abs(fft(out2));
figure(5);
plot(fmsignal2, fftshift(fmmsignal2));
xlabel("frequency");
ylabel("magnitude");
title("modulated signal 2 ");

% plotting modulated signal 3 in frequency domain
fmsignal3 = (resample_freq/mslen3)  * (-mslen3/2: mslen3/2 - 1);
fmmsignal3 = abs(fft(out3));
figure(6);
plot(fmsignal3, fftshift(fmmsignal3));
xlabel("frequency");
ylabel("magnitude");
title("modulated signal 3 ");

% forcing three signals to have the same length (maximum length of these signals) to be able to add
 % them
maximumlength = max(mslen1, max(mslen2, mslen3));
signal1 = [out1;zeros(maximumlength-mslen1, 1)];
signal2 = [out2;zeros(maximumlength-mslen2, 1)];
signal3 = [out3;zeros(maximumlength-mslen3, 1)];
% getting s(t)
mod_signal = signal1 + signal2 + signal3;

% plotting modulated signal in time domain
tms = (0: maximumlength - 1) * (1 / resample_freq);
fms = (-maximumlength/2 : maximumlength/2 - 1) * (resample_freq / maximumlength);
figure(7);
subplot(2,1,1);
plot(tms, mod_signal);
xlabel("time");
ylabel("amplitude");
title("Modulated signal");

% plotting modulated signal in frequency domain
fftmodsignal = abs(fft(mod_signal));
subplot(2,1,2);
plot(fms, fftshift(fftmodsignal));
xlabel("frequency");
ylabel("magnitude");
title("Modulated signal");

% demodulation  carriers
time_demod = (0:maximumlength - 1) * (1/resample_freq);
c1 = cos(2*pi*(f1) * time_demod);
c2 = cos(2*pi*(f2) * time_demod);
c3 = sin(2*pi*(f2) * time_demod);


% demodulation & filtering demodulation output
fres = (-maximumlength/2 : maximumlength/2 - 1) * (resample_freq / maximumlength);

demodulatesignal(mod_signal, c1, length(audio1), sample_freq1, resample_freq, 24000, fres, 'out1');
demodulatesignal(mod_signal, c2, length(audio2), sample_freq2, resample_freq, 24000, fres, 'out2');
demodulatesignal(mod_signal, c3, length(audio3), sample_freq3, resample_freq, 24000, fres, 'out3');
 
% phase shift 10 degrees
[c1_10, c2_10, c3_10] = getcarrier(f1, f2, time_demod, 10);
demodulatesignal(mod_signal, c1_10, length(audio1), sample_freq1, resample_freq, 24000, fres, 'out1_10');
demodulatesignal(mod_signal, c2_10, length(audio2), sample_freq2, resample_freq, 24000, fres, 'out2_10');
demodulatesignal(mod_signal, c3_10, length(audio3), sample_freq3, resample_freq, 24000, fres, 'out3_10');
 
% phase shift 30 degrees
[c1_30, c2_30, c3_30] = getcarrier(f1, f2, time_demod, 30);
demodulatesignal(mod_signal, c1_30, length(audio1), sample_freq1, resample_freq, 24000, fres, 'out1_30');
demodulatesignal(mod_signal, c2_30, length(audio2), sample_freq2, resample_freq, 24000, fres, 'out2_30');
demodulatesignal(mod_signal, c3_30, length(audio3), sample_freq3, resample_freq, 24000, fres, 'out3_30');

% phase shift 90 degrees
[c1_90, c2_90, c3_90] = getcarrier(f1, f2, time_demod, 90);
demodulatesignal(mod_signal, c1_90, length(audio1), sample_freq1, resample_freq, 24000, fres, 'out1_90');
demodulatesignal(mod_signal, c2_90, length(audio2), sample_freq2, resample_freq, 24000, fres, 'out2_90');
demodulatesignal(mod_signal, c3_90, length(audio3), sample_freq3, resample_freq, 24000, fres, 'out3_90');

% shifting signal 1 carrier 2 hz and demodulating it again
carrier1_shifted2hz=cos(2*pi*(f1-2)*time_demod);
demodulatesignal(mod_signal, carrier1_shifted2hz, length(audio1), sample_freq1, resample_freq, 24000, fres, 'out1-2hz');

% shifting signal 1 carrier 10 hz and demodulating it again
carrier1_shifted10hz=cos(2*pi*(f1-10)*time_demod);
demodulatesignal(mod_signal, carrier1_shifted10hz, length(audio1), sample_freq1, resample_freq, 24000, fres, 'out1-10hz');


% demodulation function
function  demodulatesignal(s, carr, oldlen, fsold, fsnew, fpass, fres, filename)
    %demodulating signal
    ress = s .* carr';

    %filtering demodulated signal
    res = lowpass(ress, fpass , fsnew);

    % plotting modulated signal in frequency domain
    fftress = abs(fft(res));
    figure();
    plot(fres, fftshift(fftress));
    title(strcat(filename, " demodulated magnitude"));

    % resampling signal to its original sampling frequency
    [P, Q] = rat(fsold/fsnew);
    res = resample(res, P, Q);
    res = res(1: oldlen);
     % writing output signal in output .wav file 
    audiowrite(strcat(filename, '.wav'), res, fsold)
end

% function to get carrier with specific frequency
function [c1, c2, c3] = getcarrier(f1, f2, tdem, phaseshiftdeg)
     phaseshiftrad = (phaseshiftdeg * pi) / 180; 
     c1 = cos(2*pi*f1*tdem + phaseshiftrad);
     c2 = cos(2*pi*f2*tdem + phaseshiftrad);
     c3 = sin(2*pi*f2*tdem + phaseshiftrad);
end
