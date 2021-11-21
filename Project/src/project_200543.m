% TODO
% 1. Otázkou je jestli můžu použít buffer i na
% t bo indexace


%% INFO / Task 1

clc;
clear;
% close all;

% print basic info header
fprintf('---------------------------------\n');
fprintf('[INFO] VUT ID: 200543\n');
fprintf('[INFO] Author: Martin Juricek\n');
fprintf('[INFO] IACS FME BUT @2021\n');
fprintf('---------------------------------\n');

% read specified audio file from assignment
[s, Fs] = audioread('../audio/200543');

% get signal as row vector
s = s';

% length of samples
len_samples = length(s);

% length in seconds
len_seconds = len_samples / Fs;

% maximum of signal
max_s = max(s);

% minimum of signal
min_s = min(s);

% print
fprintf('[INFO] Length of samples: %i\n',len_samples);
fprintf('[INFO] Length in sec: %.2f\n',len_seconds);
fprintf('[INFO] Max: %.2f\n',max_s);
fprintf('[INFO] Min: %.2f\n',min_s);

% plot
% time to sample
t = (0:len_samples -1) / Fs;

figure(1)
plot(t, s)
title('Input signal');
xlabel('t [s]');
xlim([-0.1 len_seconds+0.1]);
ylim([min_s-0.01 max_s+0.01]);

%% Task 2 
% http://noel.feld.cvut.cz/vyu/ucz/cv12/index.htm

% centering signal
s_cent = s - mean(s);

% normaliazing signal
s_norm = s_cent / abs(max_s);

% split norm. signal to frames to 1024 with oversampling 512
frames = buffer(s_norm, 1024, 512, 'nodelay')';
t_frames = buffer(t, 1024, 512, 'nodelay')';

% choose "nice" frames: 5, 11

% plot
figure(2)
plot(t_frames(5,:), frames(5,:) , 'r');
title('"Nice" sounded signal');
xlabel('t [s]');
xlim([t_frames(5,1)-0.005 t_frames(5,end)+0.005]);

%% Task 3

% defined N 
N = 1024;

% get sounded signal frame
sig = frames(5,:);

% compute dft
trans_s = dft(sig,N);

% compute fft for compare
comp_s = fft(sig);

% get frequencies
f = (0:length(trans_s)-1)/length(trans_s) * Fs;

% plot DFT
figure(3)
subplot(2,1,1);
stem(f, abs(trans_s), 'b');
grid on;
title('DFT module');
xlabel('frequency [Hz]');
xlim([0 Fs/2]);

% visual compare FFT
subplot(2,1,2);
stem(f, abs(trans_s), 'b');
hold on;
stem(f, abs(comp_s), 'r');

grid on;
title('FFT module');
xlabel('frequency [Hz]');
xlim([0 Fs/2]);

legend({'DFT','FFT'});
hold off;

%% Task 4
N = 1024;
S = 512;

[sgr,f_sgr,t_sgr] = spectrogram(s_norm, N, S,[], Fs/2);
P = 10 * log10(1 / N * abs(sgr).^2);

figure(4)
pcolor(t_sgr,f_sgr,P);
shading flat;
z = colorbar;
title('Signal spectogram');
ylabel(z,'power spectral density [dB]');
xlabel('t [s]');
ylabel('frequency [Hz]');
axis tight;
colormap(jet);

%% Task 5

f1 = 350;
f2 = 690;
f3 = 1030;
f4 = 1370;

%% Task 6

t = (1:1/Fs:len_seconds);
cos_out = cos(2 * pi * f1 * t) + cos(2*pi*f2*t) + cos(2*pi*f3*t) + cos(2*pi*f4*t);
% audiowrite('../audio/4cos.wav',cos_out,Fs,'BitsPerSample',16);
% sound(cos_out,Fs)

[sgr,f_sgr,t_sgr] = spectrogram(cos_out, N, S,[], Fs);
P = 10 * log10(1 / N * abs(sgr).^2);

figure(6)
pcolor(t_sgr,f_sgr,P);
shading flat;
z = colorbar;
title('Disturbing cosines spectogram');
ylabel(z,'power spectral density [dB]');
xlabel('t [s]');
ylabel('frequency [Hz]');
axis tight;
colormap(jet);

%% Task 7

freq = [f1,f2,f3,f4];
[s_filtred, z, p, b, a] = four_filter(freq, Fs/2, s);

%% Task 8

figure(8)
for i=1:length(freq)
    subplot(length(freq),1,i);
    zplane(z(:,i),p(:,i));
    title(['Zeros and poles: filter n. ', num2str(i)]);
end

%% Task 9

figure(9)

% freqz(b(1,:),a(1,:),[],Fs/2); 

for i=1:length(freq)
    [h,w] = freqz(b(i,:),a(i,:));
    subplot(length(freq),1,i);
    plot(w / 2 / pi * Fs/2, abs(h), 'LineWidth', 2)
    title(['Filter n.', num2str(i),' freq. char. module |H(e^{j\omega})|'])
    xlim([0 + i*150 350 + i*300]);
    ylim([-0.1 1]);
    xlabel('Frequency [Hz]');
end

%% Task 10

% centering signal
s_f_cent = s_filtred - mean(s_filtred);

% normaliazing signal
s_f_norm = s_f_cent / abs(max(s_filtred));

sound(s_f_norm, Fs)
audiowrite('../audio/clean_bandstop.wav', s_f_norm, Fs,'BitsPerSample',16);

% [s,f,t] = spectrogram(s,1024,512,[],Fs/2);
% P = 10 * log10(1 / N * abs(s).^2);
% pcolor(t,f,P),shading flat; colorbar;

%% Functions

function [trans_s] = dft(s,N)
    mat_B = dftmtx(N);
    trans_s = s * mat_B;
end

function [s, z_s, p_s, b_s, a_s] = four_filter(f, Fs, y)
    % https://documentation.sas.com/doc/en/pgmsascdc/9.4_3.5/casforecast/casforecast_dfil_details03.htm
    % Nyquist Frequency    
    Fn = Fs/2;                                                          
    % Passband Ripple
    Rp =  3;                                                           
    % Passband Ripple (Attenuation)
    Rs = 40;

    z_s = [];
    p_s = [];

    a_s = [];
    b_s = [];

    for i=1:length(f)
        % Stopband Frequency (Normalised)
        Wp = [f(i)-20 f(i)+20]/Fn;                                                  
        
        % Passband Frequency (Normalised)
        Ws = [f(i)-50 f(i)+50]/Fn;                                                  
        
        % Elliptic Order Calculation
        [n,Wp] = buttord(Wp,Ws,Rp,Rs);                                     
        
        % Elliptic Filter Design: Zero-Pole-Gain 
        [z,p,k] = butter(n,Wp,'stop');
        
        % Elliptic Filter Design: transfer function coefficients
        [b,a] = butter(n,Wp,'stop');
        
        b_s = [b_s; b];
        a_s = [a_s; a]; 

        % [z, p, k] = tf2zpk(b, a);
        
        z_s = [z_s z]; 
        p_s = [p_s p];

        figure(7)
        subplot(length(f),1,i);

        [h, t] = impz(b,a,50);
        plot(t, h, 'LineWidth',2);
        title(['Impulse response frequency ', num2str(f(i)),' [Hz]']);
        xlabel('samples');

        fprintf("[INFO] Coeficient of filter for frequency %i [HZ]\n",f(i));
        fprintf("[INFO] [b] = [");
        fprintf(' %.2g ', b);
        fprintf(']\n');
        fprintf("[INFO] [a] = [");
        fprintf(' %.2g ', a);
        fprintf(']\n');
        fprintf('\n');
        
        % Second-Order Section For Stability
        [sos,g] = zp2sos(z,p,k);
        
        if i == 1
            % Filter signal
            s = filtfilt(sos, g, y);
        else
            s = filtfilt(sos, g, s);
        end
    end
end