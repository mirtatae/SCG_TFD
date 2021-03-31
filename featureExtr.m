function [all_dom_freq2,T_spectro1] = featureExtr(Sig, Fs, t_IF, PCT_input_1, PCT_input_2, STFT_input_1, STFT_input_2, STFT_input_3,SPWV_input_1, SPWV_input_2, SPWV_input_3,Rem_out_peak)
%% Code Description
%
%
% Please cite the paper if you used this code:
%
% A. Taebi and H. A. Mansy, "Analysis of seismocardiographic signals using
% polynomial chirplet transform and smoothed pseudo Wigner-Ville
% distribution," 2017 IEEE Signal Processing in Medicine and Biology
% Symposium (SPMB), Philadelphia, PA, 2017, pp. 1-6, doi:
% 10.1109/SPMB.2017.8257022.
%
% Version History
% 
%
% + Inputs:
%        Sig:   the signal to be analyzed
%         Fs:   sampling frequency
%       t_IF:   a specified time that we are interested to know the IF of
%               the signal at it.
%     
% 
% + Outputs:
%     
%
% + Note:
%
% + NEXT STEP TO IMPROVE THE CODE (temporary note):
%        
% Author: Amirtaha Taebi
% Biomedical Acoustics Research Laboratory
% University of Central Florida

%% assigning default values

t = 0:1/Fs:(length(Sig)-1)/Fs;
y = Sig';
y = y - mean(y);
length(y)

if nargin < 3
    t_IF = [];
end

if isempty(t_IF)
    t_IF = min(t);
end

if t_IF > max(t)
    error('argValue:notFound', 'Error! Requested time is out of the signal length range. Please choose a proper value!');
elseif t_IF < 0
    error('argValue:notFound', 'Error! Time should always be positive. Please choose a proper value!');
end

%% Switches
Plot_on = 0;
TFD_amp_square_on = 1;              % switch between amp and amp square for the spectrum of different TFD techniques
TFD_amp_normalizing_on = 1;         % normalizes the spectrum amp of different TFD techniques
VCG_high_freq_off = 0;

%% defining required parameters
% figure position parameters
left = 50; 
bottom = 50; 
width = 500; 
height = 400;

%% chirplet transform of VCG
iSig = y;

if Plot_on == 1
    hFig1 = figure(1);
    set(hFig1, 'Position', [left, bottom, width, height]);
    plot(t,abs(iSig))
    set(gcf,'Color','w');
end

[Spec1,f1,t1] = FPCT(iSig',Fs,0,PCT_input_1,PCT_input_2);

if VCG_high_freq_off == 1
    Spec1(125:length(f1),:) = 0;
end

if TFD_amp_square_on == 1
    Spec1 = Spec1.^2;
end

if Plot_on == 1
    hFig2 = figure(2);
    set(hFig2, 'Position', [left, bottom, width, height]);
    surf(t1,f1,abs(Spec1),'EdgeColor','none');
    axis([min(t1) max(t1) min(f1) max(f1)]);
    xlabel('Time [s]');
    ylabel('Frequency [Hz]');
    title('TFD of synthetic VCG by using PCT');
    colormap hot
    colorbar('eastoutside','FontSize',8);
end

% finding instantaneous freq. (using global maximum ?energy)
[Mx_Spec1,i_mx_Spec1] = max(Spec1);
f_dom_Spec1 = f1(i_mx_Spec1);

% finding the spectrum: first steps
[Max_Spec1,i_max_Spec1] = max(Spec1,[],2);
[Max_Max_Spec1,i_max_max_Spec1] = max(Max_Spec1);


if TFD_amp_normalizing_on == 1
    Spec1 = Spec1 / Max_Max_Spec1;    % normalizing the spectreum
    Mx_Spec1 = Mx_Spec1 / Max_Max_Spec1;
end

%% spectrogram

[Spectro1,F_spectro1,T_spectro1, PSD_spectro1] = spectrogram(y,STFT_input_1, STFT_input_2, STFT_input_3,Fs); % Display the spectrogram

if VCG_high_freq_off == 1
    PSD_spectro1(25:length(F_spectro1),:) = 0;
end

if TFD_amp_square_on == 0
    PSD_spectro1 = sqrt(PSD_spectro1);
end

if Plot_on == 1
    hFig4 = figure(4);
    set(hFig4, 'Position', [left, bottom, width, height]);
    surf(T_spectro1,F_spectro1,PSD_spectro1,'EdgeColor','none');


    xlabel('Time [s]');   
    ylabel('Frequency [Hz]');
    title('TFD of synthetic VCG by using STFT');
    colormap hot;   
    colorbar('eastoutside','FontSize',8);
end

[Mx_PSD_spectro1,i_mx_PSD_spectro1] = max(PSD_spectro1);
f_dom_PSD_spectro1 = F_spectro1(i_mx_PSD_spectro1);

[Max_PSD_spectro1,i_max_PSD_spectro1] = max(PSD_spectro1,[],2);
[Max_Max_PSD_spectro1,i_max_max_PSD_spectro1] = max(Max_PSD_spectro1);  % finding the max amplitude of the spectrum

if TFD_amp_normalizing_on == 1
    PSD_spectro1 = PSD_spectro1 / Max_Max_PSD_spectro1;    % normalizing the spectreum
    Mx_PSD_spectro1 = Mx_PSD_spectro1 / Max_Max_PSD_spectro1;
end

%% Wigner-Ville Distribution
[tfr,t_wv,f_wv] = wv(y,t);

if VCG_high_freq_off == 1
    tfr(:,245:length(f_wv)) = 0;
end

tfr = max(tfr,0);       % to remove the negative energy values

if TFD_amp_square_on == 0
    tfr = sqrt(abs(tfr));
end

[Mx_tfr,i_mx_tfr] = max(tfr');
f_dom_tfr = f_wv(i_mx_tfr);

[Max_tfr,i_max_tfr] = max(tfr,[],1);
[Max_Max_tfr,i_max_max_tfr] = max(Max_tfr);  % finding the max amplitude of the spectrum

if TFD_amp_normalizing_on == 1
    tfr = tfr / Max_Max_tfr * Rem_out_peak;    % normalizing the spectreum
    Mx_tfr = Mx_tfr / Max_Max_tfr;
end

if Plot_on == 1
    figure(28)
    subplot(311);
    plot(y);
    [F_wv, T_wv] = meshgrid(f_wv, t_wv);
    subplot(312);
    mesh(F_wv, T_wv, abs(tfr));
    subplot(313);
    surf(f_wv, t_wv, abs(tfr),'EdgeColor','none');
    axis([0 f_wv(length(f_wv)) 0 t_wv(length(t_wv))]);
    xlabel('frequency');
    ylabel('time');
    view(90, -90);
    colormap('hot');
end

%% Smoothed pseudo Wigner-Ville distribution

[tfr_spwv,t_spwv,f_spwv] = tfrspwv(y', 1:length(y), SPWV_input_1, SPWV_input_2, SPWV_input_3, 0);

if VCG_high_freq_off == 1
    tfr_spwv(245:length(f_spwv),:) = 0;
end


if TFD_amp_square_on == 0
    tfr_spwv = sqrt(abs(tfr_spwv));
end

[Mx_spwv,i_mx_spwv] = max(tfr_spwv);
f_dom_spwv = Fs*f_spwv(i_mx_spwv);

[Max_spwv,i_max_spwv] = max(tfr_spwv,[],2);
[Max_Max_spwv,i_max_max_spwv] = max(Max_spwv);  % finding the max amplitude of the spectrum

if TFD_amp_normalizing_on == 1
    tfr_spwv = tfr_spwv / Max_Max_spwv;    % normalizing the spectreum
    Mx_spwv = Mx_spwv / Max_Max_spwv;
end

all_dom_freq2 = {f_dom_PSD_spectro1', f_dom_Spec1, f_dom_tfr, f_dom_spwv};

end