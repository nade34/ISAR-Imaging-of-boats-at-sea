%% 11 November 2020
% Institute: University of Cape Town 
% Course: EEE4022s
% Name: Nadir Mahomed
% Student Number: ABBNAD002
% Supervisor: Dr M.Y. Abdul Gaffar
%% ------------------------------------------------------------------EMC-ATWS------------------------------------------------------------------------------
clear all;
close all;

c = 299792458;
% % Esp_CAP77-58INB_1_P728_B1_STC76_HRR
% % Tigress_INB_1_P727_B1_STC88_HRR : Checked (a few usable ones)
% % Tigress_INB_1_P728_B1_STC88_HRR : Checked
% % Tigress_OUT_2_P727_B1_STC88_HRR : Corrupt
% % Tigress_OUT_2_P728_B1_STC88_HRR : Checked
% % Tigress_TURN_1_P728_B1_STC76_HRR : Checked
% % Tigress_TURN_2_P728_B1_STC76_HRR : Checked
%
EsperanceDataset = load('Esp_CAP77-58INB_1_P728_B1_STC76_HRR ');
%EsperanceDataset = load('Esp_CAP77-51TURN_1min_rev_1_P728_B1_STC76_HRR');
% Obtain the relevant data from the Esperance Dataset
RadarDataset.HRR_profiles = EsperanceDataset.Settings.HRR.HRR_calib_velcomppeak_win.';
[NumberofHRRprofiles NumberofPulsesinISARBurst] = size(RadarDataset.HRR_profiles);
frequence_step_MHz = EsperanceDataset.Settings.Pattern.FStep(2)*1e6;
RangeResolution_m = c/(2*NumberofPulsesinISARBurst*frequence_step_MHz);
RadarDataset.Range_axis = (0:1:(NumberofPulsesinISARBurst-1))*RangeResolution_m;
RadarDataset.BRF_Hz = 154;


WindowMatrix = EsperanceDataset.Settings.HRR.Win_Func;
sb_Calib_Matrix =  EsperanceDataset.Settings.HRR.corr_filter_vector;

HRR_profiles = ifft(EsperanceDataset.Settings.HRR.TgtBin_velcomp_peak.*WindowMatrix.*sb_Calib_Matrix).';
HRR_profiles = HRR_profiles(185:587,:);
xaxis = RadarDataset.Range_axis;
yaxis = (1:1:size(RadarDataset.HRR_profiles,1));
% Dataset = load('Altera_CWLFM_Dataset.mat');
% HRR_profiles = Dataset.HRR_profiles(:,:);
% ProfileRepetitionFreq = Dataset.ProfileRepitionFreq;
% xaxis = Dataset.Range_axis;

N = size(HRR_profiles,2);
M = size(HRR_profiles,1);
window_range = [0.3,0.8];         %[0.2,0.45,0.75,1.2]
for range = 1:2
    window = window_range(range);
    fs = 154;
    PRF = fs;
    N = size(HRR_profiles,2);
    M = size(HRR_profiles,1);
    CPTWL = roundn(window*fs,-2);
    nfft = CPTWL;               % Numeber off fft points
    win = kaiser(CPTWL);        % Window Function
    factor = 0.75;               % Overlap factor
    overlap = factor*CPTWL;     % Overlap samples
    polyorder = 2;
    n = 6;
    wlen = length(win);
    f = (wlen/2:1:(wlen/2-1))*fs/wlen;
    [Autofocus,indexes,yes] = ISAR_IMAGE_GENERATIONV4(HRR_profiles,win,overlap,nfft,fs,xaxis,factor, polyorder);
    if indexes == 0
        disp ("window length does not give contrast peaks");
        continue;
    else
        if range == 1
            overlap_matrix = yes;
        else
            overlap_matrix = horzcat(overlap_matrix,yes);
        end
    end
end
overlap_matrix
counter = 1;
non_overlap(1:2,counter) = 0;
int_length = size(overlap_matrix,2);
for i = 1:int_length
    non_overlap
    index = overlap_matrix(1,1);
    wlen = overlap_matrix(5,1);
    non_overlap(1,counter) = index;
    non_overlap(2,counter) = wlen;
    StartProfile = overlap_matrix(3,1);
    StopProfile = overlap_matrix(4,1);
    overlap_matrix(:,1) = [];
    A = StartProfile:1:StopProfile;
    for check = 1:size(overlap_matrix,2)
        index = overlap_matrix(1,check);
        wlen = overlap_matrix(5,check);
        StartProfile2 = overlap_matrix(3,check);
        StopProfile2 = overlap_matrix(4,check);
        B = StartProfile2:1:StopProfile2;
        tf = ismember(A,B);
        if sum(tf)>0
            non_overlap(1,counter) = index;
            non_overlap(2,counter) = wlen;
        end
    end
    counter = counter +1;
end
non_overlap
counter = 1;
final_indexes(1,counter) = non_overlap(1,1);
final_indexes(2,counter) = non_overlap(2,1);
for redund = 1:size(non_overlap,2)
    if ismember(non_overlap(1,redund),final_indexes(1,:))
        continue
    else
        counter = counter +1;
        final_indexes(1,counter) = non_overlap(1,redund);
        final_indexes(2,counter) = non_overlap(2,redund);
    end
    
end
for i = 1:size(final_indexes,2)
    
    CPTWL = final_indexes(2,i);
    nfft = CPTWL;               % Numeber off fft points
    win = kaiser(CPTWL);        % Window Function
    overlap = factor*CPTWL;     % Overlap samples
    wlen = length(win);
    f = (wlen/2:1:(wlen/2-1))*fs/wlen;
    optimum_matrix = CPTWL_optV3(Autofocus,win, overlap, nfft, PRF,xaxis,final_indexes(1,i),factor,n);
    if i == 1
        newoptimum_matrix = optimum_matrix
    else
        newoptimum_matrix = horzcat(newoptimum_matrix,optimum_matrix)
    end
end
redf = 0;
finalf = 0;
optimum_matrix =  newoptimum_matrix;
if size(optimum_matrix,2)==1
    average_contrast = 0;
else
    average_contrast = mean(optimum_matrix(3,:));
end
for final = 1:size(optimum_matrix,2)
    %if optimum_matrix(3,final)> average_contrast
    redf = redf+1;
    nlen = optimum_matrix(2,final);
    index = optimum_matrix(1,final);
    overlap = factor*nlen;
    hop = nlen-overlap;
    [new_contrast,StopProfile,ISAR,Entropy,CentreProfile_time,ocontrast,oentropy,ISAR_new] = image_generation (nlen, hop, 'n',Autofocus,index,PRF);
    invert = 'n';
    finalf(1,redf) = index
    finalf(2,redf) = nlen
    finalf(3,redf) = ocontrast
    finalf(4,redf) = oentropy
    finalf(5,redf) = optimum_matrix(5,final)
    % end
end
finalf
finalISARImages (Autofocus,win, overlap, nfft, PRF,xaxis,factor,finalf)
function finalISARImages (y,win, overlap, nfft, PRF,Range_axis,factor,Optimum_matrix)
counter = 1;
for i = 1:size(Optimum_matrix,2)
    [max_contrast,index_1] = min(Optimum_matrix(4,:))
    ylen = size(y,1);                                                                            % Signal Length
    wlen = Optimum_matrix(2,index_1)                                                            % window function length should be 1024
    CPTWL = roundn(wlen*1/PRF,-2);
    %% Calculate the number of frames to be taken, given the signal size and amount of overlap
    overlap = wlen*factor;
    hop = wlen-overlap;
    frames = 1+floor((ylen-wlen)/(hop));
    %% Executing the STFT to generate multiple ISAR Images
    %% STFT
    StartProfile = 1+(Optimum_matrix(1,index_1)-1)*hop;                                                 % Abdul Gaffar
    StopProfile = wlen+(Optimum_matrix(1,index_1)-1)*hop;                                               % Abdul Gaffar
    CentreProfile = (StartProfile + StopProfile)/2;                                                     % Abdul Gaffar
    CentreProfile_time = roundn(CentreProfile*1/PRF,-2);                                                % Abdul Gaffar
    Ess = y(StartProfile:StopProfile,:).*kaiser(wlen).*repmat(kaiser(wlen),1,size(y,2));                % windowing of the sampled data that moves 'overlap' samples for respective frame
    Enew = Ess;
    %% Compute ISAR image
    ISAR = fftshift(fft(Enew,[],1),1);
    %% Contrast
    answer = max(ISAR);
    answer2 = max(answer);                                                                              % Getting the peak value
    ISAR_new = ISAR./answer2;
    limit = -30;                                                                                        % Limiting SNR
    for row = 1:size(ISAR_new,1)
        for column = 1:size(ISAR_new,2)
            if 20*log10(abs(ISAR_new(row,column)))<limit
                ISAR_new(row,column) = ((10^(limit/20))/(abs(ISAR_new(row,column))))*ISAR_new(row,column);
            end
        end
    end
    B = abs(ISAR).^2;
    C = mean(mean(B(),1));
    D = B - C;
    E = sqrt((mean(mean(D.^2,1))));
    Contrast = E/C;
    contrast_matrix(1,i+1) = Contrast;
    N = B/(sum(sum(B,1)));
    oentropy = -sum(sum(N.*log(N),1));
    [max_down_range_energy, energy_cell_index] = max(sum(20*log10(B)));
    B_zdt = round(30*(size(B,1)/PRF))
    zero_doppler = size(B,1)/2;
    negative_freqs = sum(20*log10(B(zero_doppler-B_zdt:zero_doppler,energy_cell_index)));
    positive_freqs = sum(20*log10(B(zero_doppler:zero_doppler+B_zdt,energy_cell_index)));
    E_threshold = 350; %dB;
    Energy_test = abs(negative_freqs - positive_freqs)
    if Energy_test < E_threshold
        Optimum_matrix(4,index_1) = 0;
    else
        Optimum_matrix(4,index_1) = oentropy;
    end
    if Optimum_matrix(4,index_1)==0
        Optimum_matrix(4,index_1) = 100000;
        continue;
    end
    %% Entropy
    N = B/(sum(sum(B,1)));
    Entropy = -sum(sum(N.*log(N),1));
    %% Figure
    f = (-wlen/2:1:(wlen/2-1))*PRF/wlen;
    figure;
    imagesc(Range_axis,f,20*log10(abs(ISAR_new)));
     title(sprintf("Image Rank = %0.5g, CenterTime = %0.5g s,  \n Entropy = %0.5g ,Contrast = %0.5g, \n Optimum CPTWL = %0.5g, Initial CPTWL = %0.5g ",counter,CentreProfile_time, oentropy, Contrast, CPTWL,roundn(Optimum_matrix(5,index_1)*1/PRF,-2))); % Abdul Gaffar
    colormap(flipud(gray));
    ylabel('Doppler Frequency (Hz)');
    xlabel('Range (m)');
    colorbar;
    axis xy;
    Optimum_matrix(4,index_1) = 100000;
    counter = counter +1;
end
end
function [Autofocus, indexes,overlap_matrix] =  ISAR_IMAGE_GENERATIONV4 (HRRP,win, overlap, nfft, fs, xaxis,factor, polyorder)
Alignment = alignment(HRRP,polyorder);
[Autofocus,Power_vector] = autofocus(Alignment);
[indexes,overlap_matrix] = GenISARImages(Autofocus,win, overlap, nfft, fs,xaxis,factor);
end
%% ------------------------------------------------------ RANGE ALIGNMENT-------------------------------------------------------------------------------------------------------------------------------------------------------------
function alignedmatrix = alignment(hrrp,polyorder)
M = size(hrrp,2);                                                           % Range bins, Fast time data
N = size(hrrp,1);                                                           % Slow time data
vector = (0:1:(M-1));
%% Step 1 : Determine the reference profile. For this case we are picking the first profile as the reference profile
ref_prof = hrrp(1,:);                                                       % First profile, M=1, is chosen as the reference profile
%% Step 2 : Perform cross correlation to determine the time delay of respective bin for alignment purposes (Working with magnitudes only)
mag_ref = abs(ref_prof);                                                    % Taking the magnitude of the reference profile
corr_matrix = conv2(fliplr(abs(hrrp)),mag_ref);                             % Calculates the correlation between ref and the range profiles using the convolution operation
%% Step 3 : Locate the peak and required bin shifts
[MaxP,Index] = max(corr_matrix.');                                          % Computing max and it's location of each bin
plot(MaxP);
binshifts = Index - M ;               
% Calculate the pos of peak, remember with corr index -3,-2,-1,0,1,2,3...
%% Step 4 : Smoothen Delay curve.
x = 0:1:N-1;                                                                % Represents the number of range profiles
calc_coeff = polyfit(x,binshifts,polyorder);                                % Fitting the delay curve 'D' with a low order polynomial
smooth_data = polyval(calc_coeff,x);                                        % Retreiving co-efficients for smoothen curve 's'
%% Step 5 : Element wise multiplication with phase vector
% Abdul Gaffar: simpler way of computing phi without needing a for loop
phi = exp(-1i*2*pi*(smooth_data.').*vector/M);                              % Phase shift vector
alignedmatrix = ifft(phi.*(fft(hrrp,[],2)),[],2);                           % Aligning matrix
end

%% ----------------------------------------------------------Autofocus-----------------------------------------------------------------------------------------------------------------------------------------------------------
function [phase_compensated,Power_vector] = autofocus(Aligned_HRRP)
hrrp = Aligned_HRRP;
M = size(hrrp,2);
N = size(hrrp,1);
%% Step 1 : Calculate the range cell with the minimum variance for each HRRP
V = var(hrrp,[],1);                                                        % Here we are calculating the variance for each range bin - 1 refers to the dimension (range bin) and 0 is syntax for M-1
%% Step 2 : Using average power to ensure dominant scatterer criteria  and locating range cell with minimum variance
mag_hrrp = abs(hrrp);                                                      % We are working explicitly with the magnitude
Power = mag_hrrp.^2;                                                       % Calculating the power across each range cell and profile
V_temp = max(V);                                                           % Using a reference Dominant Scatterer
Power_vector  = sum(Power,1);                                              % Calculating the totol power for each range cell
multi = 1;

for k = 1:1:length(V)                                                      % Checking to see if minimum variance cell agrees with dominant scatterer criteria PrangeCell > Pavg
    P_new = Power_vector;
    P_new(:,k) = [];                                                       % Removes power at range cells at index k
    avg_power = (1/(M-1))*(sum(P_new));                                    % Calculating average power across all range cells where k!=m
    if avg_power < 1                                                       % Check to see if signal is contaminated with noise - noise produces a low average power
        multi = 6;                                                         % Scales average power threshold by a calue of 6 if signal is contaminated with noise
    end
    if Power_vector(k)> multi*avg_power
        if V(k)< V_temp
            V_temp = V(k);
            DS = k;                                                        % DS = Dominant Scatterer
        end
    end
end

%% Step 3 : Calculate the phase difference
Phase_shift = angle(hrrp(1,DS));                                           % Calculating phase for reference range cell
phase_array = angle(hrrp(:,DS)) - Phase_shift;                             % Calculating phase difference across all range cells in reference range.
%% Step 4 : Apply Phase shifts to all range profiles
phase_vector = exp(-1i*phase_array);
phase_compensated = hrrp.*repmat(phase_vector,1,M);                        % Apply phase difference to each range profile
end

%% --------------------------------------------------------------BEST ISAR FORMATION----------------------------------------------------------------------------------------------------------------------------------------------
%% STFT Function definition
function [indexes,overlap_matrix] = GenISARImages (y,win, overlap, nfft, PRF,Range_axis,factor)
ylen = size(y,1);           % Signal Length
wlen = length(win);         % window function length should be 1024
CPTWL = roundn(wlen*1/PRF,-2);
%% Calculate the number of frames to be taken, given the signal size and amount of overlap
overlap = wlen*factor;
hop = wlen-overlap;
frames = 1+floor((ylen-wlen)/(hop));
%% Contrast matrix
contrast_matrix = zeros(1,frames);
centre_time_matrix = zeros(1,frames);
%% Executing the STFT to generate multiple ISAR Images
for i = 0:frames-1
    %% STFT
    StartProfile = 1+i*hop;                                 % Abdul Gaffar
    StopProfile = wlen+i*hop;                               % Abdul Gaffar
    CentreProfile = (StartProfile + StopProfile)/2;         % Abdul Gaffar
    CentreProfile_time = roundn(CentreProfile*1/PRF,-2);    % Abdul Gaffar
    centre_time_matrix(1,i+1) = CentreProfile_time;
    Ess = y(StartProfile:StopProfile,:).*kaiser(wlen).*repmat(kaiser(wlen),1,size(y,2));       % windowing of the sampled data that moves 'overlap' samples for respective frame
    Enew = Ess;
    %% Compute ISAR image
    ISAR = fftshift(fft(Enew,[],1),1);
    %% Contrast
    %% ---------------- Removal of 0 Doppler -----------------------------
    B_zdm = round(40*(size(ISAR,1)/PRF))    ;
        No_ZERO_DOPPLER = ISAR;
        % Removal of zero Doppler 
        No_ZERO_DOPPLER((size(ISAR,1)/2)-B_zdm:(size(ISAR,1)/2)+B_zdm,:) = [];
        % Recalculate Contrast 
        B = abs(No_ZERO_DOPPLER).^2;
        C = mean(mean(B,1));
        D = B - C;
        E = sqrt((mean(mean(D.^2,1))));
        Contrast = E/C;

    contrast_matrix(1,i+1) = Contrast;
    center_instance(1,i+1) = CentreProfile_time;
    %% Entropy
    N = B/(sum(sum(B,1)));
    Entropy = -sum(sum(N.*log(N),1));
end
%% Locating peaks of contrast matrix
counter = 3;
increment = 1;
peak_matrix(1:2,:) = 0;
plot(center_instance,contrast_matrix);
ylabel('Contrast of each image instance')
xlabel('Image instances, \tau, (s)');
title('Peak contrast analysis criteria') 
while(1)
    if counter > length(contrast_matrix)-2
        break;
    end
    if (contrast_matrix(1,counter)>contrast_matrix(1,counter-1) && contrast_matrix(1,counter)>contrast_matrix(1,counter+1))
        if (contrast_matrix(1,counter-1)>contrast_matrix(1,counter-2) && contrast_matrix(1,counter+1)>contrast_matrix(1,counter+2))
            peak_matrix(1,increment) = contrast_matrix(1,counter);
            peak_matrix(2,increment) = counter;
            increment = increment +1;
        end
    end
    counter = counter+1;
end
indexes = 0 ;
overlap_matrix = 0;
if peak_matrix(2,1) == 0
    indexes = 0;
    overlap_matrix =0;
else
    cakes = peak_matrix(2,:)
    lmao = 1;
    %% STFT
    for i = 1:length(cakes)
        StartProfile = 1+(cakes(i)-1)*hop;
        StopProfile = wlen+(cakes(i)-1)*hop;
        CentreProfile = (StartProfile + StopProfile)/2;
        CentreProfile_time = roundn(CentreProfile*1/PRF,-2);
        subset = y(StartProfile:StopProfile,:);
        Ess = subset.*kaiser(wlen).*repmat(kaiser(wlen),1,size(y,2));         % windowing of the sampled data that moves 'overlap' samples for respective frame
        ISAR = fftshift(fft(Ess,[],1),1);
        %% --------------------Contrast-----------------------------
        %% Step 1 : Calculate Power
        B = abs(ISAR).^2;
        %% Step 2 : Calculate image intensity
        C = mean(mean(B,1));
        %% Step 3 : Calculate standard deviation of image intensity
        D = B - C;
        E = sqrt((mean(mean(D.^2,1))));
        %% Step 4 : Normalise image intensity and calculate contrast
        new_contrast = E/C;
        %% Entropy
        N = B/(sum(sum(B,1)));
        Entropy = - sum(sum(N.*log(N),1));
        %% --------------------------------IT-------------------------------------------
        %Locate bin with max energy
        [max_down_range_energy, energy_cell_index] = max(sum(20*log10(B)));
        %Select bandwidth B_ZDT
        B_zdt = round(30*(size(B,1)/PRF));
        zero_doppler = size(B,1)/2;
        %Sum negative and positive frequencies
        negative_freqs = sum(20*log10(B(zero_doppler-B_zdt:zero_doppler,energy_cell_index)));
        positive_freqs = sum(20*log10(B(zero_doppler:zero_doppler+B_zdt,energy_cell_index)));
        %Select energy threshold
        E_threshold = 350; %dB;
        %Check to see if image contains both negative and positve frequencies 
        Energy_test = abs(negative_freqs - positive_freqs);
        if Energy_test < E_threshold
            IT = 0;
        else
            IT = 1;
        end
        if IT == 1
            indexes(lmao) = cakes(i)
            overlap_matrix(1,lmao) = cakes(i);
            overlap_matrix(2,lmao) = CentreProfile_time;
            overlap_matrix(3,lmao) = StartProfile;
            overlap_matrix(4,lmao) = StopProfile;
            overlap_matrix(5,lmao) = wlen;
            
            lmao =lmao+1;
            %cakes
        end
    end
end

end
function optimum_matrix = CPTWL_optV3(y,win, overlap, nfft, PRF,Range_axis,indexes,factor,n)
flag = 1;
index = indexes;
ylen = size(y,1);                                                       % Signal Length
Actual_CPTWL = length(win);
wlen = length(win)
overlap = factor*wlen;                                                  % window function length should be 1024
hop = wlen-overlap;
[new_contrast,StopProfile,ISAR,Entropy,CentreProfile_time,original,oentropy,ISAR_new]= image_generation (wlen, hop, 'n',y,index,PRF);
Entropy
%% Check to see if contrast increases or decreases or increases on first iteration
counter = 0;
h_counter = counter;
h = 1;
nlen = wlen + (2*(n-counter)*h);                                        % Increase window size by 2n
overlap = factor*nlen;                                                  % Calculate new overlap to factor the size of the new window length
hop = nlen-overlap;                                                     % Abdul Gaffar
Test = nlen+(index-1)*hop;
if Test > ylen
    wlen = length(win);
    overlap = factor*wlen;
    hop = wlen-overlap;
    [new_contrast,StopProfile,ISAR,Entropy,CentreProfile_time,ocontrast,oentropy,ISAR_new] = image_generation (wlen, hop, 'n',y,index,PRF);
    optimum_matrix(1,flag) = index;
    optimum_matrix(2,flag) = wlen;
    optimum_matrix(3,flag) = ocontrast;
    optimum_matrix(4,flag) = Entropy;
    optimum_matrix(5,flag) = wlen;
    invert = 'n';
end
%Image_print1(PRF,wlen,ISAR,'n',Entropy,new_contrast,Range_axis,CentreProfile_time,lol)
if Test < ylen
    [new_contrast,StopProfile,ISAR,Entropy,CentreProfile_time,second,oentropy,ISAR_new] = image_generation (nlen, hop, 'n',y,index,PRF);
    Entropy
    StopProfile = nlen+(index-1)*hop;
    new_cmatrix = zeros(2,n+1);
    new_cmatrix(1,2) = second;
    new_cmatrix(1,1) = original;
    new_cmatrix(2,1) = wlen;
    new_cmatrix(2,2) = nlen;
    %% Determin if contrast increases or decreases with the addition of 2n to the window length.
    if original < second
        while 1
            if h_counter == n
                break;
            end
            new_cmatrix
            while (new_contrast>new_cmatrix(1,counter+1))
                h = h+1;
                counter = counter +1
                nlen = nlen + (2*(n - h_counter)*h)
                overlap = factor*nlen;
                hop = nlen-overlap;
                Tester = nlen+(index-1)*hop;
                if Tester > ylen
                    break;
                end
                [new_contrast,StopProfile,ISAR,Entropy,CentreProfile_time,ocontrast,oentropy,ISAR_new] = image_generation (nlen, hop, 'n',y,index,PRF);
                new_cmatrix(1,counter+2) = ocontrast
                new_cmatrix(2,counter+2) = nlen
            end
            counter = counter+1;
            h_counter = h_counter+1;
            h = 1;
            nlen = nlen + (2*(n - h_counter)*h);
            overlap = factor*nlen;
            hop = nlen-overlap;
            Tester = nlen+(index-1)*hop;
            if Tester > ylen
                break;
            end
            %% Storing all the contrasts as window length increases
            [new_contrast,StopProfile,ISAR,Entropy,CentreProfile_time,ocontrast,oentropy,ISAR_new] = image_generation (nlen, hop, 'n',y,index,PRF);
            %Image_print1(PRF,wlen,ISAR,'n',Entropy,new_contrast,Range_axis,CentreProfile_time,index);
            new_cmatrix(1,counter+2) = ocontrast
            new_cmatrix(2,counter+2) = nlen
            if new_contrast < new_cmatrix(1,counter+1)
                break
            end           
        end
    else
        h_counter
        nlen = wlen - (2*(n-h_counter)*h);                                            % Increase window size by 2n
        overlap = factor*nlen;                                                  % Calculate new overlap to factor the size of the new window length
        hop = nlen-overlap;                                                     % Abdul Gaffar
        while 1
            if h_counter == n
                break;
            end
            if StopProfile > ylen
                continue;
            end
            if hop<0
                break;
            end
            %% Storing all the contrasts as window length decreases
            [new_contrast,StopProfile,ISAR,Entropy,CentreProfile_time,ocontrast,oentropy,ISAR_new] = image_generation (nlen, hop, 'n',y,index,PRF);
            new_cmatrix(1,counter+3) = ocontrast
            new_cmatrix(2,counter+3) = nlen
            Image_print1(PRF,wlen,ISAR,'n',Entropy,new_contrast,Range_axis,CentreProfile_time,index);
            if new_contrast < new_cmatrix(1,counter+2)
                break
            else
                while (new_contrast>new_cmatrix(1,counter+2))
                    h = h +1;
                    counter = counter +1;
                    nlen = nlen - (2*(n - h_counter)*h);
                    overlap = factor*nlen;
                    hop = nlen-overlap;
                    if StopProfile > ylen
                        continue;
                    end
                    if hop<0
                        break;
                    end
                    [new_contrast,StopProfile,ISAR,Entropy,CentreProfile_time,ocontrast,oentropy,ISAR_new] = image_generation (nlen, hop, 'n',y,index,PRF);
                    new_cmatrix(1,counter+3) = ocontrast
                    new_cmatrix(2,counter+3) = nlen
                    
                end
                counter = counter+1;
                h_counter = h_counter +1;
                h=1;
                nlen = nlen - (2*(n - h_counter)*h);
                overlap = factor*nlen;
                hop = nlen-overlap
            end
        end
    end
    %% Comparison between the original image and the Optimised window version
    wlen = length(win);
    overlap = factor*wlen;
    hop = wlen-overlap;
    [new_contrast,StopProfile,ISAR,Entropy,CentreProfile_time,ocontrast,oentropy,ISAR_new] = image_generation (wlen, hop, 'n',y,index,PRF);
    %Image_print1(PRF,wlen,ISAR,'n',Entropy,new_contrast,Range_axis,CentreProfile_time,index);
    [max_contrast,apple] = max(new_cmatrix(1,:));
    nlen = new_cmatrix(2,apple)
    overlap = factor*nlen
    hop = nlen-overlap
    [new_contrast,StopProfile,ISAR,Entropy,CentreProfile_time,ocontrast,oentropy,ISAR_new] = image_generation (nlen, hop, 'n',y,index,PRF);
    optimum_matrix(1,flag) = index;
    optimum_matrix(2,flag) = nlen;
    optimum_matrix(3,flag) = ocontrast
    optimum_matrix(4,flag) = Entropy
    optimum_matrix(5,flag) = wlen;
    invert = 'n';
    optimum_matrix
end
%Image_print(PRF,nlen,ISAR_new,invert,Entropy,ocontrast,Range_axis,CentreProfile_time,index);


    function Image_print(PRF,wlen,ISAR,invert,Entropy,new_contrast,Range_axis,CentreProfile_time,lol)
        CPTWL = roundn(wlen*1/PRF,-2);
        Contrast = new_contrast;
        f = (-wlen/2:1:(wlen/2-1))*PRF/wlen;
        figure;
        imagesc(Range_axis,f,20*log10(abs(ISAR)));
        title(sprintf("Image %d CPTWL Optimized,Optimum CPTWL = %0.5g s, CenterTime = %0.5g s \n Entropy = %0.5g, Contrast = %0.5g",lol,CPTWL,CentreProfile_time,Entropy, Contrast)); % Abdul Gaffar
        colormap(flipud(gray));
        ylabel('Doppler Frequency (Hz)');
        xlabel('Range (m)');
        colorbar;
        axis xy;
    end
    function Image_print2(PRF,wlen,ISAR,invert,Entropy,new_contrast,Range_axis,CentreProfile_time,lol)
        CPTWL = roundn(wlen*1/PRF,-2);
        Contrast = new_contrast;
        f = (-wlen/2:1:(wlen/2-1))*PRF/wlen;
        figure;
        imagesc(Range_axis,f,20*log10(abs(ISAR)));
        title(sprintf("Image %d CPTWL Optimized,Optimum CPTWL = %0.5g s, CenterTime = %0.5g s \n Entropy = %0.5g, Contrast = %0.5g",lol,CPTWL,CentreProfile_time,Entropy, Contrast)); % Abdul Gaffar
        colormap(flipud(gray));
        ylabel('Doppler Frequency (Hz)');
        xlabel('Range (m)');
        colorbar;
        axis xy;
    end

    function Image_print1(PRF,wlen,ISAR,invert,Entropy,new_contrast,Range_axis,CentreProfile_time,lol)
        CPTWL = roundn(wlen*1/PRF,-2);
        Contrast = new_contrast;
        f = (-wlen/2:1:(wlen/2-1))*PRF/wlen;
        figure;
        imagesc(Range_axis,f,20*log10(abs(ISAR)));
        title(sprintf("Original Image %d,Optimum CPTWL = %0.5g s, CenterTime = %0.5g s \n Entropy = %0.5g, Contrast = %0.5g",lol,CPTWL,CentreProfile_time,Entropy, Contrast)); % Abdul Gaffar
        colormap('jet');
        ylabel('Doppler Frequency (Hz)');
        xlabel('Range (m)');
        colorbar;
        axis xy;
    end

    function [new_contrast,StopProfile,ISAR,Entropy,CentreProfile_time,ocontrast,oentropy,ISAR_new] = image_generation (nlen, hop, fix,y,index,PRF)
        %% STFT
        StartProfile = 1+(index-1)*hop;
        StopProfile = nlen+(index-1)*hop;
        CentreProfile = (StartProfile + StopProfile)/2;
        CentreProfile_time = roundn(CentreProfile*1/PRF,-2);
        subset = y(StartProfile:StopProfile,:);
        Ess = subset.*kaiser(nlen).*repmat(kaiser(nlen),1,size(y,2));         % windowing of the sampled data that moves 'overlap' samples for respective frame
        ISAR = fftshift(fft(Ess,[],1),1);
        %% --------------------Contrast-----------------------------
        %% Step 1 : Calculate Power
        B = abs(ISAR).^2;
        %% Step 2 : Calculate image intensity
        C = mean(mean(B,1));
        %% Step 3 : Calculate standard deviation of image intensity
        D = B - C;
        E = sqrt((mean(mean(D.^2,1))));
        %% Step 4 : Normalise image intensity and calculate contrast
        new_contrast = E/C;
        %% Entropy
        N = B/(sum(sum(B,1)));
        Entropy = - sum(sum(N.*log(N),1));
        %% --------------------------------- Normalisation ----------------------------------------
        % Locate the peak Doppler
        Peak = max(ISAR);
        Peak_Doppler = max(Peak);
        % Normalisation of each Doppler component
        ISAR_new = ISAR./Peak_Doppler;
        % Limiting the dB magnitude to supress spectral ringing
        limit = -30;
        for row = 1:size(ISAR_new,1)
            for column = 1:size(ISAR_new,2)
                if 20*log10(abs(ISAR_new(row,column)))<limit
                    ISAR_new(row,column) = ((10^(limit/20))/(abs(ISAR_new(row,column))))*ISAR_new(row,column);
                end
            end
        end
        %% Contrast
        B = abs(ISAR).^2;
        C = mean(mean(B,1));
        D = B - C;
        E = sqrt((mean(mean(D.^2,1))));
        ocontrast = E/C;
        %% Entropy
        N = B/(sum(sum(B,1)));
        oentropy = -sum(sum(N.*log(N),1));
    end
end
function [new_contrast,StopProfile,ISAR,Entropy,CentreProfile_time,ocontrast,oentropy,ISAR_new] = image_generation (nlen, hop, fix,y,index,PRF)
%% STFT
StartProfile = 1+(index-1)*hop;
StopProfile = nlen+(index-1)*hop;
CentreProfile = (StartProfile + StopProfile)/2;
CentreProfile_time = roundn(CentreProfile*1/PRF,-2);
subset = y(StartProfile:StopProfile,:);
Ess = subset.*kaiser(nlen).*repmat(kaiser(nlen),1,size(y,2));         % windowing of the sampled data that moves 'overlap' samples for respective frame
ISAR = fftshift(fft(Ess,[],1),1);
%% --------------------Contrast-----------------------------
%% Step 1 : Calculate Power
B = abs(ISAR).^2;
%% Step 2 : Calculate image intensity
C = mean(mean(B,1));
%% Step 3 : Calculate standard deviation of image intensity
D = B - C;
E = sqrt((mean(mean(D.^2,1))));
%% Step 4 : Normalise image intensity and calculate contrast
new_contrast = E/C;
%% Entropy
N = B/(sum(sum(B,1)));
Entropy = - sum(sum(N.*log(N),1));
%% --------------------------------- Normalisation ----------------------------------------
% Locate the peak Doppler
Peak = max(ISAR);
Peak_Doppler = max(Peak);
% Normalisation of each Doppler component
ISAR_new = ISAR./Peak_Doppler;
% Limiting the dB magnitude to supress spectral ringing
limit = -30;
for row = 1:size(ISAR_new,1)
    for column = 1:size(ISAR_new,2)
        if 20*log10(abs(ISAR_new(row,column)))<limit
            ISAR_new(row,column) = ((10^(limit/20))/(abs(ISAR_new(row,column))))*ISAR_new(row,column);
        end
    end
end
%% ---------------- Removal of 0 Doppler -----------------------------
%% Zero Doppler Removal
size(ISAR,1);
zero_doppler = size(ISAR,1)/2;
B = abs(ISAR).^2;
% Bandwidth to ignore
B_zdm = round(40*(size(B,1)/PRF))    ;

        B = abs(B).^2;
        C = mean(mean(B,1));
        D = B - C;
        E = sqrt((mean(mean(D.^2,1))));
        ocontrast = E/C;

%% Entropy
N = B/(sum(sum(B,1)));
oentropy = -sum(sum(N.*log(N),1));

end