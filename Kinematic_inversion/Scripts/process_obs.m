clearvars
close all
clc



% Path to functions:
addpath Matlab_functions/
% Path to save results:
save_path='../DATA/';



%% User-defined parameters:
new_sampling=4;   % New sampling rate for traces (Hz)
window=512;       % New trace length in samples
fmin=0.02;        % Minimum frequency (Hz) for pass-band filter
fmax=0.1;         % Maximum frequency (Hz) for pass-band filter



%% Processing:
% Path to directory with SAC files (better to also specify a common pattern 
% in their filenames if there are other files or folders that are not part
% of the analysis). To work well, there should be three channels per 
% station, one SAC file per channel:
files = dir('../DATA/RAW/*.sac');
if length(files) < 3
    disp('Not enough files found. Stopping!')
    return
end
nst = length(files)/3;
disp(['Number of stations: ', num2str(nst)])


% Numerical processing. Loop over files (channels), which is done in
% alphabetival order (E, N, Z):
data = zeros(window,length(files));
k = 1;
for i=1:length(files)
    filename = append(files(i).folder,'/',files(i).name);

    % Retrieve waveform data and headers for each channel:
    [WAV, Hdr, ~, ~]=read_sac(filename);
    
    % Display station name on console:
    if  contains(filename,'E.sac')
        disp(Hdr.KSTNM)
    end
    
    % Remove mean and linear trend from channel:
    st = WAV-mean(WAV);
    st = detrend(st);
    
    % Filter and integrate channel:
    ch = Hdr.KCMPNM;
    if strcmp(ch(2),'H') == 1
        % For broadband instruments:
        [st,~,~] = pssa3(st,1/(Hdr.DELTA),2,fmin,fmax,0,0,0,'Titulo');
    else
        % For acelerometers:
        [~,st,~] = pssa3(st,1/(Hdr.DELTA),2,fmin,fmax,0,0,0,'Titulo');
    end

    % Resample channel:
    st = resample(st,new_sampling,round(1/Hdr.DELTA));
    
    % Trim channel:
    st = st(1:round(window));
    
    % Save all channels in a single matrix:
    data(:,k) = st;
    k = k+1;
end

% Separate channels in different matrices:
dE = zeros(window,nst);
dN = zeros(window,nst);  
dZ = zeros(window,nst);
j = 1;
for i=1:3:length(data(1,:))
    East = data(:,i);
    North = data(:,i+1);
    Vert = data(:,i+2);
    
    dE(:,j) = East;
    dN(:,j) = North;
    dZ(:,j) = Vert;
    j = j+1;
end

% Reshape and concatenate channels:
E = zeros(window*nst,1); 
N = zeros(window*nst,1); 
Z = zeros(window*nst,1);
for i=1:length(dE(1,:))
    a1 = dE(:,i);
    a2 = dN(:,i);
    a3 = dZ(:,i);
    
    E(window*(i-1)+1:window*i) = a1;
    N(window*(i-1)+1:window*i) = a2;
    Z(window*(i-1)+1:window*i) = a3;
end

% Save results is ASCII text files:
save(fullfile(save_path,'real_disp_x'), 'N', '-ascii')
save(fullfile(save_path,'real_disp_y'), 'E', '-ascii')
save(fullfile(save_path,'real_disp_z'), 'Z', '-ascii')



%% Plotting:
scrsz = [1 1 1366 768];
fig=figure('Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/1.5]);
subplot(3,1,1)
set(subplot(3,1,1), 'Position', [0.05, 0.71, 0.92, 0.24])
plot(N); xlim([0 length(N)]) 
title('Concatenation for inversion'); set(gca,'xticklabel',[])
legend('North'); grid on
subplot(3,1,2)
set(subplot(3,1,2), 'Position', [0.05, 0.40, 0.92, 0.24])
plot(E); xlim([0 length(E)]) 
legend('East'); grid on; set(gca,'xticklabel',[])
subplot(3,1,3)
set(subplot(3,1,3), 'Position', [0.05, 0.09, 0.92, 0.24])
plot(Z); xlim([0 length(Z)]) 
legend('Vertical'); grid on
xlabel('Samples')
