% MATLAB function that reads seismic data in SAC format from a particular
% directory (pathfiles), and arranges them for processing by the master 
% script.
function [wavfull, hdrfull]=load_obs_data(pathfiles)



% Path to directory with SAC files (better to also specify a common pattern 
% in their filenames if there are other files or folders that are not part
% of the analysis). To work well, there should be three channels per 
% station, one SAC file per channel:
files = dir(pathfiles);
addpath(files(1).folder)
nst = length(files)/3;

% Organize waveforms and headers (three channels per station):
fullN = cell(nst,2); fullE = cell(nst,2); fullZ = cell(nst,2);
cn = 1; ce = 1; cz = 1;
for i=1:length(files)
    filename = files(i).name;
    
    % Retrieve waveform data and headers for each channel:
    if contains(filename,'N.sac') || contains(filename,'1.sac')
        [wavN, HdrN, ~, ~]=read_sac(filename);
        fullN(cn,:) = {wavN, HdrN};
        cn = cn+1;
    elseif contains(filename,'E.sac') || contains(filename,'2.sac')
        [wavE, HdrE, ~, ~]=read_sac(filename);
        fullE(ce,:) = {wavE, HdrE};
        ce = ce+1;
    elseif contains(filename,'Z.sac')
        [wavZ, HdrZ, ~, ~]=read_sac(filename);
        fullZ(cz,:) = {wavZ, HdrZ};
        cz = cz+1;
    end
end

% Create final arrays for waveforms and headers:
wavfull = cell(nst,1);
hdrfull = cell(nst,1);
for i=1:nst
    wavfull(i) = {[fullN{i,1}, fullE{i,1}, fullZ{i,1}]};
    hdrfull(i) = {[fullN{i,2}, fullE{i,2}, fullZ{i,2}]};
end