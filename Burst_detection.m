clear all; close all;

wrkpath = 'I:\MATLAB';
datapath = 'I:\data';
opath    = 'I:\Burst';


f1 = 13; % start frequency
f2 = 30; % stop frequency
%%

cd(wrkpath)
addpath(genpath(fullfile(wrkpath,'matlab_code')))
addpath(genpath(fullfile(wrkpath,'fieldtrip-20180619')));
ft_defaults;



%--------------------------------------------------------------------------
% go to folder where data is
%--------------------------------------------------------------------------

files_STN = dir(fullfile(datapath,'*STN.mat'));
if ~exist(opath,'dir'), mkdir(opath); end

%%
for l=1:length(files_STN)
    
    % generate path
    files_CFA = fullfile(datapath, extractBefore(files_STN(l).name,'_STN'), '_CFA.mat');
    files_RFA = fullfile(datapath, extractBefore(files_STN(l).name,'_STN'),  '_RFA.mat');

    CFA_dat=load(files_CFA);
    RFA_dat=load(files_RFA);
    STN_dat=load([datapath,'\' files_STN(l).name]); 

    rec_name = extractBefore(files_STN(l).name,'_STN');

    %% downsample data to 300 Hz
    cfg = []; 
    cfg.resamplefs      =300;
    cfg.detrend         ='yes';
    cfg.trials          = 'all'
     [RFA_data_300]  = ft_resampledata(cfg, RFA_dat);   
     [CFA_data_300]  = ft_resampledata(cfg, CFA_dat);
     [STN_data_300]  = ft_resampledata(cfg, STN_dat);
     
     

     clear STN_dat CFA_dat RFA_dat

%--------------------------------------------------------------------------
%% Signal decomposing using Wavelet transformation
%--------------------------------------------------------------------------
    pectrum_STN=cell(length(RFA_data_300.trial),1);
    spectrum_RFA=cell(length(RFA_data_300.trial),1);
    spectrum_CFA=cell(length(RFA_data_300.trial),1);


    for t=1:length(RFA_data_300.trial)
    [spectrum_STN_t,freqoi,timeoi] = ft_specest_wavelet(cell2mat(STN_data_300.trial(t)),cell2mat(STN_data_300.time(t)),'freqoi',f1:1:f2,'timeoi', 0:0.01:5,'width', 10, 'gwidth',5);
    [spectrum_CFA_t,freqoi,timeoi] = ft_specest_wavelet(cell2mat(CFA_data_300.trial(t)),cell2mat(CFA_data_300.time(t)),'freqoi',f1:1:f2,'timeoi', 0:0.01:5,'width', 10, 'gwidth',5);
    [spectrum_RFA_t,freqoi,timeoi] = ft_specest_wavelet(cell2mat(RFA_data_300.trial(t)),cell2mat(RFA_data_300.time(t)),'freqoi',f1:1:f2,'timeoi', 0:0.01:5,'width', 10, 'gwidth',5);

    spectrum_STN{t,1}=spectrum_STN_t;
    spectrum_RFA{t,1}=spectrum_RFA_t;
    spectrum_CFA{t,1}=spectrum_CFA_t;
    end

    
%--------------------------------------------------------------------------
%% finde peek in frequency profil
%--------------------------------------------------------------------------

    for j=1:32
        for t=1:length(spectrum_CFA)
             RFA_1ch_trials(t,:,:)=squeeze(real(spectrum_RFA{t,1}(j,:,:)));
             CFA_1ch_trials(t,:,:)=squeeze(real(spectrum_CFA{t,1}(j,:,:)));
             STN_1ch_trials(t,:,:)=squeeze(real(spectrum_STN{t,1}(j,:,:)));
        end


    freq_mean_RFA_r(j,:)=nanmean(squeeze(nanmean(RFA_1ch_trials(:,:,:),1)),2)';
    freq_mean_CFA_r(j,:)=nanmean(squeeze(nanmean(CFA_1ch_trials(:,:,:),1)),2)';
    freq_mean_STN_r(j,:)=nanmean(squeeze(nanmean(STN_1ch_trials(:,:,:),1)),2)';


    [Max_STN,pmax_STN(j)]=find(freq_mean_STN_r(j,:)==max(freq_mean_STN_r(j,:)));
    f_max_STN(j)=freqoi(pmax_STN(j));

    [Max_RFA,pmax_RFA(j)]=find(freq_mean_RFA_r(j,:)==max(freq_mean_RFA_r(j,:)));
    f_max_RFA(j)=freqoi(pmax_RFA(j));

    [Max_CFA,pmax_CFA(j)]=find(freq_mean_CFA_r(j,:)==max(freq_mean_CFA_r(j,:)));
    f_max_CFA(j)=freqoi(pmax_CFA(j));

    end

    save([opath 'f_max' '.mat'],'-v7.3','f_max_STN','f_max_RFA','f_max_CFA');

%--------------------------------------------------------------------------
%% rectify and smooth data
%--------------------------------------------------------------------------
    envHigh_STN=NaN(32,18,length(spectrum_CFA_t));
    envHigh_RFA=NaN(32,18,length(spectrum_CFA_t));
    envHigh_CFA=NaN(32,18,length(spectrum_CFA_t));

    envLow_CFA=NaN(1,length(spectrum_CFA_t));
    envLow_RFA=NaN(1,length(spectrum_CFA_t));
    envLow_STN=NaN(1,length(spectrum_CFA_t));

    for i=1:32
        for t=1:length(spectrum_CFA)
            STN_dat_nan=real(squeeze(spectrum_STN{t,1}(i,pmax_STN(i),:)));
            STN_dat_t=STN_dat_nan(~isnan(STN_dat_nan));

            RFA_dat_nan=real(squeeze(spectrum_RFA{t,1}(i,pmax_RFA(i),:)));
            RFA_dat_t=RFA_dat_nan(~isnan(RFA_dat_nan));

            CFA_dat_nan=real(squeeze(spectrum_CFA{t,1}(i,pmax_CFA(i),:)));
            CFA_dat_t=CFA_dat_nan(~isnan(CFA_dat_nan));

            [envHigh_STN(i,t,1:length(STN_dat_t)), envLow_STN(1:length(STN_dat_t))] = envelope(STN_dat_t,20,'peak');
            [envHigh_RFA(i,t,1:length(RFA_dat_t)), envLow_RFA(1:length(RFA_dat_t))] = envelope(RFA_dat_t,20,'peak');
            [envHigh_CFA(i,t,1:length(CFA_dat_t)), envLow_CFA(1:length(CFA_dat_t))] = envelope(CFA_dat_t,20,'peak');

            STN_75_percentile(i,t)=prctile(envHigh_STN(i,t,:),75);
            RFA_75_percentile(i,t)=prctile(envHigh_RFA(i,t,:),75);
            CFA_75_percentile(i,t)=prctile(envHigh_CFA(i,t,:),75);
        end
    end
    save([opath 'envelope' '.mat'],'-v7.3','envHigh_STN','envHigh_RFA','envHigh_CFA');
    save([opath 'percentile' '.mat'],'-v7.3','STN_75_percentile','RFA_75_percentile','CFA_75_percentile');


%--------------------------------------------------------------------------
%% detect sections over threshold
%--------------------------------------------------------------------------
    STN_burst=[];
    RFA_burst=[];
    CFA_burst=[];
    
    STN_burst_amplitude=[];
    RFA_burst_amplitude=[];
    CFA_burst_amplitude=[];

    for j=1:32
        for t=1:length(spectrum_CFA)
            for i=1: length(envHigh_STN(1,1,:))
                if envHigh_STN(j,t,i)>= STN_75_percentile(j,t)
                    STN_burst(j,t,i)=10;  % if envelope of signal is higher than 75 percentile, save 10 (time in ms)
                    STN_burst_amplitude(j,t,i)=envHigh_STN(j,t,i);
                else
                    STN_burst(j,t,i)=0;
                    STN_burst_amplitude(j,t,i)=0;
                end

                if envHigh_RFA(j,t,i)>= RFA_75_percentile(j,t)
                    RFA_burst(j,t,i)=10;
                    RFA_burst_amplitude(j,t,i)=envHigh_RFA(j,t,i);
                else
                    RFA_burst(j,t,i)=0;
                    RFA_burst_amplitude(j,t,i)=0;
                end

                if envHigh_CFA(j,t,i)>= CFA_75_percentile(j,t)
                    CFA_burst(j,t,i)=10;
                    CFA_burst_amplitude(j,t,i)=envHigh_CFA(j,t,i);
                else
                    CFA_burst(j,t,i)=0;
                    CFA_burst_amplitude(j,t,i)=0;
                end
            end
        end
    end

    save([opath 'burst_0_1' '.mat'],'-v7.3','CFA_burst','RFA_burst','STN_burst');
    save([opath 'burst_amplitude_raw' '.mat'],'-v7.3','CFA_burst_amplitude','RFA_burst_amplitude','STN_burst_amplitude');

%--------------------------------------------------------------------------
%% analyse length/mean amplitude of beta bursts
%--------------------------------------------------------------------------

    %calculate length of STN bursts
    STN_burst_length=[];
    STN_burst_length_n=[];
    STN_burst_amp_total=[];
    for j=1:32
        for t=1:length(spectrum_STN)
            k=1;
            STN_burst_length(j,t,k)=0; 
            STN_burst_amp_total(j,t,k)=0;
            STN_burst_length_n(j,t,k)=0;
            for i=1:(length(STN_burst(j,t,:))-1)
                if STN_burst(j,t,i)>0
                    STN_burst_length(j,t,k)=STN_burst_length(j,t,k)+STN_burst(j,t,i);
                    STN_burst_amp_total(j,t,k)=STN_burst_amp_total(j,t,k)+STN_burst_amplitude(j,t,i);
                    STN_burst_length_n(j,t,k)=STN_burst_length_n(j,t,k)+1;
                elseif STN_burst(j,t,i)== 0 && STN_burst(j,t,i+1)> 0;
                    k=k+1;
                    STN_burst_length(j,t,k)=0;
                    STN_burst_amp_total(j,t,k)=0;
                    STN_burst_length_n(j,t,k)=0;
                end
            end
        end
    end
  
STN_burst_amp_mean=STN_burst_amp_total./STN_burst_length_n;




 %%   
   %calculate length of CFA bursts
    CFA_burst_length=[];
    CFA_burst_length_n=[];
    CFA_burst_amp_total=[];
    for j=1:32
        for t=1:length(spectrum_CFA)
            k=1;
            CFA_burst_length(j,t,k)=0; 
            CFA_burst_amp_total(j,t,k)=0;
            CFA_burst_length_n(j,t,k)=0;
            for i=1:(length(CFA_burst(j,t,:))-1)
                if CFA_burst(j,t,i)>0
                    CFA_burst_length(j,t,k)=CFA_burst_length(j,t,k)+CFA_burst(j,t,i);
                    CFA_burst_amp_total(j,t,k)=CFA_burst_amp_total(j,t,k)+CFA_burst_amplitude(j,t,i);
                    CFA_burst_length_n(j,t,k)=CFA_burst_length_n(j,t,k)+1;
                elseif CFA_burst(j,t,i)== 0 && CFA_burst(j,t,i+1)> 0;
                    k=k+1;
                    CFA_burst_length(j,t,k)=0;
                    CFA_burst_amp_total(j,t,k)=0;
                    CFA_burst_length_n(j,t,k)=0;
                end
            end
        end
    end
  
    CFA_burst_amp_mean=CFA_burst_amp_total./CFA_burst_length_n;

   %calculate length of RFA bursts
    RFA_burst_length=[];
    RFA_burst_length_n=[];
    RFA_burst_amp_total=[];
    
    for j=1:32
        for t=1:length(spectrum_RFA)
            k=1;
            RFA_burst_length(j,t,k)=0; 
            RFA_burst_amp_total(j,t,k)=0;
            RFA_burst_length_n(j,t,k)=0;
            for i=1:(length(RFA_burst(j,t,:))-1)
                if RFA_burst(j,t,i)>0
                    RFA_burst_length(j,t,k)=RFA_burst_length(j,t,k)+RFA_burst(j,t,i);
                    RFA_burst_amp_total(j,t,k)=RFA_burst_amp_total(j,t,k)+RFA_burst_amplitude(j,t,i);
                    RFA_burst_length_n(j,t,k)=RFA_burst_length_n(j,t,k)+1;
                elseif RFA_burst(j,t,i)== 0 && RFA_burst(j,t,i+1)> 0;
                    k=k+1;
                    RFA_burst_length(j,t,k)=0;
                    RFA_burst_amp_total(j,t,k)=0;
                    RFA_burst_length_n(j,t,k)=0;
                end
            end
        end
    end
  
    RFA_burst_amp_mean=RFA_burst_amp_total./RFA_burst_length_n;

    RFA_burst_length(RFA_burst_length==0) = NaN;
    CFA_burst_length(CFA_burst_length==0) = NaN;
    STN_burst_length(STN_burst_length==0) = NaN;
    
    
    save([opath 'burst_length_CFA' '.mat'],'-v7.3','CFA_burst_length');
    save([opath 'burst_length_RFA' '.mat'],'-v7.3','RFA_burst_length');
    save([opath 'burst_length_STN' '.mat'],'-v7.3','STN_burst_length');

    save([opath 'burst_amplitude_CFA' '.mat'],'-v7.3','CFA_burst_amp_mean');
    save([opath 'burst_amplitude_RFA' '.mat'],'-v7.3','RFA_burst_amp_mean');
    save([opath 'burst_amplitude_STN' '.mat'],'-v7.3','STN_burst_amp_mean');


end 




