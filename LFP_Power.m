clear all; close all;

wrkpath = 'I:\MATLAB';
datapath = 'I:\data';
opath    = 'I:\LFP_power';


%% This analysis is based on the burst detection protocoll described by Tinkhauser et al.: The modulatory effect of adaptive deep brain stimulation on beta bursts in Parkinson's disease (DOI: 10.1093/brain/awx010)

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

    load(files_CFA);
    load(files_RFA);
    load([datapath,'\' files_STN(l).name]); 

    rec_name = extractBefore(files_STN(l).name,'_STN');


    % 
    % %% Time-frequency analysis I.Hanning taper, fixed window length0.5
    % RFA data
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmfft';
    cfg.taper        = 'hanning';
    cfg.foi          = 0:0.5:100;                      
    cfg.t_ftimwin    = 3./cfg.foi;
    cfg.toi          =  0:0.05:1;
    cfg.pad            ='nextpow2'
    cfg.trials      = 'all' ;
    cfg.keeptrials  = 'yes';
    TFRhann_RFA = ft_freqanalysis(cfg,RFA_data);

    save([folder 'TFRhann_RFA_', char(rec_name),'_', '.mat'],'TFRhann_RFA' );

    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmfft';
    cfg.taper        = 'hanning';
    cfg.foi          = 0:0.5:100;                      
    cfg.t_ftimwin    = 3./cfg.foi;
    cfg.toi          =  0:0.05:1;
    cfg.pad            ='nextpow2'
    cfg.trials      = 'all' ;
    cfg.keeptrials  = 'yes';
    TFRhann_STN = ft_freqanalysis(cfg,STN_data);

    save([folder 'TFRhann_STN_', char(rec_name), '.mat'],'TFRhann_STN' );


    % CFA data
    cfg              = [];
    cfg.output       = 'pow';
    cfg.method       = 'mtmfft';
    cfg.taper        = 'hanning';
    cfg.foi          = 0:0.5:100;                      
    cfg.t_ftimwin    = 3./cfg.foi;
    cfg.toi          =  0:0.05:1;
    cfg.pad            ='nextpow2'
    cfg.trials      = 'all' ;
    cfg.keeptrials  = 'yes';
    TFRhann_CFA = ft_freqanalysis(cfg,CFA_data);

    save([folder 'TFRhann_CFA_', char(rec_name), '.mat'],'TFRhann_CFA' );


    n_trials =length(RFA_data.trial);


    %% Extract theta band
    %----------------------------------------------------------------------------------------------
    CFA_pow=TFRhann_CFA.powspctrm(:,1:32,:);
    RFA_pow=TFRhann_RFA.powspctrm(:,1:32,:);
    STN_pow=TFRhann_STN.powspctrm(:,1:32,:);
    freq= TFRhann_STN.freq;


    CFA_pow_beta          =  nan(n_trials,32);
    RFA_pow_beta          =  nan(n_trials,32);
    CFA_pow_gamma_high    =  nan(n_trials,32);
    RFA_pow_gamma_high    =  nan(n_trials,32);
    STN_pow_beta          =  nan(n_trials,32);
    STN_pow_gamma_high    =  nan(n_trials,32);
    
    for i = 1:32
        CFA_pow_beta(:,i)       = mean(CFA_pow(:,i,26:50),3);
        RFA_pow_beta(:,i)       = mean(RFA_pow(:,i,26:50),3);
        CFA_pow_gamma_high(:,i) = mean(CFA_pow(:,i,125:165),3);
        RFA_pow_gamma_high(:,i) = mean(RFA_pow(:,i,125:165),3);

        STN_pow_beta(:,i)       = mean(STN_pow(:,i,26:50),3);
        STN_pow_gamma_high(:,i) = mean(STN_pow(:,i,125:165),3);
    end

    save([folder_band   'RFA_pow_' rec_name '.mat'],'freq', 'RFA_pow_beta','RFA_pow_gamma_high')
    save([folder_band   'CFA_pow_' rec_name  '.mat'],'freq', 'CFA_pow_beta','CFA_pow_gamma_high')
    save([folder_band   'STN_pow_' rec_name  '.mat'],'freq', 'STN_pow_beta','STN_pow_gamma_high')

end











