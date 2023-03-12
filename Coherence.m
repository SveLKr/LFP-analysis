clear all; close all;

wrkpath = 'I:\MATLAB';
datapath = 'I:\data';
opath    = 'I:\Coherence';


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


for l=1:length(files_STN)
    
    % generate path
    files_CFA = fullfile(datapath, extractBefore(files_STN(l).name,'_STN'), '_CFA.mat');
    files_RFA = fullfile(datapath, extractBefore(files_STN(l).name,'_STN'),  '_RFA.mat');
    
    sf=20000; % Sampling frequency
    winlen = 2*sf; % windowlength for calculation of coherence spectra
    
    rec=[extractBefore(files_STN(l).name,'_STN')];
    disp(rec)
    Subject_CFA='_00009CFA';
    Subject_RFA='_00026RFA';
    
        
        folder           = fullfile(opath,rec);
        mkdir(folder);
        folder_RFA_CFA           = fullfile(folder,'RFA_CFA');
        mkdir(folder_RFA_CFA);
        folder_STN_CFA           = fullfile(folder,'STN_CFA');
        mkdir(folder_STN_CFA);
        folder_STN_RFA           = fullfile(folder,'STN_RFA');
        mkdir(folder_STN_RFA);
        

        
        load(files_CFA);
        load(files_RFA);
        load([datapath,'\' files_STN(l).name]); 
        
        for k=1:length(RFA_data_t.trial)
            data_RFA_spon=cell2mat(RFA_data_t.trial(k));
            data_CFA_spon=cell2mat(CFA_data_t.trial(k));
            data_STN_spon=cell2mat(STN_data_t.trial(k));
            winlen=10000;
            fprinf('%s\n',num2str(winlen))
            for i= 1:32
                for  j= 1:32
                    V_CFA_RFA=nc([data_CFA_spon(i,:)' data_RFA_spon(j,:)'],sf,winlen,0,100);
                    V_CFA_RFA_cor(k,:,i,j)      = V_CFA_RFA(:,4);
                    
                    V_STN_RFA=nc([data_STN_spon(i,:)' data_RFA_spon(j,:)'],sf,winlen,0,100);
                    V_STN_RFA_cor(k,:,i,j)      = V_STN_RFA(:,4);
                    
                    V_STN_CFA=nc([data_STN_spon(i,:)' data_CFA_spon(j,:)'],sf,winlen,0,100);
                    V_STN_CFA_cor(k,:,i,j)      = V_STN_CFA(:,4);
                    F_vec=V_STN_RFA(:,1);
                    
                end
                fprinf('%s\n',num2str(i))
            end
        end
        %%
        F_vec=V_STN_RFA(:,1);
        
        %%
        
        V_CFA_RFA_beta  = nan(length(RFA_data_t.trial),32,32);
        V_CFA_RFA_gamma_high = nan(length(RFA_data_t.trial),32,32);
        
        V_STN_RFA_beta  = nan(length(RFA_data_t.trial),32,32);
        V_STN_RFA_gamma_high = nan(length(RFA_data_t.trial),32,32);
        
        V_STN_CFA_beta  = nan(length(RFA_data_t.trial),32,32);
        V_STN_CFA_gamma_high = nan(length(RFA_data_t.trial),32,32);
        
        for k=1:length(RFA_data_t.trial)
            for i = 1:32
                for j=1:32
                    V_CFA_RFA_beta(k,i,j)         = mean(V_CFA_RFA_cor(k,9:16,i,j));         % 16-31 Hz (beta)
                    V_CFA_RFA_gamma_high(k,i,j)   = mean(V_CFA_RFA_cor(k,39:51,i,j));       % 76-100 Hz (high gamma)
                    
                    
                    V_STN_RFA_beta(k,i,j)         = mean(V_STN_RFA_cor(k,9:16,i,j));         % 16-31 Hz (beta)
                    V_STN_RFA_gamma_high(k,i,j)   = mean(V_STN_RFA_cor(k,39:51,i,j));       % 76-100 Hz (high gamma)
                    
                    V_STN_CFA_beta(k,i,j)         = mean(V_STN_CFA_cor(k,9:16,i,j));         % 16-31 Hz (beta)
                    V_STN_CFA_gamma_high(k,i,j)   = mean(V_STN_CFA_cor(k,39:51,i,j));       % 76-100 Hz (high gamma)
                end
            end
        end
        
        % Save data in both freq bands
        save(fullfile(folder, 'F_vec.mat'),'F_vec');
        save(fullfile(folder_RFA_CFA, 'V_CFA_RFA_beta.mat'),'V_CFA_RFA_beta');
        save(fullfile(folder_RFA_CFA, 'V_CFA_RFA_gamma_high.mat'),'V_CFA_RFA_gamma_high');
        
        save(fullfile(folder_STN_RFA, 'V_STN_RFA_beta.mat'),'V_STN_RFA_beta');
        save(fullfile(folder_STN_RFA, 'V_STN_RFA_gamma_high.mat'),'V_STN_RFA_gamma_high');
        
        save(fullfile(folder_STN_CFA, 'V_STN_CFA_beta.mat'),'V_STN_CFA_beta');
        save(fullfile(folder_STN_CFA, 'V_STN_CFA_gamma_high.mat'),'V_STN_CFA_gamma_high');


    
    close all;
    clearvars -except wrkpath datapath opath  files k files_CFA files_RFA files_STN l
    
end
