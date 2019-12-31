function demo_osses2016a
% function demo_osses2016a
%
% 1. Description:
%       It generates the Fluctuation Strength figures from the paper Osses, 
%       Garcia, Kohlrausch (2016).
%   
% 2. Stand-alone example:
%       demo_osses2016a;
%
% 3. Additional info:
%       Tested cross-platform: Yes
%       Original file name: r20160613_FS_PAK.m
%
% 4. Reference:
%       Osses, A., García, R., & Kohlrausch, A. (2016). Modelling the sensation 
%       of fluctuation strength. Proc Mtgs Acoust, 28(050005), 1–8. 
%       doi: http://dx.doi.org/10.1121/2.0000410
%
% Programmed by Alejandro Osses (ale.a.osses@gmail.com), HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 13/06/2016
% Last update on: 13/06/2016 
% Last use on   : 12/11/2019 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dir_where = [fs_tue_basepath 'auxdata' filesep 'osses2016a' filesep]; 
dir_stim      = [dir_where 'Stimuli'      filesep];
dir_stim_real = [dir_where 'Stimuli-real' filesep];
    
% addpath(dir_model);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;

bDo_AM         = 1;
bDo_FM         = 1;
bDo_AMBBN      = 1;
bDo_realsounds = 1;

%%%
bPlot_fi = 0; % Plots specific fluctuation strength
%%%

FontSize = 18;

h     = [];
hname = [];
N     = 2*44100; % 2 s at 44100 Hz
% dur   = 4; % 4 s, real sounds should be 4-s long (after truncating the sounds)
durAM = 2; % it does not have any sense to keep it in 4 s, since FS is constant
durFM = 2;

Bark_ = 0.5*[1:47];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[filenames, FS_theo] = il_ICA_FS_Create; % gets filenames and theoretical values
FS   = nan(size(FS_theo));
FSmm = nan(size(FS_theo),2);

%         filenames                  FS value from PAK
files = {'Spr_1_1-mono-excerpt.wav', 1.11; ... 
         'Spr_1_2-mono-excerpt.wav', 1.21; ... 
         'Spr_2-mono-excerpt.wav'  , 0.38; ... 
         'M6_2-mono-excerpt.wav'   , 0.21; ... 
         'M5-mono-excerpt.wav'     , 0.56; ... 
         'Tier1-mono-excerpt.wav'  , 1.77; ... 
         'RR2-mono-excerpt.wav'    , 0.02}; ...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx_col = 1;
h = [];

% % 1. Reference again but with loaded cal value (idx = 1):
idx = 1;
[insig,fs]  = audioread([dir_stim filenames{idx}]);
dBSPL(idx)  = rmsdb(insig)+100;

Nlength     = round(durAM*fs);

Nstart = 1;
insig1 = insig(Nstart:Nstart+Nlength-1);

[FS_,fi,outs]  = FluctuationStrength_TUe(insig1,fs,N);

FS(idx,idx_col) = median(FS_);
[xx,idx_max] = max(FS_);
[xx,idx_min] = min(FS_);
FSmm(idx,1:2) = [FS_(idx_min) FS_(idx_max)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

idx = 2:7;
if bDo_AM

    for i = 1:length(idx)
        [insig,fs]  = audioread([dir_stim filenames{idx(i)}]);
        dBSPL(idx(i))  = rmsdb(insig)+100;

        Nlength     = round(durAM*fs);

        Nstart = 1;
        insig1 = insig(Nstart:Nstart+Nlength-1);

        [FS_,fi,outs]  = FluctuationStrength_TUe(insig1,fs,N);
        FS(idx(i),idx_col) = median(FS_);
        if bPlot_fi
            [xx,idx_max] = max(FS_);
            [xx,idx_min] = min(FS_);
            FSmm(idx(i),1:2) = [FS_(idx_min) FS_(idx_max)];

            figure;
            plot(Bark_,fi(idx_max,:)); grid on, hold on
            plot(Bark_,fi(idx_min,:),'r--'); 
            xlabel('Frequency [Bark]')
            xlim([0 max(Bark_)])
            title( filenames{idx(i)} );
        end
    end

    xoffset = 0.05;
    Xvar = [1:length(idx)];
    XvarLabel = [1 2 4 8 16 32];
    figure;
    hl1 = plot( Xvar-xoffset,FS(idx,idx_col) ,'s-','LineWidth',2); grid on, hold on;
    hl2 = plot( Xvar+xoffset,FS_theo(idx) ,'r<-');
    ha = gca;
    set(ha,'XTick',Xvar);
    set(ha,'XTickLabel',XvarLabel);
    set(ha,'FontSize',FontSize);
    xlabel('f_m_o_d [Hz]');
    xlabel('Fluctuation strength [vacil]');
    title('AM-tones')
    h(end+1) = gcf;
    hname{end+1} = 'res-AM-tones';
    xlim([min(Xvar)-0.4 max(Xvar)+0.4])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx = 14:19;
% idx = 17; % fmod = 8 Hz
if bDo_FM
    for i = 1:length(idx)
        [insig,fs]  = audioread([dir_stim filenames{idx(i)}]);
        dBSPL(idx(i))  = rmsdb(insig)+100;

        Nlength     = round(durFM*fs);
        % Nlength = min(Nlength,length(insig));
        if Nlength > length(insig)
            insig = [insig; insig; insig; insig];
        end

        Nstart = 1;
        insig1 = insig(Nstart:Nstart+Nlength-1);

        [FS_,fi,outs]  = FluctuationStrength_TUe(insig1,fs,N);
        FS(idx(i),idx_col) = median(FS_);
        if bPlot_fi
            [xx,idx_max] = max(FS_);
            [xx,idx_min] = min(FS_);
            FSmm(idx(i),1:2) = [FS_(idx_min) FS_(idx_max)];

            figure;
            plot(Bark_,fi(idx_max,:)); grid on, hold on
            plot(Bark_,fi(idx_min,:),'r--'); 
            xlabel('Frequency [Bark]')
            xlim([0 max(Bark_)])
            plot(Bark_,outs.mdept,Bark_,outs.kp,Bark_,outs.gzi);
            title( filenames{idx(i)} )
        end
    end

    xoffset = 0.05;
    Xvar = [1:length(idx)];
    XvarLabel = [1 2 4 8 16 32];
    figure;
    plot( Xvar-xoffset,FS(idx,idx_col) ,'s-','LineWidth',2), grid on, hold on;
    plot( Xvar+xoffset,FS_theo(idx) ,'r<-')
    ha = gca;
    set(ha,'XTick',Xvar);
    set(ha,'XTickLabel',XvarLabel);
    set(ha,'FontSize',FontSize);
    xlabel('f_m_o_d [Hz]')
    ylabel('Fluctuation strength [vacil]')
    title('FM-tones')
    h(end+1) = gcf;
    hname{end+1} = 'res-FM-tones';
    xlim([min(Xvar)-0.4 max(Xvar)+0.4])
end

idx = 8:13;
if bDo_AMBBN
    for i = 1:length(idx)
        [insig,fs]= audioread([dir_stim filenames{idx(i)}]);
        dBSPL(idx(i)) = rmsdb(insig)+100;

        Nlength     = round(durAM*fs);

        Nstart = 1;
        insig1 = insig(Nstart:Nstart+Nlength-1);

        [FS_,fi,outs]  = FluctuationStrength_TUe(insig1,fs,N);
        FS(idx(i),idx_col) = median(FS_);
        [xx,idx_max] = max(FS_);
        [xx,idx_min] = min(FS_);
        FSmm(idx(i),1:2) = [FS_(idx_min) FS_(idx_max)];

        if idx(i) == 11 & bPlot_fi == 1
            figure;
            plot(Bark_,fi(idx_max,:)/0.5,'k-','LineWidth',2); grid on, hold on
            % plot(Bark_,fi(idx_min,:),'r--'); 
            xlabel('Frequency [Bark]')
            xlim([0 max(Bark_)])
            title('AM BBN, f_m_o_d=8 [Hz]');
            fact = 0.05;
            plot(Bark_,fact*outs.mdept(idx_max,:).^outs.p_m,'m-.','LineWidth',2); 
            plot(Bark_,fact*outs.kp(idx_max,:).^outs.p_k ,'b--','LineWidth',2);
            plot(Bark_,fact*outs.gzi(idx_max,:).^outs.p_g,'rx','LineWidth',2);
            legend({'f_i','m_i','k_i','gz_i'},'Location','NorthWest')
            ylabel('Specific fluctuation strength [vacil/Bark]')
            set(gca,'FontSize',FontSize)
            h(end+1) = gcf;
            hname{end+1} = 'AM-BBN-case';
        end

    end
    %%%

    xoffset = 0.05;
    Xvar = [1:length(idx)];
    XvarLabel = [1 2 4 8 16 32];

    figure;
    plot( Xvar-xoffset,FS(idx,idx_col) ,'s-','LineWidth',2), grid on, hold on;
    plot( Xvar+xoffset,FS_theo(idx) ,'r<-')
    ha = gca;
    set(ha,'XTick',Xvar);
    set(ha,'XTickLabel',XvarLabel);
    set(ha,'FontSize',FontSize);
    xlabel('f_m_o_d [Hz]')
    ylabel('Fluctuation strength [vacil]')
    title('AM BBN')
    legend({'Our model','Literature'},'FontSize',FontSize)
    h(end+1) = gcf;
    hname{end+1} = 'res-AMBBN';
    xlim([min(Xvar)-0.4 max(Xvar)+0.8])
end

idxi = 19; % value from last AM-BBN
idx  = idxi;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if bDo_realsounds

    for i = 1:length(files)
        [insig,fs]  = audioread([dir_stim_real files{i,1}]);
        insig       = il_from_dB(-6)*insig;
        SPL(i)      = rmsdb(insig)+100;

        [FS_,fi_,outs]  = FluctuationStrength_TUe(insig,fs,N);
        FS(idx+i,idx_col) = median(FS_);
        if length(FS_)>1
            FS_(end) = []; % excludes the last sample
        end
        [xx,idx_max] = max(FS_);
        [xx,idx_min] = min(FS_);
        FSmm(idx+i,1:2) = [FS_(idx_min) FS_(idx_max)];

        FS_theo(idx+i,1) = files{i,2};

        if i == 6 & bPlot_fi == 1
            figure;
            plot(Bark_,fi_(idx_max,:)/0.5,'k-','LineWidth',2); grid on, hold on
            % plot(Bark_,fi(idx_min,:),'r--'); 
            xlabel('Frequency [Bark]')
            xlim([0 max(Bark_)])
            title('Duck''s quacking')
            plot(Bark_,outs.mdept(idx_max,:).^outs.p_m,'m-.','LineWidth',2); 
            plot(Bark_,outs.kp(idx_max,:).^outs.p_k ,'b--','LineWidth',2);
            plot(Bark_,outs.gzi(idx_max,:).^outs.p_g,'rx','LineWidth',2);
            legend({'f_i','m_i','k_i','gz_i'},'Location','NorthWest')
            ylabel('Specific fluctuation strength [vacil/Bark]');
            set(gca,'FontSize',FontSize)
        end

    end

    idx = [idxi+1:idxi+length(files)];

    BarsL = FS(idx,idx_col)-FSmm(idx,1);
    BarsU = FSmm(idx,2)-FS(idx,idx_col);
    xoffset = 0.11;
    Xvar = [1:length(idx)];
    XvarLabel = [1 2 23 29 31 34 61];
    figure;
    errorbar( Xvar-xoffset,FS(idx,idx_col),BarsL,BarsU,'s','LineWidth',2), grid on, hold on;
    plot( Xvar+xoffset,FS_theo(idx) ,'r<')
    plot(6-xoffset  ,FS(idx(6),idx_col),'s','MarkerFaceColor','b') % asterisk
    plot(6-7*xoffset,FS(idx(6),idx_col),'k*','MarkerSize',12,'LineWidth',2) % asterisk
    ha = gca;
    set(ha,'XTick',Xvar);
    set(ha,'XTickLabel',XvarLabel);
    set(ha,'FontSize',FontSize);
    xlabel('Track Nr.')
    ylabel('Fluctuation strength [vacil]')
    title('Everyday sounds + pink noise')
    legend({'Our model','Literature'},'FontSize',FontSize)
    h(end+1) = gcf;
    hname{end+1} = 'res-everyday-sounds';
    ylim([-0.05 2.45]);
    xlim([min(Xvar)-1 max(Xvar)+1])

    figure;
    hl2 = plot( Xvar+xoffset,FS_theo(idx) ,'r<'); grid on, hold on;
    ha = gca;
    set(ha,'XTick',Xvar);
    set(ha,'XTickLabel',XvarLabel);
    set(ha,'FontSize',FontSize);
    xlabel('Track Nr.')
    ylabel('Fluctuation strength [vacil]')
    title('Everyday sounds + pink noise')
    legend(hl2,'Literature','FontSize',FontSize)
    h(end+1) = gcf;
    hname{end+1} = 'res-everyday-sounds-literature';
    ylim([-0.05 2.45]);
    xlim([min(Xvar)-1 max(Xvar)+1])

end
% idx_col = idx_col+1;

FS
[FS_theo FS FSmm]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EOF

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Inline functions:

function gain = il_from_dB(gain_dB)
% function gain = il_from_dB(gain_dB)

gain = 10 .^ (gain_dB / 20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [filename, FS_theory_all] = il_ICA_FS_Create
% function [filename, FS_theory_all] = il_ICA_FS_Create

FS_theory_all = [];
filename = [];

%%% 1. Reference sound:
filename{end+1} = 'AM-tone-fc-1000_fmod-4_mdept-100-SPL-60-dB.wav';
FS_theory = 1 ;
% [fluct(i),fi,outs] = FluctuationStrength_Garcia(insig(:,i), fs, N);
%%%
FS_theory_all = [FS_theory_all; FS_theory'];

%%% 2. AM tones at 70 dB
fmod      = [1     2     4     8    16     32]; % Hz
FS_theory = [0.39  0.84  1.25  1.3   0.36  0.06] ;
FS_theory_all = [FS_theory_all; FS_theory'];

filename(end+1:end+6) = {   'AM-tone-fc-1000_fmod-1_mdept-100-SPL-70-dB.wav', ...
                            'AM-tone-fc-1000_fmod-2_mdept-100-SPL-70-dB.wav', ...
                            'AM-tone-fc-1000_fmod-4_mdept-100-SPL-70-dB.wav', ...
                            'AM-tone-fc-1000_fmod-8_mdept-100-SPL-70-dB.wav', ...
                            'AM-tone-fc-1000_fmod-16_mdept-100-SPL-70-dB.wav', ...
                            'AM-tone-fc-1000_fmod-32_mdept-100-SPL-70-dB.wav'};

%%% 3. AM BBN at 60 dB
fmod      = [1     2     4     8    16     32]; % Hz
FS_theory = [1.1205  1.5840 1.7978 1.5660 0.4793 0.1418] ;
FS_theory_all = [FS_theory_all; FS_theory'];

filename(end+1:end+6) = {   'randomnoise-Fc-8010_BW-15980_Fmod-1_Mdept-100_SPL-60.wav', ...
                            'randomnoise-Fc-8010_BW-15980_Fmod-2_Mdept-100_SPL-60.wav', ...
                            'randomnoise-Fc-8010_BW-15980_Fmod-4_Mdept-100_SPL-60.wav', ...
                            'randomnoise-Fc-8010_BW-15980_Fmod-8_Mdept-100_SPL-60.wav', ...
                            'randomnoise-Fc-8010_BW-15980_Fmod-16_Mdept-100_SPL-60.wav', ...
                            'randomnoise-Fc-8010_BW-15980_Fmod-32_Mdept-100_SPL-60.wav'};

%%% 4. FM tones at 70 dB
f = 1500; % Hz
deltaf    = 700; % Hz
SPL       = 70; % Hz
fmod      = [1     2     4     8    16     32]; % Hz
FS_theory = [0.85  1.17  2     0.7   0.27  0.02] ;
FS_theory_all = [FS_theory_all; FS_theory'];

filename(end+1:end+6) = {   'FM-tone-fc-1500_fmod-1_deltaf-700-SPL-70-dB.wav', ...
                            'FM-tone-fc-1500_fmod-2_deltaf-700-SPL-70-dB.wav', ...
                            'FM-tone-fc-1500_fmod-4_deltaf-700-SPL-70-dB.wav', ...
                            'FM-tone-fc-1500_fmod-8_deltaf-700-SPL-70-dB.wav', ...
                            'FM-tone-fc-1500_fmod-16_deltaf-700-SPL-70-dB.wav', ...
                            'FM-tone-fc-1500_fmod-32_deltaf-700-SPL-70-dB.wav'};