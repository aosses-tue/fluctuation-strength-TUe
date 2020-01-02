function Hweight = Get_Hweight_fluctuation(fs)
% function Hweight = Get_Hweight_fluctuation(fs)
% 
% Returns the Hweight filter.
% 
% Inputs:
% params: Struct specifying filter characteristics.
% fs: Sampling frequency.
% 
% Outputs:
% Hweight: The digital filter.
% 
% Author: Alejandro Osses/Rodrigo Garcia
% Original file name: Get_Hweight_fluctuation2014
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

try
    load(sprintf('Hweight-%.0f-Hz-LP.mat',fs));
    load(sprintf('Hweight-%.0f-Hz-HP.mat',fs));
    Hweight = [Hweight_HP; Hweight_LP];

catch
    
    warning('Make sure you are using MATLAB 2014a or later...')
    
	% Design parameters of band-pass filter
    sf1 = 0.5; % 0.5
    pf1 = 3.1; % 2
    pf2 = 12; % Hz % 8
    sf2 = 20;
    passAtt1 = 17.5;
    passAtt2 = 14;

    Hweight_lp = designfilt(   'lowpassiir', ...
                            'PassbandFrequency'  , pf2, ...
                            'StopbandFrequency'  , sf2, ...
                            'PassbandRipple'      , 3, ...
                            'StopbandAttenuation', passAtt2, ... % 100
                            'SampleRate'          , fs);

    Hweight_hp = designfilt(   'highpassiir', ...
                            'StopbandFrequency'  , sf1, ...
                            'PassbandFrequency'  , pf1, ...
                            'StopbandAttenuation', passAtt1, ...
                            'PassbandRipple'      , 3, ...
                            'SampleRate'          , fs);

    Hweight_HP = Hweight_hp.Coefficients; % second-order sections
    Hweight_LP = Hweight_lp.Coefficients; % second-order sections
    dirout = [Get_TUe_paths('MATLAB') 'Psychoacoustics' delim 'FluctuationStrength_TUe' delim 'private' delim];
    save(sprintf('%sHweight-%.0f-Hz-LP.mat',dirout,fs),'Hweight_LP');
    save(sprintf('%sHweight-%.0f-Hz-HP.mat',dirout,fs),'Hweight_HP');
    
    Hweight = [Hweight_HP; Hweight_LP];

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargout == 0
    N = fs;
    K = N/2;
    hResponse = []; 
    for i = 1:size(Hweight,1); 
        hResponse(:,i) = freqz(Hweight(i,1:3),Hweight(i,4:6),K);
    end
    
    hTot = hResponse(:,1);
    for i = 2:size(Hweight,1)
        hTot = hTot.*hResponse(:,i);
    end
    
    figure;
    plot(f,20*log10(abs(hTot))); 
    xlim([0 50]);
    xlabel('Frequency [Hz]')
    ylabel('Relative amplitude [dB]')
end
