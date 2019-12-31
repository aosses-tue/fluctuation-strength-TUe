function [fluct,fi,outs] = FluctuationStrength_TUe(insig_b, fs, N)
% function [fluct,fi,outs] = FluctuationStrength_TUe(insig_b, fs, N)
%
% 1. Description:
%       Frame-based, off-line implementation of the Fluctuation Strength 
%       algorithm based. The algorithm was adapted from the Roughness model.
% 
% 2. Stand-alone example:
% 
%       % 2.1. Fluctuation strength of an input signal (use audioread to get
%       %      insig and fs.
%           N = round(2*fs); % this is the default
%           [fluct,fi,outs] = FluctuationStrength_TUe(insig, fs, N);
% 
%       2.2. Generating NBN and BBN and then applying the algorithm:
%         Fmodeff = 4;
%         BW = Fmodeff/0.64;
%         Fmod = 0;
%         dur = 2;
%         fs = 44100;
%         fc = 1000;
%         Mdept = 0;
%         SPL = 70;
%         insig = AM_random_noise_BW(fc,BW,SPL,dur,fs,Fmod,Mdept);
%         N = length(insig);
%         % sound(insig,fs); 
%         [fluct fi outs] = FluctuationStrength_TUe(insig, fs, N);
% 
%         BW = 2000;
%         fc = 2000; % then sound from 1000 to 3000 Hz
%         insig2 = AM_random_noise_BW(fc,BW,SPL,dur,fs,Fmod,Mdept);
%         N = length(insig2);
%         % sound(insig2,fs); 
%         [fluct2,fi,outs] = FluctuationStrength_TUe(insig2, fs, N);
% 
% 3. Additional info:
%       Tested cross-platform: Yes
%
% Programmed by Alejandro Osses V./Rodrigo Garcia L., HTI, TU/e, the Netherlands, 2014-2016
% Created on    : 03/02/2016
% Last update on: 25/09/2016 
% Last use on   : 03/11/2016 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model_par = il_Get_fluctuation_strength_params(N,fs);
model_par.debug = 'none';
Chno = model_par.Chno;

% model_par = ef(model_par,'window_type','cosine');

Dz  = 0.5; % Bark resolution. This is not fully automate, so automate it within the inline functions in case you use a value different from 0.5
t_b = ( 1:length(insig_b) )/fs;

overlap = round(0.9*N);
insig_b = buffer(insig_b,N,overlap,'nodelay');
t_b     = buffer(t_b    ,N,overlap,'nodelay');
nFrames = size(insig_b,2);
fluct   = zeros(1,nFrames); % Memory allocation

%% ei = peripheral_stage(insig,fs,N);
% 1. Cosine window:
window = ones(N,1);
attackrelease = 50;
window = il_Do_cos_ramp(window,fs,attackrelease,attackrelease);

for iFrame = 1:nFrames
    
    insig = insig_b(:,iFrame);
    t(iFrame,1) = t_b(1,iFrame);
    % Apply window to frame
    insig = transpose(window .* insig);

    %% 2. Peripheral stages
    % 2.1 Peripheral hearing system (transmission factor a0)
    %     (see 'model_par.a0_in_time' == 1, in _debug version):
    %
    % 4096th order FIR filter:
    insig = il_PeripheralHearingSystem_t(insig,fs); % since 14/05
    
    % 2.2 Excitation patterns
    %     (see model_par.filterbank == 'terhardt', in _debug version):
    dBFS = 100; % unit amplitude corresponds to 100 dB (AMT Toolbox convention)
    ei   = il_TerhardtExcitationPatterns_v3(insig,fs,Chno,dBFS);
    z    = 0.5:Dz:23.5; % Bark
    fc   = bark2hz(z);
    flow = bark2hz(z-.5); flow(1) = 0.01;
    fup  = bark2hz(z+.5);
    BWHz = fup - flow;
        
    %% 3. Modulation depth (estimation)
    [mdept,hBPi] = il_modulation_depths(ei,model_par.Hweight);
    
    bExtraCompensation = 0;
    if bExtraCompensation == 1
        method_vr = 2;
        lvls = rmsdb(ei')+dBFS;
        switch method_vr
            case 0
                val = ones(size(mdept)); % nothing done..
                extraweight = ones(size(mdept));
            case 1
                val = From_dB(lvls);
                val = val/max(val);
                % model_par.cal = (1/0.4374)*model_par.cal; % 0.5693
                model_par.cal = (1/0.5094)*model_par.cal; 
                extraweight = val;
            case 2
                lvls = max(lvls,0);
                lvls(isinf(lvls)) = 0;
                val = phon2sone(lvls);
                val = val/max(val);
                % model_par.cal = (1/0.5002)*model_par.cal; % 0.4978
                model_par.cal = (1/0.5916)*model_par.cal; 
                extraweight = val;
            case 3
                % This does not work
                lvls = max(lvls,0);
                lvls(isinf(lvls)) = 0;
                val = lvls/max(lvls);
                model_par.cal = (1/0.6883)*model_par.cal; % 0.3618
                extraweight = val;
        end
        
        % mdept = mdept.*val;
        % disp('Applying extra compensation as suggested by Duisters 2005, see his equation 5.13...')
    end
    
    %% 4. Cross-correlation coefficient:
    %     (see model_par.dataset == 0, in _debug version)
                    
    % % here cross-correlation is computed before band-pass filtering:
    % Ki = il_cross_correlation(inoutsig); % with hBPi Ki goes down but not as much as 'it should'
    Ki = il_cross_correlation(hBPi);
    [fi_,mdept,kp,gzi] = il_specific_fluctuation(mdept,Ki,model_par);

    kp_fr(iFrame,:)= kp;
    gzi_fr(iFrame,:) = gzi;
    md_fr(iFrame,:) = mdept;
    % fi_ = fi_ .* extraweight; warning('fi is changed')
    fi(iFrame,:)  = (1/Dz)*model_par.cal * fi_;
    fluct(iFrame) = Dz*sum(fi(iFrame,:));
    
end

outs.mdept = md_fr;
outs.kp    = kp_fr;
outs.gzi   = gzi_fr;
outs.cal   = model_par.cal;
outs.p_g   = model_par.p_g;
outs.p_k   = model_par.p_k;
outs.p_m   = model_par.p_m;
outs.t     = t; % time stamp
outs.z     = transpose(z);

disp('')
%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mdept,hBPi,ei] = il_modulation_depths(ei,Hweight)
    
[Chno,Nc] = size(ei);
mdept   = zeros(1,Chno);

ei      = transpose(abs(ei));
h0      = mean(ei);
ei      = ei - repmat(h0,Nc,1);

if ~isnumeric( Hweight )
    % older versions of MATLAB
    hBPi = filter(Hweight,ei); % getting the envelopes
else
    hBPi = sosfilt(Hweight,ei);
end

try 
    % In case LTFAT toolbox is installed (overloads rms from signal processing toolbox)
    hBPrms = rms(hBPi,'dim',1); 
catch
    % uses the default rms calculation from the signal processing toolbox
    hBPrms = rms(hBPi,1);
end
hBPi = transpose(hBPi);

idx = find(h0>0);
mdept(idx) = hBPrms(idx)./h0(idx);

idx = find(h0==0);
mdept(idx) = 0;

idx = find(h0<0);
if length(idx) ~= 0
    error('There is an error in the algorithm')
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ki = il_cross_correlation(hBPi)
[Chno,~] = size(hBPi);

ki = zeros(2,Chno);
for k=1:Chno-2
    try
        cfac = cov(hBPi(k,:),hBPi(k+2,:));
    catch
        error('You do not have the function cov (stats toolbox). Contact me at ale.a.osses@gmail.com to solve this problem');
    end
    den  = diag(cfac);
    den  = sqrt(den*den');

    if den(2,1) > 0 % Pearson correlation
        ki(1,k) = cfac(2,1)/den(2,1);
    elseif den(2,1) == 0
        ki(1,k) = 0;
    else
        warning('Cross correlation factor less than 1')
        ki(1,k) = 0;
    end
end

try
    ki(1,Chno-1) = interp1([0 0.5], ki(1,Chno-3:Chno-2),1,'spline');
    ki(1,Chno  ) = interp1([0 0.5], ki(1,Chno-2:Chno-1),1,'spline');
    ki(2,2     ) = interp1([0.5 1], ki(1,3:4),0,'spline');
    ki(2,1     ) = interp1([0.5 1], ki(1,2:3),0,'spline');
catch
    ki(1,Chno-1) = ki(1,Chno-2);
    ki(1,Chno  ) = ki(1,Chno-2);
    ki(2,1) = ki(1,3);
    ki(2,2) = ki(1,3);
end

ki(2,3:Chno) = ki(1,1:Chno-2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fi,mdept,kp,gzi] = il_specific_fluctuation(mdept,Ki,model_par,dataset)

if nargin <4
    dataset = 0;
end

gzi = model_par.gzi; 
p_g = model_par.p_g;
p_m = model_par.p_m;
p_k = model_par.p_k;

Chno = length(gzi);

fi = zeros(1,Chno);

switch dataset
    case {0,90,99}

        % Version 3: % Improves approximation for FM tones
        thres = 0.7;
        
        idx = find(mdept>thres);
        exceed = mdept(idx)-thres;
        mdept(idx) = thres+(1-thres)*exceed;
        md    = min(mdept,ones(size(mdept)));

    case 1
        md    = min(mdept,ones(size(mdept)));
        md    = mdept-0.1*ones(size(mdept));
        md    = max(mdept,zeros(size(mdept)));
end

kp     = Ki(1,:).*Ki(2,:);
kpsign = (sign(kp));
kp     = abs(kp);

switch dataset
    case {0,90,99}
        fi = gzi.^p_g .* md.^p_m .* (kp.^p_k).*kpsign;
    case 1
        fi = gzi.^p_g*md.^p_m*kp.^p_k;
    otherwise
        error('Dataset does not include the calculation of fi')
end

mdept = md;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function params = il_Get_fluctuation_strength_params(N,fs)
% function params = il_Get_fluctuation_strength_params(N,fs)

dataset = 0; % 0 = Approved version

if nargin < 1
    N = 2*fs; % 2 seconds
end

params         = struct;
params.fs      = fs;
params.N       = N;
params.Chno    = 47;
params.debug   = 'none';

% dataset = 0; % 0 = Approved version
params.window_type = 'cosine';
params.filterbank = 'terhardt'; 
params.p_g     = 1; 
params.p_m     = 1.7;  
params.p_k     = 1.7; % warning('Temporal value') 
params.a0_in_time = 1;
params.a0_in_freq = ~params.a0_in_time;

params.cal     = 0.2490; % on 15/06/2016
params.bIdle   = 1; % v5
%%%        

params.Hweight = Get_Hweight_fluctuation(fs);
params.gzi     = il_Get_gzi_fluctuation(params.Chno);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outsig = il_Do_cos_ramp(insig,fs,attack_ms,release_ms)
%  Applies a cosine ramp with attack and release times given in [ms]

sig_len = length(insig);
r =  cos_ramp(sig_len,fs,attack_ms,release_ms);
try
    outsig = transpose(r).*insig;
catch
    outsig = r.*insig;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function outsig = il_PeripheralHearingSystem_t(insig,fs)
% function outsig = il_PeripheralHearingSystem_t(insig,fs)
% 
% Applies the effect of transmission from free field to the cochlea to a
% given signal. Time domain version.
% 
% Inputs:
% insig: The signal to process. insig has to be a row vector.
% fs: Sampling frequency,
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

K = 2^12; % FIR filter order 

% B = il_calculate_a0(fs,K);
B = il_calculate_a0_idle(fs,K);

outsig = filter(B,1,[insig zeros(1,K/2)]);
outsig = outsig(K/2+1:end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B = il_calculate_a0(fs,N)
% Compensation of the transmission factor from Free-field, taken from
% Fastl2007, Fig. 8.18, page 226

df    = fs/N;
N0    = round(20/df)+1;
Ntop  = round(20e3/df)+1;
qb    = N0:Ntop;
freqs = (qb-1)*df;

Bark = il_Get_Bark(N,qb,freqs);

a0tab = [ % lower slope from middle ear (fig_a0.c, see Figure_Psychoacoustics_tex)
    0       -999
    0.5     -34.7
    1       -23
    1.5     -17
    2       -12.8
    2.5     -10.1
    3       -8
    3.5     -6.4
    4       -5.1
    4.5     -4.2
    5       -3.5
    5.5     -2.9
    6       -2.4
    6.5     -1.9
    7       -1.5
    7.5     -1.1 % 850 Hz
    8       -0.8
    8.5     0
    10      0     % 1.2 kHz
    12      1.15
    13      2.31
    14      3.85
    15      5.62
    16      6.92
    16.5    7.38
    17      6.92  % 3.5 kHz
    18      4.23
    18.5    2.31
    19      0     % 5.4 kHz
    20      -1.43
    21		-2.59
    21.5	-3.57
    22		-5.19
    22.5	-7.41
    23		-11.3
    23.5	-20
    24		-40
    25		-130
    26		-999
];

a0            = zeros(1,N);
a0(qb)        = il_From_dB(interp1(a0tab(:,1),a0tab(:,2),Bark(qb)));
a0(isnan(a0)) = 0;

B = il_create_a0_FIR(freqs,a0(qb),N,fs);

if nargout == 0
    il_create_a0_FIR(freqs,a0(qb),N,fs);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B = il_calculate_a0_idle(fs,N)
% No resonance of the ear canal accounted for.

df    = fs/N;
N0    = round(20/df)+1;
Ntop  = round(20e3/df)+1;
qb    = N0:Ntop;
freqs = (qb-1)*df;

Bark = il_Get_Bark(N,qb,freqs);

a0tab = [
    0       0
    10      0
    19      0
    20      -1.43
    21		-2.59
    21.5	-3.57
    22		-5.19
    22.5	-7.41
    23		-11.3
    23.5	-20
    24		-40
    25		-130
    26		-999
];

a0            = zeros(1,N);
a0(qb)        = il_From_dB(interp1(a0tab(:,1),a0tab(:,2),Bark(qb)));
a0(isnan(a0)) = 0;

B = il_create_a0_FIR(freqs,a0(qb),N,fs);

if nargout == 0
    il_create_a0_FIR(freqs,a0(qb),N,fs);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function B = il_create_a0_FIR(f,a0,N,fs)

f = [0 f fs/2];
a0 = [a0(1) a0 a0(end)];

B = fir2(N,f/(fs/2),a0);

if nargout == 0
    [H1,Fn]=freqz(B,1,N/2);
    
    figure;
    plot(fs/2*Fn/pi, 20*log10(abs([H1])));
    xlabel('Frequency [Hz]')
    legend([num2str(N) ' taps']);
    title('FIR filter to be used as approximation to isolation curve')
    xlim([0 fs/2])
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ei,ei_f,freq] = il_TerhardtExcitationPatterns_v3(insig,fs,Chno,dBFS)
% function [ei,ei_f,freq] = il_TerhardtExcitationPatterns_v3(insig,fs,Chno,dBFS)

corr = dBFS + 3;

dB2calibrate = rmsdb(insig)+dBFS;

% General parameters
params = il_calculate_params(insig,fs);
N01   = params.N01;
freqs = params.freqs;

dfreq = fs/params.N;
freq = dfreq*(1:params.N); % freqs and freq are the same array, but freqs starts at bin N0

% Transforms input signal to frequency domain
insig = il_From_dB(corr)*fft(insig)/params.N; % 3 dB added to adjust the SPL values to be put into slope equations

% Use only samples that fall into the audible range
Lg  = abs(insig(params.qb));
LdB = il_To_dB(Lg);

% Use only components that are above the hearing threshold
whichL = find(LdB > params.MinExcdB);
nL     = length(whichL);

% Steepness of slopes
S1 = -27;			
S2 = zeros(1,nL);
for w = 1:nL;
    steep = -24 - (230 / freqs(whichL(w))) + (0.2 * LdB(whichL(w))); 
    if steep < 0
        S2(w) = steep;
    end
end

whichZ      = zeros(2,nL);
whichZ(1,:)	= floor(2 * params.Barkno(whichL+N01));
whichZ(2,:)	=  ceil(2 * params.Barkno(whichL+N01));

% Calculate slopes from steep values
Slopes = zeros(nL,Chno);
Slopes_dB = nan(nL,Chno);

for l = 1:nL
    Li = LdB(whichL(l));
    zi = params.Barkno(whichL(l)+N01);
    
    for k = 1:whichZ(1,l)
        zk = k * 0.5;
        delta_z = zi - zk;
        Stemp =	(S1 * delta_z) + Li;
        if Stemp > params.MinBf(k)
            Slopes(l,k) = il_From_dB(Stemp);
            Slopes_dB(l,k) = Stemp;
        end
    end

    for k = whichZ(2,l):Chno
        zk = k * 0.5;
        delta_z = zk - zi;
        Stemp = S2(l) * delta_z + Li;
        if Stemp > params.MinBf(k)
            Slopes(l,k) = il_From_dB(Stemp);
            Slopes_dB(l,k) = Stemp;
        end
    end 
end

% Excitation patterns:
%   Each frequency having a level above the absolute threshold is looked at.
%   The contribution of that level (and frequency) onto the other critical
%   band levels is computed and then assigned.
ExcAmp  = zeros(nL,Chno);
ei      = zeros(Chno,params.N);
for i = 1:Chno
    etmp = zeros(1,params.N);
    for l = 1:nL
        N1tmp = whichL(l);
        N2tmp = N1tmp + N01;

        if whichZ(1,l) == i
            ExcAmp(N1tmp,i)	= 1;
        elseif whichZ(2,l) == i
            ExcAmp(N1tmp,i)	= 1;
        elseif whichZ(2,l) > i
            ExcAmp(N1tmp,i) = Slopes(l,i+1)/Lg(N1tmp);
        else % whichZ(1,l) < k
            ExcAmp(N1tmp,i) = Slopes(l,i-1)/Lg(N1tmp);
        end

        etmp(N2tmp) = ExcAmp(N1tmp,i)*insig(N2tmp); % for each level, the level is projected to that in the respective critical band i
    end

    if nargout >= 2
        ei_f(i,:) = il_To_dB(abs(etmp));
    end
    ei(i,:) = 2*params.N*real(ifft(etmp)); % figure; plot(To_dB(abs(etmp)),'x','LineWidth',4); xlim([1950 2050])
end

outsig = sum(ei,1);
gain = dB2calibrate - (rmsdb(outsig)+dBFS);

ei = il_From_dB(gain)*ei;

disp('');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function params = il_calculate_params(x,fs)

params      = struct;
params.N    = length(x);

% Defines audible range indexes and frequencies
df           = fs/params.N;
N0           = round(20/df)+1; % start at 20 Hz
Ntop         = round(20e3/df)+1; % start at 20 kHz
params.N01   = N0-1;
params.qb    = N0:Ntop;
params.freqs = (params.qb-1)*df;

[params.Barkno,Bark_raw] = il_Get_Bark(params.N,params.qb,params.freqs);
% Loudness threshold related parameters
params.MinExcdB = il_calculate_MinExcdB(params.N01,params.qb,params.Barkno);
params.MinBf    = il_calculate_MinBf(params.N01,df,Bark_raw,params.MinExcdB);
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MinExcdB = il_calculate_MinExcdB(N01,qb,Barkno)

HTres = [
    0		130
    0.01    70
    0.17    60
    0.8     30
    1       25
    1.5     20
    2		15
    3.3     10
    4		8.1
    5		6.3
    6		5
    8		3.5
    10		2.5
    12		1.7
    13.3	0
    15		-2.5
    16		-4
    17		-3.7
    18		-1.5
    19		1.4
    20		3.8
    21		5
    22		7.5
    23      15
    24      48
    24.5 	60
    25		130
];

MinExcdB            = zeros(1,length(qb));
MinExcdB(qb-N01)    = interp1(HTres(:,1),HTres(:,2),Barkno(qb));
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MinBf = il_calculate_MinBf(N01,df,Bark,MinExcdB)
    
Cf = round(Bark(2:25,2)'/df)-N01+1;
Bf = round(Bark(1:25,3)'/df)-N01+1;  

zb      = sort([Bf Cf]);
MinBf   = MinExcdB(zb);

function [Bark,Bark_raw] = il_Get_Bark(N,qb,freqs)

Bark_raw = [
    0   0       50      0.5
    1   100     150     1.5
    2   200     250     2.5
    3   300     350     3.5
    4   400     450     4.5
    5   510     570     5.5
    6   630     700     6.5
    7   770     840     7.5
    8   920     1000	8.5
    9   1080	1170	9.5
    10  1270	1370	10.5
    11  1480	1600	11.5
    12  1720	1850	12.5
    13  2000	2150	13.5
    14  2320	2500	14.5
    15  2700	2900	15.5
    16  3150	3400	16.5
    17  3700	4000	17.5
    18  4400	4800	18.5
    19  5300	5800	19.5
    20  6400	7000	20.5
    21  7700	8500	21.5
    22  9500	10500	22.5
    23  12000	13500	23.5
    24  15500   20000   24.5
]; 


Bark_sorted = [  sort([Bark_raw(:,2);Bark_raw(:,3)]),... % frequencies
                 sort([Bark_raw(:,1);Bark_raw(:,4)])];   % Bark

Bark     = zeros(1,N/2+1);
Bark(qb) = interp1(Bark_sorted(:,1),Bark_sorted(:,2),freqs);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gzi = il_Get_gzi_fluctuation(Chno)
% function gzi = il_Get_gzi_fluctuation(Chno)
%
% Returns gzi parameters using the specified number of channels.

Chstep = 0.5;
    
% Hz:   100 250   519   717 926 1084 1255 1465 1571   1972 2730 4189   15550
g0 = [0,  1,  2.5,  4.9,6.5,  8,   9,  10,  11,  11.5,  13,  15,  17.5,   24;
      1,  1,  1  ,  1  ,1  ,  1,   1,   1,   1,   1  ,   1, 0.9,   0.7, 0.5];
g0 = transpose(g0);

gzi = interp1(g0(:,1),g0(:,2),(1:Chno)*Chstep);
gzi(isnan(gzi)) = g0(end,2); % 0
   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gain_dB = il_To_dB(gain)
% function gain_dB = il_To_dB(gain)
% 
% 1. Description:
%   To_dB: Convert voltage gain to decibels.
%   gain_dB = To_dB(gain)

gain_dB = 20 * log10(gain);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gain = il_From_dB(gain_dB,div)
% function gain = il_From_dB(gain_dB,div)
%
% 1. Description:
%       From_dB: Convert decibels to voltage gain (if div = 20, default).
%       gain = From_dB(gain_dB)

if nargin < 2
    div = 20;
end

gain = 10 .^ (gain_dB / div);