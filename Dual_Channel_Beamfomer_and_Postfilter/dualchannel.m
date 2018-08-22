% This file performs the dual channel speech enhancement based on GSC beamformer and postfilter.
% Refer to "Analysis of Two-Channel Generalized Sidelobe Canceller (GSC) With Post-Filtering".

clear
clc

%% Load audio file

audioname = 'carmix';
[x, fs0 ] = audioread([audioname,'.wav']);
fs      = 16e3;
x       = resample(x,fs,fs0);
Ntime   = size(x,1);
Nmic    = size(x,2);
if Nmic ~= 2
    error('The number of channels does not equal 2!')
end

%% Set parameters

NFFT        = 256;
f           = 0:fs/NFFT:fs/2;
full        = 0:fs/NFFT:fs-fs/NFFT;
len         = 17e-2;
fragsize    = 256;
overlap     = 0.5;
c           = 340;
lamdaB0     = 1.54;
omegalow    = 1;
omegahigh   = 3;
gamma0      = 4.6;
alpha       = 0.92;
alphas      = 0.80;
alphadconst = 0.85;
KL          = 9;  
KR          = 50;  
phi0        = 0.25;
% b           = [1 3 4 3 1] / 12;
b           = hanning(5).';
b           = b/sum(b);
Gmin        = 10^(-20/20); % -20dB
vmax        = 5;
EPIS        = 1E-10;
angle1      = 0;
angle2      = 0.1;
xia         = 0.98;
ximin       = 10^(-15/10);

%% set window

win_1       = hanning(fragsize).';
win_2       = hanning(NFFT).';
win_pre     = sqrt(win_1);
win_post    = sqrt(win_2);
win_pre     = sqrt(win_pre);
win_post    = sqrt(win_post);
scale_fac   = sqrt(NFFT/sum(win_1));
scale_postfac = 1.0 / sqrt(NFFT/sum(win_2));
win_pre     = win_pre * scale_fac;
win_post    = win_post * scale_postfac;

%% Cutframes

% ptout = zeros(Nmic, Nframe, fragsize);
for m = 1:Nmic
    [ptout(m,:,:), Nframe] = cutframe(x(:,m), fragsize, overlap);
end

%% Process

Z       = zeros(Nmic,NFFT/2+1);
pt      = zeros(Nmic,fragsize);
pxx     = zeros(NFFT,Nmic);
q       = ones(Nframe,NFFT/2+1);
p       = zeros(Nframe,NFFT/2+1);
yxx     = zeros(Nframe,NFFT/2+1);
yt      = zeros(Nframe,NFFT);
alphad  = zeros(Nframe,NFFT/2+1);
v       = zeros(Nframe,NFFT/2+1);
gamma   = ones(Nframe,NFFT/2+1);
gammas  = ones(Nframe,NFFT/2+1);
xi      = zeros(Nframe,NFFT/2+1);
phi     = zeros(Nframe,NFFT/2+1);
SY      = zeros(Nframe,NFFT);
MY      = zeros(Nframe,NFFT);
SU      = zeros(Nframe,NFFT);
MU      = zeros(Nframe,NFFT);
lamdad  = zeros(Nframe,NFFT/2+1);
GH1     = zeros(Nframe,NFFT/2+1);
G       = ones(Nframe,NFFT/2+1);
lamdaBY = zeros(Nframe,NFFT/2+1);
lamdaBU = zeros(Nframe,NFFT/2+1);
omega   = zeros(Nframe,NFFT/2+1);
phi_tilt = zeros(1,Nframe);
YMAT    = zeros(Nframe,NFFT);
UMAT    = zeros(Nframe,NFFT);

deltak  = 2*pi*full * len * sin(angle1*pi/180) / c + angle2*pi/180;
W       = 0.5 * conj([exp(1i*deltak/2); exp(-1i*deltak/2)]);
B       = 0.5 * conj([exp(1i*deltak/2); -exp(-1i*deltak/2)]);
coh     = sinc(2*pi*full*len/c);       % diffused noise field
H       = 1i*imag(exp(1i*deltak).*coh) ./ (1-real(exp(1i*deltak).*coh));

% compute OMEGA_HIGH
% sincwl = min(sinc(2*pi*f*len/c), 0.998);
% OMEGA_HIGH = 0.57 * cot(deltak(1:NFFT/2+1)/2).^2 .* (1-sincwl).^2 ./ max((1-sincwl.*cos(deltak(1:NFFT/2+1))).^2, EPIS);
% OMEGA_HIGH = 0.57 * cot(deltak(1:NFFT/2+1)/2).^2;
% OMEGA_HIGH = 0.57 * sin(2*pi*f*len*sin(0)/c).^2 ./ (sin(deltak(1:NFFT/2+1)/2).^2 .* sin(2*pi*f*len*sin(0)/c-deltak(1:NFFT/2+1)/2).^2);

hbar = waitbar(0,'Please wait...');
for n = 1:Nframe
    waitbar(n/Nframe,hbar);
    
%     % ***************   GSC part   **************************************
%     pt1 = [zeros(1,(NFFT-fragsize)/2), reshape(ptout(1,n,:),1,fragsize), zeros(1,(NFFT-fragsize)/2)];
%     pt2 = [zeros(1,(NFFT-fragsize)/2), reshape(ptout(2,n,:),1,fragsize), zeros(1,(NFFT-fragsize)/2)];
%     pt  = [ pt1 .* win_pre.' ; pt2 .* win_pre.'];
    pt1 = [reshape(ptout(1,n,:),1,fragsize) .* win_pre, zeros(1,NFFT-fragsize)];
    pt2 = [reshape(ptout(2,n,:),1,fragsize) .* win_pre, zeros(1,NFFT-fragsize)];
    pt  = [ pt1 .' pt2.' ];
    pxx = fft(pt,NFFT);
    Z   = pxx.';
    Y   = sum(conj(W) .* Z) - conj(H) .* sum(conj(B) .* Z);
    U   = sum(conj(B) .* Z);

    % **********   Postfilter part   *************************************
    if n == 1 
        SY(n,:)     = abs(Y).^2;
        MY(n,:)     = abs(Y).^2;
        lamdad(n,:) = abs(Y(1:NFFT/2+1)).^2;
        SU(n,:)     = abs(U).^2;
        MU(n,:)     = abs(U).^2;
        GH1(n,:)    = ones(1,NFFT/2+1);
        gamma(n,:)  = ones(1,NFFT/2+1);
        
        param_Y = initialise_parameters(MY(n,:).',fs,'imcra');
        param_U = initialise_parameters(MU(n,:).',fs,'imcra');
        
        alphad(n,:)     = alphadconst + (1-alphadconst) * p(n,:);
        lamdad(n+1,:)   = alphad(n,:) .* lamdad(n,:) + (1-alphad(n,:)) .* abs(Y(1:NFFT/2+1)).^2;
        
        xi(n,:) = xia + (1-xia) * max(gamma(n,:)-1,0);     % initialization of a priori SNR 
    else
        Ymat    = buildmat(Y,length(b));
        Umat    = buildmat(U,length(b));
        SY(n,:) = alphas * SY(n-1,:) + (1-alphas) * b * abs(Ymat).^2;
        SU(n,:) = alphas * SU(n-1,:) + (1-alphas) * b * abs(Umat).^2;
        param_Y = noise_estimation((abs(Y).^2).','imcra',param_Y);
        param_U = noise_estimation((abs(U).^2).','imcra',param_U);
        MY(n,:) = param_Y.noise_ps.';
        MU(n,:) = param_U.noise_ps.';
        
        
        lamdaBY(n,:)    = SY(n,1:NFFT/2+1) ./ max(MY(n,1:NFFT/2+1), EPIS);
        lamdaBU(n,:)    = SU(n,1:NFFT/2+1) ./ max(MU(n,1:NFFT/2+1), EPIS);
        
%         omega(n,:)      = max(SY(n,1:NFFT/2+1) - MY(n,1:NFFT/2+1), 0) ./ max(SU(n,1:NFFT/2+1) - MU(n,1:NFFT/2+1), 0.001*MY(n,1:NFFT/2+1));
        omega(n,:)      = max(SY(n,1:NFFT/2+1) - MY(n,1:NFFT/2+1), 0) ./ max(SU(n,1:NFFT/2+1) - MU(n,1:NFFT/2+1), EPIS);
        
        % decide q according to fig.4
        for k = 1:NFFT/2+1
            if lamdaBY(n,k) > lamdaB0
                if lamdaBU(n,k) > lamdaB0
                    % ***** compute phi, using eq.33 ***********
                    if omega(n,k) <= omegalow
                        phi(n,k) = 0;
                    else
                        if omegalow < omega(n,k) && omega(n,k) <= omegahigh      % OMEGA_HIGH(k)       % replace omegahigh with OMEGA_HIGH(k)
                            phi(n,k) = (omega(n,k) - omegalow) / (omegahigh - omegalow);
%                         if omegalow < omega(n,k) <= OMEGA_HIGH(k)     % OMEGA_HIGH(k)       % replace omegahigh with OMEGA_HIGH(k)
%                             phi(n,k) = (omega(n,k) - omegalow) / (OMEGA_HIGH(k) - omegalow);
                        else
                            phi(n,k) = 1;
                        end
                    end
                    % ******************************************
                else
                    phi(n,k) = 1;
                end
            else
                phi(n,k) = 0;
            end
        end
        
        % compute eq.34
        phi_tilt(n) = sum(phi(n,KL:KR)) / (KR-KL+1);
        
        gammas(n,:) = abs(Y(1:NFFT/2+1)).^2 ./ MY(n,1:NFFT/2+1);
        
        % compute q, using eq.35
        if phi_tilt(n) > phi0
            for k = 1:NFFT/2+1
                if (gammas(n,k) <= 1)  || (phi_tilt(n) <= phi0)
                    q(n,k) = 1;
                else 
                    q(n,k) = max( (gamma0-gammas(n,k))/(gamma0-1) , 1-phi(n,k) );
                end
            end
        else
            q(n,:) = ones(1,NFFT/2+1);
        end
        
        q(n,:) = min(q(n,:),1);
        q(n,:) = max(q(n,:),0);

        % compute xi, using eq.37
        gamma(n,:) = abs(Y(1:NFFT/2+1)).^2 ./ max(lamdad(n,1:NFFT/2+1), EPIS);
 
        xi(n,:) = alpha * GH1(n-1,:).^2 .* gamma(n-1,:) + (1-alpha) * max(gamma(n,:)-1,0);
        xi(n,:) = max(xi(n,:), ximin);
        v(n,:)  = gamma(n,:) .* xi(n,:) ./ (1+xi(n,:));
        
        % compute GH1, using eq.38
        GH1(n,:)        = ones(1,NFFT/2+1);
        index           = find(v(n,:) > vmax);
        GH1(n,index)    = xi(n,index) ./ (1+xi(n,index));
        index           = find(0<v(n,:) & v(n,:)<=vmax);
        GH1(n,index)    = xi(n,index) .* exp(0.5*expint(v(n,index)))./ (1+xi(n,index));        
        
        % compute p, using eq.36
        p(n,:) = 1 ./ (1 + q(n,:) .* (1+xi(n,:)).*exp(-v(n,:)) ./ (1-q(n,:)+EPIS));
        p(n,:) = min(p(n,:), 1);
        p(n,:) = max(p(n,:), 0);
        
        % compute eq.39 eq.40
        alphad(n,:) = alphadconst + (1-alphadconst) * p(n,:);
        if n < Nframe
            lamdad(n+1,:) = alphad(n,:) .* lamdad(n,:) + (1-alphad(n,:)) .* abs(Y(1:NFFT/2+1)).^2;
        end
        
    end
    
    % apply filter and ifft
    G(n,:)      = (GH1(n,:).^p(n,:)) .* (Gmin.^(1-p(n,:)));
    yxx(n,:)    = G(n,:) .* Y(1:NFFT/2+1);
    yt(n,:)     = real(ifft([yxx(n,:),conj(fliplr(yxx(n,2:end-1)))],NFFT));
    yt(n,:)     = yt(n,:) .* win_post; 
    
end
close(hbar)

%% OLA

yout = ola(yt,fragsize,overlap);


%% .wav output

audiowrite([audioname,'_tc_out_',num2str(fs),'_',num2str(NFFT),'_',num2str(fragsize),'.wav'], yout, fs);


