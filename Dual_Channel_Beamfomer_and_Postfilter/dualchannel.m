% This file performs the dual channel speech enhancement based on GSC beamformer and postfilter.
% Refer to "Analysis of Two-Channel Generalized Sidelobe Canceller (GSC) With Post-Filtering".
% If you have any problem, please contact me freely by email: destinygxx@gmail.com.

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

NFFT        = 512;
f           = 0:fs/NFFT:fs/2;
len         = 17e-2;
fragsize    = 320;
overlap     = 0.5;
c           = 340;
lamdaB0     = 1.54;
omegalow    = 1;
omegahigh   = 3;
gamma0      = 4.6;
alpha       = 0.92;
alphas      = 0.80;
alphadconst = 0.85;
KL          = 8;  
KR          = 180;  
phi0        = 0.25;
b           = [1 3 4 3 1] / 12;
% b           = hanning(9).';
% b           = b/sum(b);
Gmin        = 10^(-20/20); % -20dB
vmax        = 5;
EPIS        = 1E-10;
angle1      = 5;
angle2      = 5;
xia         = 0.98;
ximin       = 10^(-15/10);

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
SY      = zeros(Nframe,NFFT/2+1);
MY      = zeros(Nframe,NFFT/2+1);
SU      = zeros(Nframe,NFFT/2+1);
MU      = zeros(Nframe,NFFT/2+1);
lamdad  = zeros(Nframe,NFFT/2+1);
GH1     = zeros(Nframe,NFFT/2+1);
G       = ones(Nframe,NFFT/2+1);
lamdaBY = zeros(Nframe,NFFT/2+1);
lamdaBU = zeros(Nframe,NFFT/2+1);
omega   = zeros(Nframe,NFFT/2+1);
phi_tilt = zeros(1,Nframe);
win     = hanning(fragsize);

deltak  = 2*pi*f * len * sin(angle1*pi/180) / c + angle2*pi/180;
W       = 0.5 * conj([exp(1i*deltak/2); exp(-1i*deltak/2)]);
B       = 0.5 * conj([exp(1i*deltak/2); -exp(-1i*deltak/2)]);
% coh = zeros(1,NFFT/2+1);        % incoherent noise field
% coh = exp(-1i*(2*pi*f*len/c)*sin(pi/2));    % coherent noise field
coh     = sinc(2*pi*f*len/c);       % diffused noise field
H       = 1i*imag(exp(1i*deltak).*coh) ./ (1-real(exp(1i*deltak).*coh));

hbar = waitbar(0,'Please wait...');
for n = 1:Nframe
    waitbar(n/Nframe,hbar);
    
    % ***************   GSC part   **************************************
    pt  = [reshape(ptout(1,n,:),1,fragsize) .* win.'; reshape(ptout(2,n,:),1,fragsize) .* win.'];
    pxx = fft(pt.',NFFT);
    Z   = pxx(1:NFFT/2+1,:).';
    Y   = sum(conj(W) .* Z) - conj(H) .* sum(conj(B) .* Z);
    U   = sum(conj(B) .* Z);
    
    % **********   Postfilter part   *************************************
    if n == 1 
        SY(n,:)     = abs(Y).^2;
        MY(n,:)     = abs(Y).^2;
        lamdad(n,:) = abs(Y).^2;
        SU(n,:)     = abs(U).^2;
        MU(n,:)     = abs(U).^2;
        GH1(n,:)    = ones(1,NFFT/2+1);
        gamma(n,:)  = ones(1,NFFT/2+1);
        
        param_Y = initialise_parameters(MY(n,:).',fs,'imcra');
        param_U = initialise_parameters(MU(n,:).',fs,'imcra');
        
        alphad(n,:)     = alphadconst + (1-alphadconst) * p(n,:);
        lamdad(n+1,:)   = alphad(n,:) .* lamdad(n,:) + (1-alphad(n,:)) .* abs(Y).^2;
        
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
        
        lamdaBY(n,:)    = SY(n,:) ./ max(MY(n,:), EPIS);
        lamdaBU(n,:)    = SU(n,:) ./ max(MU(n,:), EPIS);
        
        omega(n,:)      = max(SY(n,:) - MY(n,:), 0) ./ max(SU(n,:) - MU(n,:), 0.001*MY(n,:));
        
        % decide q according to fig.4
        for k = 1:NFFT/2+1
            if lamdaBY(n,k) > lamdaB0
                if lamdaBU(n,k) > lamdaB0
                    % ***** compute phi, using eq.33 ***********
                    if omega(n,k) <= omegalow
                        phi(n,k) = 0;
                    else
                        if omegalow < omega(n,k) <= omegahigh
                            phi(n,k) = (omega(n,k) - omegalow) / (omegahigh - omegalow);
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
        
        gammas(n,:) = abs(Y).^2 ./ MY(n,:);
        
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
        gamma(n,:) = abs(Y).^2 ./ max(lamdad(n,:), EPIS);
 
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
            lamdad(n+1,:) = alphad(n,:) .* lamdad(n,:) + (1-alphad(n,:)) .* abs(Y).^2;
        end
        
    end
    
    % apply filter and ifft
    G(n,:)      = GH1(n,:).^p(n,:) .* Gmin.^(1-p(n,:));
    yxx(n,:)    = G(n,:) .* Y;
    yt(n,:)     = real(ifft([yxx(n,:),conj(fliplr(yxx(n,2:end-1)))],NFFT));

end
close(hbar)

%% OLA

yout = ola(yt,fragsize,overlap);


%% audio output

audiowrite([audioname,'_tc_out.wav'], yout, fs);

