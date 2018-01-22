% FLEET steady slice profiles
% Algorithm: 
%   1) Design first pulse for starting flip angle, calculate its Mxy and
%      its Mz profiles
%   2) Use residual Mz and first pulse's Mxy profile to calculate target
%      beta profile for current pulse
%   3) Design current pulse, calculate Mz profile after it
%   4) Go to (2) until last pulse. 

% How much B1 sensitivity is there, couple with B0? Is it worse than conventional VFA?
% Will want to try SMS. Blips can push echoes apart. Could set it up to be
% roughly uniform across a range of B1. Also they get a navigator and it
% could be OK to change last flip angle for example to avoid overflip, so
% long as the drop in signal is consistent across B1 amplitudes. 

addpath util/

T1g = 1500; % ms, ~T1 of gray matter at 3T
T1w = 750; % ms, ~T1 of white matter at 3T
T1 = Inf; % which T1 to use in pulse design. 
% set T1 = Inf; to ignore longitudinal relaxation
TRseg = 2000/30; % time between segment excitations

% parameters
Nseg = 3; % number of segments
tb = 4; % time-bandwidth product
N = 2000; % number of time points in the first segment pulse
zPadFact = 4;   % zero-pad factor; keep high to capture increasing spatial
                % frequencies in the slice profile with each segment
cancelAlphaPhs = true;  % cancel the alpha phase using the beta polynomial
                        % (no reason to change this)
winTruncFact = 1.5; % widening factor compared to first pulse; varying
                    % this between 1 and zPadFact will trade off increase
                    % in pulse duration versus slice profile consistency
                    

% get initial flip angle given this number of segments. 
flip = zeros(Nseg,1);flip(end) = 90;
for jj = Nseg-1:-1:1
    flip(jj) = atand(sind(flip(jj+1)));
end

b = zeros(zPadFact*N,Nseg); % beta polynomials
rf = zeros(zPadFact*N,Nseg); % rf pulses
Mz = ones(zPadFact*N,1); % Mz profile after each pulse
Mxy = zeros(zPadFact*N,Nseg); % slice profiles

% design first RF pulse; Mz is relaxed here
[~,b(zPadFact*N/2-N/2+1:zPadFact*N/2+N/2,1)] = dzrf(N,tb,'st','ls',0.0001,0.01);
% shift beta so that RF is centered, and we get no phase in beta response
B = ft(b(:,1));
B = B.*exp(-1i*2*pi/(N*zPadFact)*1/2*(-N*zPadFact/2:N*zPadFact/2-1)');
b(:,1) = ift(B);
b(:,1) = b(:,1)*sind(flip(1)/2); % scale to first flip in passband
% get RF from this centered + scaled beta polynomial
rf(:,1) = b2rf(b(:,1));

if cancelAlphaPhs
    % cancel a phase by absorbing into b
    % Note that we don't have to do this later since it is done here
    a = b2a(b(:,1));
    b(:,1) = ifft(fft(b(:,1)).*exp(1i*angle(fft((a(:))))));
    rf(:,1) = b2rf(b(:,1));
end

% get the min-phase alpha and its response
a = b2a(b(:,1));
A = ft(a);

% calculate the beta filter response
B = ft(b(:,1));

% use them to get Mxy
Mxy(:,1) = 2*conj(A(:)).*B; % this is the magnetization profile we want all
                            % other pulses to produce

% Amplitude of next pulse's Mxy profile will be |Mz*2*a*b| = |Mz*2*sqrt(1-abs(B).^2)*B|. 
% If we set this = |Mxy_1|, we can solve for |B| via solving quadratic equation
% 4*Mz^2*(1-B^2)*B^2 = |Mxy_1|^2. 
% Subsequently solve for |A|, and get phase of A via min-phase, and 
% then get phase of B by dividing phase of A from first pulse's Mxy phase.
for jj = 2:Nseg
    
    % calculate Mz profile after previous pulse
    Mz = Mz.*(1-2*abs(B).^2)*exp(-TRseg/T1)+(1-exp(-TRseg/T1));
    
    % set up quadratic equation to get |B|
    cq = -abs(Mxy(:,1)).^2;
    bq = 4*Mz.^2;
    aq = -4*Mz.^2;
    Bmag = sqrt((-bq+sqrt(bq.^2-4*aq.*cq))./(2*aq));
    %Bmag2 = (-bq-sqrt(bq.^2-4*aq.*cq))./(2*aq);
    % get A - easier to get complex A than complex B since |A| is
    % determined by |B|, and phase is gotten by min-phase relationship
    A = ft(b2a(ift(Bmag))); % Phase of B doesn't matter here since only profile mag is used by b2a
    % trick: now we can get complex B from ratio of Mxy and A
    B = Mxy(:,1)./(2*conj(A(:)).*Mz);
    Mxy(:,jj) = Mz.*(2*conj(A(:)).*B);
    
    % get polynomial
    b(:,jj) = ift(B);
    
    if winTruncFact < zPadFact % then we need to apply windowing
        winLen = (winTruncFact-1)*N;
        Npad = N*zPadFact - winTruncFact*N;
        % TODO: Replace with a tukeywin (built-in MATLAB func)? 
        window = blackman((winTruncFact-1)*N);
        % split it in half; stick N zeros in the middle
        window = [window(1:winLen/2);ones(N,1);window(winLen/2+1:end)];
        window = [zeros(Npad/2,1);window;zeros(Npad/2,1)];
        b(:,jj) = b(:,jj).*window;
        % recalculate B and Mxy
        B = ft(b(:,jj));
        A = ft(b2a(b(:,jj))); % Phase of B doesn't matter here since only profile mag is used by b2a
        Mxy(:,jj) = Mz.*(2*conj(A(:)).*B);
    end
    
    % get RF
    rf(:,jj) = b2rf(b(:,jj)).';
        
end

Mxy = Mxy.*repmat(exp(1i*2*pi/N*(zPadFact*N/2)*(-N/2:1/zPadFact:N/2-1/zPadFact)'),[1 Nseg]);

% also look at case with pulse 1 scaled to different flips
Mz = ones(zPadFact*N,1);
Mxy_sameRF = zeros(zPadFact*N,Nseg);
[A,B] = abr(-1i*[rf(:,1);0],[ones(1,zPadFact*N) -0.5]*2*pi/(zPadFact*N),...
    -N*zPadFact/2:1:N*zPadFact/2-1);
Mxy_sameRF(:,1) = 2*conj(A).*B;
for jj = 2:Nseg
    
    % calculate Mz profile after previous pulse
    Mz = Mz.*(1-2*abs(B).^2)*exp(-TRseg/T1)+(1-exp(-TRseg/T1));
    [A,B] = abr(-1i*[rf(:,1);0]*sind(flip(jj)/2)/sind(flip(1)/2),...
        [ones(1,zPadFact*N) -0.5]*2*pi/(zPadFact*N),...
        -N*zPadFact/2:1:N*zPadFact/2-1);
    Mxy_sameRF(:,jj) = 2*Mz.*conj(A).*B;
        
end
Mxy_sameRF = Mxy_sameRF.*...
    repmat(exp(1i*2*pi/N*(zPadFact*N/2)*(-N/2:1/zPadFact:N/2-1/zPadFact)'),[1 Nseg]);

% truncate the RF
if winTruncFact < zPadFact % then we need to truncate
    pulseLen = winTruncFact*N;
    Npad = N*zPadFact - pulseLen;
    rf = rf(Npad/2+1:Npad/2+pulseLen,:);
end

% compare stability of summed signals 
figure;
subplot(2,2,1)
plot(real(rf));
hold on
plot(imag(rf));
title 'RF pulses'
ylabel 'radians'
xlabel 'time index'

subplot(2,2,2)
plot(abs(sum(Mxy)));
hold on
plot(abs(sum(Mxy_sameRF)));
legend('SLR-VFA','Conventional VFA')
xlabel 'Segment'
ylabel '|\int M_{xy}| (a.u.)'
title 'Amplitude of integrated slice profile'

subplot(2,2,3)
plot(real(Mxy));
hold on
plot(imag(Mxy));
title 'SLR-VFA Mxy Profiles'
ylabel '/M_0'
xlabel 'frequency index'

subplot(2,2,4)
plot(real(Mxy_sameRF));
hold on
plot(imag(Mxy_sameRF));
title 'Conventional VFA Mxy Profiles (First pulse scaled)'
ylabel '/M_0'
xlabel 'frequency index'

