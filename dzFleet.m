function [rf, Mxy, rfRef] = dzFleet(Nseg, tb, Npts, zPadFact, winFact, ...
    aphs, plotAll, T1, TRseg, useMz, seSeq, tbRef)
% DZFLEET
% Generate RF pulses for steady slice profiles in FLEET
%
%   Input Arguments
%
% Nseg  The number of EPI segments (i.e. shots).
%   Optional Input Arguments
% tb        The time-bandwidth product. Default = 4
% Npts      The number of time points in the first segment pulse. Default = 2000
% zPadFact  Zero-pad factor. Keep high to capture increasing spatial
%           frequencies in the slice profile with each segment. Default = 4
% winFact   Widening factor compared to first pulse; varying this between 1
%           and zPadFact will trade off increase in pulse duration versus
%           slice profile consistency. Default = 1.75
% aphs      Cancel the alpha phase using the beta polynomial (no reason to
%           change this, according to Will). Default = true.
% plotAll   Generate plots of RF waveforms and slice profiles. Default = false
% T1        Longitudinal relaxation rate to use in pulse design [ms].
%           Default = Inf. Consider T1 GM = 1500 ms and T1 WM = 750 ms.
% TRseg     Time between excitations [ms]. Default = 60 ms.
% useMz     Use the residual Mz from the previous pulse. Default = true.
%           If false, then the pulses will be conventional VFA
% seSeq     false: GRE sequence (typical); true: SE sequence
% tbRef     time-bandwidth of refocusing pulse if seSeq == true
%
%   Output arguments
%
% rf    The RF waveform for each segment. An Npts*zPadFact x Nseg complex
%       matrix of the RF waveforms. Dimensions are arbitrary
% Mxy   The transverse slice profile after each segment.
% rfRef The refocusing RF waveform (if seSeq == true)
%
% Algorithm:
%   1) Design first pulse for starting flip angle, calculate its Mxy and
%      its Mz profiles
%   2) Use residual Mz and first pulse's Mxy profile to calculate target
%      beta profile for current pulse
%   3) Design current pulse, calculate Mz profile after it
%   4) Go to (2) until last pulse.
%
% Written by Will Grissom
%
% 10/24/2017:   Modified by Avery Berman (ajberman@mgh.harvard.edu) to make
%               this a function that will return the RF waveforms as
%               output.

% Check input parameters
if Nseg < 1
    error('Nseg must be positive integer');
end
if ~exist('tb', 'var') || isempty(tb)
    tb = 4;
end
if ~exist('Npts', 'var') || isempty(Npts)
    N = 2000;
else
    N = Npts;  % Change to Will's original syntax
end
if ~exist('zPadFact', 'var') || isempty(zPadFact)
    zPadFact = 4;
end
if ~exist('aphs', 'var') || isempty(aphs)
    cancelAlphaPhs = true;
else
    cancelAlphaPhs = aphs;  % Change to Will's original syntax
end
if ~exist('winFact', 'var') || isempty(winFact)
    winTruncFact = 1.75;
else
    winTruncFact = winFact;  % Change to Will's original syntax
end
if ~exist('plotAll', 'var') || isempty(plotAll)
    plotAll = false;
end
if ~exist('T1', 'var') || isempty(T1)
    T1 = Inf; % set T1 = Inf; to ignore longitudinal relaxation
end
if ~exist('TRseg', 'var') || isempty(TRseg)
    TRseg = 60; % time between segment excitations, ms
end
if ~exist('useMz', 'var') || isempty(useMz)
    useMz = true;
end
if ~exist('seSeq', 'var') || isempty(seSeq)
    seSeq = false; % false: GRE sequence (typical); true: SE sequence
end
if ~exist('tbRef', 'var') || isempty(tbRef)
    tbRef = 8; % time-bandwidth of refocusing pulse if SE sequence
end

if ( winTruncFact > zPadFact )
    warning('dzFleet:win', 'winFact is greater than zPadFact, will set winFact = zPadFact.');
    winTruncFact = zPadFact;
elseif ( winTruncFact < 1 )
    error('dzFleet:win', 'winFact must be greater than or equal to 1.');
end

if seSeq % design a refocusing pulse
    bRef = zeros(zPadFact*N,1); % refocusing beta polynomial
    [rfRef,bRef(zPadFact*N/2-N/2+1:zPadFact*N/2+N/2,1)] = dzrf(N,tbRef,'se','ls',0.01,0.01);
    Bref = ft(bRef);
    Bref = Bref./max(abs(Bref));
    BrefMag = abs(Bref); % get beta profile so we can account for it in calculating Mz
    ArefMag = abs(sqrt(1-BrefMag.^2));
    flipRef = 2*asind(BrefMag(zPadFact*N/2+1));
end

% get initial flip angle given this number of segments.
flip = zeros(Nseg,1);flip(end) = 90;
for jj = Nseg-1:-1:1
    if ~seSeq
        flip(jj) = atand(sind(flip(jj+1)));
    else
        flip(jj) = atand(cosd(flipRef)*sind(flip(jj+1)));
    end
end

b = zeros(zPadFact*N,Nseg); % beta polynomials
rf = zeros(zPadFact*N,Nseg); % rf pulses
Mz = ones(zPadFact*N,1); % Mz profile after each pulse
Mxy = zeros(zPadFact*N,Nseg); % slice profiles

% design first RF pulse; Mz is relaxed here
b(zPadFact*N/2-N/2+1:zPadFact*N/2+N/2,1) = dzrf(N,tb,'st','ls',0.01,0.01);
% shift beta so that RF is centered, and we get no phase in beta response
B = ft(b(:,1));
B = B.*exp(-1i*2*pi/(N*zPadFact)*1/2*(-N*zPadFact/2:N*zPadFact/2-1)');
b(:,1) = ift(B);
b(:,1) = b(:,1)./max(abs(B))*sind(flip(1)/2); % scale to first flip in passband.
% we have to divide by max(abs(B)) because dzrf does not guarantee that B
% has a max of 1.
% get RF from this centered + scaled beta polynomial

a = b2a(b(:,1));
if cancelAlphaPhs
    % cancel a phase by absorbing into b
    % Note that we don't have to do this later since it is done here
    b(:,1) = ifft(fft(b(:,1)).*exp(1i*angle(fft((a(:))))));
end
rf(:,1) = b2rf(b(:,1));

% get the min-phase alpha and its response
a = b2a(b(:,1));
A = ft(a);

% calculate the beta filter response
B = ft(b(:,1));

if winTruncFact < zPadFact % then we need to apply windowing
    winLen = (winTruncFact-1)*N;
    Npad = N*zPadFact - winTruncFact*N;
    % TODO: Replace with a tukeywin (built-in MATLAB func)?
    window = blackman((winTruncFact-1)*N);
    % split it in half; stick N ones in the middle
    window = [window(1:winLen/2);ones(N,1);window(winLen/2+1:end)];
    window = [zeros(Npad/2,1);window;zeros(Npad/2,1)];
    % apply windowing to first pulse for consistency
    b(:,1) = b(:,1).*window;
    rf(:,1) = b2rf(b(:,1));
    % recalculate B and Mxy
    B = ft(b(:,1));
    A = ft(b2a(b(:,1))); % Phase of B doesn't matter here since only profile mag is used by b2a
end

% use A and B to get Mxy
if ~seSeq
    Mxy(:,1) = 2*conj(A(:)).*B; % this is the magnetization profile we want all
                                % other pulses to produce
else
    Mxy(:,1) = 2*A(:).*conj(B).*Bref.^2;
end

% Amplitude of next pulse's Mxy profile will be |Mz*2*a*b| = |Mz*2*sqrt(1-abs(B).^2)*B|.
% If we set this = |Mxy_1|, we can solve for |B| via solving quadratic equation
% 4*Mz^2*(1-B^2)*B^2 = |Mxy_1|^2.
% Subsequently solve for |A|, and get phase of A via min-phase, and
% then get phase of B by dividing phase of A from first pulse's Mxy phase.
for jj = 2:Nseg

    % calculate Mz profile after previous pulse
    if ~seSeq % gre sequence
        Mz = Mz.*(1-2*abs(B).^2)*exp(-TRseg/T1)+(1-exp(-TRseg/T1));
    else % se sequence
        Mz = Mz.*(1-2*(abs(A(:).*BrefMag).^2 + abs(ArefMag.*B).^2)); % second term is about 1%
    end
    %figure;plot(Mz);

    if useMz % design the pulses accounting for the actual Mz profile (the full method)

      % set up quadratic equation to get |B|
      cq = -abs(Mxy(:,1)).^2;
      if ~seSeq
        bq = 4*Mz.^2;
        aq = -4*Mz.^2;
        %figure;plot(real(Mz));hold on;plot(abs(Mxy(:,1)));
      else
        bq = 4*(BrefMag.^4).*Mz.^2;
        aq = -4*(BrefMag.^4).*Mz.^2;
        %figure;plot(abs(real(Mz.*Bref.^2)));hold on;plot(abs(Mxy(:,1)));
      end
      Bmag = sqrt((-bq+real(sqrt(bq.^2-4*aq.*cq)))./(2*aq));
      Bmag(isnan(Bmag)) = 0;
      % get A - easier to get complex A than complex B since |A| is
      % determined by |B|, and phase is gotten by min-phase relationship
      A = ft(b2a(ift(Bmag))); % Phase of B doesn't matter here since only profile mag is used by b2a
      % trick: now we can get complex B from ratio of Mxy and A
      B = Mxy(:,1)./(2*conj(A(:)).*Mz);

    else % design assuming ideal Mz (conventional VFA)

      B = B*sind(flip(jj)/2)/sind(flip(jj-1)/2);
      A = ft(b2a(ift(B))); % Phase of B doesn't matter here since only profile mag is used by b2a

    end

    % get polynomial
    b(:,jj) = ift(B);

    if winTruncFact < zPadFact % then we need to apply windowing
        b(:,jj) = b(:,jj).*window;
        % recalculate B and Mxy
        B = ft(b(:,jj));
        A = ft(b2a(b(:,jj))); % Phase of B doesn't matter here since only profile mag is used by b2a
    end

    if ~seSeq
        Mxy(:,jj) = Mz.*(2*conj(A(:)).*B);
    else
        Mxy(:,jj) = Mz.*(2*A(:).*conj(B).*Bref.^2);
    end

    % get RF
    rf(:,jj) = b2rf(b(:,jj)).';

end

Mxy = Mxy.*repmat(exp(1i*2*pi/N*(zPadFact*N/2)*(-N/2:1/zPadFact:N/2-1/zPadFact)'),[1 Nseg]);

% also look at case with pulse 1 scaled to different flips
if ~seSeq
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
end

% truncate the RF
if winTruncFact < zPadFact % then we need to truncate
    pulseLen = winTruncFact*N;
    Npad = N*zPadFact - pulseLen;
    rf = rf(Npad/2+1:Npad/2+pulseLen,:);
end

% compare stability of summed signals
if plotAll

    figure;
    subplot(2,2,1)
    plot(real(rf));
    hold on
    plot(imag(rf));
    title 'RF pulses'
    ylabel 'radians'
    xlabel 'time index'
    legendText = {};
    for ii = 1:Nseg
        legendText{ii} = ['Seg ' num2str(ii)];
    end
    legend(legendText);

    subplot(2,2,2)
    plot(abs(sum(Mxy)));
    hold on
    if ~seSeq
        plot(abs(sum(Mxy_sameRF)));
        legend('SLR-VFA','Conventional VFA')
    else
        legend('SLR-VFA');
    end
    xlabel 'Segment'
    ylabel '|\int M_{xy}| (a.u.)'
    title 'Amplitude of integrated slice profile'

    x = -N/2:1/zPadFact:N/2-1/zPadFact;
    subplot(2,2,3)
    plot(x,real(Mxy));
    hold on
    plot(x,imag(Mxy));
    c = axis;
    axis([-4*tb 4*tb c(3:4)]);
    title 'SLR-VFA Mxy Profiles'
    ylabel '/M_0'
    xlabel 'frequency index'

    if ~seSeq
        subplot(2,2,4)
        plot(x,real(Mxy_sameRF));
        hold on
        plot(x,imag(Mxy_sameRF));
        c = axis;
        axis([-4*tb 4*tb c(3:4)]);
        title 'Conventional VFA Mxy Profiles (First pulse scaled)'
        ylabel '/M_0'
        xlabel 'frequency index'
    end

    % Plot slice profile magn and phase
    Mxy_mag = abs(Mxy);
    Mxy_pha = angle(Mxy);
    Mxy_pha(Mxy_mag < 0.02) = 0;

    figure;
    subplot(2,1,1);
    plot(x,Mxy_mag);
    c = axis;
    axis([-4*tb 4*tb c(3:4)]);
    title 'SLR-VFA Mxy Profiles'
    ylabel 'magnitude'
    xlabel 'frequency index'

    subplot(2,1,2);
    plot(x,Mxy_pha);
    c = axis;
    axis([-4*tb 4*tb c(3:4)]);
    title 'SLR-VFA Mxy Profiles'
    ylabel 'phase  [rad]'
    xlabel 'frequency index'

end  % plotAll
