% script to test the dzFleet function which designs
% same-profile RF pulses for the FLEET sequence

Nseg = 3; % number of EPI segments/RF pulses
tb = 4; % time bandwidth product of excitation pulses
Npts = 2000; % number of time points in first segment's pulse
zPadFact = 4; % Zero-pad factor. Keep high to capture increasing spatial
              % frequencies in the slice profile with each segment.
winFact = 1.75; % Widening factor compared to first pulse; varying this between 1
                % and zPadFact will trade off increase in pulse duration versus
                % slice profile consistency.
cancelAlphaPhs = true; % Cancel the alpha phase using the beta polynomial
plotAll = true; % Generate plots of RF waveforms and slice profiles
T1 = 1000; % Longitudinal relaxation rate to use in pulse design [ms].
TRseg = 100; % Time between excitations [ms].
finalFlip = 50; % 50 works well for se and other settings here; 90 can be used for GRE
TE = 50; % echo time [ms] - relevant only for spin echo
seSeq = true; % false: GRE sequence (typical); true: SE sequence
tbRef = 8; % time-bandwidth of refocusing pulse if seSeq == true
useMz = true; % Use the residual Mz from the previous pulse. Default = true.
               % If false, then the pulses will be conventional VFA

addpath util/

[rf, Mxy] = dzFleet(Nseg, tb, Npts, zPadFact, winFact, cancelAlphaPhs, plotAll, T1, ...
              TRseg, TE, finalFlip, useMz, seSeq, tbRef);
