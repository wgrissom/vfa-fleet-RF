tbExc = 4;
tbRef = 6.5;
tbRefScale = 1.5;  % This seemed quite high but was what I needed in order to get the signal consistency across shots below 1% diff.
zPad = 4;
winFact = 1.75;
Npts = 1024;
Nseg = 2;

[rf, Mxy, rfRef] = dzFleet(Nseg, tbExc, Npts, zPad, winFact, true, true, ...
    [], [], true, true, tbRef);
