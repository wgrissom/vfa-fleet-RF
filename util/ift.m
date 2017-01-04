function output=ift(in)
% Performs fftshift(ifft(fftshift( input)))
% 1D inverse FT
output = ifftshift(ifft(ifftshift(in)));
