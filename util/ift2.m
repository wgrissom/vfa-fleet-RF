function output=ift2(in)
% Performs fftshift(ifft2(fftshift( input)))
% 2D inverse FT
output = fftshift(ifft2(fftshift(in)));
