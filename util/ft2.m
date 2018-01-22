function output=ft2(in)
% Performs fftshift(fft2(fftshift(input)))
% 2D forward FT
output = fftshift(fft2(fftshift(in)));
