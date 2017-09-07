function [H, f, n_hi, n_lo] = ...
                dispnmatrix(frange,Nout,nfreq,DD,Tout,direc)
% Calculates the de-dispersion matrix and the number of leading
% (n_hi) and trailing (n_lo) elements that need to be removed during
% overlap-save.
%
% Inputs:
% --------
%   frange    - 2-element vector containing highest and lowest freqs
%   Nout      - length of FFT to be analysed
%   nfreq     - number of frequency channels
%   DD        - DM*Dconst 
%   Tout      - sampling period of data
%   direc     - Dispersion = 1; De-dispersion = -1 (default)
%
% Outputs:
% --------
%
%   H       - de-dispersion matrix
%   f       - mean frequency of each frequency channel
%   n_hi    - number of leading elements to be removed
%   n_lo    - number of trailing elements to be removed
%
% Description:
% ------------
% Calculates the de-dispersion matrix used to multiply the Fourier 
% transformed, filterbanked data. Also returns the number of leading and
% trailing elements that need to be removed.
% 
% Changes:
% --------
%
% Author           Date         Comments
% ---------------  -----------  ----------------------------------------
% D. Hicks         21-Apr-2014  Original version
%
% ----------------------------------------------------------------------
 
if ~exist('direc', 'var')
    direc = -1;
end;
        
if direc ~= -1 && direc ~= 1
    direc = -1;
    disp('Warning: direc can only be 1 or -1');
end;

%direc
 
% Vector of frequency bin assignments. Note that the highest frequency
% bin is assigned a value of frange(2)-df where frange(2) is the Nyquist
% frequency
bandwidth = (frange(2) - frange(1));
deltaf = bandwidth/(Nout*nfreq); % frequency spacing
fabs = frange(1) + (0:(Nout*nfreq-1))*deltaf;
%fabs = linspace(frange(1), frange(2), Nout*nfreq ); %Not quite correct
 
% Absolute freqs in each channel
fc = reshape(fabs, Nout, nfreq);
% Mean of each freq channel
f = mean(fc,1);
% Expand to create matrix 
f0c = repmat(f, Nout, 1);
% Phase dispersion
fphi = DD*((fc-f0c).^2)./(f0c.^2)./fc*1E6;

%Dispersion matrix: De-dispersion has direc = -1; Dispersion has direc = 1
H = exp(complex(0,direc)*2*pi*fphi);

% Calculate convolution overlap region
%fcmin = [fabs(nnmin(1)) fabs(nnmax(1))];
fcmin = [ fc(1,1), fc(Nout,1) ]; %CHECK THIS !!!!!!
fcmin0 = mean(fcmin);
 
% CHECK WHETHER n_hi and n_lo ARE FLIPPED DEPENDING ON THE DIRECTION
% OF THE DISPERSION (i.e. DISPERSING OR DE-DISPERSING)
n_hi = ceil(DD*(1/fcmin0^2 - 1/max(fcmin)^2)/Tout);
n_lo = ceil(DD*(1/min(fcmin)^2 - 1/fcmin0^2)/Tout);
 
% If overlap exceeds length of vector, force there to be no overlap
% (useful for debugging with short time series).
if (n_hi + n_lo) >= Nout,
    n_hi=0;
    n_lo=0;
    disp('Time series too short for given DM');
end;

return
end

