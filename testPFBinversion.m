fprintf('\nTest of OS-PFB Inversion via FFT\n');

%% GLOBAL PARAMETERS

% PFB parameters
N = 32;  % number of output channels - power of 2, min OS_Nu, max 256
OS_Nu = 4;  % OS numerator - should be a sub-multiple of N
OS_De = 3;  % OS denominator
tapsPerChan = 8;
quantisePFBcoeff = 16;  % choose between 0 (no quantisation), 8, 12, 13, 14, 15 or 16
displayFilterDesign = 1;  % 1 to display filter design plots, 0 otherwise

% Width of PFB channel passband in MHz = spacing of PFB output channels
fineChanPassband = 0.01;

% Length of forward FFT to process fine channels
fftLength = 2^8;

% Length of test vector blocks (spacing of impusles)
blockLength = 2*N*fftLength;

 
%% GENERATE TEST VECTOR - ONE BLOCK (input to PFB)

WaveType = 1;  % 0 for pulsar, 1 for impulse
impulseOffset = blockLength/4;  % location of impulse within each block
impulseWidth = 1;  % number of samples width of impusle
freqSampleOut = N*fineChanPassband;  % sample rate in MHz
period = 0.001;  % simulated pulsar period in seconds
noise = 0.0;  % sets SNR of simulated pulsar signal
DM = 0.0;   % sets DM of simulated pulsar signal
testVectorFilename = 'test_vec.bin';

fprintf('\nGenerating test vector...\n');
genTestVector(WaveType,impulseOffset,impulseWidth,blockLength,1,freqSampleOut,period,noise,DM,testVectorFilename);


%% DESIGN PFB PROTOTYPE FILTER

Ntaps = N*tapsPerChan + 1;  % must be odd; one more than N*tapsPerChan
    
fprintf('designing PFB prototype filter\n');
if (displayFilterDesign)
    fprintf('\nPress any key to continue...\n');
end;
designPFB(N,OS_Nu,OS_De,Ntaps-1,fftLength,quantisePFBcoeff,displayFilterDesign);  % subtract 1 from num taps because design functions adds 1


%% PFB CHANNELISE - ONE BLOCK

% minimum input size is (blockLength/OS_factor) - can be longer
fprintf('\nChannelising...\n');

PFBchannelizerComplex(N,OS_Nu,OS_De,(blockLength*OS_De)/OS_Nu,1,testVectorFilename,'fine_channel_');


%% PROCESS EACH FINE CHANNEL

inputOffset = 64;  % number of samples to drop at the start of the PFB output data, to ensure impulse within window
equaliseRipple = 1;  % 1 to equalise PFB ripple, 0 to not
fprintf('\nProcessing each channel...\n');
for chan = 1:N
    fprintf('channel %d\n', chan);
    fineChannelProc(chan,fftLength,OS_Nu,OS_De,inputOffset,strcat('fine_channel_',int2str(chan),'.dump'),strcat('chunk_',int2str(chan),'.mat'),equaliseRipple);
end;


%% COMBINE CHUNKS, BACK-TRANSFORM AND COMPARE TO ORIGINAL

fprintf('\nCombining channels and back transforming...\n');
compareOffset = (Ntaps-1)/2 + 1 - (OS_De*N/OS_Nu)*inputOffset;
invertComplex(N,OS_Nu,OS_De,blockLength/2,'chunk_',testVectorFilename,compareOffset);

fprintf('\nDone! Press any key to close plots and exit...\n\n');
pause;
close all;
clear all;
