clear;
close all;
path(pathdef);
addpath('../src/');
addpath('../src/shape');
addpath('../tools/');

par = CPMParameters();
par.L = 4;
par.m = 1;
par.hK = 1;
par.hP = 2;
par.gmskBT = 0.3;
par.frPulseStr = 'GMSK';
T = 64;

modulator = CPMModulator(par,T); 

numSyms = 2^15;
symSeq = randi([0 par.M-1],1,numSyms);
iqSamples = modulator.modulate(symSeq);

plotSpectrum( iqSamples, T, ...
              strcat('CPM L:',num2str(par.L),'-',par.frPulseStr, ', M:', num2str(par.M), ', h:', num2str(par.h())) ) ;
pctBandwidth = 99;
normBandwidth = CPMBandwidth(par, pctBandwidth);
disp(strcat('Bandwidth (',num2str(pctBandwidth),'pct.):',num2str(normBandwidth),'/T'));


function plotSpectrum( signal, normSampleRate, titleStr )
% sanity checks 
signal = transpose(signal(:));
assert( ~isempty(signal), '[plotSpectrum] input signal is empty...');
assert( normSampleRate>0, '[plotSpectrum] normSampleRate must be positive...');
if nargin==2
    titleStr=[];
end
%
fftSize = 512*normSampleRate;
if length(signal) < fftSize
    signal = [signal, zeros(1,fftSize-length(signal))];
end
window = hann(fftSize); %kaiser(fftSize,5);
ovlSize = round(7*fftSize/8);
freqRangeType = 'twosided';
[ energyPerFreqBin, normFreqBins ] = ...
    pwelch( signal, window, ovlSize, fftSize, normSampleRate, freqRangeType);
normEnergyPerFreqBin = energyPerFreqBin/max(energyPerFreqBin);

normFreqBins = ( normFreqBins - normSampleRate/2 );
normEnergyPerFreqBin = fftshift( normEnergyPerFreqBin );
normEnergyPerFreqBin = movmean(normEnergyPerFreqBin,31);
normEnergyPerFreqBindB = 10*log10( normEnergyPerFreqBin );
normEnergyPerFreqBindB = normEnergyPerFreqBindB - max(normEnergyPerFreqBindB);
h=plot(normFreqBins,normEnergyPerFreqBindB,'k-');
set(h,'LineWidth',2.5);
xlabel('Norm. frequency, fT','FontName','Times');
ylabel('Norm. energy per frequency, dB-peak','FontName','Times');
axis([-2 2 -40 0]);
%set(gca,'XTick',-2:0.5:2);
grid;
if ~isempty(titleStr)
    title(titleStr,'FontName','Times');
end

end
