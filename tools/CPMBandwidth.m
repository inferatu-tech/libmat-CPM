function [normBandwidth] = CPMBandwidth(CPMpar, pctBandwidth)
% [normBandwidth,filterTaps] = CPMBandwidth(CPMpar, pctBandwidth)
% computes the normalized bandwidth of a CPM signal by Monte-Carlo
% simulation.
% 
% --- INPUTS
%        CPMpar: CPM parameter set 
%  pctBandwidth: percent energy bandwidth 
% 
% --- OUTPUTS
% normBandwidth: bandwidth (normalized to symbol rate) at indicated percentile
% 

% fixed configuration 
T = 64; 
numSyms = 8192; 
normFiltDelay= 30;
normFiltLen=2*normFiltDelay;
if nargin==1
    pctBandwidth = 99;
end
%
modulator = CPMModulator(CPMpar,T); 
symSeq = randi([0 CPMpar.M-1],1,numSyms);
signalIn = modulator.modulate(symSeq);
powerIn = mean(abs(signalIn).^2);  
%
hypNormBandwidth = 0.75:0.01:2;
filtCutOff = hypNormBandwidth/T;
normBandwidth = [];
for i=1:length(filtCutOff) 
    filtTaps = firls(normFiltLen*T, [0 filtCutOff(i) filtCutOff(i) 1], [1 1 0 0]);
    signalOut = filter(filtTaps,1,[signalIn, zeros(1,normFiltLen*T)]);
    signalOut = signalOut((1+normFiltDelay)*T:end);
    powerOut = mean(abs(signalOut).^2);
    fracPower = powerOut/powerIn;
    disp(fracPower)
    if fracPower>=(pctBandwidth/100)
        normBandwidth = hypNormBandwidth(i);
        break;
    end
end
assert(~isempty(normBandwidth),...
    '[CPMBandwidth]:Could not find percentile bandwidth in the search range...');



