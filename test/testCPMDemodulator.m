clear;
close all;
path(pathdef);
addpath('../src/');
addpath('../src/shape');

par = CPMParameters();
par.L = 2;
par.m = 2;
par.gmskBT = 0.3;
par.hK = 1;
par.hP = 4;
par.frPulseStr = 'RC';
T = 4;

modulator = CPMModulator(par,T); 
demodulator = CPMDemodulator(par,T);
numSyms = 1024;

numTests = 10;
for iTest=1:numTests
    symSeq= randi([0 par.M-1],1,numSyms);
    initPhaseState = randi([0, modulator.numPhaseStates-1]);
    modulator.reset(initPhaseState);
    iqSamples = modulator.modulate(symSeq);
    demodulator.reset(initPhaseState);
    [symSeqEst, stSeqEst, bestPathMetric] = demodulator.demodulate(iqSamples);
    numSymErrs = sum(symSeqEst~=symSeq);
    disp(strcat('Sym error rate:',num2str(numSymErrs),'/',num2str(numSyms)));
    disp(strcat('Best path metric:',num2str(bestPathMetric)));
end