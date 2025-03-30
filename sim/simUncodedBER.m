% simulated uncoded SER 
clear;
close all;
path(pathdef);
addpath('../src/');
addpath('../src/shape');
addpath('../tools/');

% signal specification
par = CPMParameters();
par.L = 3;
par.m = 1;
par.hK = 1;
par.hP = 2;
par.frPulseStr = 'GMSK';
par.gmskBT = 0.33;
T = 4;
normBandwidth = 1; %find this using tools/GMSKBandwidth
% modem 
modulator = CPMModulator(par,T); 
demodulator = CPMDemodulator(par,1);
% receive filter specification
normFiltDelay= 20;
normFiltLen=2*normFiltDelay;
rcvFiltTaps = firls(normFiltLen*T, [0 normBandwidth/T normBandwidth/T 1], [1 1 0 0]);
% simulation specification 
numSymsPerTest = 1024;
numTests = 25;
SNRdB = -3:0.5:12;
bitMat = randi([0,1], numTests, numSymsPerTest*par.m );

% bit mapper
mapper = CPMBitMapper(par.m); 
SERbySNR = zeros(size(SNRdB));
BERbySNR = zeros(size(SNRdB));
for iSNR=1:length(SNRdB)
    numSymErrs = 0;
    numBitErrs = 0;
    for iTest=1:numTests    
        bitSeq = bitMat(iTest,:);
        symSeq = mapper.map(bitSeq);
        initPhaseState = randi([0, modulator.numPhaseStates()-1]);
        modulator.reset(initPhaseState);
        sigSamples = modulator.modulate(symSeq);
        sigPower = mean(abs(sigSamples).^2);
        noiPower = (T/normBandwidth)*sigPower*10^(-SNRdB(iSNR)/10);
        noiSamples = sqrt(noiPower/2)*(randn(size(sigSamples))+1i*randn(size(sigSamples)));
        demodulator.reset(initPhaseState); %known initial phase state 
        %demodulator.reset(0); %unknown phase state
        rcvrInput= sigSamples+noiSamples;
        filtOutput = filter(rcvFiltTaps,1,[rcvrInput, zeros(1,normFiltDelay*T)]);
        demodInput = filtOutput(1+normFiltDelay*T:T:end);
        [symSeqEst, ~, bestPathMetric] = demodulator.demodulate(demodInput);
        bitSeqEst = mapper.demap(symSeqEst);
        numSymErrs = numSymErrs + sum(symSeqEst~=symSeq);
        numBitErrs = numBitErrs + sum(bitSeqEst~=bitSeq);
    end
    SERbySNR(iSNR) = numSymErrs/(numTests*numSymsPerTest);
    BERbySNR(iSNR) = numBitErrs/(numTests*numSymsPerTest*par.m );
    %plot
    figure(1);
%     h=semilogy(SNRdB(1:iSNR),SERbySNR(1:iSNR),'k-o',...
%                SNRdB(1:iSNR),BERbySNR(1:iSNR),'k-*');
    h=semilogy(SNRdB(1:iSNR),BERbySNR(1:iSNR),'k-o');
    set(h,'LineWidth',1.5);
    set(h,'MarkerSize',3.5);
    xlabel('SNR, dB','FontName','Times');
    ylabel('Error Rate','FontName','Times');
%     legend('SER','BER');
legend('BER');
    axis([min(SNRdB) max(SNRdB) 1e-4 1]);
    set(gca,'XTick',SNRdB);
    axis square;
    grid;
end




