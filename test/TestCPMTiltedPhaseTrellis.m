clear;
close all;
path(pathdef);
addpath('../src/');
addpath('../src/shape');

par = CPMParameters();
par.L = 3;
par.m = 1;
par.hK = 1;
par.hP = 2;
par.frPulseStr = 'REC';

T = 8; 
tr = CPMTiltedPhaseTrellis(par,T); 
disp(strcat('  Number of phase states: ',num2str(tr.numStates)));
disp(strcat('Number of branches/state: ',num2str(tr.numBranchesPerState)));
 
format compact;
figure(1);
phiRad=pi*[-1:0.01:1];
h = plot(cos(phiRad), sin(phiRad),'k-');
set(gca,'XTick',-1:0.2:1);
set(gca,'YTick',-1:0.2:1);
grid;
hold on;
for St=0:tr.numStates-1
    for u=0:tr.numBranchesPerState-1
        disp(strcat('St:',num2str(St),...
                    ',Inp:', num2str(u),...
                    ',nSt:', num2str(tr.nextState(St,u)),...
                    ',phi/pi:', num2str(tr.samplesRad(St,u)/pi)));
        h = plot( cos(tr.samplesRad(St,u)), sin(tr.samplesRad(St,u)), 'bo');
        set(h,'MarkerSize',3);
        set(h,'LineWidth',2);
        axis([-1 1 -1 1]);
        axis square;
    end
end


for St=0:tr.numStates-1
    disp(strcat('Num branches to state-',num2str(St),':',...
                num2str(length(tr.branchesToState(St)))));
end

for St=0:tr.numStates-1
    disp(strcat('Num branches from state-',num2str(St),':',...
                num2str(length(tr.branchesFromState(St)))));
end