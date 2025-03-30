function [ frPulse ] = RaisedCosineFreqPulse( L, T )
% [ frPulse ] = RaisedCosineFreqPulse( L, T ) computes  
% a raised cosine frequency pulse 
%
% --- INPUTS
%       L: memory parameter
%       T: number of samples per step - needs to be >> 1 for accuracy
%
% --- OUTPUTS
% frPulse: phase pulse as a 1 x LT vector
%

% sanity checks
assert(L>=1,'[RaisedCosineFreqPulse]: Invalid L...');
assert(T>=1,'[RaisedCosineFreqPulse]: Invalid T...');
t = 0:L*T-1;
frPulse = (1 - cos(2*pi*t/(L*T)))/(2*L*T);
frPulse = 0.5*frPulse /sum(frPulse);

end %RaisedCosineFreqPulse 
            
