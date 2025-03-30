function [ frPulse ] = RecFreqPulse( L, T )
% [ frPulse ] = RecFreqPulse( L, T ) computes  
% a rectangular frequency pulse 
%
% --- INPUTS
%       L: memory parameter
%       T: number of samples per step - needs to be >> 1 for accuracy
%
% --- OUTPUTS
% frPulse: phase pulse as a 1 x LT vector
%

% sanity checks
assert(L>=1,'[RecPhasePulse]: Invalid L...');
assert(T>=1,'[RecPhasePulse]: Invalid T...');
frPulse = ones(1,L*T)/(2*L*T);
end %RecFreqPulse 
            
