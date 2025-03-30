function [ frPulse ] =HalfSineFreqPulse( L, T )
% [ frPulse ] =HalfSineFreqPulse( L, T ) computes 
% a half-sine frequency pulse 
%
% --- INPUTS
%       L: memory parameter
%       T: number of samples per step - needs to be >> 1 for accuracy
%
% --- OUTPUTS
% frPulse: phase pulse as a 1 x LT vector
%

% sanity checks
assert(L>=1,'[HalfSineFreqPulse]: Invalid L...');
assert(T>=1,'[HalfSineFreqPulse]: Invalid T...');
t = (0:L*T-1);
frPulse = pi*sin(pi*t/(L*T))/(4*L*T);
end %HalfSineFreqPulse 
            
