function [ frPulse ] = GMSKFreqPulse( gauBT, L, T )
% [ frPulse ] = GMSKFreqPulse( gauBT, L, T ) computes the frequency
% pulse for a GMSK modulation specified by the inputs
%
% --- INPUTS
%   gauBT: BT value of the underlying frequency pulse
%       L: memory parameter
%       T: number of samples per step - needs to be >> 1 for accuracy
%
% --- OUTPUTS
%  frPulse: frequency pulse as a 1 x LT vector
%

% sanity checks
assert(gauBT>0,'[GMSKFreqPulse]: Invalid gauBT...');
assert(L>=1,'[GMSKFreqPulse]: Invalid L...');
assert(T>=1,'[GMSKFreqPulse]: Invalid T...');
% gaussian stdev
sigma = sqrt(log(2))/( 2*pi*gauBT );
% computes the frequency pulse
frPulse = double.empty();
t = 1;
while length(frPulse) < L*T
    fVal = (0.5/T) * (GMSKQfunc((t/T-0.5)/sigma) - GMSKQfunc((t/T+0.5)/sigma));
    frPulse = [fVal, frPulse, fVal];
    t=t+1;
end
frPulse = frPulse - frPulse(1);
frPulse = frPulse(1:L*T);
frPulse = 0.5*frPulse / sum(frPulse);

    % Gaussian Q function
    function q = GMSKQfunc(x)
        q = 0.5*erfc(x/sqrt(2));
    end

end %GMPKPhasePulse 
            
