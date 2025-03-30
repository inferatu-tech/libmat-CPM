classdef CPMParameters  

    properties (SetAccess = 'public', GetAccess = 'public')
        L;  % memory
        m;  % number of bits per symbol (log2(M))
        hK; % numerator of modulation index
        hP; % denominator of moulation index 
        frPulseStr; % frequency pulse string
        gmskBT; %special case, only specified for frPulseStr = 'GMSK'  
    end

    methods (Access = 'public')
        % default constructor 
        function obj = CPMParameters()   
        end
        
        % provides access to modulation index
        function [ modIdx ] = h( obj )
            modIdx = obj.hK/obj.hP;
        end
        
        % provides access to alphabet size
        function [ alphSize ] = M( obj )
            alphSize = 2^obj.m;
        end

        % performs sanity checks
        function v = isValid(obj)
            v = true;
            v = v & (obj.L>=1); assert(v, '[CPMParameters]:Invalid L...');
            v = v & (obj.m>=1); assert(v, '[CPMParameters]:Invalid m...');
            v = v & (obj.hK>=1); assert(v, '[CPMParameters]:Invalid hK...');
            v = v & (gcd(obj.hP,obj.hK)==1); assert(v, '[CPMParameters]:hP and hK must be relatively prime...');
            v = v & (strcmp(obj.frPulseStr,'REC')|strcmp(obj.frPulseStr,'RC')|strcmp(obj.frPulseStr,'HS')|strcmp(obj.frPulseStr,'GMSK')); assert(v, '[CPMParameters]:Invalid frPulseStr...');
            v = v & (~strcmp(obj.frPulseStr,'GMSK')||obj.gmskBT>0); assert(v, '[CPMParameters]:invalid gmskBT parameter...');
        end
    end %public methods

end % 