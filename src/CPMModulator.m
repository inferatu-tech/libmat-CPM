classdef CPMModulator < handle

    % member variables
    properties (SetAccess = 'private', GetAccess = 'public')
        phaseTrellis;
        amrRotPerSampleRad;
    end

    % member variables
    properties (SetAccess = 'private', GetAccess = 'private')
        T; %shorthand
        M;
        phaseState;
        amrRotRad;
    end

    methods (Access = 'public')
        % constructs with CPM parameters
        function obj = CPMModulator( CPMpar, T )
            obj.phaseTrellis = CPMTiltedPhaseTrellis( CPMpar, T );
            obj.T = T;
            obj.M = CPMpar.M;
            obj.amrRotPerSampleRad = pi*CPMpar.h()*(obj.M-1)/obj.T;
            obj.reset();
        end

        % resets the phase state
        function reset( obj, initPhaseState )
            if nargin == 1
                initPhaseState = 0;
            end
            assert(initPhaseState>=0 & initPhaseState<obj.numPhaseStates(),...
                   '[CPMModulator::reset] invalid initPhaseState...' );
            obj.phaseState = initPhaseState;
            obj.amrRotRad = 0;
        end

        % provides access to number of phase states
        function [ numSt ] =  numPhaseStates( obj )
            numSt = obj.phaseTrellis.numStates;
        end

        % maps a sequence of symbols onto samples, using 
        % samplePower [optional] power in modulated samples
        function [ sampleSeq ] = modulate( obj, symSeq, samplePower )
            % perform sanity checks
            assert(~isempty(symSeq),...
                   '[CPMModulator::modulate] symSeq is empty...');  
            assert(min(symSeq)>=0 & max(symSeq)<obj.M, ...
                  '[CPMModulator::modulate] invalid symbols in symSeq ...' );
            if nargin == 2
                samplePower = 1;
            end
            assert(samplePower>0, ...
                   '[CPMModulator::modulate] samplePower needs to be > 0...' );
            sampleSeq = sqrt(samplePower) * obj.modulateImpl( symSeq ); 
        end

    end % public methods

    methods (Access = 'private')
        %
        function sampleVec = modulateImpl( obj, symSeq )
            numSyms = length(symSeq);
            samplesByStep= zeros(numSyms,obj.T);
            for symIdx=1:length(symSeq)
                thisSym = symSeq(symIdx);
                thisSymAmrRotRad = obj.amrRotRad + (0:obj.T-1)*obj.amrRotPerSampleRad;
                thisSymPhRad = obj.phaseTrellis.samplesRad(obj.phaseState,thisSym) - thisSymAmrRotRad;
                samplesByStep(symIdx,:) = exp(1i*thisSymPhRad);
                obj.phaseState = obj.phaseTrellis.nextState(obj.phaseState, thisSym);
                obj.amrRotRad = obj.amrRotRad + obj.amrRotPerSampleRad*obj.T;
            end
            obj.amrRotRad = mod(obj.amrRotRad,2*pi);
            % arrange into row vector
            numSamples = numel(samplesByStep);
            sampleVec = reshape(samplesByStep',numSamples,1)';
        end

    end

end


