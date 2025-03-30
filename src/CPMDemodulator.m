classdef CPMDemodulator < handle

    properties ( Constant )
        minMetric = 0;
        maxMetric = Inf;
    end

    % member variables
    properties (SetAccess = 'private', GetAccess = 'public')
        phaseTrellis;
        amrRotPerSampleRad;
    end

    % member variables
    properties (SetAccess = 'private', GetAccess = 'private')
        T; %shorthand 
        M;
        amrRotRad;
        survPaths;
    end

    methods (Access = 'public')
        % constructs with CPM parameters
        function obj = CPMDemodulator( CPMpar, T )
            obj.phaseTrellis = CPMTiltedPhaseTrellis( CPMpar, T );
            obj.T = T;
            obj.M = CPMpar.M;
            obj.amrRotPerSampleRad = pi*CPMpar.h()*(obj.M-1)/obj.T;
            obj.reset();
        end

        % resets to the factor settings
        function reset( obj, initPhaseState )
            initPhaseStateSpecified = (nargin==2);
            if initPhaseStateSpecified
                assert(initPhaseState>=0 & initPhaseState<obj.numPhaseStates(),...
                       '[CPMDemodulator::reset] invalid initPhaseState...' );
            end
            obj.amrRotRad = 0;
            obj.survPaths = CPMTrellisPath.empty( obj.phaseTrellis.numStates, 0 );
            % initialize survivors
            for phState = 0:obj.phaseTrellis.numStates-1
                obj.survPaths(1+phState).metric = obj.maxMetric + eps;
                obj.survPaths(1+phState).stateSeq = phState;
            end
            if initPhaseStateSpecified 
                obj.survPaths(1+initPhaseState).metric = obj.minMetric;
            end
        end

        % provides access to number of phase states
        function [ numSt ] =  numPhaseStates( obj )
            numSt = obj.phaseTrellis.numStates;
        end

        % de-maps a sample sequence to a modulating symbol sequence  
        % phState [optional] starting phase state
        % all channel (phase and amplitude) effects are assumed to be
        % have been accounted for. 
        function [ symSeq, stateSeq, bestSurvMetric ] = demodulate( obj, sampleSeq )
            % perform sanity checks
            assert(~isempty(sampleSeq),...
                   '[CPMDemodulator::demodulate] sampleSeq is empty...');
            assert(mod(length(sampleSeq),obj.T)==0,...
                   '[CPMDemodulator::demodulate] sampleSeq should represent integer number of symbols...');
           
            % process step by step 
            numSteps = length(sampleSeq)/obj.T;
            samplesByStep = transpose(reshape( sampleSeq(:), obj.T, numSteps ));
            for idxStep = 1:numSteps
                thisStepAmrRotRad = obj.amrRotRad + (0:obj.T-1)*obj.amrRotPerSampleRad;
                rotSamplesThisStep = samplesByStep(idxStep,:) .* exp(1i*thisStepAmrRotRad);
                obj.updateSurvivors(rotSamplesThisStep);
                obj.amrRotRad = obj.amrRotRad + obj.amrRotPerSampleRad*obj.T;
            end
            % extract bit and state sequence from best survivor 
            [ symSeq, stateSeq, bestSurvMetric ] = obj.finalize();
            bestSurvMetric = bestSurvMetric/numSteps;
        end

    end % public methods

    methods (Access = 'private')
        % 
        function [] = updateSurvivors( obj, samplesThisStep )
            % make a copy before updating survivors 
            survPathsCopy = obj.survPaths;
            % for all ending states 
            for phState = 0:obj.phaseTrellis.numStates-1
                brTo = obj.phaseTrellis.branchesToState(phState);
                bestBrIdx = 0; % initialize
                bestPathMetric = obj.maxMetric;
                for brIdx=1:length(brTo)
                    brMetric = sum( abs(samplesThisStep- exp(1i*brTo(brIdx).phaseRad)).^2 );
                    pathMetric = survPathsCopy(1+brTo(brIdx).begSt).metric + brMetric;
                    if pathMetric<=bestPathMetric
                        bestPathMetric=pathMetric;
                        bestBrIdx = brIdx;
                    end
                end %for all branches
                % update
                obj.survPaths(1+phState) = survPathsCopy(1+brTo(bestBrIdx).begSt); 
                obj.survPaths(1+phState).stateSeq(end+1) = phState;
                obj.survPaths(1+phState).inpSymSeq(end+1) = brTo(bestBrIdx).inpSym;
                obj.survPaths(1+phState).metric = bestPathMetric;
            end %for all end states
        end %updateSurvivors

        % 
        function [ symSeq, stateSeq, bestSurvMetric] = finalize( obj )
            bestSurvSt = 0;
            bestSurvMetric = obj.maxMetric;
            for phState = 0:obj.phaseTrellis.numStates-1
                if obj.survPaths(1+phState).metric < bestSurvMetric
                    bestSurvMetric = obj.survPaths(1+phState).metric;
                    bestSurvSt = phState;
                end
            end
            symSeq = obj.survPaths(1+bestSurvSt).inpSymSeq;
            stateSeq = obj.survPaths(1+bestSurvSt).stateSeq;
        end

    end%private methods

end% CPMDemodulator 


