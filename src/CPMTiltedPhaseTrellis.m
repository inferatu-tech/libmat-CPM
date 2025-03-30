classdef CPMTiltedPhaseTrellis < handle

    % CPMTiltedPhaseTrellis described th time-invariant (tilted) CPM phase trellis
    % per Rimoldi.
    %
    % The underlying upsample factor (Tfull) is relatively large
    % to accommodate both transmit- and receive-side modeling, via the specification
    % of an upsample factor (T), representative of DAC-rate/symbol-rate for the former,
    % and the fractional processing spacing for the latter (typ. T/2 or T).
    %
    % All external state and input indices are zero-based.

    properties (Constant)
        Tfull = 640; %underlying full upsample factor
    end

    % member variables
    properties (SetAccess = 'private', GetAccess = 'public')
        L;
        M;
        T;
        hK;
        hP;
        numStates;
    end

    % member variables
    properties (SetAccess = 'private', GetAccess = 'private')
        phPulseFull; % phase pulse defined at full upsample factor
        stateBuf;
        lutBrFrState; % indexed by input bit in the second dimension
        lutBrToState;
        h; % modulation index 
        decFacFromFull;
        decBegFromFull;
    end

    methods (Access = 'public')
        % constructs with paramater object, and upsample factor (T) 
        function obj = CPMTiltedPhaseTrellis( CPMpar, T )
            assert( CPMpar.isValid(), '[CPMTiltedPhaseTrellis] Invalid CPMpar');
            % sanity checks
            if nargin == 1
                T=obj.Tfull;
            end
            assert( mod(obj.Tfull,T)==0,...
                '[CPMTiltedPhaseTrellis]: T is incompatible with Tfull... ');
            obj.L = CPMpar.L;
            obj.M = CPMpar.M;
            obj.hK = CPMpar.hK;
            obj.hP = CPMpar.hP;
            obj.h = CPMpar.h();
            obj.T = T;
            obj.decFacFromFull = obj.Tfull/obj.T;
            obj.decBegFromFull = 1; %ceil(obj.decFacFromFull/2);
            switch CPMpar.frPulseStr
                case 'REC'
                    obj.phPulseFull = cumsum( RecFreqPulse(obj.L,obj.Tfull) );
                case 'RC'
                    obj.phPulseFull = cumsum( RaisedCosineFreqPulse(obj.L,obj.Tfull) );
                case 'HS'
                    obj.phPulseFull = cumsum( HalfSineFreqPulse(obj.L,obj.Tfull) );  
                case 'GMSK'
                    obj.phPulseFull = cumsum( GMSKFreqPulse(CPMpar.gmskBT,obj.L,obj.Tfull) );
                otherwise
                    assert(CPMTiltedPhaseTrellis,'[CPMTiltedPhaseTrellis] Unknown frequency pulse type...');
            end
            % compute preliminary quantities
            obj.numStates = obj.hP*obj.M^(obj.L-1);
            obj.stateBuf = zeros(1,obj.L); %[u(n-1) u(n-2) u(n-L+1) v]
            % computes luts
            obj.initialize();
            %
            %disp('initialized');
            %for st=0:obj.numStates-1
            %obj.writeStateToBuffer(st);
            %disp(strcat('st:',num2str(st),'fromBuf:',num2str(obj.readStateFromBuffer())));
            %end
        end

        % provides access to upsample factor
        function [nSamp] = numSamplesPerStep( obj )
            nSamp = obj.T;
        end

        % provides access to upsample factor
        function [nBr] = numBranchesPerState( obj )
            nBr = obj.M;
        end

        % looks up next state given current state and input
        function [ nextState ] = nextState( obj, curSt, inpBit )
            nextState = obj.lutBrFrState(1+curSt,1+inpBit).endSt;
        end

        % looks up phase samples given present state and input
        function [ sRad ] = samplesRad( obj, curSt, inpBit )
            sRad = obj.lutBrFrState(1+curSt,1+inpBit).phaseRad;
        end

        % looks up all branches to a state
        function [ brTo ] = branchesToState( obj, endSt )
            brTo = obj.lutBrToState(1+endSt,:);
        end

         % looks up all branches from a state
         function [ brFr ] = branchesFromState( obj, begSt )
            brFr = obj.lutBrFrState(1+begSt,:);
        end

    end % public methods

    methods (Access = 'private')

        % called just once - so be explicit
        function writeStateToBuffer( obj, state )
            % state = P * [ u(n-1) + u(n-2)*M + ... + u(n-L)*M^(L-1) ]+ v
            stIdx = state;
            v = mod(stIdx,obj.hP);
            obj.stateBuf(obj.L) = v;
            stIdx = (stIdx - v)/obj.hP;
            for l=1:obj.L-1
                u = mod(stIdx,obj.M);
                obj.stateBuf(l) = u;
                stIdx = (stIdx - u)/obj.M;
            end
        end

        % called just once - so be explicit
        function state = readStateFromBuffer( obj )
            uIdx = 0;
            for l=1:obj.L-1
                uIdx = uIdx + obj.M^(l-1) * obj.stateBuf(l);
            end
            state = obj.hP*uIdx+ obj.stateBuf(obj.L);
        end

        % called just once - so be explicit
        function initialize( obj )
        
            obj.lutBrFrState = CPMTrellisBranch.empty(obj.numStates,0);
            phaseFullRad = zeros(1,obj.Tfull);
            % go over states
            for begSt = 0:obj.numStates-1
                for u = 0:obj.M-1 % input bit
                    obj.writeStateToBuffer(begSt);
                    v = obj.stateBuf(obj.L);
                    obj.stateBuf = [u, obj.stateBuf(1:obj.L-1)];
                    for t = 0:obj.Tfull-1
                        phRad = (pi*obj.h)*(2*v+(obj.M-1)*((obj.L-1)+((t+1)/obj.Tfull)));
                        for l=0:obj.L-1
                            phRad = phRad + 2*pi*obj.h*(2*obj.stateBuf(l+1)-(obj.M-1))*obj.phPulseFull(1+l*obj.Tfull+t);
                        end
                        phaseFullRad(1+t) = phRad;
                    end % for t
                    obj.stateBuf(end) = mod(obj.stateBuf(end)+ v,obj.hP);
                    endSt = obj.readStateFromBuffer();
                    % decimate to desired T
                    phaseRad = phaseFullRad(obj.decBegFromFull:obj.decFacFromFull:end);
                    % make branch
                    br = CPMTrellisBranch;
                    br.begSt = begSt;
                    br.inpSym = u;
                    br.endSt = endSt;
                    br.phaseRad = phaseRad;
                    % add branch to fwd lut  
                    obj.lutBrFrState(1+begSt,1+u) = br;
                end %u
            end %st

            obj.lutBrToState = CPMTrellisBranch.empty(obj.numStates,0);
            brCountByToState = zeros(1,obj.numStates);
            for begSt = 0:obj.numStates-1
                for u = 0:obj.M-1
                    br = obj.lutBrFrState(1+begSt,1+u);
                    brCount = brCountByToState(1+br.endSt);
                    obj.lutBrToState(1+br.endSt,brCount+1) = br;
                    brCountByToState(1+br.endSt) = brCount+1;
                end
            end

        end %initialize
    end %private methods 

end %CPMTiltedPhaseTrellis


