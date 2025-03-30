classdef CPMBitMapper

    properties (SetAccess = 'private', GetAccess = 'public')
        numBitsPerSym;  
    end

    methods (Access = 'public')
        % construct with modulation alphabet size
        function obj = CPMBitMapper( nbitsPerSym )  
            assert(nbitsPerSym>=1,'[BitMapper]: nbitsPerSym must be at least 1...');
            obj.numBitsPerSym = nbitsPerSym;
        end
        
        % maps a bit sequence to a symbol sequence 
        function [ symSeq ] = map( obj, bitSeq )
            bitSeq = bitSeq(:);
            if ~isrow(bitSeq) 
                bitSeq = bitSeq';
            end
            numBits = length(bitSeq);
            assert(mod(numBits,obj.numBitsPerSym)==0,'[BitMapper::map] input bit sequence does not map to an integer number of symbols...');
            numSyms = numBits/obj.numBitsPerSym;
            bitsBySym = reshape(bitSeq,obj.numBitsPerSym,numSyms)';
            symSeq = sum(bitsBySym .* 2.^(0:obj.numBitsPerSym-1), 2)';
        end

        % maps a symbol sequence to a bit sequence 
        function [ bitSeq ] = demap( obj, symSeq )
            symSeq = symSeq(:);
            if ~iscolumn(symSeq) 
                symSeq = symSeq';
            end
            numSyms = length(symSeq);
            bitsBySym = zeros(numSyms, obj.numBitsPerSym );
            for ibit = 0:obj.numBitsPerSym-1
                thisBitBySym = mod(symSeq,2);
                bitsBySym(:,ibit+1) = thisBitBySym;
                symSeq = (symSeq - thisBitBySym)/2;
            end
            bitSeq = reshape(bitsBySym',1,numSyms*obj.numBitsPerSym);
        end

    end %public methods

end % BitMapper


