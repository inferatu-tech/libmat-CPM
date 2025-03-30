classdef CPMTrellisPath 
  
    properties (SetAccess = 'public', GetAccess = 'public')
        stateSeq; % state sequence 
        inpSymSeq;% input symbol sequence 
        metric; % metric for the path 
    end

    methods (Access = 'public')
        function obj = CPMTrellisPath()
        end
    end %public methods

end % CPMTrellisBranch