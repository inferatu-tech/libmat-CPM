classdef CPMTrellisBranch< handle

    properties (SetAccess = 'public', GetAccess = 'public')
        begSt; % beginning state of the branch
        inpSym;% input symbol that causes transition
        endSt; % ending state of the branch 
        phaseRad; % phase samples in radians  
    end

    methods (Access = 'public')
        function obj = CPMTrellisBranch()
        end
    end %public methods

end % CPMTrellisBranch