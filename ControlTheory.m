classdef ControlTheory
    %CONTROLTHEORY Utility class for control systems analysis and validation.
    % Provides static methods for validating Control Theory general stuff.

    methods(Access = public, Static)
        
        function ValidateStateSpaceMatrices(A, B, C, D)
            %VALIDATESTATESPACEMATRICES Validate state-space matrix dimensions.
            % Checks consistency of A, B, C and optional D matrices.

            arguments
                A (:,:) double {mustBeFinite, mustBeReal}
                B (:,:) double {mustBeFinite, mustBeReal}
                C (:,:) double {mustBeFinite, mustBeReal}
                D (:,:) double {mustBeFinite, mustBeReal} = []
            end

            nx = size(A,1);
            nu = size(B,2);
            ny = size(C,1);

            validateattributes(A, {'double'}, {'square'}, ...
            'ValidateStateSpaceMatrices', 'A');

            if size(B,1) ~= nx
                error("MATLAB:sizeDimensionsMustMatch", ...
                    "B must have same number of rows as A.");
            end

            if size(C,2) ~= nx
                error("MATLAB:sizeDimensionsMustMatch", ...
                    "C must have same number of columns as A.");
            end

            if ~isempty(D)
                if size(D,1) ~= ny
                    error("MATLAB:sizeDimensionsMustMatch", ...
                        "D must have same number of rows as C.");
                end

                if size(D,2) ~= nu
                    error("MATLAB:sizeDimensionsMustMatch", ...
                        "D must have same number of columns as B.");
                end
            end
        end

        function ValidateObserverMatrices(A, C, L)
            %VALIDATEOBSERVERMATRICES Validate observer model dimensions.
            % Checks consistency of state-space matrices and observer gain L.
            arguments
                A (:,:) double {mustBeFinite, mustBeReal}
                C (:,:) double {mustBeFinite, mustBeReal}
                L (:,:) double {mustBeFinite, mustBeReal}
            end

            nx = size(A,1);
            ny = size(C,1);
            
            ControlTheory.ValidateStateSpaceMatrices(A, zeros(nx,1), C);
            
            if size(L,1) ~= nx
                error("MATLAB:sizeDimensionsMustMatch", ...
                    "Observer gain L must have same number of rows as A.");
            end

            if size(L,2) ~= ny
                error("MATLAB:sizeDimensionsMustMatch", ...
                    "Observer gain L must have same number of columns as outputs.");
            end
        end

        function ValidatePoleCount(A, poles)
            %

            arguments
                A (:, :) double {mustBeFinite, mustBeReal}
                poles (:, 1) double {mustBeFinite}
            end
            
            nx = size(A, 1);

            validateattributes(A, {'double'}, {'square'}, ...
                'ValidatePoleCount', 'A');

            if numel(poles) ~= nx
                error("MATLAB:sizeDimensionsMustMatch", ...
                "Number of poles must be equal to number of states.");
            end
        end

        function ValidateControllability(A, B)
            %VALIDATECONTROLLABILITY Validate state-space model controllability.
            % Throws an error if the controllability matrix is rank deficient.

            arguments
                A (:, :) double {mustBeFinite, mustBeReal}
                B (:, :) double {mustBeFinite, mustBeReal}
            end

            ControlTheory.ValidateStateSpaceMatrices(A, B, eye(size(A,1)));

            if rank(ctrb(A, B)) < size(A, 1)
                error("MATLAB:rankDeficientMatrix", ...
                    "System is not fully controllable.");
            end
        end

        function ValidateObservability(A, C)
            %VALIDATEOBSERVABILITY Validate state-space model observability.
            % Throws an error if the observability matrix is rank deficient.
            
            arguments
                A (:, :) double {mustBeFinite, mustBeReal}
                C (:, :) double {mustBeFinite, mustBeReal}
            end

            nx = size(A, 1);

            ControlTheory.ValidateStateSpaceMatrices(A, zeros(nx,1), C);

            if rank(obsv(A, C)) < nx
                error("MATLAB:rankDeficientMatrix", ...
                    "System is not fully observable.");
            end
        end
    end
end