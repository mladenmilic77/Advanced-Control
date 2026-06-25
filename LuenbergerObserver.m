classdef LuenbergerObserver < handle
    %LUENBERGEROBSERVER Generic discrete-time Luenberger state observer.
    % Estimates system states from input and output measurements.

    properties (Access = public)
        A % State matrix
        B % Input matrix
        C % Output matrix
        L % Observer gain
        Ts % Sample time
    end

    properties (SetAccess = private)
        xHat % Estimated states
        yHat % Estimated output
        error % y - yHat
    end

    methods (Access = private, Static)
        function obj = LuenbergerObserver(params)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                params struct % Observer parameters
            end

            obj.A = params.A;
            obj.B = params.B;
            obj.C = params.C;
            obj.L = params.L;
            obj.Ts = params.Ts;

            nx = size(obj.A, 1);
            ny = size(obj.C, 1);

            obj.xHat = zeros(nx, 1);
            obj.yHat = zeros(ny, 1);
            obj.error = zeros(ny, 1);
        end
    end

    methods (Access = public, Static)
        function obj = Create(params)
            %
            
            arguments
                params.A (:,:) double {mustBeFinite, mustBeReal}
                params.B (:,:) double {mustBeFinite, mustBeReal}
                params.C (:,:) double {mustBeFinite, mustBeReal}
                params.ObserverPoles (:,1) double {mustBeFinite, mustBeReal}
                params.Ts (1,1) double {mustBePositive, mustBeFinite}
            end

            ControlTheory.ValidateStateSpaceMatrices(params.A, ...
            params.B, params.C);

            ControlTheory.ValidateObservability(params.A, params.C);

            ControlTheory.ValidatePoleCount(params.A, params.ObserverPoles);

            params.L = place(params.A', params.C', params.ObserverPoles)';

            ControlTheory.ValidateObserverMatrices(params.A, ...
            params.C, params.L);

            obj = LuenbergerObserver(params);
        end
    end

    methods (Access=public)
        function Reset(obj)
            %RESET Reset observer estimates.
            % Restores estimated states, output and output error to zero.
            
            arguments
                obj 
            end

            nx = size(obj.A, 1);
            ny = size(obj.C, 1);

            obj.xHat = zeros(nx, 1);
            obj.yHat = zeros(ny, 1);
            obj.error = zeros(ny, 1);
        end

        function Initialize(obj, x0)
            %INITIALIZE Initialize observer state estimates.
            % Sets the initial state estimate and updates the estimated output.

            arguments
                obj 
                x0 (:, 1) double {mustBeFinite, mustBeReal} % Initial state
            end

            if numel(x0) ~= size(obj.A,1)
                error("MATLAB:sizeDimensionsMustMatch", ...
                    "Initial state estimate must have the same number of elements as system states."); %#ok<CPROP>
            end

            obj.xHat = x0;
            obj.yHat = obj.C * obj.xHat;
            obj.error = zeros(size(obj.C, 1), 1);
        end

        function xHat = Step(obj, u, y)
            %STEP Perform one discrete-time observer update.
            % Updates state estimate using input and measured output.

            arguments
                obj 
                u (:, 1) double {mustBeFinite, mustBeReal}
                y (:, 1) double {mustBeFinite, mustBeReal}
            end

            if numel(u) ~= size(obj.B, 2)
                error("MATLAB:sizeDimensionsMustMatch", ...
                    "Input vector u must have same number of rows as input matrix B has columns."); %#ok<CPROP>
            end

            if numel(y) ~= size(obj.C, 1)
                error("MATLAB:sizeDimensionsMustMatch", ...
                    "Output vector y must have same number of rows as output matrix C has rows."); %#ok<CPROP>
            end

            obj.error = y - obj.C * obj.xHat;
            obj.xHat = obj.A * obj.xHat + obj.B * u + obj.L * obj.error;
            obj.yHat = obj.C * obj.xHat;

            xHat = obj.xHat;
        end

        function SetPoles(obj, poles)
            %SETPOLES Set new discrete observer poles.
            % Recalculates observer gain matrix L using pole placement.

            arguments
                obj 
                poles (:, 1) double {mustBeFinite, mustBeReal} % Desired poles
            end

            ControlTheory.ValidatePoleCount(obj.A, poles);
            obj.L = place(obj.A', obj.C', poles)';
        end
    end
end