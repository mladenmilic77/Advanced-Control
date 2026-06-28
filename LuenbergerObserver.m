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
        dHat % Estimated disturbance
        nx % Number of states
        ny % Number of outputs
        nd % Number of disturbance states
    end

    methods (Access = private, Static)
        function obj = LuenbergerObserver(params)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            arguments
                params struct % Observer parameters
            end

            obj.nx = size(params.A, 1);
            obj.ny = size(params.C, 1);
            obj.nd = size(params.Ad, 1);
            obj.Ts = params.Ts;

            obj.A = [params.A, params.Bd; zeros(obj.nd, obj.nx), params.Ad];
            obj.B = [params.B; zeros(obj.nd, size(params.B, 2))];
            obj.C = [params.C, zeros(obj.ny, obj.nd)];
            obj.L = params.L;

            n_total = obj.nx + obj.nd;
            obj.xHat = zeros(n_total, 1);
            obj.yHat = zeros(obj.ny, 1);
            obj.error = zeros(obj.ny, 1);
            obj.dHat = zeros(obj.nd, 1);
        end

        function [Ad, Bd_expanded] = ComputeDisturbanceDynamics(type, Bd, Ts, Freq, nx)
            %COMPUTEDISTURBANCEDYNAMICS Generates the discrete-time state-space transition
            % matrices for the selected disturbance model.
            
            arguments
                type string
                Bd (:,:) double
                Ts (1,1) double
                Freq (1,1) double
                nx (1,1) double
            end

            nd_inputs = size(Bd, 2);

            switch lower(type)
                case "constant"
                    % d(k+1)=d(k)
                    Ad=eye(nd_inputs);
                    Bd_expanded = Bd;
                case "ramp"
                    % d(k+1) = d(k) + Ts*v(k)
                    % v(k+1) = v(k)
                    Ad_single = [1, Ts; 0, 1];
                    Ad = blkdiag(Ad_single); % Fallback
                    
                    if nd_inputs > 1
                        Ad_cells = repmat({Ad_single}, 1, nd_inputs);
                        Ad = blkdiag(Ad_cells{:});
                    end

                    % Scaling Bd
                    Bd_expanded = zeros(nx,nd_inputs * 2);
                    Bd_expanded(:, 1:2:end) = Bd;
                case "sinusoidal"
                    omega = 2 * pi * Freq;
                    Ad_single = [cos(omega * Ts),  sin(omega * Ts); -sin(omega * Ts), cos(omega * Ts)];
                    Ad = blkdiag(Ad_single); % Fallback

                    if nd_inputs > 1
                        Ad_cells = repmat({Ad_single}, 1, nd_inputs);
                        Ad = blkdiag(Ad_cells{:});
                    end
                    
                    Bd_expanded = zeros(nx, nd_inputs * 2);
                    Bd_expanded(:, 1:2:end) = Bd;
            end
        end
    end

    methods (Access = public, Static)
        function obj = Create(params)
            %
            
            arguments
                params.A (:,:) double {mustBeFinite, mustBeReal}
                params.B (:,:) double {mustBeFinite, mustBeReal}
                params.C (:,:) double {mustBeFinite, mustBeReal}
                params.Bd (:,:) double {mustBeFinite, mustBeReal}
                params.ObserverPoles (:,1) double {mustBeFinite, mustBeReal}
                params.Ts (1,1) double {mustBePositive, mustBeFinite}

                % Disturbance configuration parameters
                params.DisturbanceType (1,1) string {mustBeMember(params.DisturbanceType, ["constant","ramp","sinusoidal"])} = "constant"
                params.DisturbanceFrequency (1,1) double {mustBeReal,mustBeFinite} = 0
            end

            if size(params.Bd,1) ~= size(params.A,1)
                error("MATLAB:sizeDimensionsMustMatch", ...
                    "Bd must have same number of rows as A.");
            end

            params.nx = size(params.A, 1);

            [params.Ad, params.Bd] = LuenbergerObserver.ComputeDisturbanceDynamics( ...
                params.DisturbanceType, params.Bd, params.Ts, params.DisturbanceFrequency, params.nx);

            n_dist = size(params.Ad, 1);
            A_aug = [params.A, params.Bd; zeros(n_dist, params.nx),  params.Ad];
            C_aug = [params.C, zeros(size(params.C,1), n_dist)];
            B_aug = [params.B; zeros(n_dist, size(params.B,2))];

            ControlTheory.ValidateStateSpaceMatrices(A_aug, B_aug, C_aug);

            ControlTheory.ValidateObservability(A_aug, C_aug);

            ControlTheory.ValidatePoleCount(A_aug, params.ObserverPoles);

            params.L = place(A_aug', C_aug', params.ObserverPoles)';

            ControlTheory.ValidateObserverMatrices(A_aug, C_aug, params.L);

            % Not great need some fixing here disturbance part made a mess
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

            n_total = obj.nx + obj.nd;

            obj.xHat = zeros(n_total, 1);
            obj.yHat = zeros(obj.ny, 1);
            obj.error = zeros(obj.ny, 1);
            obj.dHat = zeros(obj.nd, 1);
        end

        function Initialize(obj, x0, d0)
            %INITIALIZE Initialize observer state estimates.
            % Sets the initial state estimate and updates the estimated output.

            arguments
                obj 
                x0 (:, 1) double {mustBeFinite, mustBeReal} % Initial state
                d0 (:, 1) double {mustBeFinite, mustBeReal} = [] % Disturbance Initial state
            end

            if numel(x0) ~= obj.nx
                error("MATLAB:sizeDimensionsMustMatch", ...
                    "Initial state estimate must have the same number of elements as system states."); %#ok<CPROP>
            end

            if isempty(d0)
                d0 = zeros(obj.nd, 1);
            elseif numel(d0) ~= obj.nd
                error("MATLAB:sizeDimensionsMustMatch", ...
                    "Initial disturbance vector must match configured disturbance state size."); %#ok<CPROP>
            end

            obj.xHat = [x0; d0];
            obj.yHat = obj.C * obj.xHat;
            obj.error = zeros(obj.ny, 1);
            obj.dHat = d0;
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

            if numel(y) ~= obj.ny
                error("MATLAB:sizeDimensionsMustMatch", ...
                    "Output vector y must have same number of rows as output matrix C has rows."); %#ok<CPROP>
            end

            obj.error = y - obj.C * obj.xHat;
            obj.xHat = obj.A * obj.xHat + obj.B * u + obj.L * obj.error;
            obj.yHat = obj.C * obj.xHat;

            xHat = obj.xHat(1:obj.nx);
            obj.dHat = obj.xHat(obj.nx+1:end);
        end

        function dHat = GetDisturbance(obj)
            %

            arguments
                obj 
            end

            dHat = obj.dHat;
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