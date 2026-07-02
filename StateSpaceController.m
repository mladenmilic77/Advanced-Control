classdef StateSpaceController < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = public)
        K  % State feedback gain
        Kr % Reference/feedforward gain
    end

    properties (SetAccess = private)
        u    % Current controller output
        uMin % Minimum control signal
        uMax % Maximum control signal
        nx   % Number of states
        nr   % Number of reference inputs
        nu   % Number of control outputs
    end

    methods (Access = private, Static)
        function obj = StateSpaceController(options)
            %STATESPACECONTROLLER Internal constructor.
            %   Creates controller object using already computed controller
            %   gains and validated options.

            arguments
                options struct
            end
            
            obj.K = options.K;
            obj.Kr = options.Kr;
            obj.uMin = options.uMin;
            obj.uMax = options.uMax;
            
            obj.nu = size(obj.K, 1);
            obj.nx = size(obj.K, 2);
            obj.nr = size(obj.Kr, 2);
            obj.u = zeros(obj.nu, 1);
        end
    end

    methods (Access = public, Static)
        function obj = Create(options)
            %CREATE Create and initialize state-space controller object.
            %   Computes K from desired poles and Kr for reference tracking.
            
            arguments
                options.A (:,:) double {mustBeReal, mustBeFinite}
                options.B (:,:) double {mustBeReal, mustBeFinite}
                options.C (:,:) double {mustBeReal, mustBeFinite}
                options.desiredPoles (:,1) double {mustBeReal, mustBeFinite}
                options.uMin (:,1) double {mustBeReal, mustBeFinite} = 0
                options.uMax (:,1) double {mustBeReal, mustBeFinite} = 1
            end

            nx = size(options.A, 1);

            if size(options.A, 2) ~= nx
                error("MATLAB:sizeDimensionsMustMatch", ...
                    "A must be square.");
            end

            if length(options.desiredPoles) ~= nx
                error("MATLAB:sizeDimensionsMustMatch", ...
                    "Number of desired poles must match number of states.");
            end

            if options.uMin > options.uMax
                error("MATLAB:invalidInput", ...
                    "uMin must be less than or equal to uMax.");
            end

            options.K = place(options.A, options.B, options.desiredPoles);
            options.Kr = inv(options.C * ...
                ((eye(nx) - options.A + options.B * options.K) \ options.B));
            
            obj = StateSpaceController(options);
        end
    end

    methods (Access = public)
        function uClamped = Step(obj, x, r, d)
            %STEP Compute controller output.
            %   Calculates control signal using:
            %   u[k] = -K*xHat[k] + Kr*r[k] - dHat[k]

            arguments
                obj 
                x (:,1) double {mustBeReal, mustBeFinite}
                r (:,1) double {mustBeReal, mustBeFinite}
                d (:,1) double {mustBeReal, mustBeFinite}
            end

            u = -obj.K * x + obj.Kr * r - d; %#ok<PROP>

            obj.u = min(max(u, obj.uMin), obj.uMax); %#ok<PROP>

            uClamped = obj.u;
        end

        function Reset(obj)
            %RESET Reset controller output.
            %   Restores stored control signal to zero.

            arguments
                obj
            end

            obj.u = zeros(obj.nu, 1);
        end
    end
end