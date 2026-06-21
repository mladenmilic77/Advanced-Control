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
end