classdef ARXLSIdentifier < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here

    properties (Access = public)
        na     % Number of output regressors
        nb     % Number of input regressors
        lambda % Forgetting factor
    end

    properties (SetAccess = private)
        theta % Estimated parameters vector
        P     % Covariance matrix
        yHat  % Predicted output
        error % Predicted error
    end

    properties (Access = private)
        uHistory % Input history
        yHistory % Output history
        phi      % Regression vector
    end

    methods (Access = private, Static)
        function obj = ARXLSIdentifier(options)
            %ARXRLSIDENTIFIER Internal constructor for RLS identifier.
            %   Initializes ARX model order, forgetting factor,
            %   parameter vector, covariance matrix and data history.

            arguments
                options struct
            end

            obj.na = options.na;
            obj.nb = options.nb;
            obj.lambda = options.lambda;

            parameterCount = obj.na + obj.nb;
            obj.theta = zeros(parameterCount, 1);
            obj.P = options.initialCovariance * eye(parameterCount);

            obj.yHat = 0;
            obj.error = 0;

            obj.uHistory = zeros(obj.nb, 1);
            obj.yHistory = zeros(obj.na, 1);
            obj.phi = zeros(parameterCount, 1);
        end
    end

    methods (Access = public, Static)
        function obj = Create(options)
            %CREATE Create and initialize ARX RLS identifier.
            %   Constructs recursive least squares identifier
            %   using selected ARX order and forgetting factor.

            arguments
                options.na (1,1) double {mustBeInteger, mustBePositive} = 2
                options.nb (1,1) double {mustBeInteger, mustBePositive} = 2
                options.lambda (1,1) double {mustBeGreaterThan(options.lambda,0), mustBeLessThanOrEqual(options.lambda, 1)} = 0.98
                options.initialCovariance (1,1) double {mustBePositive} = 1000
            end

            obj = ARXLSIdentifier(options);
        end
    end
end