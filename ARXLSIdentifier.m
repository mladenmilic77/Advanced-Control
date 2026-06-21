classdef ARXLSIdentifier < handle
    %ARXLSIDENTIFIER ARX model parameter estimation using Least Squares.
    % Identifies discrete-time system models from input-output data.

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

    methods (Access = public)
        function yHat = Step(obj, u, y)
            %STEP Perform one recursive least squares update.
            %   Updates ARX parameter estimates using new input-output sample.
            
            arguments
                obj 
                u (1,1) double % System input
                y (1,1) double % Measured output
            end

            % Regression vector
            obj.phi = [-obj.yHistory; obj.uHistory];

            % Output prediction
            obj.yHat = obj.phi' * obj.theta;

            % Prediction error
            obj.error = y - obj.yHat;

            % RLS adaption gain vector
            gain = (obj.P * obj.phi)/(obj.lambda + obj.phi' * obj.P * obj.phi);

            % Parameter estimates update
            obj.theta = obj.theta + gain * obj.error;

            % Covariance matrix update
            obj.P = (1 / obj.lambda) * (obj.P - gain * (obj.phi' * obj.P));

            % Update histories
            obj.yHistory = [y; obj.yHistory(1:end-1)];
            obj.uHistory = [u; obj.uHistory(1:end-1)];

            % Return predicted output
            yHat = obj.yHat;
        end
    
        function Reset(obj)
            %RESET Reset recursive least squares estimator.
            %   Restores parameter estimates, covariance matrix,
            %   prediction variables and input-output histories.

            arguments
                obj 
            end

            parameterCount = obj.na + obj.nb;
            obj.theta = zeros(parameterCount, 1);
            obj.phi = zeros(parameterCount, 1);
            
            obj.P = 1000 * eye(parameterCount);
            
            obj.yHat = 0;
            obj.error = 0;
            
            obj.uHistory = zeros(obj.nb, 1);
            obj.yHistory = zeros(obj.na, 1);
        end
        
        function theta = GetParameters(obj)
            %GETPARAMETERS Return estimated ARX parameter vector.
            %   Returns current recursive least squares parameter estimates.

            arguments
                obj 
            end

            theta = obj.theta;
        end

        function yHat = GetPrediction(obj)
            %GETPREDICTION Return latest predicted output.
            %   Returns most recent ARX model output prediction.
            
            arguments
                obj 
            end

            yHat = obj.yHat;
        end

        function error = GetError(obj)
            %GETERROR Return latest prediction error.
            %   Returns difference between measured and predicted output.

            arguments
                obj 
            end

            error = obj.error;
        end
    end
end