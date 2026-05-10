classdef DCMotor < handle
    %DCMOTOR Discrete-time armature-controlled DC motor model.
    %   Simulates DC motor electrical and mechanical dynamics in
    %   state-space form using voltage as input and speed as output.
    %   Differential equations:
    %       Electrical: u(t)=L*d(i)/dt+R*i+Ke*ω
    %       Mechanical: Kt*i(t)=J*d(ω)/dt+b*ω+τL

    properties (Access = public)
        % Mechanical parameters
        J          % Motor inertia
        b          % Viscous friction coefficient
        Kt         % Torque constant

        % Optional load parameters
        r          % Shaft/pulley/wheel radius
        m          % Load mass
        F          % Load force
        loadTorque % External load torque

        % Electrical parameters
        R          % Armature resistance
        L          % Armature inductance
        Ke         % Back-EMF constant
    end

    properties (SetAccess = private)
        % Mechanical properties
        theta      % Angular position
        omega      % Angular velocity
        alpha      % Angular acceleration
        torque     % Torque

        % Electrical properties
        V          % Voltage
        i          % Current
        eb         % Back EMF
    end

    properties (Constant)
        g = 9.81   % Gravitational acceleration
    end

    methods (Access = private,Static)
        function obj = DCMotor(parameters)
            %DCMOTOR Internal constructor for DC motor initialization.
            %   Creates and initializes a DC motor object using the
            %   provided physical parameters structure.
            
            % Mechanical properties
            obj.J  = parameters.J;
            obj.b  = parameters.b;
            obj.Kt = parameters.Kt;

            % Electrical properties
            obj.R  = parameters.R;
            obj.L  = parameters.L;
            obj.Ke = parameters.Ke;

            % Optional load parameters
            obj.r = parameters.r;
            obj.m = parameters.m;
            obj.F = parameters.F;

            % Dynamic states
            obj.theta = 0;
            obj.omega = 0;
            obj.alpha = 0;
            obj.i = 0;
            obj.V = 0;
            obj.eb = 0;

            % Torque variables
            obj.torque = 0;
            obj.loadTorque = parameters.loadTorque;
        end
    end

    methods (Access = public, Static)
        function parameters = Parameters(options)

            arguments
                % Mechanical parameters
                options.J (1,1) double {mustBePositive} = 0.01;
                options.b (1,1) double {mustBeNonnegative} = 0.001;
                options.Kt (1,1) double {mustBePositive} = 0.01;

                % Electrical parameters
                options.R (1,1) double {mustBePositive} = 1;
                options.L (1,1) double {mustBePositive} = 0.5;
                options.Ke (1,1) double {mustBeNonnegative} = 0.01;

                % Optional load parameters
                options.r (1,1) double = NaN;
                options.m (1,1) double = NaN;
                options.F (1,1) double = NaN;
                options.loadTorque (1,1) double {mustBeNonnegative} = 0;
            end

            %Forwarded DC motor model verification
            loadTorquePresent = options.loadTorque > 0;
            rPresent = ~isnan(options.r);
            mPresent = ~isnan(options.m);
            FPresent = ~isnan(options.F);
            hasLoadDescription = rPresent || mPresent || FPresent;

            if hasLoadDescription

                if loadTorquePresent
                    error("MATLAB:badargs", ...
                        "Use either loadTorque directly or describe the load using (F;r) or (m;r).")
                elseif mPresent && FPresent
                    error("MATLAB:badargs", "Use either force F or mass m, not both.")
                elseif ~rPresent || ~(mPresent || FPresent)
                    error("MATLAB:badargs", ...
                        "Incomplete load description. Use either (F;r) or (m;r).");
                elseif options.r <= 0
                    error("DCMotor:InvalidRadius", "Radius r must be positive.")
                elseif FPresent
                    if options.F < 0
                        error("MATLAB:badargs", "Force F must be nonnegative.");
                    else
                        options.loadTorque = options.F * options.r;
                    end
                else
                    if options.m < 0
                        error("MATLAB:badargs", "Mass m must be nonnegative.");  
                    else
                        options.F = options.m * DCMotor.g;
                        options.loadTorque = options.F * options.r;
                    end
                end
                
            end

            parameters = options;
        end
        
        function obj = Create(options)

            arguments
                % Mechanical parameters
                options.J (1,1) double;
                options.b (1,1) double;
                options.Kt (1,1) double;

                % Electrical parameters
                options.R (1,1) double;
                options.L (1,1) double;
                options.Ke (1,1) double;

                % Optional load parameters
                options.r (1,1) double;
                options.m (1,1) double;
                options.F (1,1) double;
                options.loadTorque (1,1) double;
            end

            args = namedargs2cell(options);
            obj = DCMotor(DCMotor.Parameters(args{:}));
        end
    end
end