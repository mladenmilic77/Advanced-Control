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

    properties (Access = private)
        % Simulation parameters
        Ts % Sampling time

        % Continuous-time state-space matrices
        A  % State matrix
        B  % Input matrix
        E  % Disturbance matrix
        C  % Output matrix
        D  % Feedthrough matrix

        % Discrete-time state-space matrices
        Ad % State matrix
        Bd % Input matrix
        Ed % Disturbance matrix
        Cd % Output matrix
        Dd % Feedthrough matrix

        % Internal system states
        x  % State vector [omega; i]
    end

    properties (Constant)
        g = 9.81   % Gravitational acceleration
    end

    methods (Access = private, Static)
        function obj = DCMotor(parameters)
            %DCMOTOR Internal constructor for DC motor initialization.
            %   Creates and initializes a DC motor object using the
            %   provided physical parameters structure.
            arguments
                parameters struct
            end
            
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

            % Simulation parameters
            obj.Ts = 0.001;

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

            % Internal state vector
            obj.x = [obj.omega; obj.i];

            % Build state-space model
            obj.BuildStateSpaceModel();
        end

        function options = LoadOptionsChecker(options)
            %LOADOPTIONSCHECKER Validate and resolve load description.
            
            arguments
                options struct
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
        end
    end

    methods (Access = private)
        function BuildStateSpaceModel(obj)
            %BUILDSTATESPACEMODEL Build DC motor state-space model.
            %   Starting from:
            %       Electrical: u(t)=L*d(i)/dt+R*i+Ke*ω
            %       Mechanical: Kt*i(t)=J*d(ω)/dt+b*ω+τL
            %   States:
            %       x=[x1; x2]=[omega; i]
            %   Then:
            %       dx1/dt = (Kt/J)*x2-(b/J)*x1-(1/J)*τL
            %       dx2/dt = (1/L)*u-(R/L)*x2-(Ke/L)*x1
            %   Therefore:
            %       dx/dt = A*x+B*u+E*d
            %           y = C*x+D*u

            arguments
                obj 
            end

            % Matrices definition
            obj.A = [-obj.b/obj.J obj.Kt/obj.J;
                    -obj.Ke/obj.L -obj.R/obj.L];
            obj.B = [0;
                    1/obj.L];
            obj.E = [-1/obj.J;
                    0];
            obj.C = [1 0];
            obj.D = [0 0];

            % Controllability and Observability checks
            ControlTheory.ValidateControllability(obj.A, obj.B);
            ControlTheory.ValidateObservability(obj.A, obj.C);

            % Transformation
            sysc = ss(obj.A, [obj.B obj.E], obj.C, obj.D);
            sysd = c2d(sysc,obj.Ts,"zoh");

            % Discrete model
            obj.Ad = sysd.A;
            obj.Bd = sysd.B(:,1);
            obj.Ed = sysd.B(:,2);
            obj.Cd = sysd.C;
            obj.Dd = sysd.D(:,1);
        end
    end

    methods (Access = public, Static)
        function parameters = Parameters(options)
            %PARAMETERS Validate and prepare DC motor parameters.
            %   Creates a validated parameters structure and resolves
            %   optional load descriptions into equivalent load torque.

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

            parameters = DCMotor.LoadOptionsChecker(options);
        end
        
        function obj = Create(varargin)
            %CREATE Create and initialize a DC motor object.
            %   Constructs a DC motor model using validated physical
            %   parameters provided through name-value arguments.

            if nargin == 0
                parameters = DCMotor.Parameters();
            elseif nargin == 1 && isstruct(varargin{1})
                args = namedargs2cell(varargin{:});
                parameters = DCMotor.Parameters(args{:});
            else
                parameters = DCMotor.Parameters(varargin{:});
            end

            obj = DCMotor(parameters);
        end
    end

    methods (Access = public)
        function SetSampleTime(obj, Ts)
            %SETSAMPLETIME Set discrete simulation sample time.
            %   Updates sample time and rebuilds discrete 
            %   state-space model.

            arguments
                obj 
                Ts (1, 1) double {mustBePositive}
            end

            obj.Ts = Ts;
            obj.BuildStateSpaceModel;
        end

        function omega = Step(obj, V)
            %STEP Advance DC motor simulation by one discrete time step.
            %   Updates motor states using the discrete state-space model
            %   and returns current angular velocity.

            arguments
                obj 
                V (1,1) double % Input voltage
            end

            % State update: x[k+1] = Ad*x[k] + Bd*V[k] + Ed*τL[k]
            obj.x = obj.Ad * obj.x + obj.Bd * V + obj.Ed * obj.loadTorque;

            % Update object states
            obj.omega = obj.x(1);
            obj.i = obj.x(2);

            % Update derived quantities
            obj.alpha = (obj.Kt * obj.i - obj.b * obj.omega - obj.loadTorque) / obj.J;
            obj.theta = obj.theta + obj.omega * obj.Ts;
            obj.V = V;
            obj.eb = obj.Ke * obj.omega;
            obj.torque = obj.Kt * obj.i;

            % System output
            omega = obj.Cd * obj.x + obj.Dd * V;
        end

        function Reset(obj)
            %RESET Reset DC motor dynamic states.
            %   Restores all simulation states and outputs to zero.
            arguments
                obj 
            end

            % Dynamic states
            obj.theta = 0;
            obj.omega = 0;
            obj.alpha = 0;
            obj.i = 0;
            obj.V = 0;
            obj.eb = 0;

            % Internal state vector
            obj.x = [obj.omega; obj.i];
        end

        function x = GetState(obj)
            %GETSTATE Return current DC motor state vector.
            %   Returns internal state vector [omega; i].

            arguments
                obj 
            end

            x = obj.x;
        end

        function SetState(obj, omega, i)
            %SETSTATE Set DC motor internal state values.
            %   Updates angular velocity and armature current states.

            arguments
                obj 
                omega (1,1) double % Angular velocity
                i (1,1) double % Current
            end

            obj.omega = omega;
            obj.i = i;

            obj.x = [obj.omega; obj.i];

            obj.eb = obj.Ke * obj.omega;
            obj.torque = obj.Kt * obj.i;
        end
        
        function SetLoad(obj, options)
            %SETLOAD Set external load disturbance.
            %   Defines load using loadTorque, or calculates it from
            %   force/radius or mass/radius.

            arguments
                obj 
                options.loadTorque (1,1) double {mustBeNonnegative} = 0 % Load torque disturbance
                options.F (1,1) double = NaN;
                options.m (1,1) double = NaN;
                options.r (1,1) double = NaN;
            end

            options = DCMotor.LoadOptionsChecker(options);

            obj.loadTorque = options.loadTorque;
            obj.F = options.F;
            obj.m = options.m;
            obj.r = options.r;
        end
    end
end