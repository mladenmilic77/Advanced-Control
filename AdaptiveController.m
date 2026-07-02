classdef AdaptiveController < handle
    %ADAPTIVECONTROLLER Identification + observer + controller supervisor.
    %   Coordinates:
    %       ARXLSIdentifier -> LuenbergerObserver -> StateSpaceController
    %       -> PWM -> Plant

    properties (Access = public)
        plant
        pwm
        identifier
        observer
        controller
        Ts
        reference
    end

    properties (SetAccess = private)
        identificationActive
        identificationDone

        theta
        previousTheta
        thetaDifference

        uDuty
        voltage
        y
        xHat
        dHat

        sampleCounter
        identificationCounter
    end

    properties (Access = private)
        identificationSamples
        minimumIdentificationSamples
        maximumIdentificationSamples
        thetaTolerance

        excitationLow
        excitationHigh
        excitationSwitchSamples

        controllerPoles
        observerPoles

        uMin
        uMax
    end

    methods (Access = private, Static)
        function obj = AdaptiveController(options)

            arguments
                options struct
            end
            
            obj.plant = options.plant;
            obj.pwm = options.pwm;
            obj.identifier = options.identifier;

            obj.observer = [];
            obj.controller = [];

            obj.Ts = options.Ts;
            obj.reference = options.reference;

            obj.minimumIdentificationSamples = options.minimumIdentificationSamples;
            obj.maximumIdentificationSamples = options.maximumIdentificationSamples;
            obj.thetaTolerance = options.thetaTolerance;

            obj.excitationLow = options.excitationLow;
            obj.excitationHigh = options.excitationHigh;
            obj.excitationSwitchSamples = options.excitationSwitchSamples;

            obj.controllerPoles = options.controllerPoles;
            obj.observerPoles = options.observerPoles;

            obj.uMin = options.uMin;
            obj.uMax = options.uMax;

            obj.identificationActive = true;
            obj.identificationDone = false;

            obj.uDuty = 0;
            obj.voltage = 0;
            obj.y = 0;
            obj.xHat = [];
            obj.dHat = 0;

            obj.sampleCounter = 0;
            obj.identificationCounter = 0;
        end
    end

    methods (Access = private)
        function uDuty = GenerateExcitation(obj)
            %GENERATEEXCITATION Generate identification duty cycle signal.

            arguments
                obj 
            end

            blockIndex = floor(obj.identificationCounter / obj.excitationSwitchSamples);

            if mod(blockIndex, 2) == 0
                uDuty = obj.excitationHigh;
            else
                uDuty = obj.excitationLow;
            end
        end

        function BuildControlLoop(obj)
            %BUILDCONTROLLOOP Build observer and controller from identified model.

            arguments
                obj 
            end

            [A, B, C, ~] = obj.identifier.GetStateSpaceModel();

            nx = size(A, 1);

            if isempty(obj.controllerPoles)
                obj.controllerPoles = linspace(0.25, 0.45, nx).';
            end

            Bd = B;

            nd = size(Bd, 2);

            if isempty(obj.observerPoles)
                obj.observerPoles = linspace(0.05, 0.20, nx + nd).';
            end

            obj.controller = StateSpaceController.Create(A = A, B = B, C = C, ...
                desiredPoles = obj.controllerPoles, uMin = obj.uMin, uMax = obj.uMax);

            obj.observer = LuenbergerObserver.Create(A = A, B = B, C = C, ...
                Bd = Bd, ObserverPoles = obj.observerPoles, Ts = obj.Ts, DisturbanceType = "constant");

            obj.observer.Initialize(zeros(nx, 1), 0);

            obj.xHat = obj.observer.Step(obj.uDuty, obj.y);
            obj.dHat = obj.observer.GetDisturbance();
            obj.identificationActive = false;
            obj.identificationDone = true;
        end
    end

    methods (Access = public, Static)
        function obj = Create(options)
            %CREATE Create adaptive controller object.
            
            arguments
                options.plant
                options.pwm PWM
                options.identifier ARXLSIdentifier
                
                options.Ts (1,1) double {mustBePositive, mustBeFinite} = 0.001
                options.reference (1,1) double {mustBeFinite, mustBeReal} = 5
                options.minimumIdentificationSamples (1,1) double {mustBeInteger, mustBePositive} = 100
                options.maximumIdentificationSamples (1,1) double {mustBeInteger, mustBePositive} = 3000
                options.thetaTolerance (1,1) double {mustBePositive, mustBeFinite} = 1e-5

                options.excitationLow (1,1) double {mustBeGreaterThanOrEqual(options.excitationLow,0), mustBeLessThanOrEqual(options.excitationLow,1)} = 0.08
                options.excitationHigh (1,1) double {mustBeGreaterThanOrEqual(options.excitationHigh,0), mustBeLessThanOrEqual(options.excitationHigh,1)} = 0.18
                options.excitationSwitchSamples (1,1) double {mustBeInteger, mustBePositive} = 100

                options.controllerPoles (:,1) double {mustBeFinite, mustBeReal} = []
                options.observerPoles (:,1) double {mustBeFinite, mustBeReal} = []
                options.uMin (1,1) double {mustBeFinite, mustBeReal} = 0
                options.uMax (1,1) double {mustBeFinite, mustBeReal} = 1
            end

            obj = AdaptiveController(options);
        end
    end

    methods (Access = public)
        function y = Step(obj)
            %STEP Run one adaptive control sample.

            arguments
                obj
            end

            obj.sampleCounter = obj.sampleCounter + 1;

            if obj.identificationActive
                obj.uDuty = obj.GenerateExcitation();
                
                obj.voltage = obj.pwm.Step(obj.uDuty, obj.Ts);
                obj.y = obj.plant.Step(obj.voltage);
                
                obj.identifier.Step(obj.uDuty, obj.y);
                
                obj.identificationCounter = obj.identificationCounter + 1;

                obj.theta = obj.identifier.GetParameters();

                if ~isempty(obj.previousTheta)
                    obj.thetaDifference = abs(obj.theta - obj.previousTheta);

                    thetaConverged = obj.identificationCounter >= obj.minimumIdentificationSamples && all(obj.thetaDifference < obj.thetaTolerance);
                    maxReached = obj.identificationCounter >= obj.maximumIdentificationSamples;

                    if thetaConverged || maxReached
                        obj.BuildControlLoop();
                    end
                end

                obj.previousTheta = obj.theta;
            else
                % Control phase

                obj.uDuty = obj.controller.Step(obj.xHat, obj.reference, obj.dHat);
                obj.voltage = obj.pwm.Step(obj.uDuty, obj.Ts);
                obj.y = obj.plant.Step(obj.voltage);
                obj.xHat = obj.observer.Step(obj.uDuty, obj.y);
                obj.dHat = obj.observer.GetDisturbance();
            end

            y = obj.y;
        end

        function SetReference(obj, reference)
            %SETREFERENCE Set desired reference vector.

            arguments
                obj
                reference (:,1) double {mustBeReal, mustBeFinite}
            end

            obj.reference = reference;
        end

        function Reset(obj)
            %RESET Reset adaptive controller and all internal components.
            %   Restores the adaptive controller to its initial state.

            arguments
                obj
            end

            % Reset subsystem objects
            obj.plant.Reset();
            obj.pwm.Reset();
            obj.identifier.Reset();

            if ~isempty(obj.observer)
                obj.observer.Reset();
            end

            if ~isempty(obj.controller)
                obj.controller.Reset();
            end

            % Reset supervisor flags
            obj.identificationActive = true;
            obj.identificationDone = false;

            % Reset internal signals
            obj.uDuty = 0;
            obj.voltage = 0;
            obj.y = 0;

            obj.xHat = [];
            obj.dHat = 0;

            % Reset identification monitoring
            obj.theta = [];
            obj.previousTheta = [];
            obj.thetaDifference = inf;

            obj.sampleCounter = 0;
            obj.identificationCounter = 0;
        end
    end
end