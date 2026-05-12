classdef PWM < handle
    %PWM Pulse Width Modulation model for DC motor voltage control.
    %   Converts duty cycle command into motor input voltage using
    %   average or switching PWM behavior.

    properties (Access = public)
        supplyVoltage % DC voltage
        frequency     % DC switching frequency
        mode          % PWM mode:"average" or "switching"
    end

    properties (SetAccess = private)
        dutyCycle     % Current duty cycle
        outputVoltage % Voltage sent to motor
        time          % Internal PWM time
    end

    properties (Access = private)
        period        % PWM period
    end

    methods (Access = private, Static)
        function obj = PWM(options)
            %PWM Internal constructor for PWM model initialization.
            %   Creates PWM object using validated configuration options.
            arguments
                options struct
            end

            obj.supplyVoltage = options.supplyVoltage;
            obj.frequency = options.frequency;
            obj.mode = options.mode;
            
            obj.period = 1 / obj.frequency;

            obj.dutyCycle = 0;
            obj.outputVoltage = 0;
            obj.time = 0;
        end
    end

    methods (Static)
        function obj = Create(options)
            %CREATE Create and initialize PWM object.
            %   Constructs PWM model using supply voltage, switching
            %   frequency and selected operating mode.

            arguments
                options.supplyVoltage (1,1) double {mustBePositive} = 24
                options.frequency (1,1) double {mustBePositive} = 10E3
                options.mode (1,1) string {mustBeMember(options.mode, ["average", "switching"])} = "average"
            end

            obj = PWM(options);
        end
    end

    methods
        function voltage = Step(obj, dutyCycle, Ts)
            %STEP Advance PWM model by one simulation step.
            %   Converts duty cycle into output voltage according to the
            %   selected PWM mode and updates internal time.

            arguments
                obj 
                dutyCycle (1,1) double {mustBeBetween(dutyCycle,0,1)}
                Ts (1,1) double {mustBePositive}
            end

            obj.dutyCycle = dutyCycle;

            switch obj.mode
                case "average"
                    obj.outputVoltage = obj.dutyCycle * obj.supplyVoltage;
                case "switching"
                    localTime = mod(obj.time, obj.period);
                    onTime = obj.dutyCycle * obj.period;

                    if localTime < onTime
                        obj.outputVoltage = obj.supplyVoltage;
                    else
                        obj.outputVoltage = 0;
                    end   
            end

            obj.time = obj.time + Ts;
            voltage = obj.outputVoltage;
        end

        function Reset(obj)
            %RESET Reset PWM dynamic states.
            %   Restores duty cycle, output voltage and internal time to zero.
            
            arguments
                obj 
            end

            obj.dutyCycle = 0;
            obj.outputVoltage = 0;
            obj.time = 0;
        end
    end
end