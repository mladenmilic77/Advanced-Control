clear; clc; close all;

%% Simulation setup
Ts = 0.001;
N  = 100000;

%% Create plant
motor = DCMotor.Create();
motor.SetSampleTime(Ts);

%% Create PWM
pwm = PWM.Create( ...
    supplyVoltage = 24, ...
    frequency = 10000, ...
    mode = "switching");

%% Create identifier
identifier = ARXLSIdentifier.Create( ...
    na = 2, ...
    nb = 2, ...
    lambda = 0.996, ...
    initialCovariance = 1000);

%% Create adaptive controller
adaptive = AdaptiveController.Create( ...
    plant = motor, ...
    pwm = pwm, ...
    identifier = identifier, ...
    Ts = Ts, ...
    reference = 5, ...
    minimumIdentificationSamples = 300, ...
    maximumIdentificationSamples = 1000, ...
    thetaTolerance = 1e-5, ...
    excitationLow = 0.08, ...
    excitationHigh = 0.18, ...
    excitationSwitchSamples = 100, ...
    uMin = 0, ...
    uMax = 1);

%% Logs
yLog = zeros(1, N);
uLog = zeros(1, N);
voltageLog = zeros(1, N);
dHatLog = zeros(1, N);
thetaDiffLog = zeros(1, N);
identificationLog = false(1, N);

%% Run simulation
for k = 1:N
    yLog(k) = adaptive.Step();

    uLog(k) = adaptive.uDuty;
    voltageLog(k) = adaptive.voltage;
    dHatLog(k) = adaptive.dHat;
    thetaDiffLog(k) = norm(adaptive.thetaDifference, Inf);
    identificationLog(k) = adaptive.identificationActive;
end

%% Find switch sample
switchIndex = find(~identificationLog, 1, "first");

%% Plot output
figure;
plot(yLog, "LineWidth", 1.5);
hold on;
plot(adaptive.reference * ones(1, N), "--", "LineWidth", 1.5);

if ~isempty(switchIndex)
    xline(switchIndex, "--");
end

grid on;
xlabel("Sample");
ylabel("Speed");
legend("Measured speed", "Reference", "Control started");
title("Adaptive controller output");

%% Plot duty
figure;
plot(uLog, "LineWidth", 1.5);

if ~isempty(switchIndex)
    xline(switchIndex, "--");
end

grid on;
xlabel("Sample");
ylabel("Duty cycle");
title("Duty cycle command");

%% Plot voltage
figure;
plot(voltageLog, "LineWidth", 1.5);

if ~isempty(switchIndex)
    xline(switchIndex, "--");
end

grid on;
xlabel("Sample");
ylabel("Voltage [V]");
title("PWM output voltage");

%% Plot disturbance estimate
figure;
plot(dHatLog, "LineWidth", 1.5);

if ~isempty(switchIndex)
    xline(switchIndex, "--");
end

grid on;
xlabel("Sample");
ylabel("Estimated disturbance");
title("Estimated disturbance");

%% Plot theta convergence
figure;
semilogy(thetaDiffLog, "LineWidth", 1.5);

if ~isempty(switchIndex)
    xline(switchIndex, "--");
end

grid on;
xlabel("Sample");
ylabel("Theta difference");
title("Identifier parameter convergence");