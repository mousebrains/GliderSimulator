% Rip through G3 deployments and get the energy for various operations

myDir = fileparts(mfilename('fullpath'));
figDir = myDir;
dataDir = '/Users/pat/Dropbox/Gliders';

deployments = { ...
    '07-24-2018.osu683_KECK', ...
    '07-24-2018.osu684_KECK', ...
    '10-01-2018.osu684_NHline_test_deployment', ...
    '01-25-2019.osu685_NHline_test_deployment', ...
    '04-28-2019.osu683_PugetSound_deployment', ...
    '04-28-2019.osu685_PugetSound_deployment', ...
    '05-26-2019.osu683_GulfMexico_deployment', ...
    '05-26-2019.osu684_GulfMexico_deployment', ...
    '05-26-2019.osu685_GulfMexico_deployment', ...
    '05-26-2019.osu686_GulfMexico_deployment'};

if ~exist('a', 'var') || ~istable(a)
    a = table();
    for index = 1:numel(deployments)
        b = loadDeployment(index, dataDir, deployments{index});
        a = [a; b];
    end % for index
end % if

calcSurfaceEnergy(1, figDir, a); % extract the surface energy as a function of dtIridium

fitEnergy = calcDiveClimbPower(10, figDir, a); % Estimate the dive/climb and thruster power usage and return residual

x = calcOilPower(20, figDir, a, max(0, a.energy1)); % Estimate power to move one CC out and in

function x = calcSurfaceEnergy(fig, figDir, a)
% A lot is going on at the same time at the surface, so rather than sort
% out the individual components just estimate the energy as
% energy = constant + slope * dtIridium
% so constant represent the base energy
% while slope represents the energy/second of Iridium call time
tv = tic;
tbl = table();
tbl.t = a.t;
tbl.dt = a.dt;
tbl.energy = a.energy1; % Use the coulomb counter energy, it should be more accuracte
tbl.surface = a.surface;
tbl.iridium = a.c_iridium_on == 1;
tbl.index = a.index;
tbl.grp = findgroups(cumsum([0; abs(diff(a.surface))]));

x = rowfun(@calcSurfRowFun, tbl, ...
    'InputVariables', {'t', 'dt', 'energy', 'surface', 'iridium'}, ...
    'GroupingVariables', {'grp', 'index'}, ...
    'SeparateInputs', true, ...
    'OutputVariableNames', {'t0', 't1', 't2', 't3', 'dt', 'dtIridium', 'energy'}, ...
    'OutputFormat', 'table');

q = (x.dtIridium > 0) & (x.dtIridium < 4000);

indices = find(q);
x.dtPrev = nan(size(x.grp));
x.dtPrev(q) = seconds(x.t3(indices-1) - x.t2(indices-1));


x = x(q,:);
x.power = x.energy ./ x.dt;

fprintf('Surface energy fit\n');
mdl = fitlm(x, 'linear', ...
    'ResponseVar', 'energy', ...
    'PredictorVars', 'dtIridium', ...
    'RobustOpts', true)

fprintf('Fit of previous down time to Iridium call time\n');
mdl1 = fitlm(x, 'linear', ...
    'ResponseVar', 'dtIridium', ...
    'PredictorVars', 'dtPrev', ...
    'RobustOpts', true)

fprintf('Took %.2f seconds to calculate the surface energy\n', toc(tv));

y = feval(mdl, x.dtIridium);

figure(fig);
[cnts0, bins] = histcounts(x.energy, linspace(0, 20000, 50));
cnts1 = histcounts(y, bins);
abins = (bins(1:end-1) + bins(2:end)) / 2;
plot(abins, cnts0, '-', abins, cnts1, '-');
grid on;
xlabel('Energy per surfacing (Joules)');
ylabel('Counts');
legend('Measured energy/surfacing (Joules)', ...
    sprintf('energy=%.0f+%.2f*callTime', mdl.Coefficients.Estimate), ...
    'Location', 'best');

if ~isempty(figDir)
    print(fullfile(figDir, 'SurfaceEnergy.pdf'), '-dpdf', '-fillpage')
end % if
end % calcSurfaceEnergy

function [t0, t1, t2, t3, timeTotal, timeIridium, energyTotal] = calcSurfRowFun(t, dt, energy, q, iridium)
t0 = NaT;
t1 = t0;
t2 = t0;
t3 = t0;
timeTotal = 0;
timeIridium = 0;
energyTotal = 0;
    
if isempty(t), return; end
if sum(q) == 0
    t2 = min(t);
    t3 = max(t);
    return;
end % if

t0 = t(1);
t1 = t(end);
t2 = min(t(q));
t3 = max(t(q));
timeTotal = seconds(max(t(q), [], 'omitnan') - min(t(q), [], 'omitnan'));
timeIridium = sum(dt(iridium), 'omitnan');
energyTotal = sum(energy(q), 'omitnan');
end % myFun

function fitEnergy = calcDiveClimbPower(fig, figDir, a)
% Estimate the dive/climb base power and the thruster power
tbl = table();
tbl.dive = a.diving;
tbl.climb = a.climbing;
tbl.thruster = (a.c_thruster_on > 0) | (a.m_thruster_current > 0);
tbl.thrusterPower = tbl.thruster .* a.m_thruster_power;
tbl.qOil = a.m_is_de_pump_moving > 0;
tbl.power = a.power;

q = tbl.dive | tbl.climb;

fprintf('Dive/Climb Power\n');
mdl = fitlm(tbl, 'linear', ...
    'Intercept', false, ...
    'ResponseVar', 'power', ...
    'PredictorVars', {'dive', 'climb'}, ...
    'RobustOpts', true, ...
    'Exclude', ~q | tbl.qOil | tbl.thruster)

tbl.resid = tbl.power - feval(mdl, tbl.dive, tbl.climb);

mdl1 = fitlm(tbl, 'linear', ...
    'Intercept', false, ...
    'ResponseVar', 'power', ...
    'PredictorVars', {'thruster', 'thrusterPower'}, ...
    'RobustOpts', true, ...
    'Exclude', ~q | tbl.qOil)

fitPower = zeros(size(tbl.power));
fitPower(q) = feval(mdl, tbl.dive(q), tbl.climb(q)) ...
    + feval(mdl1, tbl.thruster(q), tbl.thrusterPower(q));
fitEnergy = fitPower .* a.dt;

% Plots

qDive  = tbl.dive  & ~tbl.qOil & ~tbl.thruster;
qClimb = tbl.climb & ~tbl.qOil & ~tbl.thruster;
yD = feval(mdl, tbl.dive(qDive), tbl.climb(qDive));
yC = feval(mdl, tbl.dive(qClimb), tbl.climb(qClimb));

figure(fig);
[cnts0, bins] = histcounts(tbl.power(qDive), linspace(0,4,20));
cnts1 = histcounts(tbl.power(qClimb), bins);
cntsD = histcounts(yD, bins);
cntsC = histcounts(yC, bins);

abins = (bins(1:end-1) + bins(2:end)) / 2;
plot(abins, cnts0, '-', ...
    abins, cntsD, '-', ...
    abins, cnts1, '-', ....
    abins, cntsC, '-');
grid on;
xlabel('Power (Watts)');
ylabel('Counts');
legend(sprintf('Measured dive power mean=%.2f median=%.2f', ...
    mean(tbl.power(qDive), 'omitnan'), ...
    median(tbl.power(qDive), 'omitnan')), ...
    sprintf('Fit dive power=%.2f', mdl.Coefficients.Estimate(1)), ...
    sprintf('Measured climb power mean=%.2f median=%.2f', ...
    mean(tbl.power(qClimb), 'omitnan'), ...
    median(tbl.power(qClimb), 'omitnan')), ...
    sprintf('Fit climb power=%.2f', mdl.Coefficients.Estimate(2)), ...
    'Location', 'northeast');
title('Power with not oil nor thruster');

qThruster = q & tbl.thruster & ~tbl.qOil;
yThruster = feval(mdl1, tbl.thruster(qThruster), tbl.thrusterPower(qThruster));

figure(fig+1)
[cnts0, bins] = histcounts(tbl.resid(qThruster), linspace(0,20,40));
cnts1 = histcounts(yThruster, bins);
cnts2 = histcounts(tbl.thrusterPower(qThruster), bins);
abins = (bins(1:end-1) + bins(2:end)) / 2;
plot(abins, cnts0, '-', ...
    abins, cnts1, '-', ...
    abins, cnts2, '-');
grid on;
xlabel('Power (Watts)');
ylabel('Counts');
legend('Measured Thruster Power', ...
    sprintf('power=%.2f+%.4f*m_thruster_power', mdl1.Coefficients.Estimate), ...
    'm_thruster_power', ...
    'Location', 'best', 'Interpreter', 'none');
title('Power-Dive/Climb fit when thruster is on but not the oil pump');

if ~isempty(figDir)
    print(fig, fullfile(figDir, 'DiveClimb.pdf'), '-dpdf', '-fillpage');
    print(fig+1, fullfile(figDir, 'Thruster.pdf'), '-dpdf', '-fillpage');
end % if
end % calcDiveClimbPower

function x = calcOilPower(fig, figDir, a, energy)
tbl = table();
tbl.t = a.t;
tbl.energy = energy;
tbl.depth = a.m_depth;
tbl.qOil = a.m_is_de_pump_moving > 0;
tbl.oilVol = a.m_de_oil_vol;
tbl.index = a.index;
tbl.grp = findgroups(cumsum([0; abs(diff(tbl.qOil))]));

x = rowfun(@oilRowFun, tbl, ...
    'InputVariables', {'t', 'energy', 'qOil', 'oilVol', 'depth'}, ...
    'GroupingVariables', {'grp', 'index'}, ...
    'SeparateInputs', true, ...
    'OutputVariableNames', {'t', 'dt', 'dVol', 'energy', 'depth'}, ...
    'OutputFormat', 'table');

q = (x.energy > 0) & (x.dt > duration(0,0,0));
x = x(q,:);

x.power = x.energy ./ seconds(x.dt);

qOut = x.dVol > 0;
qIn  = x.dVol < 0;

fprintf('Oil out\n');
mdlOut = fitlm(x, 'interactions', ...
    'RobustOpts', true, ...
    'ResponseVar', 'energy', ...
    'PredictorVars', {'dVol', 'depth'}, ...
    'Exclude', ~qOut)

fprintf('Oil in\n');
mdlIn = fitlm(x, 'interactions', ...
    'RobustOpts', true, ...
    'ResponseVar', 'energy', ...
    'PredictorVars', {'dVol', 'depth'}, ...
    'Exclude', ~qIn)

yOut = feval(mdlOut, x.dVol(qOut), x.depth(qOut));
yIn = feval(mdlIn, x.dVol(qIn), x.depth(qIn));

[cnts0, bins] = histcounts(x.energy(qOut), 0:50:5000);
cnts1 = histcounts(yOut, bins);
abins = (bins(1:end-1) + bins(2:end)) / 2;

figure(fig);
plot(abins, cnts0, '-', ...
    abins, cnts1, '-');
grid on;
xlabel('Energy Per Pump Out Cycle (Joules)');
ylabel('Counts');
legend('Measured pump out', ...
    sprintf('Fit pump out=%.3f*dVol+%.3g*dVol*depth', mdlOut.Coefficients.Estimate([2,4])), ...
    'Location', 'best');
title('Energy per oil pump out cycle');

[cnts0, bins] = histcounts(x.energy(qIn), 0:5:500);
cnts1 = histcounts(yIn, bins);
abins = (bins(1:end-1) + bins(2:end)) / 2;

figure(fig+1);
plot(abins, cnts0, '-', ...
    abins, cnts1, '-');
grid on;
xlabel('Energy Per Pump In Cycle (Joules)');
ylabel('Counts');
legend('Measured pump in', ...
    sprintf('Fit pump in=%.3f*dVol+%.3g*dVol*depth', mdlIn.Coefficients.Estimate([2,4])), ...
    'Location', 'best');
title('Energy per oil pump in cycle');

if ~isempty(figDir)
    print(fig, fullfile(figDir, 'OilOut.pdf'), '-dpdf', '-fillpage');
    print(fig+1,fullfile(figDir, 'OilIn.pdf'), '-dpdf', '-fillpage');
end % if
end % calcOilPower

function [t, dt, dVol, energyTotal, depthCenter] = ...
    oilRowFun(t, energy, q, oilVol, depth)

if (numel(t) < 3) || (sum(q) == 0)
    dt = 0;
    dVol = 0;
    energyTotal = 0;
    depthCenter = 0;
    t = NaT;
    return;
end % if

qq = q  | [false; q(1:end-1)]; % Extend by one slot
qq = qq | [false; false; q(1:end-2)]; % Extend by second slot

oilVol = oilVol(qq & ~isnan(oilVol));

dt = max(t(qq)) - min(t(qq));
t = min(t(qq));
energyTotal = sum(energy(qq), 'omitnan');
dVol = sign(oilVol(end) - oilVol(1)) * (max(oilVol) - min(oilVol));
depthCenter = median(depth(qq), 'omitnan');
end % oilRowFun

function data = loadDeployment(index, dataDir, deployment)
tv = tic;
fn = osgl_pass0_filename(dataDir, deployment, 'dbd');
% a = osgl_read_netCDF(fn, ...
%     'm_present_time', 'm_de_oil_vol', 'm_pitch', 'm_depth', ...
%     'm_is_battpos_moving', 'm_is_de_pump_moving', ...
%     'm_air_pump', 'c_iridium_on', ...
%     'c_thruster_on', 'm_thruster_current', 'm_thruster_power', ...
%     'm_bms_aft_current', 'm_bms_ebay_current', ...a
%     'm_bms_pitch_current', 'm_bms_main_battery_voltage', ...
%     'm_vacuum_air_bag_inflated', 'm_vacuum_air_pump_on', ...
%     'm_vacuum_change_since_air_pump_on', 'm_air_fill', ...
%     'm_coulomb_amphr');
a = osgl_read_netCDF(fn);
[~, ix] = unique(a.m_present_time);
ix = ix(~isnan(a.m_present_time(ix)));

names = fieldnames(a);
n = numel(a.m_present_time);

data = table();

for i = 1:numel(names)
    name = names{i};
    if isequal(name, 'uniqueID'), continue; end
    if size(a.(name),1) == n
        data.(name) = a.(name)(ix);
    end % if
end % for i

data.t = datetime(data.m_present_time, 'ConvertFrom', 'posixtime');
data.index = index + zeros(size(data.t));

data.dt = [nan; diff(data.m_present_time)];
data.dt(1) = median(data.dt(2:end), 'omitnan'); 

msk = ~isnan(data.m_depth);
data.m_depth(~msk) = interp1(data.t(msk), data.m_depth(msk), data.t(~msk));

data.dOilVol = [nan; diff(data.m_de_oil_vol)];

% Get a power estimate using the coulomb counter
% It has a finite precision, so interpolate from when it
% ticks over to the next time it ticks over
msk = ~isnan(data.m_coulomb_amphr_total);
t = data.t(msk);
Q = data.m_coulomb_amphr_total(msk);
Q = Q - min(Q);
dQ = [0; diff(Q)];
q = dQ ~= 0; % transitions due to finite size
Q = interp1(t(q), Q(q), data.t, 'linear', 'extrap');
dQ = [nan; diff(Q) * 3600]; % Convert from amp-hours to amps
dQ(1) = median(dQ, 'omitnan');
data.Q = Q;
data.dQ = dQ;
data.voltage = medfilt1(data.m_bms_main_battery_voltage, 3);
data.energy1 = data.voltage .* dQ;
data.power1 = data.energy1 ./ data.dt;

% Get power estimate just using BMS
data.current = (data.m_bms_aft_current + data.m_bms_ebay_current + data.m_bms_pitch_current);
data.power = data.current .* data.m_bms_main_battery_voltage;
data.power2 = medfilt1(data.current, 3) .* data.voltage;
data.energy = data.power .* data.dt;
data.energy2 = data.power2 .* data.dt;

data.airPump = calcAirPumpOn(data); % When the air pump is actually on
data.surface = calcSurface(data); % This is air pump extended to let the bladder deflat and GPS fixes

[data.diving, data.climbing] = calcDiveClimb(data); % Get when the glider is diving/climbing
fprintf('Took %.2f seconds to load %s\n', toc(tv), deployment)
end % loadDeployment

function z = calcAirPumpOn(a)
% Get air pump on and air pump off
x = a.m_vacuum_air_pump_on; % What the reference pressure is when pump is on
x(x>0) = 1; % Pump is actually on
x(x<0) = -1; % Pump is actually off
x(isnan(x)) = 0; % No information on pump state
q = x ~= 0; % Pump state known
y = x(q);
indices = find(y(1:end-1) == y(2:end)); % Identical consective states
y(indices+1) = 0; % Drop duplicate state points
x(q) = y;
z = cumsum(x) > 0; % Identify when the air pump is actually on
end % calcAirPumpOn

function q = calcSurface(a)
% Logical from when the air pump turns on until after the last valid GPS fix
q = a.airPump ...
    | (a.c_iridium_on > 0) ...
    | (a.m_appear_to_be_at_surface > 0) ...
    ; % Airpump, Iridium, or GPS are on
end % calcSurface

function [qDive, qClimb] = calcDiveClimb(a)
% Return when the glider is diving/climbing and not on the surface
[~, indices] = findpeaks(a.m_depth, ...
    'MinPeakHeight', 20, ...
    'MinPeakDistance', 75); % ~300 seconds apart

fprintf('Found %d bottom turns', numel(indices));

depth = a.m_depth;
qDive = false(size(depth));
qClimb = false(size(depth));

if isempty(indices) % Nothing found
    return; 
end % if

[~,ix] = min(depth(1:indices(1)));
qDive(ix:indices(1)) = true; % First dive before first index

for i = 2:numel(indices)
    ii = indices(i-1):indices(i);
    [~, ix] = min(depth(ii));
    ix = ix + ii(1) - 1;
    qClimb(ii(1):ix) = true;
    qDive(ix:ii(end)) = true;
end % for

[~, ix] = min(depth(indices(end):end));
ix = ix + indices(end) - 1;
qClimb(indices(end):ix) = true;

qDive  = qDive & ~a.surface;
qClimb = qClimb & ~a.surface & ~qDive;
end % calcDiveClimb


