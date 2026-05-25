%% Load chem1d output data in matrix y
% Variable names are loaded into cell array a
[y,t,a] = readchem1d('./yiend.dat');

% Assign some pointers
iTemp = find(strcmpi('temp',a));
iDensity = find(strcmpi('density',a));
iMassFlow = find(strcmpi('massflow',a));
iCH4  = find(strcmpi('CH4',a));
iO2   = find(strcmpi('O2',a));
iCO2  = find(strcmpi('CO2',a));
iH2O  = find(strcmpi('H2O',a));
iCO   = find(strcmpi('CO',a));
iNO   = find(strcmpi('NO',a));
iOH   = find(strcmpi('OH',a));

% Add more species here, if you want

% Put spatial coordinate in array x
x = y(:,strcmpi('x(i)',a));

fprintf('Flame temperature: %e K\n', y(end,iTemp));
fprintf('Burning velocity : %e cm/s\n', y(1,iMassFlow)/y(1,iDensity));
fprintf('CO emission      : %e\n', y(end,iCO));
fprintf('NO emission      : %e\n', y(end,iNO));

%% Plot T vs x
figure(1);
plot(x, y(:,iTemp), '.-');
xlabel('x [cm]');
ylabel('T [K]');

%% Plot species mass fractions vs x
figure(2);
% Select some species
ivar = [iCH4, iO2, iCO2, iH2O, iCO, iOH];
plot(x, y(:,ivar), '.-');
xlabel('x [cm]');
ylabel('Mass fraction [-]');
% Add a legend
legend(a{ivar});

% Uncomment next line for log scale
% set(gca, 'Yscale', 'log', 'Ylim', [1e-8 1]);

%% Plot NO mass fraction vs x
figure(3);
plot(x, y(:,iNO), '.-');
xlabel('x [cm]');
ylabel('NO mass fraction [-]');

%% Plot density vs x
figure(4);
plot(x, y(:,iDensity), '.-');
xlabel('x [cm]');
ylabel('Density [kg/m3]');
