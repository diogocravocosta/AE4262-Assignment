%% Chem1D series for NOx comparison
% Runs phi = 0.5, 0.8, 1.0 for CH4 (or H2, change settings.csf first)
% Stores sL, Tb, NO, CO for each phi for later plotting

val = [0.5, 0.8, 1.0];
N = length(val);

startchem1dcmd = 'chem1d';
outdir = './NOx_results';
mkdir(outdir);

% Pre-allocate results storage
sL = zeros(1,N);
Tb = zeros(1,N);
NO = zeros(1,N);
CO = zeros(1,N);

% Save original settings
copyfile('settings.csf','settings.csf.bck');

hwb = waitbar(0,'Starting chem1d series','Name','NOx series');

for n = 1:N

    % Restore settings file
    copyfile('settings.csf.bck','settings.csf');
    
    % Append new phi to settings
    fid = fopen('settings.csf','a');
    fprintf(fid,'[BOUNDARY_EQUIVALENCERATIO]\n%e %e\n',val(n),val(n));
    fprintf(fid,'[PREPROCESSING_STARTSOLUTIONFILE]\nyiend.dat\n');
    fclose(fid);

    % Run chem1d
    waitbar(n/N,hwb,sprintf('Running phi = %.2f',val(n)));
    fprintf('\nSTART RUN %i: phi = %.2f\n\n',n,val(n));
    system(startchem1dcmd);
    
    % Read output
    [y,~,a] = readchem1d('yiend.dat');
    
    % Extract variables
    iTemp     = find(strcmpi('temp',a));
    iDensity  = find(strcmpi('density',a));
    iMassFlow = find(strcmpi('massflow',a));
    iCO       = find(strcmpi('CO',a));
    iNO       = find(strcmpi('NO',a));
    
    % Store results (last grid point = burnt gas)
    Tb(n) = y(end, iTemp);
    sL(n) = y(1, iMassFlow) / y(1, iDensity);  % cm/s
    NO(n) = y(end, iNO);
    CO(n) = y(end, iCO);
    
    % Save output file
    copyfile('yiend.dat', sprintf('%s/yi_%.2f.dat', outdir, val(n)));

end

close(hwb);
fprintf('\nFINISHED\n');

% Restore original settings
movefile('settings.csf.bck','settings.csf');

%% Display results
fprintf('\nphi\t\tsL(cm/s)\tTb(K)\t\tNO\t\tCO\n');
for n = 1:N
    fprintf('%.2f\t\t%.2f\t\t%.1f\t\t%.3e\t%.3e\n', val(n),sL(n),Tb(n),NO(n),CO(n));
end

%% Plot NO and CO vs phi
figure(1);
plot(val, NO, 'ro-', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('\phi [-]');
ylabel('NO mass fraction [-]');
title('NO emissions vs equivalence ratio (1D flame)');
grid on;

figure(2);
plot(val, CO, 'bo-', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('\phi [-]');
ylabel('CO mass fraction [-]');
title('CO emissions vs equivalence ratio (1D flame)');
grid on;

figure(3);
plot(val, sL, 'ks-', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('\phi [-]');
ylabel('Burning velocity s_L [cm/s]');
title('Laminar burning velocity vs equivalence ratio');
grid on;

figure(4);
plot(val, Tb, 'ms-', 'LineWidth', 1.5, 'MarkerSize', 8);
xlabel('\phi [-]');
ylabel('T_b [K]');
title('Flame temperature vs equivalence ratio');
grid on;