%% Set the variable range

valstrt = 0.50;
valend  = 1.30;
absd    = 0.05;

val = (valstrt:absd:valend);
N = length(val);
plot(val,'.-')

% Some other settings
outdir = './test';
startchem1dcmd = 'chem1d';

%% Start the loop

% Create output directory
mkdir(outdir);
% Save the original settings file
copyfile('settings.csf','settings.csf.bck');

% Show nice wait bar
hwb = waitbar(0,'Starting chem1d series','Name','C1D series');

for n = 1:N

    % Restore settings file
    copyfile('settings.csf.bck','settings.csf');
    
    % Adjust settings file
    fid = fopen('settings.csf','a');
    fprintf(fid,'[BOUNDARY_EQUIVALENCERATIO]\n%e %e\n',val(n),val(n));

    % This script can also be used to change other parameters.
    % Uncomment one of the lines below to change temperature or pressure.
    % fprintf(fid,'[BOUNDARY_INLETTEMPERATURE]\n%e %e\n',val(n),val(n));
    % fprintf(fid,'[BOUNDARY_INLETPRESSURE]\n%e\n',val(n));
    
    fprintf(fid,'[PREPROCESSING_STARTSOLUTIONFILE]\nyiend.dat\n');
    fclose(fid);

    % Run chem1d
    waitbar(n/N,hwb,sprintf('Running %e',val(n)));
    fprintf('\nSTART RUN %i WITH VALUE = %e\n\n',n,val(n));
    system(startchem1dcmd);
    
    % Copy output file
    copyfile('yiend.dat',sprintf('%s/yi_%5.3f.dat',outdir,val(n)));
    copyfile('siend.dat',sprintf('%s/si_%5.3f.dat',outdir,val(n)));
   
end

close(hwb);
fprintf('\nFINISHED CHEM1D SERIES\n\n');

% Restore original settings file
movefile('settings.csf.bck','settings.csf');
