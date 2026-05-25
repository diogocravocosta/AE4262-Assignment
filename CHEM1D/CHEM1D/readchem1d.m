%READCHEM1D reads data from chem1d output file version 4
%
% [Y, T, C] = readchem1d(FILENAME) loads data from FILENAME into Y, T and
% C. The matrix Y contains the data with each variable in one column. The
% first column contains the spatial position. The variable names are loaded
% into the cell array C. The time is loaded into the scalar T.
%
% The variable names can be used to select the column from Y that you want.
%
% Example
%    [y,t,a] = readchem1d('yiend.dat');
%    plot(y(:,strcmpi('x(i)',a)),y(:,strcmpi('temp',a)))
%
function [y, t, a] = readchem1d(fname)
clear y t a;

fid = fopen(fname,'r');
if fid==-1
    error('File %s does not exist',fname);
end

while 1
    line = fgetl(fid);
    if ~ischar(line);   break;   end;
    if strcmp(line,'[TIME]')
        line = fgetl(fid);
        t = str2double(line);
    elseif strcmp(line,'[FILE_STRUCTURE_COLUMNS_CONTAINING]')
        line = fgetl(fid);
        b = textscan(line,'%s');
        a = b{1};
        ncol = length(a);
        y = fscanf(fid,'%e',[ncol,inf]);
        y = y';
        y(:,1) = [];
        a(1) = [];
        break
    end
end
fclose(fid);

end