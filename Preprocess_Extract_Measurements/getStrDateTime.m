function [datenumberString] = getStrDateTime(filelist,FormatIn)
%Get date and time from a file header- for example .wav file or .mat file
%
% INPUTS:
% - filelist = a structure containing file information (obtained by calling
% filelist = dir(fullfile(folder,'*.mat')); where folder is specified folder
% path
% - FormatIn = string specifying the format of the datetime string - e.g. 
% FormatIn= 'yyyyMMdd_HHmmss_SSS';
%
% OUTPUT:
% -datenumberString - a cell array, each cell containing a string of
% specified date time information.
%
%
% Pina Gruden, August 2021

M=size(filelist,1);
N=numel(FormatIn);

datenumberString =cell(M,1);

for k=1:M
datenumberString{k}=filelist(k).name(end-4-N+1:end-4);
end

end

