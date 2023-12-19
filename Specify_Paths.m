function [folder, folder2save2] = Specify_Paths
% Specify_Paths.m specifies paths for "TDOA_tracking_master" package
%
% Use Specify_Paths.m to CONFIGURE PATHS to recordings & folders to save 
% data/results to. 
%
% OUTPUT:
% - folder - a string specifying path to the folder where audio data is
%           located
% - folder2save2 - a structure specifying paths to where data is stored to.
%                   Has two fields:
%                   ~'rawcrosscorr' (path to where cross-correlogram 
%                   will be stored);
%                   ~'finalresults' (path to where final results will
%                    be stored)
%
%Pina Gruden, 2022, UH Manoa

% Path to Raw data:
myfolder = './Test_example/Data/wav/'; % This is where your .wav files are located
if not(isfolder(myfolder)) % If the folder doesnt exist - throw error since you need data to track.
error(['The folder ', myfolder, ' does not exists. There is no data ' ...
    'for tracking! Check your path and try again!'])
end
s=what(myfolder);
folder = [s.path,'/'];

% Path to folders where things will be stored:
myfolder='./Test_example/Data/Raw_CrossCorrelogram/'; % Raw cross-correlogram:
if not(isfolder(myfolder)) % If the folder doesnt exist
    mkdir(myfolder) %make a folder    
    disp(['WARNING: The specified folder ', myfolder, ' does not exists. Created a new folder.'])
end
s=what(myfolder); 
folder2save2.rawcrosscorr= [s.path,'/'];

% % Normalized cross-correlogram and measurements:
% s=what('./Test_example/Data/Measurements/');
% folder2save2.measurements= [s.path,'/'];

% Final tracking results
myfolder='./Test_example/Results/'; % This is where you want your final results to be saved to
if not(isfolder(myfolder)) % If the folder doesnt exist
    mkdir(myfolder) %make a folder    
    disp(['WARNING: The specified folder ', myfolder, ' does not exists. Created a new folder.'])
end
s=what(myfolder); 
folder2save2.finalresults= [s.path,'/'];


end
