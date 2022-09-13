function [folder, folder2save2] = Specify_Paths
%CONFIGURE PATHS to recordings & folder to save measurements to:

% Path to Raw data:
s=what('./Test_example/Data/wav/'); % This is where your .wav files are located
folder = [s.path,'/'];

% Path to folders where things will be stored:
% Raw cross-correlogram:
s=what('./Test_example/Data/Raw_CrossCorrelogram/');
folder2save2.rawcrosscorr= [s.path,'/'];

% % Normalized cross-correlogram and measurements:
% s=what('./Test_example/Data/Measurements/');
% folder2save2.measurements= [s.path,'/'];

% Final tracking results
s=what('./Test_example/Results/'); % This is where you want your final results to be saved to
folder2save2.finalresults= [s.path,'/'];


end
