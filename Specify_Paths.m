function [folder, folder2save2] = Specify_Paths
%CONFIGURE PATHS to recordings & folder to save measurements to:

% Path to Raw data:
s=what('./Test_example/Data/'); % This is where your .wav files are located
folder = [s.path,'/'];

% Path to folder where things will be stored:
s=what('./Test_example/Results/'); % This is where you want your results to be saved to
folder2save2= [s.path,'/'];


end
