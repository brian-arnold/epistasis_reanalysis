addpath('Utils', 'Calculations');

% setup the directory structure
pathFolder=pwd;
dataPath=[pathFolder,'/Data/input_log/']; %changed for corrected analysis
csvFiles=dir([dataPath,'/*.csv']);  % data folder path % Macbook

checkInconclusives=1; 
CheckSingleCon=0;
rescaleMethod=3;
calculationMethod='Ursell'; 
cutoffsupp=1.3;
CI=95; %Bootstrapping confidence interval parameter (alpha)
TEST = {};

% create folder for the output files
datestr = date;
analysis_folder=createFolder(pathFolder,['/Data/output_',datestr]);
figures_folder=createFolder(analysis_folder,'/figures_log');%changed for corrected analysis
publishData_folder=createFolder(analysis_folder,'/publish data_log');%changed for corrected analysis
results_folder=createFolder(analysis_folder,'/data analysis and results summary_log');%changed for corrected analysis

% define labels for the output file
expLabels={'combinations','X1','X2','X3','X4','X5',...
                         'X1X2','X1X3','X1X4','X1X5','X2X3','X2X4','X2X5','X3X4','X3X5','X4X5',...
                         'X1X2X3','X1X2X4','X1X2X5','X1X3X4','X1X3X5','X1X4X5 ','X2X3X4','X2X3X5','X2X4X5','X3X4X5',...
                         'X1X2X3X4','X1X2X3X5','X1X2X4X5','X1X3X4X5','X2X3X4X5',...
                         'X1X2X3X4X5'};


metricLabels={'NET_12','NET_13','NET_14','NET_15','NET_23','NET_24','NET_25','NET_34','NET_35','NET_45',...
             'NET_123','NET_124','NET_125','NET_134','NET_135','NET_145 ','NET_234','NET_235','NET_245','NET_345',...
             'EM_123','EM_124','EM_125','EM_134','EM_135','EM_145 ','EM_234','EM_235','EM_245','EM_345',...
             'NET_1234','NET_1235','NET_1245','NET_1345','NET_2345',...
             'EM_1234','EM_1235','EM_1245','EM_1345','EM_2345',...
             'NET_12345',...
             'EM_12345'};



