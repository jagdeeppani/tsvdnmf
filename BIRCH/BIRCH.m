% function [ ] = EstimateMeansFromBirch(k,m ,d,positiveS,negativeS,y)
function [ M_n, m_n, sigma_n] = EstimateNegativeMeansFromBirch2(Kn,d,K_birch_n,br_birch)
%COMPAREACTUALESTIMATED Summary of this function goes here
%   Detailed explanation goes here

%      load('positiveS.mat');
%      load('negitiveS.mat');

%LATER pass Kp and Kn from previous function becaous number of points may
%not be equal in all clusters.
% Kp = length(positiveS) / m;
% Kn = length(negativeS) / m;
% K = Kp + Kn;

%Preparing input files for clustering positive samples
% callString =strcat('bash -x prepInput.sh',' ',num2str(Kp));
% callString = strcat(callString,' ',num2str(Kp*2000));

%%%    TUNING PARAMETER   %%%%%
  
callString =['bash prepInput.sh',' ',num2str(Kn) , ' ',num2str(Kn*K_birch_n),' ',num2str(d),' ',num2str(br_birch)];

changedir = 'cd ../birch2 ;';
%Changes done to maintain directory structure
 callString = [changedir ,callString];
  
system(callString);
  
findClusterCentersCmnd = './birch sample.para sample.scheme sample.proj samplen_closeM.data';

%Changes done to maintain directory structure
findClusterCentersCmnd = [changedir , findClusterCentersCmnd];

system(findClusterCentersCmnd);
negativeEstimates = dlmread('../birch2/sample.para+sample.scheme+sample.proj+samplen_closeM.data-0-cluster');

M_n = negativeEstimates(:,2:d+1);
m_n = negativeEstimates(:,1);
sigma_n = negativeEstimates(:,d+2);