function [ M_n, m_n, sigma_n] = BIRCH(Kn,d,K_birch_n,br_birch)
 
callString =['bash prepInput.sh',' ',num2str(Kn) , ' ',num2str(Kn*K_birch_n),' ',num2str(d),' ',num2str(br_birch)];

changedir = 'cd ../BIRCH ;';
%Changes done to maintain directory structure
callString = [changedir ,callString];
  
system(callString);
  
findClusterCentersCmnd = './birch sample.para sample.scheme sample.proj B_k';

%Changes done to maintain directory structure
findClusterCentersCmnd = [changedir , findClusterCentersCmnd];

system(findClusterCentersCmnd);
negativeEstimates = dlmread('../BIRCH/sample.para+sample.scheme+sample.proj+B_k-0-cluster');

M_n = negativeEstimates(:,2:d+1); % each row of M_n is a center
m_n = negativeEstimates(:,1);
sigma_n = negativeEstimates(:,d+2);

% system('cd tsvd;');
