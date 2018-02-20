function [D, zeta1] = threshold_tsvdnmf(A, eps1, alpha1, beta1,tolerance)
% Assuming input will be (d*w), that is transposed, so A will not be copied
% Does not require count vector, works on real matrix A
% Columns: features, Rows: data points

addpath ../sparsesubaccess/sparsesubaccess
outpath  = 'output';
subpath=sprintf('/thresholding_%f_%f_%f_%d.mat',eps1, alpha1, beta1,tolerance);
opath=strcat(outpath,subpath);
rtime=tic;
[n,d] = size(A);
B = sort(A,1,'descend');
nu1 = B (ceil(eps1*n),: );

nu_idx = find(nu1==0);

for i=1:length(nu_idx)
    tmp = B(:,nu_idx(i));
    
    min_tmp = min( nonzeros(tmp) );     % We assume each feature is nonzero in atleast one data point, so min_temp will not be empty.
    nu1( nu_idx(i) ) = min_tmp;
end

zeta1=  alpha1*nu1; % elements of nu1 can be zero.
clear B;
fprintf('Time taken to find threshold is %f \n',toc(rtime));
nzA=nnz(A);
nzD = 0;
fprintf('NNZ of A is %d \n',nzA);

id_cols = zeros(nzA, 1);
id_rows = zeros(nzA, 1);
values = zeros(nzA, 1);

for i=1:d

        ridx = find(A(:,i) >= zeta1(i));
        lg=length(ridx);
        
        id_rows(nzD+1:nzD+lg) = ridx;
        id_cols(nzD+1:nzD+lg) = i*ones(lg,1);
        values(nzD+1:nzD+lg) = sqrt(zeta1(i))*ones(lg,1);
        nzD = nzD + lg;
        
end

D = sparse(id_rows(1:nzD),id_cols(1:nzD),values(1:nzD),n,d);
clear id_rows id_cols values;
fprintf('Time taken to find initial D is %f \n',toc(rtime));
initial_nnzD = nnz(D);
fprintf('NNZ of initial D is %d \n',initial_nnzD);

W = cell(1,d);
dpw = zeros(d,1);
for i=1:d
    W{i} = find(D(:,i));       % D is asuumed to be dataPoint x feature matrix
    dpw(i) = length(W{i});   % # of dataPoints per feature i
end
fprintf('Time taken to find W is %f \n',toc(rtime));
[dpw_val,dpw_idx] = sort(dpw,'ascend');
dpw_id = find(dpw_val,1);
R = dpw_idx(dpw_id:end);

l1=1;
lc=0;
Li_ln = zeros(d,1);
lR=length(R);
ub_t=beta1*n;
id_cols = zeros(initial_nnzD, 1);
id_rows = zeros(initial_nnzD, 1);
nzD =0;

fprintf('while loop started with lR = %d \n',lR);

while l1 <= lR-1
    i = R(l1);
    Li = W{i};
    Li_ln(i) = length(Li);
    l2=l1+1;
    
    ub = dpw(i)+ub_t;
    if (Li_ln(i) > 2*tolerance)
        while (l2 <= lR)
            id = R(l2);
            
            if  (dpw(id) >= ub)
                Wid = W{id};
                if sum(~ismember(Li, Wid))<=tolerance
                    
                    rest_records = setdiff(Wid,Li);
                    lg=length(rest_records);
                    
                    id_rows(nzD+1:nzD+lg) = rest_records;
                    id_cols(nzD+1:nzD+lg) = id*ones(lg,1);
                    nzD = nzD + lg;
                    
                    R(l2) = [];
                    lR=lR-1;
                    l2=l2-1; % As the size of R decreases auto increment by for loop may miss the next feature. so a decrement has been done
                    lc=lc+1;
                    fprintf('lc is %d  and l1 = %d, Li_ln(i) = %d \n',lc,l1,Li_ln(i))
                end
                
            end
            
            l2=l2+1;
        end
        
    end
    l1=l1+1;
end
fprintf('Finished at : l1 is %d Li = %d \n',l1-1,Li_ln(i));
D = setsparse(D, id_rows(1:nzD), id_cols(1:nzD), zeros(nzD,1));

final_nnzD = nnz(D);
fprintf('NNZ of final D is %d \n',final_nnzD);
fprintf('lc is %d\n',lc);
time_taken=toc(rtime);
fprintf('Time taken by threshold_tsvdnmf is %f secs \n',time_taken);
save(opath,'time_taken','W','Li_ln','initial_nnzD','final_nnzD','zeta1','D');
D=D';
end
