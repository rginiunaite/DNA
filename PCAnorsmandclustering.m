%% PCA and clustering, norms (distances between points)
clear all
% 
% % distances
% %% create dominant eigenvectors
% alldistances = load ('AllDistances.mat');
% 
% allshapevalues = alldistances.distances;
% 
% allshapevalues(:,all(allshapevalues == 0))=[];
% 
% covmatrix = cov(allshapevalues); % covariance matrix of all shape vectors
% 
% [eigvec, eigval] = eig(covmatrix);
% 
% [d,ind] = sort(-diag(eigval)); % sort eigenvalues, descending order
% 
% 
% v1= eigvec(:,ind(1));
% v2 =eigvec(:,ind(2));
% 
% X = [v1,v2];
% 
% save('DominEigenDistances.mat','X')
% 
% %% end of create dominant eigenvectors
% 
% 


%% after I already saved dominant eigenvectors
DominE = load('DominEigenDistances.mat'); % DominEigen.mat

X = DominE.X;

paramset = load('ParameterSets/cgDNA+ps1_posdef.mat');
data = load('Nucleosomes.mat');


%% select specific groups/organisms
vec = [data.Seq.group];

str = 'Drosophila';
idx = strcmp(str,{data.Seq.group});%,'Mouse',{data.Seq.group}) % Yeast, Drosophila, Virus, Human
% idx2 = strcmp('Mouse',{data.Seq.group})
%idxboth = or(idx,idx2);
datagroup = data.Seq(idx); %idxboth



%% for all
seqnum = length(datagroup);
lenseq =73;

%Dim2point = zeros (seqnum,2);
k=1; %count;
for j = 1:seqnum

    seq = datagroup(j).S; % from nucleosomes data
   
    if length(seq)==400
        
        [shapes, stiff] = constructSeqParms(seq, paramset);

        abs_coord = frames(shapes); % relative to absolute coordinates

        for i=1:400-lenseq 
            distances(i) = norm(abs_coord(i).rc -abs_coord(i+lenseq).rc);

        end
        
        Dim2point(k,:)  = X'*distances';
       
        k=k+1;
   end
    
    
    
end

hold on
scatter(Dim2point(:,1),Dim2point(:,2), 'g','filled')

%save('Virus.mat','Dim2point')



