load Training

Bands = {'1 Delta Band:1Hz to 4Hz' '2 Theta Band: 4 to 8Hz' '3 Alpha Band: 8Hz to 13Hz' '4 Beta Band: 13Hz to 30Hz' '5 Gamma1 Band: 30Hz to 45 Hz'}; % The frequency bands used
num_bins = 10; % Number of bins for the Coherence histogram
% alpha = 0.005; % Significance level (p < alpha) for the ttest
% OverlapTh = 0.6; % Overlap threshod (at least 100*OverlapTH% of the subjects have connections)
% CthInd = 7; % Coherence threshold ((CthInd-1)/num_bins) (eg. 4 indicates Coherence above 0.3)
OLth = 0.1:0.1:1; % Different overalp thresholds

%% Scouts and Connections
n = 68; % Number of scouts
Nc = 2346;% brain scout interconnections

matrix = zeros(n);
count = 1;
for i = 1:n
    for j = 1:i
        matrix(j, i) = count;
        count = count + 1;
    end
end

region_mapping = {
    'bankssts L',...
    'bankssts R',...
    'caudalanteriorcingulate L',...
    'caudalanteriorcingulate R',...
    'caudalmiddlefrontal L',...
    'caudalmiddlefrontal R',...
    'cuneus L',...
    'cuneus R',...
    'entorhinal L',...
    'entorhinal R',...
    'frontalpole L',...
    'frontalpole R',...
    'fusiform L',...
    'fusiform R',...
    'inferiorparietal L',...
    'inferiorparietal R',...
    'inferiortemporal L',...
    'inferiortemporal R',...
    'insula L',...
    'insula R',...
    'isthmuscingulate L',...
    'isthmuscingulate R',...
    'lateraloccipital L',...
    'lateraloccipital R',...
    'lateralorbitofrontal L',...
    'lateralorbitofrontal R',...
    'lingual L',...
    'lingual R',...
    'medialorbitofrontal L',...
    'medialorbitofrontal R',...
    'middletemporal L',...
    'middletemporal R',...
    'paracentral L',...
    'paracentral R',...
    'parahippocampal L',...
    'parahippocampal R',...
    'parsopercularis L',...
    'parsopercularis R',...
    'parsorbitalis L',...
    'parsorbitalis R',...
    'parstriangularis L',...
    'parstriangularis R',...
    'pericalcarine L',...
    'pericalcarine R',...
    'postcentral L',...
    'postcentral R',...
    'posteriorcingulate L',...
    'posteriorcingulate R',...
    'precentral L',...
    'precentral R',...
    'precuneus L',...
    'precuneus R',...
    'rostralanteriorcingulate L',...
    'rostralanteriorcingulate R',...
    'rostralmiddlefrontal L',...
    'rostralmiddlefrontal R',...
    'superiorfrontal L',...
    'superiorfrontal R',...
    'superiorparietal L',...
    'superiorparietal R',...
    'superiortemporal L',...
    'superiortemporal R',...
    'supramarginal L',...
    'supramarginal R',...
    'temporalpole L',...
    'temporalpole R',...
    'transversetemporal L',...
    'transversetemporal R'...
    };

region_map = containers.Map(region_mapping, 1:length(region_mapping));

%% Define Connections
coordinates = zeros(Nc, 2);
for m = 1:Nc
    [row, col] = find(matrix == m);
    coordinates(m, :) = [row, col];
end
selfConnections = find((coordinates(:,1)-coordinates(:,2))==0);
nonselfConnections = findAnotinB(1:Nc,selfConnections);
region_names = cell(size(coordinates, 1), size(coordinates, 2));

for i = 1:size(coordinates, 1)
    for j = 1:size(coordinates, 2)
        region_names{i, j} = region_mapping{coordinates(i, j)};
    end
end
Connections = {Nc}; % The list of all the possible connections between brain regions
for i = 1:size(coordinates, 1)
    Connections{i} = sprintf([ num2str(i) '-(' num2str(coordinates(i, 1)) ') ' region_names{i, 1} ' - (' num2str(coordinates(i, 2)) ') ' region_names{i, 2}]);
end
%%
warning('off','all')
tic
for FF =1:length(Training)
    display(sprintf('%d out of %d',FF,length(Training)))
    Coherence = Training(FF); 
    
    % Connections with a coherence value above a given threshold is choosen as
    % connected with the connection strenth of coherence
    for i = 1:num_bins
        Cthreshold = (i-1)/num_bins;
        Connected.HC (:,:,:,i)= Coherence.HC>Cthreshold;
        Connected.FTD (:,:,:,i)= Coherence.FTD>Cthreshold;
        Connected.AZ(:,:,:,i) = Coherence.AZ>Cthreshold;
    end
    % The fraction of poppulation showing the connections
    FracConnected.AZ = squeeze(sum(Connected.AZ,3)/size(Connected.AZ,3));
    FracConnected.HC = squeeze(sum(Connected.HC,3)/size(Connected.HC,3));
    FracConnected.FTD = squeeze(sum(Connected.FTD,3)/size(Connected.FTD,3));
        
    % t-test
    KKK = 1:Nc;
    for band = 1:5
        for OLthInd = 1:length(OLth)
            for CthInd = 1:10
                %%
                X = (squeeze(FracConnected.AZ(KKK,band,:)));
                Y = (squeeze(FracConnected.HC(KKK,band,:)));
                Z = (squeeze(FracConnected.FTD(KKK,band,:)));
                
                OverlapTh = OLth(OLthInd);
                KK1 = unique([find(Y(:,CthInd)>OverlapTh) ;find(X(:,CthInd)>OverlapTh)]);
                KK1 =KK1(findAnotinB(KK1,selfConnections)); % remove self connections
                KK2 = unique([find(Z(:,CthInd)>OverlapTh) ;find(X(:,CthInd)>OverlapTh)]);
                KK2 =KK2(findAnotinB(KK2,selfConnections));
                KK3 = unique([find(Z(:,CthInd)>OverlapTh) ;find(Y(:,CthInd)>OverlapTh)]);
                KK3 = KK3(findAnotinB(KK3,selfConnections));
                KK= unique([KK1;KK2;KK3]);
                KALL(band,CthInd,OLthInd).KK = KK;
                KALL(band,CthInd,OLthInd).number = length(KK);
            end
        end
    end

    Kall = [];
    CthInd = 4;
    OLthInd = 4;
    alpha = 0.008;

    for band = 1:5
        [H1, P1, ~, st1]  = ttest2(squeeze(Coherence.FTD(KKK,band,:))',squeeze(Coherence.HC(KKK,band,:))','alpha',alpha);
        [H2, P2, ~, st2]  = ttest2(squeeze(Coherence.AZ(KKK,band,:))',squeeze(Coherence.HC(KKK,band,:))','alpha',alpha);
        [H3, P3, ~, st3]  = ttest2(squeeze(Coherence.AZ(KKK,band,:))',squeeze(Coherence.FTD(KKK,band,:))','alpha',alpha);
        
        % Select only the connections that show significant difference (H=1),
        % that satisfy the overlap and thresholding criteria (connections in the KK list)
        KK = KALL(band,CthInd,OLthInd).KK;
        
        K1 = KK(findAinB(KKK(H1==1),KK)); % FTD vs HC [ 2 vs 1]
        K2 = KK(findAinB(KKK(H2==1),KK)); % AZ vs HC [3 vs 1]
        K3 = KK(findAinB(KKK(H3==1),KK)); % AZ vs FTD [3 vs 2]
        K = unique([K1;K2; K3]);
        k.k1= findAinB(K1,K);
        k.k2 = findAinB(K2,K);
        k.k3 = findAinB(K3,K);
        k.k = unique([k.k1 k.k2 k.k3]);
        k.p1 = P1(K(k.k1)); % FTD vs HC
        k.p2 = P2(K(k.k2)); % AZ vs HC
        k.p3 = P3(K(k.k3)); % AZ vs FTD
        k.K = K;
        k.st1 = st1;
        k.st2 = st2;
        k.st3 = st3;
        k.alpha = alpha;
        k.Bonf_fac = 1;% length(KK);
        allk(band) = k;              
        Kall = [Kall;allk(band).K];
    end
    Kall = unique(Kall);
    Kall_training(FF).Kall = Kall;
    Kall_training(FF).allk = allk;
    toc
end
%%
KALL_tr=[];
for i=1:100;
    KALL_tr(i,1:length(Kall_training(i).Kall)) =Kall_training(i).Kall;
end
tmp3 = KALL_tr; tmp3(KALL_tr==0)=nan;

Np = 50; % Show only the connections that are common to Np number of trials
figure
H1 = histogram(tmp3,'BinWidth',1);
Edges = H1.BinEdges(H1.Values>0);
Values = H1.Values(H1.Values>0);

% bar(1:length(Values(Values>Np)),Values(Values>Np))
% set(gca,'Xtick', 1:length(Values(Values>Np)),'XtickLabels',Connections(Edges(Values>Np)),'XTickLabelRotation',90)

%
KALL_Training=[];
for i=1:100;
    for band = 1:5;
        KALL_Training(1:length(Kall_training(i).allk(band).K),band,i)=Kall_training(i).allk(band).K;
    end;
end
Kallvalues = unique(KALL_Training(~isnan(KALL_Training)));
Kallvalues = Kallvalues(Kallvalues~=0);
Kmat = zeros(length(Kallvalues),length(Bands),length(Kall_training));
for i=1:100;
    for band = 1:5;
        Kmat((findAinB(KALL_Training(:,band,i),Kallvalues)),band,i)=1;
    end;
end

Ksum = sum(Kmat,3);
Ksum(Ksum<Np)=0;


% Initialize the coherence_mask matrix with zeros
coherence_mask = zeros(2346, 5, 1);

% Iterate over each frequency band

for band = 1:5

    % Extract significant connections for the current band
    significant_connections = Ksum(:,band) > Np;

    significant_indices = Kallvalues(significant_connections);
    
    % Mark the significant connections with 1, considering only the Kallvalues indexes
    coherence_mask(significant_indices, band, 1) = 1;
end

% Save the matrix and additional information in 'coherence_mask.mat'
save('coherence_mask_4_3_4.mat', 'coherence_mask', 'CthInd', 'OLthInd', 'alpha', 'Np');


Ipick = find(sum(Ksum,2)>60);
Ksum2 = Ksum(Ipick,:);
% figure;
% Ksum2(isnan(Ksum2)) = 0;
imagesc(Ksum2)
set(gca,'Xtick', 1:size(Bands,2),'XtickLabels',Bands,'XTickLabelRotation',90,'Ytick', 1:size(Ksum2,1),'YtickLabels',Connections(Edges(Ipick)),'YTickLabelRotation',0)
colorbar
title(sprintf('Coherance Connections showing significant differences between conditions in at least %d out of 100 trials',Np))

