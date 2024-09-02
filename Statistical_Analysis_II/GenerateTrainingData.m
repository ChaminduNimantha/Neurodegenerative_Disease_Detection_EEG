%% 
load Coherence

% Generate training and testing data
TrP = 0.7; % fraction of data to be chosen for Training
for i = 1:100
    % Healty Controls
    Tr = sort(randperm(size(Coherence.HC,3),round(size(Coherence.HC,3)*TrP)));
    Ts = findAnotinB(1:size(Coherence.HC,3),Tr);
    
    Training(i).HC = Coherence.HC(:,:,Tr);
    Training(i).TrHC = Tr;
    Training(i).TsHC = Ts;
    Training(i).Connections = Coherence.Connections;
    Training(i).bands = Coherence.bands;
    Training(i).Conditions = Coherence.Conditions;
    Training(i).mean(:,:,1) = mean(Training(i).HC,3);
    Training(i).ste(:,:,1) = std(Training(i).HC,0,3)/sqrt(size(Training(i).HC,3));
     
    Testing(i).HC = Coherence.HC(:,:,Ts);
    Testing(i).TsHC = Ts;
    Testing(i).TrHC = Tr;
    Testing(i).Connections = Coherence.Connections;
    Testing(i).bands = Coherence.bands;
    Testing(i).Conditions = Coherence.Conditions;
    Testing(i).mean(:,:,1) = mean(Testing(i).HC,3);
    Testing(i).ste(:,:,1) = std(Testing(i).HC,0,3)/sqrt(size(Testing(i).HC,3));

    % Frontotemporal Dementia Subjects
    Tr = sort(randperm(size(Coherence.FTD,3),round(size(Coherence.FTD,3)*TrP)));
    Ts = findAnotinB(1:size(Coherence.FTD,3),Tr);
    
    Training(i).FTD = Coherence.FTD(:,:,Tr);
    Training(i).TrFTD = Tr;
    Training(i).TsFTD = Ts;`ning(i).FTD,0,3)/sqrt(size(Training(i).FTD,3));
     
    Testing(i).FTD = Coherence.FTD(:,:,Ts);
    Testing(i).TsFTD = Ts;
    Testing(i).TrFTD = Tr;
    Testing(i).mean(:,:,2) = mean(Testing(i).FTD,3);
    Testing(i).ste(:,:,2) = std(Testing(i).FTD,0,3)/sqrt(size(Testing(i).FTD,3));

    
    % Alzheimer's Subjects  
    Tr = sort(randperm(size(Coherence.AZ,3),round(size(Coherence.AZ,3)*TrP)));
    Ts = findAnotinB(1:size(Coherence.AZ,3),Tr);
    
    Training(i).AZ = Coherence.AZ(:,:,Tr);
    Training(i).TrAz = Tr;
    Training(i).TsAz = Ts;
    Training(i).mean(:,:,3) = mean(Training(i).AZ,3);
    Training(i).ste(:,:,3) = std(Training(i).AZ,0,3)/sqrt(size(Training(i).AZ,3));
     
    Testing(i).AZ = Coherence.AZ(:,:,Ts);
    Testing(i).TsAz = Ts;
    Testing(i).TrAz = Tr;
    Testing(i).mean(:,:,3) = mean(Testing(i).AZ,3);
    Testing(i).ste(:,:,3) = std(Testing(i).AZ,0,3)/sqrt(size(Testing(i).AZ,3));
    
end
% save Training Training
% save Testing Testing