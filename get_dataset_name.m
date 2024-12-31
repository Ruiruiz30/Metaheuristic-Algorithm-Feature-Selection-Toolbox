function [name,feat,label] = get_dataset_name(A)
switch A
    case 1 % UCI
        load ./Dataset/arrhythmia.mat;
        name = 'arrhythmia';
    case 2 % UCI 
        load ./Dataset/ionosphere.mat;
        name = 'ionosphere';
    case 3 % Lymphoma 
        load ./Dataset/Lymphoma.mat;
        name = 'Lymphoma';
    case 4 
        load ./Dataset/CLL-SUB-111.mat;
        name = 'CLL-SUB-111';
    case 5 
        load ./Dataset/TOX-171.mat;
        name = 'TOX-171';
    case 6 
        load ./Dataset/Prostate-GE.mat;
        name = 'Prostate-GE';
    case 7 
        load ./Dataset/SMK-CAN-187.mat;
        name = 'SMK-CAN-187';
    case 8              
        load ./Dataset/nci9.mat;
        name = 'nci9';
    case 9 
        load ./Dataset/colon.mat;
        name = 'colon';
    case 10 
        load ./Dataset/ALLAML.mat;
        name = 'ALLAML';
    case 11 
        load ./Dataset/GLI-85.mat;
        name = 'GLI-85';
    case 12 
        load ./Dataset/orlraws10P.mat;
        name = 'orlraws10P';
    case 13 
        load ./Dataset/pixraw10P.mat;
        name = 'pixraw10P';
    case 14 % Yale
        load ./Dataset/yale.mat;
        name = 'Yale';
    case 15 
        load ./Dataset/warpAR10P.mat;
        name = 'warpAR10P';
    case 16 
        load ./Dataset/warpPIE10P.mat;
        name = 'warpPIE10P';
    case 17 % UCI
        load ./Dataset/gastroenterology.mat;
        name = 'gastroenterology';
    case 18 
        load ./Dataset/ArceneDataset.mat;
        name = 'ArceneDataset';
    case 19 
        load ./Dataset/CNAE9.mat;
        name = 'CNAE9';
    case 20 
        load ./Dataset/Maternal_Health_Risk.mat;
        name = 'Maternal-Health-Risk';
    case 21 %Waveform Database Generator /Physics and Chemistry   
        load ./Dataset/Waveform.mat;
        name = 'Waveform-Database';

end






end
