function [sens] = critEEG_sensorselection()

load('critEEG_chlabel.mat')

% motor based on buttonpress - nobuttonpress, collapsed over stim and cond
% press stimlocked.
% leftmotor = {'FC5', 'FC3', 'FC1', 'C5', 'C3', 'C1', 'CP3', 'CP1'}; %'CP5', 
% rightmotor = {'FC2', 'FC4', 'FC6', 'C2', 'C4', 'C6', 'CP2', 'CP4'};
leftmotor = { 'FC1' 'C3' }; %'CP5', 
rightmotor = {'FC2' 'C4'};
motor = [leftmotor rightmotor]

leftmotorind = match_str(chlabel,ft_channelselection(leftmotor,chlabel));
rightmotorind = match_str(chlabel,ft_channelselection(rightmotor,chlabel));
motorind =  match_str(chlabel,ft_channelselection(motor,chlabel));

leftocc = {'PO7', 'PO5', 'PO3', 'PO1', 'O1', 'I1'};
rightocc = {'PO4', 'PO6', 'PO8', 'PO10', 'O2', 'I2'};

% triangle without O1 O2 but with parietal P1-6
occ = { 'POz', 'Oz', 'Pz', 'PO4', 'PO3', 'P2', 'P1' , 'P3', 'P4', 'P5', 'P6'  }; 

leftoccind = match_str(chlabel,ft_channelselection(leftocc,chlabel));
rightoccind = match_str(chlabel,ft_channelselection(rightocc,chlabel));
occind =  match_str(chlabel,ft_channelselection(occ,chlabel));

occpar = {    'PO7'    'PO8'    'P6'    'P8'    'PO3'    'P7'    'O1'    'P5'    'P4'    'O2'    'PO4'}
occparind =  match_str(chlabel,ft_channelselection(occpar,chlabel));

frontal = 'F*'; % 'Fp1'    'F7'    'F3'    'FC1'    'FC5'    'FC6'    'FC2'    'F4'    'F8'    'Fp2'    'Fz'
frontalind = match_str(chlabel,ft_channelselection(frontal,chlabel));

occlatr = { 'PO7' 'P7' 'P5'  'PO8' 'P6' 'P8' }; % used also by Arazi et al
occlatrind = match_str(chlabel,ft_channelselection(occlatr,chlabel));

posterior = {'P7','P3','Pz','PO3','O1','Oz','O2','PO4','P4','P8','P1','P2','PO7','P9','P5','I1','Iz','POz','I2','P10','PO8','P6'};
posteriorind = match_str(chlabel,ft_channelselection(posterior,chlabel));

% CPP signal
cpp = {'Pz'};
cppind = match_str(chlabel,ft_channelselection(cpp, chlabel));

% POz where gamma is strongest
POz = {'POz'};
POzind = match_str(chlabel,ft_channelselection(POz, chlabel));

channel = {'Pz', 'P6', 'PO3', 'POz', 'PO4', 'O1', 'Oz', 'O2'};
quenchind = match_str(chlabel,ft_channelselection(channel, chlabel));

% handpicked lateral positive cluster showing entropy increase
channel = {'F7', 'F8', 'FC5', 'FC6', 'C3', 'C4', 'CP5', 'CP3', 'CP4'};
entropyind = match_str(chlabel,ft_channelselection(channel, chlabel));

channel = {'T7', 'T8', 'TP7', 'CP5', 'CP6', 'TP8', 'P5', 'P3', 'P6', 'P8', 'I1', 'Iz', 'I2'}; %lib-cons raw
entropy2ind = match_str(chlabel,ft_channelselection(channel, chlabel));

sens=[];
allsens = 1:48;
sens.ind = {occind; motorind; allsens'; occparind; frontalind; occlatrind; posteriorind; cppind; POzind; quenchind; entropyind; entropy2ind};
sens.leg = {'occipital'; 'motor'; 'allsens'; 'occpar'; 'frontal'; 'occlatr'; 'posterior'; 'Pz'; 'POz'; 'occ_quench'; 'front_entr'; 'occlatr_entr'};
sens.label = chlabel;



