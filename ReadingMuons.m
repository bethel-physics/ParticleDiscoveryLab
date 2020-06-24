% This program reads in the files from MIT stuff
clear all

% HOW MANY GLOBAL MUONS PER EVENT?
globalMuons=4;

% reading the file into MATLAB as a cell
files={'000D143E-9535-E311-B88B-002618943934.mod'
       '0020AF81-A835-E311-97DC-00261894398C.mod'
       '0062BD79-9235-E311-9320-003048FFCC2C.mod'
       '008B0182-5E35-E311-B039-003048FFD71E.mod'
       '009A5FA1-0D36-E311-B603-00261894388A.mod'
       '00ADCF5D-7435-E311-9347-00261894397A.mod'
       '00C4AE70-4335-E311-9743-0026189438D9.mod'
       '00EF5A25-6235-E311-8D64-002618B27F8A.mod'
       '0200DC35-9535-E311-8AE6-003048678B34.mod'
       '02495C0E-3A35-E311-A75D-003048678BE6.mod'
       '02A4FDA7-7135-E311-AEAB-003048FFCBFC.mod'
       '02A5EA3E-6F35-E311-9232-003048678F8A.mod'
       '02CC8DC7-9435-E311-8DDD-0026189438F3.mod'
       '02DC4EAD-9235-E311-B0C3-002618943985.mod'
       '02F46E80-7335-E311-B046-003048FFD732.mod'
       '041E647E-5435-E311-8A79-002618943856.mod'
       '04302453-1136-E311-95F0-00261894393E.mod'
       '0467D790-6C35-E311-B1DA-003048678FEA.mod'
       '04695FC2-5B35-E311-AADF-003048FF9AA6.mod'
       '046B867E-A835-E311-B2E8-003048FFD728.mod'
       '0490CFDC-4A35-E311-85B6-002618943920.mod'
       '049A90E4-6835-E311-A384-003048678BEA.mod'
       '04F1C9E5-A935-E311-9F92-002618943945.mod'
       '061185A7-6A35-E311-A8D1-003048678FEA.mod'
       '062E45B5-1036-E311-9CA0-0025905964BA.mod'
       '06872350-1136-E311-B9B2-002618943972.mod'
       '068A7559-9535-E311-998F-003048678B18.mod'
       '0698FAA3-A835-E311-9077-003048FFD730.mod'
       '06C90E65-1036-E311-9B29-003048FFCB74.mod'
       '06D8874F-1136-E311-99F6-003048FFCC0A.mod'
       '06F9095B-6335-E311-901D-003048D15DDA.mod'
       '081F9C81-7335-E311-9E3E-00304867904E.mod'
       '082FEFBB-9435-E311-94C2-00261894398A.mod'
       '0882F73E-5735-E311-95DE-00261894382A.mod'
       '08A7C079-9235-E311-80BB-002590593878.mod'
       '08D7F7AF-0D36-E311-9FFE-002590593902.mod'
       '0A1F51A8-4F35-E311-B159-003048679164.mod'
       '0A3AF0BB-0D36-E311-9223-0025905964A6.mod'
       '0A3B9FF6-3735-E311-8F7D-003048FFCBA8.mod'
       '0A5F6975-7035-E311-B873-0026189438EB.mod'
       '0A859F13-7135-E311-8B28-0030486792AC.mod'
       '0A8B6FAB-4935-E311-A1A2-002590593920.mod'
       '0A8C3DF6-5235-E311-B3DF-00261894398C.mod'
       '0AE1FB8D-4D35-E311-85C6-003048FFCB84.mod'
       '0AE85014-6635-E311-987C-0026189438AF.mod'
       '0AF86F54-6A35-E311-BECC-002618943918.mod'
       '0C11E993-6435-E311-ABF1-00261894389C.mod'
       '0CE5A9CE-9435-E311-A2A2-00304867C04E.mod'
       '0CFEA493-6135-E311-B1A6-003048FF9AC6.mod'
       '0E213E73-4B35-E311-AAAD-002618FDA211.mod'
       '0E2DD265-9535-E311-AF68-003048FF9AA6.mod'
       '0E74BA28-6235-E311-95F5-0030486792A8.mod'
       '0EB2D03A-9235-E311-B508-002618943930.mod'
       '100B9179-6535-E311-80B9-002618FDA210.mod'
       '10314681-AA35-E311-A2F8-002590596490.mod'
       '10627D58-5C35-E311-9BDD-003048D42D92.mod'};
lengths=[686274
         532263
         574221
         595399
         466060
         553492
         492543
         613393
         524717
         433781
         400820
         443210
         510226
         558630
         512877
         456001
         602985
         398975
         533236
         437063
         86188
         494129
         427713
         504092
         515707
         569446
         719445
         577657
         561459
         589605
         725237
         354758
         508090
         539712
         559267
         427771
         573591
         329923
         604388
         588445
         468906
         322132
         358300
         578955
         492258
         564462
         340422
         679820
         486640
         595049
         689566
         370733
         449940
         347587
         595661
         578189];
tooBig=0;
badFormat=0;
PX={};
PY={};
PZ={};
ENERGY={};
for k=1:length(files)
    fileID = fopen(files{k},'r');
    WholeFile = textscan(fileID,'%s',lengths(k),'Delimiter','\n');

    % find the indices where each new event begins
    Indices = 1;
    for i = 2:length(WholeFile{1})
        if strcmp(WholeFile{1}{i},'BeginEvent Version 0 CMS_2011 DoubleMu_Primary_Dataset')
            Indices = [Indices, i];
        end
    end

    % Separate each even into its own cell in Events (so Event1 =
    % Events{1}, Event2 = Events{2},...)
    NumEvents = length(Indices);
    for i = 1:NumEvents-1
        Events{i} = {WholeFile{1}{Indices(i):Indices(i+1)-1}};
    end
    Events{NumEvents}={WholeFile{1}{Indices(NumEvents-1):lengths(k)}};

    % These for loops find the indices where the momentum data is stored, the
    % indices where it starts and ends
    for j = 1:NumEvents
        EventGood=0;
        for i = 1:length(Events{j})

            if strcmp(Events{j}{i}, '#    Muon              px              py              pz          energy           isGlobal')
                if j == 1
                    MomentumBegin = i;
                else
                    MomentumBegin = [MomentumBegin, i];    
                end
                EventGood=1;
            end
            if strcmp(Events{j}{i}, 'EndEvent')
                if j == 1
                    MomentumEnd = i;
                else
                    MomentumEnd = [MomentumEnd, i];
                end
                if EventGood==0
                    MomentumBegin=[MomentumBegin,-42];
                    badFormat=badFormat+1;
                end
            end
        end
    end

    Check = 0;
    % THIS ASSUMES THE NUMBERS START AFTER INDEX 4
    % This creates cells for px, py, pz, and energy. px{1} stores the
    % x-momentum values from Event1, and so on.
    for i = 1:NumEvents
        areGlobal{i}=0;
        if MomentumBegin(i)>0
            for j = (MomentumBegin(i)+1):(MomentumEnd(i)-1)
                VectorOfValues = str2num(Events{i}{j}(6:length(Events{i}{j})));
                if Check == 0
                    if length(VectorOfValues)==5
                        px{i,k}{j - MomentumBegin(i)} = {VectorOfValues(1)};
                        py{i,k}{j - MomentumBegin(i)} = {VectorOfValues(2)};
                        pz{i,k}{j - MomentumBegin(i)} = {VectorOfValues(3)};
                        Energy{i,k}{j - MomentumBegin(i)} = {VectorOfValues(4)};
                        areGlobal{i,k}=areGlobal{i}+VectorOfValues(5);
                    else
                        px{i,k}{j-MomentumBegin(i)}={0};
                        py{i,k}{j-MomentumBegin(i)}={0};
                        pz{i,k}{j-MomentumBegin(i)}={0};
                        Energy{i,k}{j-MomentumBegin(i)}={0};
                        areGlobal{i,k}=-12345;
                    end
                else
                    px{i}{j - MomentumBegin(i)} = {px, VectorOfValues(1)}; %What in the world is this? The dimensions of this matrix would be so weird
                    py{i}{j - MomentumBegin(i)} = {py, VectorOfValues(2)};
                    pz{i}{j - MomentumBegin(i)} = {pz, VectorOfValues(3)};
                    Energy{i}{j - MomentumBegin(i)} = {Energy, VectorOfValues(4)};
                end
            end
        end
    end
    
    NumToKeep=0;
    CurrentKeptEvents=[];
    for i=1:NumEvents
        NumToKeep=NumToKeep+(areGlobal{i}==globalMuons);
        if areGlobal{i}==globalMuons
            CurrentKeptEvents=[CurrentKeptEvents,i];
        end
        tooBig=tooBig+(areGlobal{i}>globalMuons);
        badFormat=badFormat+(areGlobal{i}<0);
    end
    KeptEvents{k}=CurrentKeptEvents;

    % find which muon has the greatest momentum in each event
    CurrentGreatestMomentum(1:globalMuons,1:NumToKeep) = 0;
    CurrentMomMarkers(1:globalMuons,1:NumToKeep) = 0;
    for i = 1:NumToKeep
        for n=1:globalMuons
            for j = 1:length(px{KeptEvents{k}(i),k})
                p = sqrt(px{KeptEvents{k}(i),k}{j}{1}^2 + py{KeptEvents{k}(i),k}{j}{1}^2 + pz{KeptEvents{k}(i),k}{j}{1}^2);
                if p > CurrentGreatestMomentum(n,i) && (n==1 || p<CurrentGreatestMomentum(n-1,i))
                    CurrentGreatestMomentum(n,i) = p;
                    CurrentMomMarkers(n,i) = j;
                end
            end
        end
    end
    if k==1
        GreatestMomentum=CurrentGreatestMomentum;
        MomMarkers=CurrentMomMarkers;
        PX{1}=px{KeptEvents{k}(i),k}{j}{1};
        PY{1}=py{KeptEvents{k}(i),k}{j}{1};
        PZ{1}=py{KeptEvents{k}(i),k}{j}{1};
        ENERGY{1}=Energy{KeptEvents{k}(i),k}{j}{1};
        start=2;
    else
        GreatestMomentum=[GreatestMomentum,CurrentGreatestMomentum];
        MomMarkers=[MomMarkers,CurrentMomMarkers];
        start=1;
    end
    for i=start:NumToKeep
        PX={PX,px{KeptEvents{k}(i),k}{j}{1}};
        PY={PY,py{KeptEvents{k}(i),k}{j}{1}};
        PZ={PZ,pz{KeptEvents{k}(i),k}{j}{1}};
        ENERGY={ENERGY,Energy{KeptEvents{k}(i),k}{j}{1}};
    end

    % Calculate the invariant masses
%     for i = 1:NumToKeep
%         TotSqrEnergy=0;
%         TotSqrMom=0;
%         for j = 1:globalMuons
%             TotSqrEnergy=TotSqrEnergy+Energy{KeptEvents(i)}{int16(MomMarkers(j,i))}{1}^2;
%             TotSqrMom=TotSqrMom+GreatestMomentum(j,i)^2;
%         end
%         InvariantMasses{k}{i}=(1/((3*10^8)^2)) * ((1/((3*10^8)^2))*TotSqrEnergy - TotSqrMom);
%     end
     disp(strcat('Successfully read _',num2str(k),'_ files'))
end
