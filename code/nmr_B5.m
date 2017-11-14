close all
clear

filePath = char("C:\Users\Benjamin's PC\Dropbox\School Work\PHYS460A\Experiment 20\460_E20_NMR\rawdata\");

sampleDetails = [57,64,50, 100;...
                 66,71,30, 100;...
                 72,82,20, 75];
colorMap = {'YlOrRd','YlGnBu','YlGn'};
sampleName = {'Doped Water', 'Rubber', 'Ethanol'};   

for sample = 1:3
    endNum = sampleDetails(sample,1);
    startNum = sampleDetails(sample,2);
    minWidth = sampleDetails(sample,3);
    winSize = sampleDetails(sample,4);
    
    C = flipud(cbrewer('seq', colorMap{sample}, 2*abs(endNum-startNum)+1));
    C = C(1:abs(endNum-startNum)+1,:);
    
    counter = 1;
    tMax = 1;
for fileNum = startNum:-1:endNum 
    fileName = [filePath, 'T',sprintf('%04.0f',fileNum),'CH1.CSV'];
    if exist(fileName) ~= 2
        continue
    end        
    data = csvread(fileName,16,0);

    t = data(:,1)*1e3;
    vRaw = data(:,2);
    
    %% Detect background and subtract
    bgStart = 1;
    bgEnd = 300;
    bg = mean(vRaw(bgStart:bgEnd));
    vBG = vRaw - bg;
    vBG(abs(vBG) > 12) = 0;
    
    %% Filter the voltage signal
    vFilt = medfilt1(vBG,winSize,'truncate');

    %% Find time difference
    [pks,loc] = findpeaks((vFilt),'SortStr','descend','MinPeakWidth',minWidth);
    numPeaks = 2;
    pkLocs = loc(1:numPeaks);
    dT = abs(t(pkLocs(2)) - t(pkLocs(1)));

    %% Plot final data 
    figure(1); hold on;
    subplot(1,3,sample); hold on;
    plotBG = plot(t,vBG,'Color', [0.9 0.9 0.9]);
    plotFilt = plot(t,vFilt,'Color',C(counter,:));
    plotPeak = scatter(t(loc(1:2)),pks(1:2),[], C(counter,:));
    
    plotBG.Color(4) = 0.01 ;
    
    tMax = max(tMax,t(end));
    xlabel('Time, ms')
    ylabel('Voltage, V')
    title(sampleName(sample));
    axis([t(1),tMax,-2,10])
    
    counter = counter + 1;
end
end