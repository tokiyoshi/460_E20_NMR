close all
clear

filePath = char("C:\Users\Benjamin's PC\Dropbox\School Work\PHYS460A\Experiment 20\460_E20_NMR\rawdata\");
sampleDetails = [20, 33, 1, 0.1, 1;...
                 34, 43, 1, 0.1, 1;...
                 47, 56, 1, 0.1, 15];

colorMap = {'YlGnBu', 'YlGn','YlOrRd'};
sampleName = {'Doped Water', 'Rubber', 'Ethanol'};  

irData = [sampleName;cell(2,3)];

for sample = 1:3
    endNum = sampleDetails(sample,1);
    startNum = sampleDetails(sample,2);
    minWidth = sampleDetails(sample,3);
%     winSize = sampleDetails(sample,4);
    cutOff = sampleDetails(sample,5);
    
    C = flipud(cbrewer('seq', colorMap{sample}, round(1.5*abs(endNum-startNum)+1)));
    C = C(1:abs(endNum-startNum)+1,:);
    
    counter = 1;
    tMax = 1;
    tMin = 0;
    
    dTVals = [];
    mTMaxVals = [];
    
    for fileNum = startNum:-1:endNum 
        fileName = [filePath, 'T',sprintf('%04.0f',fileNum),'CH1.CSV'];
        if exist(fileName) ~= 2
            continue
        end        
        data = csvread(fileName,16,0);

        t = data(:,1)*1e3;
        vRaw = data(:,2);

        winSize = sampleDetails(sample,4)/(t(2)-t(1));

        %% Detect background and subtract
        bgStart = 1;
        bgEnd = 300;
        bg = mean(vRaw(bgStart:bgEnd));
        vBG = vRaw - bg;
        vBG(abs(vBG) > 12) = 0;

        %% Filter the voltage signal
    %     vFilt = medfilt1(vBG,winSize,'truncate');
    %     vFilt = sgolayfilt(vBG,3,winSize);
        vFilt = smoothdata(vBG,'gaussian',winSize);%,'lowess');

        %% Find time difference

        [pks,loc] = findpeaks(abs(vFilt),'SortStr','descend','MinPeakWidth',minWidth);
        numPeaks = 1;
        while t(loc(numPeaks)) < cutOff
            numPeaks = numPeaks + 1;
        end

        pkLocs = loc(numPeaks);

        dT = abs(t(pkLocs));

        dTVals(counter) = dT;
        mTMaxVals(counter) = vFilt(pkLocs);
        
        %% Plot final data 
        figure(1); hold on;
        subplot(3,1,sample); hold on;
        plotBG = plot(t,vBG,'Color', [0.9 0.9 0.9]);
        plotFilt = plot(t,vFilt,'Color',C(counter,:));
        plotPeak = scatter(t(pkLocs),vFilt(pkLocs),50, C(counter,:));

        plotBG.Color(4) = 0.2;

        tMax = max(tMax,t(end));
        tMin = min(tMin,t(1));

        xlabel('Time, ms')
        ylabel('Voltage, V')
        title(sampleName(sample));
        axis([tMin,tMax,-15,10])

        counter = counter + 1;
    end
    
    irData{2,sample} = dTVals;
    irData{3,sample} = mTMaxVals;
    
end


%% Fit exponential to calculated values
C = cbrewer('qual', 'Paired', 6);
refFileNums = [3,5,4];

for sample = 1:3
    
    refName = [filePath, 'T',sprintf('%04.0f',refFileNums(sample)),'CH1.CSV'];
    data = csvread(refName,16,0);

    tRef = data(:,1)*1e3;
    refRaw = data(:,2);
    bgStart = 1;
    bgEnd = 300;
    bg = mean(refRaw(bgStart:bgEnd));
    refBG = refRaw - bg;
    
    winSize = sampleDetails(sample,4)/(tRef(2)-tRef(1));
    refFilt = smoothdata(refBG,'gaussian',winSize);
    
    [pks,loc] = findpeaks((refFilt),'SortStr','descend','MinPeakWidth',minWidth);
    numPeaks = 1;
    pkLocs = loc(numPeaks);
    
    M0 = refFilt(pkLocs);

    dTVals = irData{2,sample};
    mTMaxVals = irData{3,sample};
    tInt = linspace(0,max(dTVals),1000);

    
    y = log(1/2-mTMaxVals./(2*M0));
    x = dTVals;
    
    fitTypeLn = fittype('-x/T1',...
                        'dependent',{'y'},'independent',{'x'},...
                        'coefficients',{'T1'});
    vFitObjLn = fit(x',y',fitTypeLn)
    
    fprintf('Sample: %s, T1: %0.9f \n\n', char(sampleName(sample)), vFitObjLn.T1)
    vFitLn = M0.*(1-2.*exp(feval(vFitObjLn,tInt)));

%     vFitUpper = M0.*
    
    [temp,locMin] = min(abs(vFitLn));
    tX = tInt(locMin)/log(2)
    
    figure(2); hold on;
    subplot(1,3,sample); hold on;
    scatter(dTVals,mTMaxVals,[], C(sample*2-1,:),'filled');
%     scatter(tInt(locMin),temp)
    fitPlot(sample) = plot(tInt,vFitLn', 'Color', C(sample*2,:));
    xlabel('Time, ms')
    ylabel('Magnetization Amplitude')
    title(sampleName(sample))
end
