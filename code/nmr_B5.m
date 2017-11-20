close all
clear
warning off

filePath = char("C:\Users\Benjamin's PC\Dropbox\School Work\PHYS460A\Experiment 20\460_E20_NMR\rawdata\");

sampleDetails = [57,64,50, 100;...
                 66,71,30, 100;...
                 72,82,20, 75];
             
colorMap = {'YlGnBu', 'YlGn','YlOrRd'};
sampleName = {'Doped Water', 'Rubber', 'Ethanol'};   

echoData = [sampleName;cell(2,3)];

for sample = 1:3
    endNum = sampleDetails(sample,1);
    startNum = sampleDetails(sample,2);
    minWidth = sampleDetails(sample,3);
    winSize = sampleDetails(sample,4);
    
    dTVals = [];
    mTMaxVals = [];
    
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

        echoLoc = max(pkLocs);

        dTVals(counter) = dT;
        mTMaxVals(counter) = vFilt(echoLoc);

        %% Plot final data 
        figure(1); hold on;
        subplot(3,1,sample); hold on;
        plotBG = plot(t,vBG,'Color', [0.9 0.9 0.9]);
        plotFilt = plot(t,vFilt,'Color',C(counter,:));
        plotPeak = scatter(t(echoLoc),vFilt(echoLoc),[], C(counter,:));

        plotBG.Color(4) = 0.01 ;

        tMax = max(tMax,t(end));
        xlabel('Time, ms')
        ylabel('Voltage, V')
        title(sampleName(sample));
        axis([t(1),tMax,-2,10])

        counter = counter + 1;
    end

    echoData{2,sample} = dTVals;
    echoData{3,sample} = mTMaxVals;

end

%% Fit exponential to calculated values
C = cbrewer('qual', 'Paired', 6);
for sample = 1:3
    dTVals = echoData{2,sample};
    mTMaxVals = echoData{3,sample};
    tInt = linspace(0,max(dTVals),1000);

    
    fitTypeLn = fittype('-x/T2 + B',...
                        'dependent',{'y'},'independent',{'x'},...
                        'coefficients',{'T2','B'});
    vFitObjLn = fit(dTVals',log(mTMaxVals)',fitTypeLn)
    
    fprintf('Sample: %s, T2: %0.9f \n\n', char(sampleName(sample)), vFitObjLn.T2)
    
    vFitLn = exp(feval(vFitObjLn,tInt));
    
    
%     fitType = fittype('M0.*exp(-(x-x0)/T2)',...
%                         'dependent',{'y'},'independent',{'x'},...
%                         'coefficients',{'M0','T2','x0'});
%     vFitObj = fit(dTVals',mTMaxVals',fitType)
%     vFit = feval(vFitObj,tInt);
    
    figure(2); hold on;
    scatter(dTVals,mTMaxVals,[], C(sample*2-1,:));
    fitPlot(sample) = plot(tInt,vFitLn', 'Color', C(sample*2,:));
%     fitPlot(sample) = plot(tInt,vFit,'Color',C(sample*2,:));%,'DisplayName',char(sampleName(sample)));
    axis([0,100,0,10])
    xlabel('Time, ms')
    ylabel('Echo Amplitude')
    legend(fitPlot,'Doped Water', 'Rubber', 'Ethanol')
end


