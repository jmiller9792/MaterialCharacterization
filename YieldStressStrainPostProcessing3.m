clc
clear
close all
%% Inputs
dataPath = 'Z:\material\Injection Molding Materials\DSM_Akulon_K-X08581\Characterization-DSM_KX-condtionedV2\'; % end in \
dataList = {...% File list of DIC output (in csv format)
    'd12out';...
    'd05out';...
    's01out';...
    's02out';...
    's03out';...
    's07out';...
    's06out';...
    };
area = 13.02*3.31;% mm^2
trimY = 3;%MPa % No stress below this value is considered for slope calculations
yieldOffset = 0.002; % Value of 0.002 is 0.2% offset strain, last value used?
increaseOnly = 1; % Trim plasticity results after maximum stress

%% read all files to collect stress and strain data
for i = 1:size(dataList,1)
    [data] = getDataStruc(dataPath,dataList(i,1));
    data.stress = data.force/area;
    [trimIndex] =  trimToe(data.eyy,data.stress,trimY);
    data.trimIndex = trimIndex;
    data.eyy2 = data.eyy(trimIndex:end);
    data.stress2 = data.stress(trimIndex:end);
    data.exx2 = data.exx(trimIndex:end);
    B(i) = data;
end
%% get slopes
fitRangeLength = 10;
for i = 1:length(B)
    [C(i)] = processData(B(i),fitRangeLength,yieldOffset)
end

E11 = mean([C.E11])
yieldStress = mean([C.yieldstress])
%% Get Abaqus Plasticity Table based on all data fit
% [p]=polyfit(C(1).abqPlastE,C(1).abqPlastS,4)
% xs = 0:max(C(1).abqPlastE)/30:max(C(1).abqPlastE)

% For all data;
stressFit = vertcat(C.abqPlastS);
strainFit = vertcat(C.abqPlastE);
[p]=polyfit(strainFit,stressFit,6);
xs = 0:max(strainFit)/30:max(strainFit);
ys = polyval(p,xs);
ys(1) = yieldStress; % Correct first yield stress point

% trim decreasing
[~,maxSind]=max(ys);

if increaseOnly 
    xs = xs(1:maxSind);
    ys = ys(1:maxSind);
end
hold on
% Plot plasticity results
plot(xs+C(1).yieldstrain,ys,'^','markerSize',12)
%% Output CSV file (must be closed)
writematrix([ys',xs'],strcat(dataPath,'abqPlastic.csv'))
disp(strcat('Abaqus Data exported to: ',dataPath,'abqPlastic.csv'))

%% Plot All
for i = 1:length(C)
    hold on
    plot(C(i).eyy2,C(i).stress2)
    x = C(i).eyy2(C(i).fitRangeResult);
    y = x*C(i).E11;
    plot(x,y,'k','linewidth',2)
    plot(C(i).yieldstrain,C(i).yieldstress,'<')
end

%% Functions

function [S] = processData(S,fitRangeLength,yieldOffset)
    % Calculate Modulus
    [S.index,S.slope,S.yint,S.fitRangeResult] = optimalSlope(S.eyy2,S.stress2,fitRangeLength);
    S.E11 = S.slope;

    % Calculate Poisson Ratio
    %     hold on
%     plot(S.exx2,S.eyy2)
    [S.nu12v1] = poissonCalc(S.exx2,S.eyy2,S.fitRangeResult);

    % Toe compensation
    [S.eyy2,S.xoffset] = toeCompensation(S.eyy2,S.slope,S.yint);
    
    % Calculate Yield Point
    [S.yieldIndex,S.yieldstress,S.yieldstrain] = yieldOffsetStrain(S.eyy2,S.stress2,S.E11,yieldOffset);

    % Plastic Stress Strain data for ABAQUS
    [S.abqPlastE,S.abqPlastS] = abqPlastData(S.eyy2, S.stress2, S.yieldIndex)
end

function [B] = getDataStruc(datapath,samplename)
    B.name = char(samplename);
    A = readmatrix(strcat(datapath,B.name,'.csv'));
    B.eyy = A(3:end,9);
    B.exx = A(3:end,8);
    B.force = A(3:end,22);
    
end

function [slope,yint,stdquality] = checkSampleRange(DataX,DataY,fitRange)
    DataX(fitRange);
    DataY(fitRange);
    slopetemp=polyfit(DataX(fitRange),DataY(fitRange),1);
    [~,~,mu]= polyfit(DataX(fitRange),DataY(fitRange),1);
    slope = slopetemp(1);
    yint = slopetemp(2);
    stdquality = mu(2);
end

function[index,p1slope,p2yint,fitRange] = optimalSlope(DataX,DataY,fitRangeLength)
% Returns index of start of fit range, best fit slope and y intercept
    k=1; % Counter to allow skipping  data indices for speed
    for j = 1:5:length(DataX)-fitRangeLength
        fitRange  = j:j+fitRangeLength;
        [slope(k),yint(k),stdquality(k)] = checkSampleRange(DataX,DataY,fitRange);
        k = k+1;
    end
%     [~,index] = min(stdquality);
    [p1slope,index] = max(slope);
%     p1slope = slope(index);
    p2yint = yint(index);
    fitRange = index:index+fitRangeLength;
end

function[DataXnew,xoffset] = toeCompensation(DataX,slope,yint)

 % Toe compensation
%     p = polyfit(exper(i).disp(fitrange),exper(i).force(fitrange),1);
    xoffset = -yint/slope;
    DataXnew = DataX-xoffset;
end

function [nu12] = poissonCalc(exx,eyy,fitInd)
    % Poisson calculation
    nu12 = mean(-exx(fitInd)./eyy(fitInd));
end

function [yieldIndex,yieldY,yieldX] = yieldOffsetStrain(DataX,DataY,slope,offset)
    % use toe compensated results
    % offset is a fraction (0.002 = 0.2%)

%     count = 0;
%     ref = DataX*(1+offset)*slope
    ref =DataX*slope-offset*slope;
    for i = 1:length(DataX)
        if ref(i)<DataY(i)  % If data drops lower than reference slope
            yieldIndex = i;
%             if count<3 
%                 yieldIndex = i; 
%             end
%             count = count+1
%         else % if data goes back up, reset counter
%             count = 0
        end
    
    end
    yieldX = DataX(yieldIndex);
    yieldY = DataY(yieldIndex);

end

function [abqPlastE,abqPlastS] = abqPlastData(strain,stress,yieldIndex)
    abqPlastE = strain(yieldIndex:end)-strain(yieldIndex);
    abqPlastS = stress(yieldIndex:end);
%     abqPlast = [abqPlastS,abqPlastE]
end

% function [DataXtrim,DataYtrim] = trimToe(DataX,DataY,trimY)
function [trimIndex] = trimToe(DataX,DataY,trimY)
%     count = 0;
    for i = 1:length(DataX)
        if trimY>DataY(i)% && count<3 % if data is less than the trim value
            trimIndex = i;
%             count = count+1
        end
    end
    
%     DataXtrim = DataX(firstIndex:end);
%     DataYtrim = DataY(firstIndex:end);
end
