function [Korrespondenzen] = punkt_korrespondenzen2(I1,I2,Mpt1,Mpt2,varargin)
% In dieser Funktion sollen die extrahierten Merkmalspunkte aus einer
% Stereo-Aufnahme mittels NCC verglichen werden um Korrespondenzpunktpaare
% zu ermitteln.
P=inputParser;
P.addOptional('window_length', 15, @isnumeric);
P.addOptional('min_corr', 0.7, @isnumeric);
P.addOptional('do_plot', false, @islogical);
P.parse(varargin{:});

winLen=P.Results.window_length;
minCorr=P.Results.min_corr;
doPlot=P.Results.do_plot;


NCorrMax=min(size(Mpt1,2),size(Mpt2,2));        %maximum number of correlation points
Nw=winLen*winLen;                           %number of points in a window
Wn1=zeros(size(Mpt1,2),Nw+1);
Wn2=zeros(Nw+1,size(Mpt2,2));           %the last row indicates how many points are in the window
img1Size=size(I1);
img2Size=size(I2);
Korrespondenzen=zeros(6,NCorrMax);

%fetch local points from image 1
for Midx=1:size(Mpt1,2)

    yu=max(Mpt1(2,Midx)-floor(winLen/2),1);
    yd=min(Mpt1(2,Midx)+floor(winLen/2),img1Size(1));
    xl=max(Mpt1(1,Midx)-floor(winLen/2),1);
    xr=min(Mpt1(1,Midx)+floor(winLen/2),img1Size(2));
    Wn1(Midx,Nw+1)=(xr-xl+1)*(yd-yu+1);
    Wn1(Midx,1:Wn1(Midx,Nw+1))=reshape(I1(yu:yd , xl:xr), 1, Wn1(Midx,Nw+1));
    W_mean=mean(Wn1(Midx , 1:Wn1(Midx,Nw+1)));
    W_sigma=std(Wn1(Midx , 1:Wn1(Midx,Nw+1)));
    Wn1(Midx , 1:Wn1(Midx,Nw+1))=(Wn1(Midx , 1:Wn1(Midx,Nw+1))-W_mean)/W_sigma;
end

%fetch local pixels from image 2
for Midx=1:size(Mpt2,2)

    yu=max(Mpt2(2,Midx)-floor(winLen/2),1);
    yd=min(Mpt2(2,Midx)+floor(winLen/2),img2Size(1));
    xl=max(Mpt2(1,Midx)-floor(winLen/2),1);
    xr=min(Mpt2(1,Midx)+floor(winLen/2),img2Size(2));
    Wn2(Nw+1,Midx)=(xr-xl+1)*(yd-yu+1);
    Wn2(1:Wn2(Nw+1,Midx),Midx)=reshape(I2(yu:yd , xl:xr),Wn2(Nw+1,Midx),1);
    W_mean=mean(Wn2(1:Wn2(Nw+1,Midx),Midx));
    W_sigma=std(Wn2(1:Wn2(Nw+1,Midx),Midx));
    Wn2(1:Wn2(Nw+1,Midx),Midx)=(Wn2(1:Wn2(Nw+1,Midx),Midx)-W_mean)/W_sigma;
end

%calculate NCC
NCC=0;NCCmax=0;
corrCnt=1;
idx2Corr=1;
for idx1=1:size(Mpt1,2)
    for idx2=1:size(Mpt2,2)
        if Wn2(Nw+1,idx2)~=Wn1(idx1,Nw+1)
            continue;
        end
        NCC=(Wn1(idx1,1:Nw)*Wn2(1:Nw,idx2))/Wn1(idx1,Nw+1);
        if NCC>NCCmax && NCC>minCorr                            %found the largest NCC, which needs to be > minCorr
            NCCmax=NCC;
            Korrespondenzen(1:5,corrCnt)=[Mpt1(1,idx1);Mpt1(2,idx1);Mpt2(1,idx2);Mpt2(2,idx2);NCCmax];
            idx2Corr=idx2;
        end
        
    end
    if NCCmax~=0
        Wn2(Nw+1,idx2Corr)=0;                   %if the corresponding point in image is found, this point can be skipped.
        Korrespondenzen(6,corrCnt)=corrCnt;
        corrCnt=corrCnt+1;
        NCCmax=0;
    end
end
Korrespondenzen=Korrespondenzen(:,1:corrCnt-1);
if doPlot
    showMatchedFeatures(I1,I2,Korrespondenzen(1:2,:)',Korrespondenzen(3:4,:)')
end

end


