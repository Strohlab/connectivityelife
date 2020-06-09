%
% YOU NEED THE FUNCTION fdr_bh FOR RUNNING THIS CODE!
%


clear all
close all

%Load the atlas
cd 'E:\ATLAS RATS';
myvoi = xff('ATLAS_ROIS_CORTEX.voi');
covarvoi= xff('callosum_ventricles.voi');
%Load the functional data (VTC)
[FileName,PathName] = uigetfile('*.mdm','Select the mdm with the MRI functional data');
allpaths=PathName;
allfiles=FileName;


myfolder = uigetdir('E:\','select the FIBER folder');
cd (myfolder)
myfiles=dir('*.mat');

%%%% Which fiber
choice=questdlg('Which Fiber contains the signal?','FIber Signal',...
'1','2','Cancel','1');
switch choice
case '1'
   fibernr=1;
case '2'
   ffibernr=0;
case 'Cancel'
    disp('User canceled');
    return
end

prompt = {'Enter amount of volumes:','Enter the TR:','Enter fiber sampling rate','remove the first ... vols'};
title = 'fMRI Acquisition';
dims = [1 35];
definput = {'600','1.5','2000','5'};
answer = inputdlg(prompt,title,dims,definput);
fMRIvals = cellfun(@(x)str2double(x), answer);
endrec=fMRIvals(1)*fMRIvals(2)*fMRIvals(3);

%variables to identify events (taken from Andreas script)
prompt2 = {'min_length (ms)','min_refrac (ms)','min_highA (%std)','min_highB (%std)',...
       'addon (ms)','sampl_freq (kHz)','Window_Base (ms)','Window_ExpAvg (ms)','downsampling'};
dlg_title2 ='Variables for identifying onsets:(see Seamari et al 2013)';
num_lines2 = 1;
defaultans2 = {'500','1000','30','10','0','2','2500','25','1'};
answer2 = inputdlg(prompt2,dlg_title2,[1, length(dlg_title2)+30],defaultans2);    
min_length=str2double(answer2{1,1}); %(event in ms)
min_refrac=str2double(answer2{2,1}); %(event in ms)
min_highA=str2double(answer2{3,1}); 
min_highB=str2double(answer2{4,1});
addon=str2double(answer2{5,1});
sampl_freq=str2double(answer2{6,1}); %(khz)
Window_Base=str2double(answer2{7,1});
Window_ExpAvg=str2double(answer2{8,1});
downsampling=str2double(answer2{9,1});
%trsnform above from ms to samplingrate
min_length=min_length*sampl_freq;
min_refrac=min_refrac*sampl_freq;
Window_Base=Window_Base*sampl_freq;
Window_ExpAvg=Window_ExpAvg*sampl_freq;
addon=addon*sampl_freq;

cd (allpaths);
mymdm = xff(allfiles);
root_vtc=mymdm.XTC_RTC;
rois_pos=zeros(myvoi.NrOfVOIs,3);

all_corr1=zeros(myvoi.NrOfVOIs,myvoi.NrOfVOIs,length(root_vtc));
all_pval1=zeros(myvoi.NrOfVOIs,myvoi.NrOfVOIs,length(root_vtc));
all_sign1=[];
datframe=[]; 
for vv=1:length(root_vtc)
    myvtc=xff(root_vtc{vv,1});
    %%%Getting rat number%%%
    myrat=strsplit(root_vtc{vv,1},'/');
    myrat=myrat{end};
    myrat=strsplit(myrat,'_');
    mycond=cell2mat(myrat(2));
    mycond=mycond(1:end-4);
    myrat=cell2mat(myrat(1));
    myrat=myrat(end-1:end);
    %%%%%%%%%%%
    factor=myvtc.Resolution;
    xstart=myvtc.XStart;
    ystart=myvtc.YStart;
    zstart=myvtc.ZStart;
    vtcdata = myvtc.VTCData;
    clear myvtc;
    timecourse=size(vtcdata,1);
    %%%%%%%%%% collect the averages for each ROI %%%%%%%%
    Rois_mean=zeros(myvoi.NrOfVOIs,timecourse);
    for nn=1:myvoi.NrOfVOIs;
        roi=myvoi.VOI(nn).Voxels;
        rois_pos(nn,:)=mean(roi,1);
        roi_timecourse=zeros(size(roi,1),timecourse);
        if min(roi(:,1))<=xstart  | min(roi(:,1))<=xstart  | min(roi(:,1))<=xstart
            warndlg(['Roi dimension does not fit the data dimension for the ROI Nr: ',num2str(nn)],'script will crash');
            return
        end
        for mm=1:size(roi,1)
            roi_timecourse(mm,:)=vtcdata(:,floor((roi(mm,1)-(xstart-1))/factor),floor((roi(mm,2)-(ystart-1))/factor),floor((roi(mm,3)-(zstart-1))/factor));
        end
        Rois_mean(nn,:)=mean(roi_timecourse,1);
    end
    all_sign1=cat(3, all_sign1, Rois_mean);
    %%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%% Collect the COVARS ROIs (CALLOSUM+VENTRICLES)
    Covar_mean=zeros(covarvoi.NrOfVOIs,timecourse);
    for nn=1:covarvoi.NrOfVOIs;
        roi=covarvoi.VOI(nn).Voxels;
        %rois_pos(nn,:)=mean(roi,1);
        roi_timecourse=zeros(size(roi,1),timecourse);
        if min(roi(:,1))<=xstart  | min(roi(:,1))<=xstart  | min(roi(:,1))<=xstart
            warndlg(['Roi dimension does not fit the data dimension for the ROI Nr: ',num2str(nn)],'script will crash');
            return
        end
        for mm=1:size(roi,1)
            roi_timecourse(mm,:)=vtcdata(:,floor((roi(mm,1)-(xstart-1))/factor),floor((roi(mm,2)-(ystart-1))/factor),floor((roi(mm,3)-(zstart-1))/factor));
        end
        Covar_mean(nn,:)=mean(roi_timecourse,1);
    end
    %%%%%%%%%%%%% LOAD THE BREATHING AND USE IT AS COVAR %%%%%%
    cd (folderbreath)
    if ww==1
        mybreath=load(['R',num2str(myrat),'_ISO.mat']);
        channels=fields(mybreath);
    elseif ww==2
         mybreath=load(['R',num2str(myrat),'_MED.mat']);
         channels=fields(mybreath);
    end
    for rr=1:length(channels)
        lala=strcat('mybreath.',channels(rr));
        cur_file=eval(lala{1});
        if strcmp(cur_file.title,'breath')  
            breath=cur_file.values;
        elseif strcmp(cur_file.title,'ScanTrig') 
            trig=cur_file.times(end);
        end
    end
    breath=breath(trig:[1000*60*30+trig]); % collect data from the 1st scan until the end
    breath=downsample(breath,1500);% subsample 
    breath=breath(6:end-1); %remove first 5 vols
    %%%%%%%%% z-score data %%%%%%%%
    all_covar=vertcat(Covar_mean,breath');
    zcovar=zscore(all_covar,0,2);
    zrois=zscore(Rois_mean,0,2);
    %%%%%%%%%%%% Connectivity to ROI %%%%%%%%%%%
    [R,P] = partialcorr(zrois',zcovar');
    all_corr1(:,:,vv)=R;  
    all_pval1(:,:,vv)=P; 
    figure(vv);
    imagesc(R); % plot the matrix
    colormap('hot'); % set the colorscheme
    colorbar; % enable colorbar
    caxis([-1 1])
    title(['RAT ',myrat,' CORREL ISOFL'])
    
    %%% ANALYZE FIBER %%%%
    if strcmp(mycond,'visstim')
    myfilename=(['RAT_',myrat,'_epis_visual.mat']);
    elseif strcmp(mycond,'resting')
    myfilename=(['RAT_',myrat,'_epis_resting.mat']);
    end
    
    for aa=1:length(myfiles)
    if strcmp(myfiles(aa).name,myfilename)
        cd(myfiles(aa).folder)
        subjruns=load(myfiles(aa).name);
        subjnameruns= fieldnames(subjruns);
        for rr=1:length(subjnameruns) %loop on the different channels
            eval(['myname=subjruns.',subjnameruns{rr},';']);
            if strcmp(myname.title,'CaRec1')    
                carec1=myname.values;
            elseif strcmp(myname.title,'CaRec2')
                carec2=myname.values;
            elseif strcmp(myname.title,'scantrig')
                trigger=myname.times;
            elseif strcmp(myname.title,'fiber2')
                stim=myname.values;
            end
        end
        clear subjruns subjnameruns;
    else
    end
    end
    if fibernr==2
            carec=carec2(floor(trigger(3)*2000):floor(trigger(3)*2000)+endrec); 
        else    
             carec=carec1(floor(trigger(3)*2000):floor(trigger(3)*2000)+endrec);  
    end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        % Detrend and invert sign of the signal (YOU NEED THE FUNCTION
        % BASELINE) AND U NEED TO KNOW THE SIGN OF YOUR SIGNAL
        % BECAUSE SOME PEOPLE SOMETIMES USE THE  2 FIBERS WITH INVERTED
        % SIGNS
        Signal1=Baseline(-carec,Window_Base);

        % predefining some variable-space
        thresh_start=Signal1;
        thresh_start1(:)=0; % will contain the threshold to define the onset of an up-state

        box1=Signal1; %binary vector, 0= down, 1 = up, preliminary
        box1(:)=0;
        onset1=box1; % binary vector containing only the onsets of up-states
        box2=box1;
        onset2=box2;


        % Calculating the final thresholds and short-averaged data
        % Determination of the thresholds for on- and offsets of an up-state
        % using a histogram SIGNAL1
        [n,xout] = hist(Signal1(2500:numel(Signal1)-2500),ceil((numel(Signal1)-5000)/1000));
        sum_n=cumsum(n);
        off_hist=find(sum_n >= max(sum_n)*(1-min_highA/100),1,'first');
        min_high1=xout(off_hist); % new definition of min_high1 based on the calculated results
        off_hist=find(sum_n >= max(sum_n)*(1-min_highB/100),1,'first');
        min_high2=xout(off_hist); % new definition of min_high2 based on the calculated results

        %Moving average with exponential weighting (s. Seamari, PlosOne, 2007)
        thresh_start1(:)=min_high1; %min_high1
        Signal_filt = tsmovavg(Signal1, 'e', Window_ExpAvg, 1);

        % calculate a preliminary on-off-vactor
        box1(Signal_filt >= thresh_start1)=0.001;

        % remove too short gaps between upstates SIGNAL1
        % here the variable min_refrac is used
        count=0;
        for i=2:numel(Signal1)

            %start to count when upstate ends
            if box1(i)==0
                count=count+1;
            end

            %if the gap is shorter than min_refrac
            %set the gap to 1
            if (count < min_refrac) && (box1(i)==0.001) && (box1(i-1)==0) 
                box1(i-count:i-1)=0.001;
                count=0;
            end

            %if the gap is longer than min_refrac
            %reset the counter
            if (count >= min_refrac) && (box1(i)==0.001) && (box1(i-1)==0)
               count=0;
            end

        end


        % correct offsets SIGNAL1 
        for z=2:numel(Signal1)
            %search for an offset
            if box1(z)==0 && box1(z-1)==0.001    
                j=z;
                %keep the vector on 1 as long as the signal is high enough
                if j+z <= numel(Signal1) && Signal1(j+1) >= min_high1*0.15 
                    j=j+1;
                    box1(z)=0.001;
                end
            end
            z=j+1;
        end

        % remove suspicious upstates Signal1
        count=0;

        for u=2:numel(Signal1)
            %start counting when an upstate starts
            if box1(u)==0.001
                count=count+1;
            end
            %if the upstate is too short (shorter than min_length)
            %set it to 0
            if (box1(u)==0) && (box1(u-1)>0) && (count < min_length)
                box1(u-count:u-1)=0;
                count=0;
            end
            %if the upstate is long enough, test if it is high enough
            %(higher than min_high2)
            %otherwise set it to 0
            if (box1(u)==0) && (box1(u-1)>0) && (count > min_length)
                if max(Signal1(u-count:u-1)) < min_high2
                    box1(u-count:u-1)=0;
                end
                count=0;
            end   
        end

        % find the S
       if strcmp(mycond,'visual') 
           stimon=[];
           for nn=2:length(box1)
               if exist('stim') && (stim(nn)>0.5) && (stim(nn-1)<0.5)
                  stimon=vertcat(stimon,nn);
               else  
               end
           end
       end
       onsets=[];
       for mm=2:length(box1)
           if box1(mm) == 0.001 && box1(mm-1) == 0
               onsets=vertcat(onsets,mm);
           else  
           end
       end
       offsets=[];
       for pp=1:(length(box1)-1)
           if box1(pp) == 0.001 && box1(pp+1) == 0
               offsets=vertcat(offsets,pp); 
           else  
           end
       end
       predon=onsets;
       predon=(floor(predon/(fMRIvals(3)*fMRIvals(2)))+1)-fMRIvals(4);
       predon(predon<=0)=[];
       
       
        %%%% CORRELATE EVENTS
        tbin=180/1;  %(time of each bin in secs/TR)
        jump=30/1;   %(time of jumping bin in secs/TR)
        vols=1195; % Volumes 
        nrpts=floor((vols-tbin)/jump)+1;
        
        for hh=1:nrpts
            zcurdata=zrois(:,(1+jump*(hh-1)):(jump*(hh-1)+tbin));
            zcurcovar=zcovar(:,(1+jump*(hh-1)):(jump*(hh-1)+tbin));
            upstates= sum(predon>=(1+jump*(hh-1)) & predon<=(jump*(hh-1)+tbin));
            %%%%%%%%%%%% Connectivity values%%%%%%%%%%%
            [R,P] = partialcorr(zcurdata',zcurcovar');
            mydat=R(:);
            mydat(mydat==1)=[];
            mydat=unique(mydat);
            signif=sum(mydat>=0.35);
            pd_norm = fitdist(mydat, 'normal');
            datframe=vertcat(datframe,[pd_norm.mu,upstates,signif,vv]); 
        end      
        
end

vv=max(datframe(:,4));
corrvals=zeros(vv,2);
permnum=50000;
for jj=1:vv
    mydaset=datframe(datframe(:,4)==jj,:);
    figure(jj);
    plot(mydaset(:,2),mydaset(:,3),'ko', 'MarkerSize',19,'MarkerEdgeColor','k',...
    'MarkerFaceColor',[0.7 0.7 0.7]);
    l = lsline ;
    set(l,'LineWidth', 4);
    [R,P]=corr(mydaset(:,2),mydaset(:,3));
    disp([R,P])
    corrvals(jj,:)=[R P];
    xlim([(min(mydaset(:,2))-2) (max(mydaset(:,2)))+2]);
    ylim([(min(mydaset(:,3))-15) (max(mydaset(:,3)))+15]);
    set(gcf,'units','points','position',[10,10,500,500])
    %%%%% Start permutation test%%%%
    allR=[];
    allP=[];
    for zz=1:permnum
        dataA=mydaset(:,2);
        dataB=mydaset(:,3);
        dataA=dataA(randperm(length(dataA)));
        dataB=dataB(randperm(length(dataB)));
        [R1,P1]=corr(dataA,dataB);
        allR=vertcat(allR,R1);
        allP=vertcat(allP,P1);
    end
    pd = fitdist(allR,'Normal');
    x_values = -1:0.001:1;
    y = pdf(pd,x_values);
    figure(100+jj)
    area( x_values,y,'FaceColor',[0.8 0.8 0.8],'LineStyle','none')
    hold on
    plot(R,0,'kd','MarkerSize',25,'MarkerFaceColor','k')
    ylim([-0.25 max(y)+0.25])   
end
    csvwrite('correlations.csv',daqtframe)








    