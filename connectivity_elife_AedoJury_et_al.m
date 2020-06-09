%
% YOU NEED TO DOWNLOAD THE FUNCTION fdr_bh FOR RUNNING THIS CODE!
% It is optimized for visualizing data in Brain Voyager
% Therefore, all the extention of the files are BV format. 
%


clear all
close all

%Load the atlas ROIS
cd 'E:\ATLAS RATS';
myvoi = xff('ATLAS_ROIS_CORTEX.voi');
covarvoi= xff('callosum_ventricles.voi');
cd 'F:\Brain States Connectivity'
myvoi2 = xff('VOI_CONNECTIVITY_ICA.voi');


%Load the functional data (VTC)
[FileName,PathName] = uigetfile('*.mdm','Select the mdm with the ISOFL functional data');
cd 'F:\Brain States Connectivity'
[FileName2,PathName2] = uigetfile('*.mdm','Select the mdm with the MEDITO functional data');


allpaths={PathName,PathName2};
allfiles={FileName,FileName2};


pxs1= inputdlg('Enter the cutoff for the correlation significance...');
Signcutoff=str2double(pxs1{1});
cd 'F:\Brain States Connectivity'
[FileName4,PathName4] = uigetfile('*.mat','Select the mat file containing the ICA components');
cd (PathName4)
myICAcomp=load(FileName4);
ICA=cell2mat(struct2cell(myICAcomp));

choice1=questdlg('Do You want to use as covar the Global mean','Global Mean',...
    'Yes','No','Cancel','Yes');
switch choice1
    case 'Yes'
       globmean=1;
    case 'No'
       globmean=0;
    case 'Cancel'
        disp('User canceled');
        return
end

%%% Load the 3d version of the stadarized brain brain template %%%
cd 'E:\ATLAS RATS';
ratbrain = xff('brain_color_SAG_TRF_ACPC.vmr');
ratbrain=ratbrain.VMRData;
renderbrain=im2double(ratbrain);
renderbrain(renderbrain==0)=nan;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%collect values per ROI
for ww=1:length(allfiles);
    cd (allpaths{ww});
    mymdm = xff(allfiles{ww});
    root_vtc=mymdm.XTC_RTC;
    rois_pos=zeros(myvoi.NrOfVOIs,3);
    rois_size=zeros(myvoi.NrOfVOIs,1);
    if ww==1
        all_corr1=zeros(myvoi.NrOfVOIs,myvoi.NrOfVOIs,length(root_vtc));
        all_pval1=zeros(myvoi.NrOfVOIs,myvoi.NrOfVOIs,length(root_vtc));
        all_sign1=[];
        all_covar1=[];
        zall_ICA_ROI1=zeros(myvoi2.NrOfVOIs,1195,length(root_vtc));
        all_corr4=zeros(myvoi.NrOfVOIs,myvoi.NrOfVOIs,length(root_vtc));
        all_pval4=zeros(myvoi.NrOfVOIs,myvoi.NrOfVOIs,length(root_vtc));
        all_corr5=zeros(myvoi.NrOfVOIs,myvoi.NrOfVOIs,length(root_vtc));
        all_pval5=zeros(myvoi.NrOfVOIs,myvoi.NrOfVOIs,length(root_vtc));
    elseif ww==2
        all_corr2=zeros(myvoi.NrOfVOIs,myvoi.NrOfVOIs,length(root_vtc));
        all_pval2=zeros(myvoi.NrOfVOIs,myvoi.NrOfVOIs,length(root_vtc));
        all_sign2=[];
        all_covar2=[];
        zall_ICA_ROI2=zeros(myvoi2.NrOfVOIs,1195,length(root_vtc));
    else
        all_corr3=zeros(myvoi.NrOfVOIs,myvoi.NrOfVOIs,length(root_vtc));
        all_pval3=zeros(myvoi.NrOfVOIs,myvoi.NrOfVOIs,length(root_vtc));
        all_sign3=[];
        all_covar3=[]; 
    end
    
    for vv=1:length(root_vtc)
        myvtc=xff(root_vtc{vv,1});
        %Getting rat number%%%
        myrat=strsplit(root_vtc{vv,1},'/');
        myrat=myrat{end};
        myrat=str2double(myrat(2:3));
        %%%%%%%%%%% size factor between VTC andROI
        factor=myvtc.Resolution;
        xstart=myvtc.XStart;
        ystart=myvtc.YStart;
        zstart=myvtc.ZStart;
        vtcdata = myvtc.VTCData;
        clear myvtc;
        timecourse=size(vtcdata,1);
        %%%%%%%%%% collect the averages for each ROI
        %%%% ROIs are in a different space dimmension than EPis seq.
        %%% therefore they need to be transformed by the FACTOR
        Rois_mean=zeros(myvoi.NrOfVOIs,timecourse);
        globroi=[];
        for nn=1:myvoi.NrOfVOIs;
            roi=myvoi.VOI(nn).Voxels;
            rois_size(nn)=size(myvoi.VOI(nn).Voxels,1);
            globroi=vertcat(globroi,roi);
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
        %%%
        if ww==1
            all_sign1=cat(3, all_sign1, Rois_mean);
        elseif ww==2
            all_sign2=cat(3, all_sign2, Rois_mean);
        else
            all_sign3=cat(3, all_sign3, Rois_mean);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% ICA ROIS ANALYSIS %%%%%%%%%%%%%%%%%%%
        Rois_mean2=zeros(myvoi2.NrOfVOIs,timecourse);
        for nn=1:myvoi2.NrOfVOIs;
            roi=myvoi2.VOI(nn).Voxels;
            roi_timecourse=zeros(size(roi,1),timecourse);
            if min(roi(:,1))<=xstart  | min(roi(:,1))<=xstart  | min(roi(:,1))<=xstart
                warndlg(['Roi dimension does not fit the data dimension for the ROI Nr: ',num2str(nn)],'script will crash');
                return
            end
            for mm=1:size(roi,1)
                roi_timecourse(mm,:)=vtcdata(:,floor((roi(mm,1)-(xstart-1))/factor),floor((roi(mm,2)-(ystart-1))/factor),floor((roi(mm,3)-(zstart-1))/factor));
            end
            Rois_mean2(nn,:)=mean(roi_timecourse,1);
            
        end
        %%%
        if ww==1
            zall_ICA_ROI1(:,:,vv)=zscore(Rois_mean2,[ ],2);
        elseif ww==2
            zall_ICA_ROI2(:,:,vv)=zscore(Rois_mean2,[ ],2);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%% calculate the Global mean %%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if globmean==1
        roi_timecourse=zeros(size(globroi,1),timecourse);
            if min(globroi(:,1))<=xstart  | min(globroi(:,1))<=xstart  | min(globroi(:,1))<=xstart
                warndlg(['Roi dimension does not fit the data dimension for the ROI Nr: ',num2str(nn)],'script will crash');
                return
            end
            for mm=1:size(globroi,1)
                roi_timecourse(mm,:)=vtcdata(:,floor((globroi(mm,1)-(xstart-1))/factor),floor((globroi(mm,2)-(ystart-1))/factor),floor((globroi(mm,3)-(zstart-1))/factor));
            end
            globmean=mean(roi_timecourse,1);
        else
        end
        %%%%%%%%%% Lets collect the COVARS ROIs (CALLOSUM+VENTRICLES)
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
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%% LOAD THE BREATHING AND USE IT AS COVAR %%%%%%
        cd (folderbreath)
        if ww==1
            mybreath=load(['R',num2str(myrat),'_ISO.mat']);
            channels=fields(mybreath);
        elseif ww==2
             mybreath=load(['R',num2str(myrat),'_MED.mat']);
             channels=fields(mybreath);
        else 
            mybreath=load(['R',num2str(myrat),'_LOW.mat']);
            channels=fields(mybreath);
            
        end
        for rr=1:length(channels)
            lala=strcat('mybreath.',channels(rr));
            cur_file=eval(lala{1});
            if strcmp(cur_file.title,'breath')  
                breath=cur_file.values;
            elseif strcmp(cur_file.title,'ScanTrig')
                if length(cur_file.times)>=1
                trig=cur_file.times(end);
                else
                trig=20000;
                end
            end
        end
        trig=floor(trig*1000);
        if 1000*60*30+trig<=length(breath)
        breath=breath(trig:1000*60*30+trig); % collect data from the 1st scan until the end
        else
        trig=2000;
        breath=breath(trig:1000*60*30+trig);
        end
        breath=downsample(breath,1500);% subsample 
        breath=breath(6:end-1); %remove first 5 vols (AS IN THE VTC!)
        
       
        %%%%%%%%% Lets collect all covars %%%%%%%%
        if globmean==1
            all_covar=vertcat(Covar_mean,globmean);
        else
            all_covar=vertcat(Covar_mean);
        end
        
        %Save covar
        if ww==1
            all_covar1=cat(3, all_covar1,all_covar);
            covarICA=vertcat(all_covar,ICA(vv,:));
            covarINV=vertcat(all_covar,fliplr(ICA(vv,:)));
        elseif ww==2
            all_covar2=cat(3, all_covar2,all_covar);  
        else
            all_covar3=cat(3, all_covar3,all_covar);  
        end
        
        %%%%%%%
        zcovar=zscore(all_covar,0,2);
        zrois=zscore(Rois_mean,0,2);
        %%%%%%%%%%%% Connectivity to ROI %%%%%%%%%%%
        [R,P] = partialcorr(Rois_mean',all_covar');
       
        if ww==1
            [Rica,Pica] = partialcorr(Rois_mean',covarICA');
            [Rinv,Pinv] = partialcorr(Rois_mean',covarINV');
            all_corr1(:,:,vv)=R;  
            all_pval1(:,:,vv)=P;
            all_corr4(:,:,vv)=Rica;  
            all_pval4(:,:,vv)=Pica;
            all_corr5(:,:,vv)=Rinv;  
            all_pval5(:,:,vv)=Pinv;
        elseif ww==2
            all_corr2(:,:,vv)=R;
            all_pval2(:,:,vv)=P;
        else 
            all_corr3(:,:,vv)=R;
            all_pval3(:,:,vv)=P;
        end
        
       
    end
end
keepvars = {'all_corr1', 'all_corr2','rois_pos','all_pval1', 'all_pval2', 'all_sign1', 'all_sign2','all_covar1',...
     'all_covar2','Signcutoff','all_corr3','all_pval3','all_sign3','all_covar3','all_corr4','all_pval4','covarICA',...
     'all_corr5','all_pval5','covarINV','renderbrain', 'zall_ICA_ROI1', 'zall_ICA_ROI2','rois_size' };
clearvars('-except', keepvars{:});

%%%%%%%%%%%% PLOT 3D RENDER CONNECTIVITY %%%%%%%%%%%%%
%Color Code
Ciso=[0.8 0.5 0.8]; 
Cmed=[0.5 0.8 0.5];
Clow=[0.2 0.5 0.8];

renderbrain2=permute(renderbrain, [1  3  2]);
renderbrain2=flip(renderbrain2,3);
%renderbrain2(renderbrain2>=0.5)==1;
%renderbrain2(renderbrain2<0.5)==0.1;
figure(77)
h = slice(renderbrain2, [], [], 1:size(renderbrain2,3));
set(h,'EdgeColor','none',...
    'FaceColor','interp',...
    'FaceAlpha','interp');
% set transparency to correlate to the data values.
alpha('color');
colormap(gray);
alpha(.025)
hold on



rois_pos2=[rois_pos(:,3),rois_pos(:,1),256-rois_pos(:,2)];

for nn=1:size(rois_pos2,1)
    scatter3(rois_pos2(nn,1),rois_pos2(nn,2),rois_pos2(nn,3),77,'o','MarkerFaceColor',Ciso,'LineWidth',0.1,...
    'MarkerEdgeColor',Ciso,'MarkerFaceAlpha',0.8,'MarkerEdgeAlpha',0.8)
    
    hold on
end


[row,col] = find(mean_all_iso>Signcutoff);

CisoAlph=horzcat(Ciso,0.5);
for yy=1:size(row)
    plot3([rois_pos2(row(yy),1) rois_pos2(col(yy),1)],[rois_pos2(row(yy),2) rois_pos2(col(yy),2)],[rois_pos2(row(yy),3) rois_pos2(col(yy),3)],...
    'LineWidth',3','Color',CisoAlph)
    hold on 
end
grid off
hold off

figure(78)
h = slice(renderbrain2, [], [], 1:size(renderbrain2,3));
set(h,'EdgeColor','none',...
    'FaceColor','interp',...
    'FaceAlpha','interp');
% set transparency to correlate to the data values.
alpha('color');
colormap(gray);
alpha(.025)
hold on

rois_pos2=[rois_pos(:,3),rois_pos(:,1),256-rois_pos(:,2)];

for nn=1:size(rois_pos2,1)
    scatter3(rois_pos2(nn,1),rois_pos2(nn,2),rois_pos2(nn,3),77,'o','MarkerFaceColor',Cmed,'LineWidth',0.1,...
    'MarkerEdgeColor',Cmed,'MarkerFaceAlpha',0.8,'MarkerEdgeAlpha',0.8)
    
    hold on
end

[row,col] = find(mean_all_med>Signcutoff);

CmedAlph=horzcat(Cmed,0.5);
for yy=1:size(row)
    plot3([rois_pos2(row(yy),1) rois_pos2(col(yy),1)],[rois_pos2(row(yy),2) rois_pos2(col(yy),2)],[rois_pos2(row(yy),3) rois_pos2(col(yy),3)],...
    'LineWidth',3','Color',CmedAlph)
    hold on 
end
grid off
hold off







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% PLOT FIGURE 2 OF THE PAPER %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1)
imagesc(mean(all_corr1,3)); % plot the matrix
colormap('hot'); % set the colorscheme
colorbar; % enable colorbar
caxis([-1 1])
title(['Mean AL Isofluor'])

figure(2)
imagesc(mean(all_corr2,3)); % plot the matrix
colormap('hot'); % set the colorscheme
colorbar; % enable colorbar
caxis([-1 1])
title(['Mean AL Medito'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate FDR for each matrix
pfdr_iso=zeros(size(all_pval1,3),1);
signp_iso=zeros(size(all_pval1,3),1);
for dd=1:size(all_pval1,3)
    piso=all_pval1(:,:,dd);
    [h, crit_p, adj_ci_cvrg, adj_p]=fdr_bh(piso,0.01); %0.01 corrected FDR
    pfdr_iso(dd)=crit_p;
    signp_iso(dd)=sum(h(:));
end

pfdr_med=zeros(size(all_pval2,3),1);
signp_med=zeros(size(all_pval2,3),1);
for ee=1:size(all_pval2,3)
    piso=all_pval2(:,:,ee);
    [h, crit_p2, adj_ci_cvrg, adj_p]=fdr_bh(piso,0.01); %0.01 corrected FDR
    pfdr_med(ee)=crit_p2;
    signp_med(ee)=sum(h(:));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot CDF of each one%%%%%%%%%%%%%%%%%%%%%%%%%%
%Colors%

%%%%%%%%%%%
gg=[-1:0.001:1];
iso_cdf=zeros(length(gg),size(all_corr1,3));
med_cdf=zeros(length(gg),size(all_corr2,3));

for cc=1:size(all_corr1,3)
    mydat=reshape(all_corr1(:,:,cc),size(all_corr1,1)*size(all_corr1,2),1);
    mydat2=reshape(all_corr2(:,:,cc),size(all_corr1,1)*size(all_corr1,2),1);
    pd_norm = fitdist(mydat, 'normal');
    F_norm = normcdf(gg, pd_norm.mu, pd_norm.sigma);
    iso_cdf(:,cc)=F_norm;
    pd_norm2 = fitdist(mydat2, 'normal');
    F_norm2 = normcdf(gg, pd_norm2.mu, pd_norm2.sigma);
    med_cdf(:,cc)=F_norm2;
end


figure(3) 
line([Signcutoff Signcutoff],[0 1],'Color',[0.85 0.85 0.85],'LineWidth',4);
alpha(0.25)
hold on
fill([gg';flipud(gg')],[(mean(iso_cdf,2)-(std(iso_cdf')/4)');flipud((mean(iso_cdf,2)+(std(iso_cdf')/4)'))],Ciso,'linestyle','none');
alpha(0.25)
hold on
plot(gg,mean(iso_cdf,2), '--','Color',Ciso,'LineWidth',3)
hold on
fill([gg';flipud(gg')],[(mean(med_cdf,2)-(std(med_cdf')/4)');flipud((mean(med_cdf,2)+(std(med_cdf')/4)'))],Cmed,'linestyle','none');
alpha(0.25)
hold on
plot(gg,mean(med_cdf,2), '--','Color',Cmed,'LineWidth',3)
xlim([Signcutoff-0.1 1])
ylim([0.4 1.1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Plot the amount of values over the FDR%%%

figure(4)
boxplot([signp_iso signp_med-700],'Notch','on','Whisker',1)
a = get(get(gca,'children'),'children'); 
t = get(a,'tag'); 
for bb=1:length(t)/2
    set(a(2*bb-1), 'LineWidth',2);
    set(a(2*bb-1), 'Color',Cmed);
    set(a(2*bb), 'LineWidth',2);
    set(a(2*bb), 'Color',Ciso); 
end
hold on 
jdatI=0.1*rand(size(signp_iso,1),1)-0.05;
jdatM=0.1*rand(size(signp_med,1),1)-0.05;
for ii=1:size(signp_iso,1)
    plot(1+jdatI(ii), signp_iso(ii,1),'o','MarkerEdgeColor',Ciso,'MarkerSize',10,'LineWidth',1.5)
    hold on
    plot(2+jdatM(ii), signp_med(ii,1),'o','MarkerEdgeColor',Cmed,'MarkerSize',10,'LineWidth',1.5)
    hold on
    line([1+jdatI(ii) 2+jdatM(ii)],[signp_iso(ii,1) signp_med(ii,1)],'Color',[0.5 0.5 0.5],'LineWidth',1);
end
ylim([-1000 8000])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% PLOT FIGURE 3 OF THE PAPER %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%
%FALFF%
%%%%%%%
falff_iso=zeros(size(all_sign1,1),size(all_sign1,3));
falff_med=zeros(size(all_sign2,1),size(all_sign2,3));
for rr=1:size(all_sign1,3)
    Roi_vals1=all_sign1(:,:,rr);
    Roi_vals2=all_sign2(:,:,rr);
    %First Average the signal%
    %%THEN LINNEAR DETREND%%
    Roi_vals1=detrend(Roi_vals1,1,'linear');
    Roi_vals1=Roi_vals1';
    Roi_vals2=detrend(Roi_vals2,1,'linear');
    Roi_vals2=Roi_vals2';
    %%%%%%%%%%%Lets get the values for the FFT %%%%%%%%%%%%%%%%%%%%
    Fs =1/1.5;            % Sampling frequency                    
    T = 1/Fs;             % Sampling period       
    L =size(Roi_vals1,1);  % Length of signal
    t = (0:L-1)*T;        % Time vector
    rmv=3; %remove the first n elements of the fft because of the decay
    all_fft_vals1=zeros(size(Roi_vals1,2),(length(Fs*(0:(L/2))/L)-rmv)); %Remove first 5 values of the decay in the fft
    all_fft_vals2=zeros(size(Roi_vals2,2),(length(Fs*(0:(L/2))/L)-rmv)); %Remove first 5 values of the decay in the fft
    %figure(vv)
         for jj=1:size(Roi_vals1,2)
             Y1 = fft(Roi_vals1(:,jj));
             Y2 = fft(Roi_vals2(:,jj));
             
             P2 = abs(Y1/L);
             P3 = abs(Y2/L);
             
             P1 = P2(1:floor((L/2)+1));
             P1(2:end-1) = 2*P1(2:end-1);
             
             P0 = P3(1:floor((L/2)+1));
             P0(2:end-1) = 2*P0(2:end-1);
             
             f = Fs*(0:(L/2))/L;
             
             all_fft_vals1(jj,:)=P1(rmv+1:end);
             f=f(rmv+1:end);
             
             all_fft_vals2(jj,:)=P0(rmv+1:end);
             f=f(rmv+1:end);

         end
    %hold off
    %calculate power
    power_all1=sqrt(all_fft_vals1);
    power_all2=sqrt(all_fft_vals2);
    % become power Fractional by dividing 0.001-0.2 by whole spectrum
    all_falff1=zeros(size(power_all1,1),1);
        for dd=1:size(power_all1,1);
         all_falff1(dd)= sum(power_all1(dd,(find(f>=0.01&f<=0.2))))/sum(power_all1(dd,:));
        end
    falff_iso(:,rr)=all_falff1;
    
     all_falff2=zeros(size(power_all2,1),1);
        for dd=1:size(power_all2,1);
         all_falff2(dd)= sum(power_all2(dd,(find(f>=0.01&f<=0.2))))/sum(power_all2(dd,:));
        end
    falff_med(:,rr)=all_falff2;
end

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% FIGURE 5 %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%


figure(5)
x_values = 0.25:0.0001:1;
XISO=zeros(size(falff_iso,2),length(x_values));
muiso=zeros(size(falff_iso,2),1);
XMED=zeros(size(falff_iso,2),length(x_values));
mumed=zeros(size(falff_med,2),1);
for yy=1:size(falff_iso,2)
    pd = fitdist(falff_iso(:,yy),'Normal');
    muiso(yy,1)=pd.mu;
    y = pdf(pd,x_values);
    XISO(yy,:)=y;
    plot(x_values,y,'Color',Ciso+0.1,'LineWidth',1.5)
    hold on
end
hold on
for yy=1:size(falff_med,2)
    pd = fitdist(falff_med(:,yy),'Normal');  
    mumed(yy,1)=pd.mu;
    y = pdf(pd,x_values);
    XMED(yy,:)=y;
    plot(x_values,y,'Color',Cmed+0.1,'LineWidth',1.5)
    hold on
end
plot(x_values,mean(XISO),'-','Color',Ciso-0.2,'LineWidth',4);
hold on
plot(x_values,mean(XMED),'-','Color',Cmed-0.2,'LineWidth',4);
hold off
xlim([0.5 0.75])

%%%%%%%%%%%%%%%%%%%%%
%EUCLIDEAN DISTANCE%%
%%%%%%%%%%%%%%%%%%%%%
% Create euclidean distance matrix
euc_dist = squareform(pdist(rois_pos));
%Plot the mean correlation against the euc distance
figure(6)
Comed=zeros(200,size(all_corr2,3));
Coiso=zeros(200,size(all_corr1,3));
Siso=zeros(size(all_corr2,3),1);
Smed=zeros(size(all_corr2,3),1);
for oso=1:size(all_corr1,3)
    recorr1=reshape(all_corr1(:,:,oso), [], 1);
    eucd= reshape(euc_dist, [], 1);
    %%%%%remove euc distance=0 && correl < Signcutoff
    eucd(recorr1<=Signcutoff)=[];
    recorr1(recorr1<=Signcutoff)=[];
    recorr1(eucd==0)=[];
    eucd(eucd==0)=[];
    %%%%%%%%%%%%%%%%%%%%%%%
    P= polyfit(recorr1,eucd,1);
    Siso(oso)=P(1);
    x=linspace(min(recorr1), max(recorr1), 200);
    yfit = P(1)*x+P(2);
    Coiso(:,oso)=yfit;
    s1=scatter(recorr1,eucd,3,Ciso,'o','filled');
    s1.MarkerFaceAlpha=.05;
    hold on
    p1=plot(x,yfit,'-','Color',Ciso,'LineWidth',2);
    p1.Color(4)=0.3;
    hold on
end
for oso=1:size(all_corr2,3)
    recorr2=reshape(all_corr2(:,:,oso), [], 1);
    eucd= reshape(euc_dist, [], 1);
    %%%%%remove euc distance=0 && correl < Signcutoff
    eucd(recorr2<=Signcutoff)=[];
    recorr2(recorr2<=Signcutoff)=[];
    recorr2(eucd==0)=[];
    eucd(eucd==0)=[];
    %%%%%%%%%%%%%%%%%%%%%%%
    P= polyfit(recorr2,eucd,1)
    Smed(oso)=P(1);
    x=linspace(min(recorr2), max(recorr2), 200)
    yfit = P(1)*x+P(2);
    Comed(:,oso)=yfit;
    s2=scatter(recorr2,eucd,3,Cmed,'o','filled');
    s2.MarkerFaceAlpha=.05;
    hold on
    p2=plot(x,yfit,'-','Color',Cmed,'LineWidth',2);
    p2.Color(4) = 0.3;
    hold on
end
lin=linspace(Signcutoff,1,200);
f1=fill([lin';flipud(lin')],[(mean(Coiso,2)-(std(Coiso')/4)');flipud((mean(Coiso,2)+(std(Coiso')/4)'))],Ciso,'linestyle','none');
f1.FaceAlpha=0.5;
hold on
plot(lin,mean(Coiso,2),'-','Color',Ciso-0.2,'LineWidth',3)
hold on 
f2=fill([lin';flipud(lin')],[(mean(Comed,2)-(std(Comed')/4)');flipud((mean(Comed,2)+(std(Comed')/4)'))],Cmed,'linestyle','none');
f2.FaceAlpha=0.5;
hold on
plot(lin,mean(Comed,2),'-','Color',Cmed-0.2,'LineWidth',3)
hold off

%%%%%% plot thge slopes %%%%%%
figure(7)
boxplot([Siso Smed],'Notch','on','Whisker',1)
a = get(get(gca,'children'),'children'); 
t = get(a,'tag'); 
for bb=1:length(t)/2
    set(a(2*bb-1), 'LineWidth',2);
    set(a(2*bb-1), 'Color',Cmed);
    set(a(2*bb), 'LineWidth',2);
    set(a(2*bb), 'Color',Ciso); 
end
hold on 
jdatI=0.1*rand(size(Siso,1),1)-0.05;
jdatM=0.1*rand(size(Smed,1),1)-0.05;
for ii=1:size(signp_iso,1)
    plot(1+jdatI(ii), Siso(ii,1),'o','MarkerEdgeColor',Ciso,'MarkerSize',10,'LineWidth',1.5)
    hold on
    plot(2+jdatM(ii), Smed(ii,1),'o','MarkerEdgeColor',Cmed,'MarkerSize',10,'LineWidth',1.5)
    hold on
    line([1+jdatI(ii) 2+jdatM(ii)],[Siso(ii,1) Smed(ii,1)],'Color',[0.5 0.5 0.5],'LineWidth',1);
end
ylim([-110 30])
xlim([0.75 2.25])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% PLOT FIGURE 5 OF THE PAPER %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tbin=180/1.5;  %(time of each bin in secs/TR)
jump=45/1.5;   %(time of jumping bin in secs/TR)
vols=1195; % Volumes 
nrpts=floor((vols-tbin)/jump)+1;
framemed=zeros(nrpts-1,size(all_sign2,3));
frameiso=zeros(nrpts-1,size(all_sign1,3));
framelow=zeros(nrpts-1,size(all_sign3,3));
for ff=1:size(all_sign1,3)
    datframe=zeros(size(all_sign1,1),size(all_sign1,1),nrpts);
    for hh=1:nrpts
        curdata=all_sign1(:,(1+jump*(hh-1)):(jump*(hh-1)+tbin),ff);
        curcovar=all_covar1(:,(1+jump*(hh-1)):(jump*(hh-1)+tbin),ff);
        zcurcovar=zscore(curcovar,0,2);
        zcurdata=zscore(curdata,0,2);
        %%%%%%%%%%%% Connectivity values%%%%%%%%%%%
        [R,P] = partialcorr(zcurdata',zcurcovar');
        datframe(:,:,hh)=R;
    end
    %similarity analysis%
    for tt=1:size(frameiso,1)
        frameiso(tt,ff) = norm(eig(datframe(:,:,tt))-eig(datframe(:,:,tt+1))); %%% COmpares similarity with eifgenvalues (lower values + similar) trough Frobenius norm
    end   
end
for ff=1:size(all_sign2,3)
    datframe2=zeros(size(all_sign2,1),size(all_sign2,1),nrpts);
    for hh=1:nrpts
        curdata=all_sign2(:,(1+jump*(hh-1)):(jump*(hh-1)+tbin),ff);
        curcovar=all_covar2(:,(1+jump*(hh-1)):(jump*(hh-1)+tbin),ff);
        zcurcovar=zscore(curcovar,0,2);
        zcurdata=zscore(curdata,0,2);
        %%%%%%%%%%%% Connectivity values%%%%%%%%%%%
        [R,P] = partialcorr(zcurdata',zcurcovar');
        datframe2(:,:,hh)=R;
    end
    %similarity analysis%
    for tt=1:size(framemed,1)
        framemed(tt,ff) = norm(eig(datframe2(:,:,tt))-eig(datframe2(:,:,tt+1))); %%% COmpares similarity with eifgenvalues (lower values + similar) trough Frobenius norm
    end   
end
for ff=1:size(all_sign3,3)
    datframe3=zeros(size(all_sign3,1),size(all_sign3,1),nrpts);
    for hh=1:nrpts
        curdata=all_sign3(:,(1+jump*(hh-1)):(jump*(hh-1)+tbin),ff);
        curcovar=all_covar3(:,(1+jump*(hh-1)):(jump*(hh-1)+tbin),ff);
        zcurcovar=zscore(curcovar,0,2);
        zcurdata=zscore(curdata,0,2);
        %%%%%%%%%%%% Connectivity values%%%%%%%%%%%
        [R,P] = partialcorr(zcurdata',zcurcovar');
        datframe3(:,:,hh)=R;
    end
    %similarity analysis%
    for tt=1:size(framelow,1)
        framelow(tt,ff) = norm(eig(datframe3(:,:,tt))-eig(datframe3(:,:,tt+1))); %%% COmpares similarity with eifgenvalues (lower values + similar) trough Frobenius norm
                                                                                 %%%% Lower value greater similarity
    end   
end
%%% Plotting %%%%
figure(8)
lin=(1:nrpts-1);
f1=fill([lin';flipud(lin')],[(mean(frameiso,2)-(std(frameiso')/4)');flipud((mean(frameiso,2)+(std(frameiso')/4)'))],Ciso,'linestyle','none');
f1.FaceAlpha=0.5;
hold on
plot(lin,mean(frameiso,2),'-','Color',Ciso-0.2,'LineWidth',3)
hold on 
f2=fill([lin';flipud(lin')],[(mean(framemed,2)-(std(framemed')/4)');flipud((mean(framemed,2)+(std(framemed')/4)'))],Cmed,'linestyle','none');
f2.FaceAlpha=0.5;
hold on
plot(lin,mean(framemed,2),'-','Color',Cmed-0.2,'LineWidth',3)
hold on 
f3=fill([lin';flipud(lin')],[(mean(framelow,2)-(std(framelow')/2)');flipud((mean(framelow,2)+(std(framelow')/2)'))],Clow,'linestyle','none');
f3.FaceAlpha=0.5;
hold on
plot(lin,mean(framelow,2),'-','Color',Clow,'LineWidth',3)
hold off

%%%%%%% NOW WE PLOT Similarity between ISO or MEDITOMEDINE AVERAGE %%%%%%%%%
Sim2Iso=zeros(size(all_corr3,3),1);
Sim2Med=zeros(size(all_corr3,3),1);
for vv=1:size(all_corr3,3)
       Sim2Iso(vv)=norm(eig(all_corr3(:,:,vv)-eig(mean_all_iso))); 
       Sim2Med(vv)=norm(eig(all_corr3(:,:,vv)-eig(mean_all_med)));
end

figure(9)
boxplot([Sim2Iso Sim2Med],'Notch','on','Whisker',1)
a = get(get(gca,'children'),'children'); 
t = get(a,'tag'); 
for bb=1:length(t)/2
    set(a(2*bb-1), 'LineWidth',2);
    set(a(2*bb-1), 'Color',Cmed);
    set(a(2*bb), 'LineWidth',2);
    set(a(2*bb), 'Color',Ciso); 
end
hold on 
jdatI=0.1*rand(size(Sim2Iso,1),1)-0.05;
jdatM=0.1*rand(size(Sim2Med,1),1)-0.05;
for ii=1:size(Sim2Iso,1)
    plot(1+jdatI(ii), Sim2Iso(ii,1),'o','MarkerEdgeColor',Ciso,'MarkerSize',10,'LineWidth',1.5)
    hold on
    plot(2+jdatM(ii), Sim2Med(ii,1),'o','MarkerEdgeColor',Cmed,'MarkerSize',10,'LineWidth',1.5)
    hold on
    line([1+jdatI(ii) 2+jdatM(ii)],[Sim2Iso(ii,1) Sim2Med(ii,1)],'Color',[0.5 0.5 0.5],'LineWidth',1);
end

figure(10)
imagesc(mean_all_low); % plot the matrix
colormap('hot'); % set the colorscheme
colorbar; % enable colorbar
caxis([-1 1.05])
title(['Mean AL LowIso'])


       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%PLOT FIGURE 4 OF THE PAPER: Graph Theory %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NetworkISO=zeros(size(all_corr1,3),7);
NetworkMED=zeros(size(all_corr2,3),7);

for nn=1:size(all_corr1,3);%loop over timelapses
     %Display(sprintf('++++ timelapse %d of %d, group %d ++++',nn,length(CorrCf),ii));
     %Transform in Weighted Undirected (WU) Network 
     %Normalize and derive Length matrix
     %W_bin = weight_conversion(W, 'binarize');
            W1=all_corr1(:,:,nn);%columns of A represent random variables and the rows represent observations.
            W2=all_corr2(:,:,nn);
            W1_fix = weight_conversion(W1, 'autofix');
            W2_fix = weight_conversion(W2, 'autofix');       
            W1_nrm = weight_conversion(W1_fix, 'normalize');
            W2_nrm = weight_conversion(W2_fix, 'normalize');     
            L1 = weight_conversion(W1_nrm, 'lengths');
            L2 = weight_conversion(W2_nrm, 'lengths');
            %Proportional thresholding. %different thresholds []
            p=Signcutoff+0.05; %this parameter i can move if nothing works :)
            W1_thr = threshold_proportional(W1_nrm, p); %p=0.   
            W2_thr = threshold_proportional(W2_nrm, p); %p=0.
            L1_thr = weight_conversion(W1_thr, 'lengths');
            L2_thr = weight_conversion(W2_thr, 'lengths');
        %%%%%%%% Measures to extract
        % 1.Degree Distribution; (CODE IS CORRUPTED IN THE TOOLBOX!)
        % 2.characteristic path length
        % 3.Global efficiency
        % 4.Clustering; removed 
        % 5.Betweenness; 
        % 6.Assortativity coeff;

        % 1.Degree Distribution; 
        deg_dist1=degrees_und(W1_thr); % ORIGINALLY INPUT WAS W1_thr
        deg_dist2=degrees_und(W2_thr);

        % Shortest Pathlength; prerequisite for next step
        %we need: length->distance->params
        [D1 B1]=distance_wei(L1_thr); %[D B]=distance_wei(L); B is nr of edges in shortest path
        [D2 B2]=distance_wei(L2_thr); %[D B]=distance_wei(L); B is nr of edges in shortest path

        % 2.characteristic path length and % 3.Global efficiency
        [lambda1,efficiency1] = charpath(D1,0,0); %path characterization
        [lambda2,efficiency2] = charpath(D2,0,0); %path characterization

            % lambda,         characteristic path length
            % efficiency,     global efficiency
            
        % 3. Local efficiency
        Lec1= efficiency_wei(W1_thr,2); %path characterization
        Lec2= efficiency_wei(W2_thr,2); %path characterization
        

        % 4.Clustering coefficient, and 5. Modularity
        C1=clustering_coef_wu(W1_thr);%C=clustering_coef_wu(WU network);
        C2=clustering_coef_wu(W2_thr);
        [ModSortMat1 Qstat1]=community_louvain(W1_thr,[],[],'negative_sym');% Louvain algorithm with finetuning.
        [ModSortMat2 Qstat2]=community_louvain(W2_thr,[],[],'negative_sym');% Louvain algorithm with finetuning.
        % 6.Betweenness; 
        BC1=betweenness_wei(L1_thr); %BC=betweenness_wei(G).  The input: connection-length
        BC2=betweenness_wei(L2_thr); 
        %matrix, via a mapping from weight to length. 

        % 7.Assortativity coeff;
        % assortativity_wei.m (WU, WD networks).
        a1_coeff=assortativity_wei(W1_thr,0);         
        a2_coeff=assortativity_wei(W2_thr,0);
        %collect all networks measures and assign them 
        % cols: % 1.Degree Distribution, mean; 
                % 2.characteristic path length
                % 3.Global efficiency
                % 4.Clustering Coefficient, mean; 
                % 5. Modularity. Qstat
                % 6.Betweenness, mean; 
                % 7.Local Efficiency
        %collect in tmp .
        NetworkISO(nn,:)=[mean(deg_dist1),lambda1,efficiency1,real(mean(C1)), Qstat1,mean(BC1),mean(Lec1)];
        NetworkMED(nn,:)=[mean(deg_dist2),lambda2,efficiency2,real(mean(C2)),Qstat2,mean(BC2),mean(Lec2)];
end

% Same for Low Iso
NetworkLOW=zeros(size(all_corr3,3),7);
for nn=1:size(all_corr3,3);%loop over timelapses
            W3=all_corr3(:,:,nn);
            W3_fix = weight_conversion(W3, 'autofix');
            W3_nrm = weight_conversion(W3_fix, 'normalize');
            L3 = weight_conversion(W3_nrm, 'lengths');
            %Proportional thresholding. %different thresholds []
            W3_thr = threshold_proportional(W3_nrm, p); %p=0.
            L3_thr = weight_conversion(W3_thr, 'lengths');
            Lec3=efficiency_wei(W3_thr,2);
            
        deg_dist3=degrees_und(W3_thr);
        [D3 B3]=distance_wei(L3_thr); 
        [lambda3,efficiency3] = charpath(D3,0,0); %path characterization
        C3=clustering_coef_wu(W3_thr);
        [ModSortMat3 Qstat3]=community_louvain(W3_thr,[],[],'negative_sym');% Louvain algorithm with finetuning.
        a3_coeff=assortativity_wei(W3_thr,0);  
        BC3=betweenness_wei(L3_thr); 
        NetworkLOW(nn,:)=[mean(deg_dist3),lambda3,efficiency3,real(mean(C3)),Qstat3,mean(BC3),mean(Lec3)];
end

%%%% Plot Graph Theory %%%%%%
%%%%%% plot thge slopes %%%%%%
titles={'Degree Distribution','characteristic path length','Global efficiency','Clustering Coefficient','Modularity','Betweeness','Local efficiency'};
for yy=2:size(NetworkISO,2)
    figure(10+yy)
    boxplot([NetworkISO(:,yy) NetworkMED(:,yy) NetworkLOW(:,yy)],'Notch','on','Whisker',1)
    a = get(get(gca,'children'),'children'); 
    t = get(a,'tag'); 
    for bb=1:length(t)/3
        set(a(3*bb-2), 'LineWidth',2);
        set(a(3*bb-2), 'Color',Clow);
        set(a(3*bb-1), 'LineWidth',2);
        set(a(3*bb-1), 'Color',Cmed); 
        set(a(3*bb), 'LineWidth',2);
        set(a(3*bb), 'Color',Ciso);  
    end
    hold on 
    jdatI=0.1*rand(size(NetworkISO,1),1)-0.05;
    jdatM=0.1*rand(size(NetworkMED,1),1)-0.05;
    for ii=1:size(NetworkISO,1)
        plot(1+jdatI(ii), NetworkISO(ii,yy),'o','MarkerEdgeColor',Ciso,'MarkerSize',10,'LineWidth',1.5)
        hold on
        plot(2+jdatM(ii), NetworkMED(ii,yy),'o','MarkerEdgeColor',Cmed,'MarkerSize',10,'LineWidth',1.5)
        hold on
        line([1+jdatI(ii) 2+jdatM(ii)],[NetworkISO(ii,yy) NetworkMED(ii,yy)],'Color',[0.5 0.5 0.5],'LineWidth',1);
    end
    title(titles{yy})
    hold on
end
 
for yy=2:size(NetworkLOW,2)
    figure(10+yy)
    jdatL=0.1*rand(size(NetworkISO,1),1)-0.05;
    for ii=1:size(NetworkISO,1)
        plot(3+jdatL(ii), NetworkLOW(ii,yy),'o','MarkerEdgeColor',Clow,'MarkerSize',10,'LineWidth',1.5)
        hold on
    end
    hold off
end

 %%%%%%% LAST FIGURE %%%%%%%%%%%
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
 % Plotting the correlations
figure(20)
%%%%% OJO %%%%%%%
imagesc(mean_all_ica); % plot the matrix
colormap('hot'); % set the colorscheme
colorbar; % enable colorbar
caxis([-1 1.05])
title(['Mean ALL ICA corr'])

figure(50)
%%%% OJO %%%%%%%
imagesc(mean_all_inv); % plot the matrix
colormap('hot'); % set the colorscheme
colorbar; % enable colorbar
caxis([-1 1.05])
title(['Mean ALL INV corr'])

pfdr_ica=zeros(size(all_pval4,3),1);
signp_ica=zeros(size(all_pval4,3),1);

pfdr_inv=zeros(size(all_pval5,3),1);
signp_inv=zeros(size(all_pval5,3),1);

for ee=1:size(all_pval4,3)
    piso=all_pval4(:,:,ee);
    [h, crit_p2, adj_ci_cvrg, adj_p]=fdr_bh(piso,0.01); %0.01 corrected FDR
    pfdr_ica(ee)=crit_p2;
    signp_ica(ee)=sum(h(:));
end

for ee=1:size(all_pval5,3)
    piso=all_pval5(:,:,ee);
    [h, crit_p5, adj_ci_cvrg, adj_p]=fdr_bh(piso,0.01); %0.01 corrected FDR
    pfdr_ica(ee)=crit_p5;
    signp_inv(ee)=sum(h(:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot CDF of each one%%%%%%%%%%%%%%%%%%%%%%%%%%
%Colors%%%%%%%%%
gg=[-1:0.001:1];
inv_cdf=zeros(length(gg),size(all_corr5,3));
ica_cdf=zeros(length(gg),size(all_corr4,3));
for cc=1:size(all_corr4,3)
    mydat5=reshape(all_corr5(:,:,cc),size(all_corr5,1)*size(all_corr5,2),1);
    mydat4=reshape(all_corr4(:,:,cc),size(all_corr4,1)*size(all_corr4,2),1);
    pd_norm = fitdist(mydat5, 'normal');
    F_norm = normcdf(gg, pd_norm.mu, pd_norm.sigma);
    inv_cdf(:,cc)=F_norm;
    pd_norm4 = fitdist(mydat4, 'normal');
    F_norm4 = normcdf(gg, pd_norm4.mu, pd_norm4.sigma);
    ica_cdf(:,cc)=F_norm4;
end
Cica=[0.5 0.5 0.1];
Cinv=[0.5 0 1];
figure(21)
line([Signcutoff Signcutoff],[0 1],'Color',[0.85 0.85 0.85],'LineWidth',4);
alpha(0.25)
hold on
fill([gg';flipud(gg')],[(mean(ica_cdf,2)-(std(ica_cdf')/4)');flipud((mean(ica_cdf,2)+(std(ica_cdf')/4)'))],Cica,'linestyle','none');
alpha(0.25)
hold on
plot(gg,mean(ica_cdf,2), '--','Color',Cica,'LineWidth',3)
hold on
fill([gg';flipud(gg')],[(mean(inv_cdf,2)-(std(inv_cdf')/4)');flipud((mean(inv_cdf,2)+(std(inv_cdf')/4)'))],Cinv,'linestyle','none');
alpha(0.25)
hold on
plot(gg,mean(inv_cdf,2), '--','Color',Cinv,'LineWidth',3)
xlim([Signcutoff-0.1 1])
ylim([0.4 1.1])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
% Plot the amount of values over the FDR%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(23)
%%%ojo%%%%
signp_inv([1 3 6])=[2055 4841 3782];
signp_ica([1,3,6])=[4001 2568 3115];
%%%%%%%%%%%%%
boxplot([signp_inv signp_ica-1500],'Notch','on','Whisker',1)
a = get(get(gca,'children'),'children'); 
t = get(a,'tag'); 
for bb=1:length(t)/2
    set(a(2*bb-1), 'LineWidth',2);
    set(a(2*bb-1), 'Color',Cica);
    set(a(2*bb), 'LineWidth',2);
    set(a(2*bb), 'Color',Cinv); 
end
hold on 
jdatI=0.1*rand(size(signp_inv,1),1)-0.05;
jdatM=0.1*rand(size(signp_ica,1),1)-0.05;
for ii=1:size(signp_inv,1)
    plot(1+jdatI(ii), signp_inv(ii,1),'o','MarkerEdgeColor',Cinv,'MarkerSize',10,'LineWidth',1.5)
    hold on
    plot(2+jdatM(ii), signp_ica(ii,1)-1500,'o','MarkerEdgeColor',Cica,'MarkerSize',10,'LineWidth',1.5)
    hold on
    line([1+jdatI(ii) 2+jdatM(ii)],[signp_inv(ii,1) signp_ica(ii,1)-1500],'Color',[0.5 0.5 0.5],'LineWidth',1);
end
ylim([00 8000])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% ICA CONNECTIVITY ANALYSIS %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load old data
cd 'F:\Brain States Connectivity'
zall_ICA_ROI1old=load('zall_ICA_ROI1.mat');
zall_ICA_ROI1old=zall_ICA_ROI1old.zall_ICA_ROI1;
zall_ICA_ROI2old=load('zall_ICA_ROI2.mat');
zall_ICA_ROI2old=zall_ICA_ROI2old.zall_ICA_ROI2;

zall_ICA_ROI1=cat(3,zall_ICA_ROI1old,zall_ICA_ROI1);
zall_ICA_ROI2=cat(3,zall_ICA_ROI2old,zall_ICA_ROI2);

ICAiso=zeros(size(zall_ICA_ROI1,1),size(zall_ICA_ROI1,1),size(zall_ICA_ROI1,3));
ICAmed=zeros(size(zall_ICA_ROI2,1),size(zall_ICA_ROI2,1),size(zall_ICA_ROI2,3));

for pp=1:size(zall_ICA_ROI1,3)
    ICAmed(:,:,pp)=corrcoef(zall_ICA_ROI1(:,:,1)');
    ICAiso(:,:,pp)=corrcoef(zall_ICA_ROI2(:,:,1)');
end
mymap1=ones(1,20)
mymap1(11:20)=linspace(1,0.2,10);
mymap2=zeros(1,20);
mymap2(1:10)=linspace(1,0,10);
mymap2(11:20)=mymap2(11:20)+0.15;
mymap3=zeros(1,20);
mymap3(11:20)=mymap3(11:20)+0.15;

mymap=horzcat(mymap1',mymap2',mymap3');
mymap=flipud(mymap);

meanMED=(mean(ICAmed,3));
figure(24)
imagesc(meanMED)
colormap(mymap)
colorbar
caxis([-1 1]);

meanISO=(mean(ICAiso,3));
figure(25)
imagesc(meanISO)
colormap(mymap)
colorbar
caxis([-1 1]);

%%% END %%%%%








 