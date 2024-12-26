%% relevant paths & data files

%run script to get behav and eyetracking data file pairs
map_asc2mat_AttDeployment

%use only participants with both kinds of data
%MasterTable_Files=MasterTable_Files(MasterTable_Files.overall,:);


% specify data directory
behavData_path='';

% specify directory to save figures
path2saveFigs='';


%%
nParticip=sum(MasterTable_Files.overall);
nTrials=20;

partic_IDs=cell(1,1);
Stim_Mat=nan(nTrials,nParticip);

for j=1:nParticip
    
    load([behavData_path filesep MasterTable_Files.mat{j}])
    Stim=BehavData.vars.Stim_seq;

    Stim_Mat(:,j)=Stim;
    partic_IDs{j,1}=MasterTable_Files.asc{j}(1:end-4);  %remove last 4 chars ('.asc')
    
end

% load large file with eye movs data
if exist('eyeMovs_byImg_UnpleasantImgs','var')
else
     load('eyeMovs_AttDepl_byImg_UnpleasantImgs.mat')
end

% load single trial data
dataFolder='';
participID='GX_38_79';
load([dataFolder filesep 'eyeMovs_AttDepl_byParticip_' participID '.mat'])

%use special folder containing copies of original images for trial #14
% these images have the circle highlighted for a better visualization
specialFolder='images_For_Figure3';
trials=fieldnames(eyeMovs.data);
trialNum=14;

imgs=eyeMovs.data.(trials{trialNum}).imgs(:,1);
GazeXYtrial=eyeMovs.data.(trials{trialNum}).GazeXY;

for thisImg=1:numel(imgs)
    eval(['img' num2str(thisImg) '=imread([specialFolder filesep imgs{thisImg}]);'])
end

imgSize=size(img1);

%load eye movs summary stats
if exist('eyeMovs_Stats','var')
else
    load('eyeMovs_AttDepl_Stats.mat', 'eyeMovs_Stats')

end


%% %%%%% CALCULATE STATS (DWELL TIME & CIRCLE BOUNDARY CROSSINGS)

MeanDwellTime=eyeMovs_Stats.data.Mean_TimeSpentInCircle;
MeanCircleBoundCrossings=eyeMovs_Stats.data.Mean_CircleBoundCrossing;
isDataOK=eyeMovs_Stats.data.isDataOK;
%remove bad data
for thisRow=1:size(eyeMovs_Stats.data.isDataOK,1)
    for thisCol=1:size(eyeMovs_Stats.data.isDataOK,2)
        if isDataOK(thisRow,thisCol)
        else
            MeanDwellTime(thisRow,thisCol)=nan;
            MeanCircleBoundCrossings(thisRow,thisCol)=nan;
        end
    end
end

%build logicals to parse trials by type and attentional condition
Logi_neutralImgs=true(size(MeanDwellTime));
Logi_FocusArous=true(size(MeanDwellTime));
Logi_FocusNonArous=true(size(MeanDwellTime));
for thisParticip=1:size(Stim_Mat,2)
    stim=Stim_Mat(:,thisParticip);
    
    Logi_neutralImgs(:,thisParticip)=stim<3;
    Logi_FocusArous(:,thisParticip)=stim==5;
    Logi_FocusNonArous(:,thisParticip)= (stim==4 | stim==2);
end
    


%calculate dwell time parsed by img type and attentional condit
MeDwTi_neut=nan(nParticip,1);
MeDwTi_unpl=nan(nParticip,1);

MeDwTi_focA=nan(nParticip,1);
MeDwTi_focNA=nan(nParticip,1);

%calculate circle boundary crossings parsed by img type and attentional condit
MeCBC_neut=nan(nParticip,1);
MeCBC_unpl=nan(nParticip,1);
MeCBC_focA=nan(nParticip,1);
MeCBC_focNA=nan(nParticip,1);



for j=1:nParticip
    %DWELL TIME
    datavector=MeanDwellTime(:,j);
    datavector_n=datavector(Logi_neutralImgs(:,j));
    datavector_u=datavector(~Logi_neutralImgs(:,j));
    
    MeDwTi_neut(j)=mean(datavector_n,'omitnan');
    MeDwTi_unpl(j)=mean(datavector_u,'omitnan');
        
    datavector_focA=datavector(Logi_FocusArous(:,j));
    datavector_focNA=datavector(Logi_FocusNonArous(:,j));
    
    MeDwTi_focA(j)=mean(datavector_focA,'omitnan');
    MeDwTi_focNA(j)=mean(datavector_focNA,'omitnan');
    
    %CIRCLE BOUNDARY CROSSINGS
    datavector=MeanCircleBoundCrossings(:,j);
    datavector_n=datavector(Logi_neutralImgs(:,j));
    datavector_u=datavector(~Logi_neutralImgs(:,j));
    
    MeCBC_neut(j)=mean(datavector_n,'omitnan');
    MeCBC_unpl(j)=mean(datavector_u,'omitnan');
        
    datavector_focA=datavector(Logi_FocusArous(:,j));
    datavector_focNA=datavector(Logi_FocusNonArous(:,j));
    
    MeCBC_focA(j)=mean(datavector_focA,'omitnan');
    MeCBC_focNA(j)=mean(datavector_focNA,'omitnan');
    
    
end


%count how many subjects were eliminated due to artifactual data

temp=logical(sum(isDataOK));

nRemoved=sum(temp==0);
logi_Removed=temp==0;


disp('done!')

%% FIGURE

fsize=12;
fname='VErdana';


fig1=figure;

%%%%%%%%%%%% PLOT SINGLE-TRIAL

w=0.16;
h=0.2;
L=0.166;
bott1=0.75;

%positions of five images
posits1=[L bott1 w h;
    L+1.03*w bott1 w h;
    L+2*(1.03*w) bott1 w h;
    L+3*(1.03*w) bott1 w h;
    L+4*(1.03*w) bott1 w h];

circleColor=[0.4 0.4 1];

nameGazesXY=fieldnames(GazeXYtrial);

for thisImg=1:numel(imgs)
    ax=axes('parent',fig1,'position', posits1(thisImg,:));
    hold(ax,'on')
    
    eval(['imagesc(flipud(img' num2str(thisImg) '))'])
    axis xy
    axis off
    
    Gaze_xy=GazeXYtrial.(nameGazesXY{thisImg});
    hold on
        
    cm=colormap('parula');
    
    plot(Gaze_xy(:,1),(Gaze_xy(:,2).*-1)+imgSize(1),'linewidth',1.6,'color','y')
    plot(Gaze_xy(1,1),(Gaze_xy(1,2).*-1)+imgSize(1),'marker','.','markersize',24,'linewidth',3,'color',cm(15,:))
    plot(Gaze_xy(end,1),(Gaze_xy(end,2).*-1)+imgSize(1),'marker','.','markersize',24,'linewidth',3,'color',cm(end,:))
    
%     text(Gaze_xy(1,1),(Gaze_xy(1,2).*-1)+imgSize(1),'S','color',[1 1 1],'fontsize',fsize-2)
%     text(Gaze_xy(end,1),(Gaze_xy(end,2).*-1)+imgSize(1),'E','fontsize',fsize-2,'color','y')
    
    
    xlim([1 imgSize(2)])
    ylim([1 imgSize(1)])


end

annotation(fig1,'textbox',...
    [0.16 0.92 0.38 0.07],'fontname',fname,...
    'String',{'Single trial example, one participant'},'edgecolor','none',...
    'FitBoxToText','off','horizontalalignment','left','fontsize',fsize-2);

annotation(fig1,'textbox',...
    [0.002 0.83 0.18 0.07],...
    'String',{'Stimulus:', 'Unpleasant'},'edgecolor','none','fontname',fname,...
    'FitBoxToText','off','horizontalalignment','left','fontsize',fsize-2);
annotation(fig1,'textbox',...
    [0.002 0.75 0.18 0.07],...
    'String',{'Focus:', 'Non-arousing'},'edgecolor','none','fontname',fname,...
    'FitBoxToText','off','horizontalalignment','left','fontsize',fsize-2);


%%%%%%%%%% PLOT ALL-PARTICIPANTS HEAT MAPS

imNumber='6370';
suffixes={'.jpg','-calm.jpg','-arouse.jpg'};

w=0.27;
h=0.3;
bott=0.37;
L=0.16;
%positions of three images
posits3=[L bott w h;
    L+1.03*w bott w h;
    L+2*(1.03*w) bott w h];


for thisFig=1:3
    image=[imNumber suffixes{thisFig}];
    imageMat=imread([specialFolder filesep image ]);

    if thisFig==1
        PVM=sum(eyeMovs_byImg_UnpleasantImgs.data.(['img_' imNumber]).PixVisitedMat_folder_03_UNPLEAS_NO_FOCUS,3);
    elseif thisFig==2
        PVM=sum(eyeMovs_byImg_UnpleasantImgs.data.(['img_' imNumber]).PixVisitedMat_folder_04_UNPLEAS_FOCUS_CALM,3);
    elseif thisFig==3
        PVM=sum(eyeMovs_byImg_UnpleasantImgs.data.(['img_' imNumber]).PixVisitedMat_folder_05_UNPLEAS_FOCUS_AROUS,3);        
    end
    
    PVM(44,:)=0;
    
    H=fspecial('gaussian',30,5);
    M=filter2(H,PVM);
    % convert small numbers into nans to remove background
        for r=1:size(M,1)
            for c=1:size(M,2)
                if M(r,c)<0.1
                    M(r,c)=nan;
                end
            end
        end
    
    ax=axes('parent',fig1,'position',posits3(thisFig,:));
    hold(ax,'on')
    imagesc(flipud(imageMat));
    [~,c]=contourf(flipud(M),3);
    c.LineColor='k';
    c.LineWidth=0.5;
    colormap autumn;
    axis off
    
end

annotation(fig1,'textbox',...
    [0.025 0.52 0.1 0.07],...
    'String',{'All','participants'},'edgecolor','none','fontname',fname,...
    'FitBoxToText','off','horizontalalignment','center','fontsize',fsize-2);

bott=0.615;
w=0.22;
h=0.1;

annotation(fig1,'textbox',...
    [0.11 bott w h],...
    'String',{'Focus-free'},'edgecolor','none','fontname',fname,...
    'FitBoxToText','off','horizontalalignment','center','fontsize',fsize-2);
annotation(fig1,'textbox',...
    [0.42 bott w h],...
    'String',{'Focus Non-Arousing'},'edgecolor','none','fontname',fname,...
    'FitBoxToText','off','horizontalalignment','center','fontsize',fsize-2);
annotation(fig1,'textbox',...
    [0.68 bott w h],...
    'String',{'Focus Arousing'},'edgecolor','none','fontname',fname,...
    'FitBoxToText','off','horizontalalignment','center','fontsize',fsize-2);

    

%%%%%%%% PLOT STATS

ms=22;

posits4=[0.17 0.108 0.6 0.21;
    0.88 0.13 0.1 0.21];


ax=axes('parent',fig1,'tickdir','out','position',posits4(1,:),...
    'fontname',fname,'fontsize',fsize,'xtick',10:10:50);
hold(ax,'on')

X=1:size(MeanDwellTime,2);
Y=100*mean(MeanDwellTime,'omitnan')./4;
E=100*std(MeanDwellTime,'omitnan')./4;

% Y=mean(MeanCircleBoundCrossings,'omitnan');
% E=std(MeanCircleBoundCrossings,'omitnan');

errorbar(X,Y,E,'linestyle','none','marker','.','markersize',ms,...
    'linewidth',1.8,'color',[0 0 0])

ylabel({'Dwell time in', 'focus circle (%)'})
xlabel('Participant')

ylim([-3.5 100])
xlim([0.5 55.5])

ax=axes('parent',fig1,'tickdir','out','position',posits4(2,:),'TickLength', [0.05 0.035],...
    'xtick',[0 1.5],'xticklabel',{'Arous-NonArous','Neutr-Unpleas'},'fontname',fname,'fontsize',fsize);
hold(ax,'on')

xtickangle(20)


Y1=100*(MeDwTi_neut-MeDwTi_unpl)./4;
X1=randn(size(Y1))./15;

Y2=100*(MeDwTi_focA-MeDwTi_focNA)./4;




plot([-0.5 2],[0 0],'k')

plot([-0.35 0.35]+1.5,[mean(Y1,'omitnan') mean(Y1,'omitnan')],'color','k','linewidth',4)
plot([-0.35 0.35],[mean(Y2,'omitnan') mean(Y2,'omitnan')],'color','k','linewidth',4)

plot(X1+1.5,Y1,'linestyle','none','marker','.','markersize',ms-9,'color',[0.55 0.55 0.55])
plot(X1,Y2,'linestyle','none','marker','.','markersize',ms-9,'color',[0.55 0.55 0.55])

%text(-0.2,20,'ns','fontsize',fsize,'fontname',fname)


[~,p,ci,st]=ttest(Y1);
if p<0.05
    text(1.3,21,'*','fontsize',fsize-2,'fontname',fname)
end

[~,p,ci,st]=ttest(Y2);

xlim([-0.5 2])
ylim([-20 21])

ylabel({'Dwell time','difference (%)'})

lettW=0.05;
lettH=0.05;
annotation(fig1,'textbox',[0.002 0.97 lettW lettH],'string','A','edgecolor','none',...
    'fontweight','bold','fontsize',fsize+14,'fontname',fname)

annotation(fig1,'textbox',[0.002 0.68 lettW lettH],'string','B','edgecolor','none',...
    'fontweight','bold','fontsize',fsize+14,'fontname',fname)

annotation(fig1,'textbox',[0.002 0.33 lettW lettH],'string','C','edgecolor','none',...
    'fontweight','bold','fontsize',fsize+14,'fontname',fname)



print(fig1,[path2saveFigs filesep 'FIG_03'],'-dpng','-r450')
%% stats

disp('%%%%%%%%%%%%%%%')
disp('DWELL TIME (%)')
MDT_PERCENT=mean(100*MeanDwellTime/4,1,'omitnan');
nEffective=sum(~isnan(MDT_PERCENT));
disp(['Median = ' num2str(median(MDT_PERCENT,'omitnan'))])
disp(['Mean = ' num2str(mean(MDT_PERCENT,'omitnan'))])
disp(['SD = ' num2str(std(MDT_PERCENT,'omitnan'))])
disp(['n = ' num2str(nEffective)])

disp('%%%%%%%%%%%%%%%')
disp('DWELL TIME DIFFERENCES')
disp('1) Arousing - NonArousing')
disp(['Arousing: Median = ' num2str(100*median(MeDwTi_focA,'omitnan')/4) '%, Mean = ' num2str(100*mean(MeDwTi_focA,'omitnan')/4) '% , SD = ' num2str(100*std(MeDwTi_focA,'omitnan')/4) '%'])
disp(['Non-Arousing: Median = ' num2str(100*median(MeDwTi_focNA,'omitnan')/4) '%, Mean = ' num2str(100*mean(MeDwTi_focNA,'omitnan')/4) '%, SD = ' num2str(100*std(MeDwTi_focNA,'omitnan')/4) '%'])

[pval, stat_no_perm, dof]=PermutePairedTest(MeDwTi_focA, MeDwTi_focNA, 1500, 'tstat', 'both');
stat_no_perm=round(stat_no_perm*1000)/1000;
pval=round(pval*1000)/1000;
disp(['t(' num2str(dof) ') = ' num2str(stat_no_perm) ', p = ' num2str(pval)])
disp('    ')

disp('2) Neutral - Unpleasant')
disp(['Neutral: Median = ' num2str(100*median(MeDwTi_neut,'omitnan')/4) '%, Mean = ' num2str(100*mean(MeDwTi_neut,'omitnan')/4) '% , SD = ' num2str(100*std(MeDwTi_neut,'omitnan')/4) '%'])
disp(['Unpleasant: Median = ' num2str(100*median(MeDwTi_unpl,'omitnan')/4) '%, Mean = ' num2str(100*mean(MeDwTi_unpl,'omitnan')/4) '%, SD = ' num2str(100*std(MeDwTi_unpl,'omitnan')/4) '%'])

[pval, stat_no_perm, dof]=PermutePairedTest(MeDwTi_neut, MeDwTi_unpl, 1500, 'tstat', 'both');
stat_no_perm=round(stat_no_perm*1000)/1000;
pval=round(pval*1000)/1000;
disp(['t(' num2str(dof) ') = ' num2str(stat_no_perm) ', p = ' num2str(pval)])
