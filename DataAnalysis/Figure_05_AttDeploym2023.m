%% LOAD ANT DATA
clearvars
AntData1=readtable('DataANT_AttDeploy01.xlsx');
AntData1=AntData1(logical(AntData1.Include),:);

% specify directory where the behavioral data is located
behavData_path='';

% specify directory to save figures
path2saveFigs='';

%% calculate

% calculate intensity and valence means

nParticip=size(AntData1,1);

%preallocate
Int_FocusFree=nan(nParticip,1);
Int_FocusNonArous=nan(nParticip,1);
Int_FocusArous=nan(nParticip,1);

Val_FocusFree=nan(nParticip,1);
Val_FocusNonArous=nan(nParticip,1);
Val_FocusArous=nan(nParticip,1);

for thisFile=1:nParticip
    
        
        load([behavData_path filesep AntData1.Mat_Filename{thisFile}])
        
        I=BehavData.vars.ResponseIntensity_seq;
        V=BehavData.vars.ResponseValence_seq;
        Stim=BehavData.vars.Stim_seq;
        
        % 1 = IMG NEUTRAL, NO FOCUS
        % 2 = IMG NEUTRAL, FOCUS NON-AROUSING
        % 3 = IMG UNPLEASANT, NO FOCUS
        % 4 = IMG UNPLEASANT, FOCUS NON-AROUSING
        % 5 = IMG UNPLEASANT, FOCUS AROUSING
        
        
        
        Int_FocusFree(thisFile)=mean(I(Stim==3),'omitnan');
        Int_FocusNonArous(thisFile)=mean(I(Stim==4),'omitnan');
        Int_FocusArous(thisFile)=mean(I(Stim==5),'omitnan');
        
        Val_FocusFree(thisFile)=mean(V(Stim==3),'omitnan');
        Val_FocusNonArous(thisFile)=mean(V(Stim==4),'omitnan');
        Val_FocusArous(thisFile)=mean(V(Stim==5),'omitnan');
        
end


% calculate deltas (AD estimate)

%NA: NonArousing
%A: Arousing
%F: Free

deltaI_NA_minus_F=Int_FocusNonArous - Int_FocusFree;
deltaI_NA_minus_A=Int_FocusNonArous - Int_FocusArous;

deltaV_NA_minus_F=Val_FocusNonArous - Val_FocusFree;
deltaV_NA_minus_A=Val_FocusNonArous - Val_FocusArous;





%% plot v1 (including ANT scores)

h=0.36;
w=0.37;
left=0.13;
bott0=0.1;
bott1=bott0+1.4*h;

posits=[left bott1 w h;
    left+1.21*w bott1 w h;
    left bott0 w h;
    left+1.21*w bott0 w h];

fsize=13;
fname='Verdana';

xTL={'Alert','Orientation','Executive'};

ms=25;

mult_factor=1;

% cols=[0.24    0.15    0.66;
%     0.0433    0.7    0.7;
%     0.7    0.54    0];
cols=[0.24    0.15    0.66;
    0.24    0.15    0.66;
    0.24    0.15    0.66];


fig1=figure;

ax1=axes('parent', fig1, 'TickDir','out', 'position', posits(1,:),...
    'xtick',[0 1 2],'xticklabel',xTL,...
    'ytick',mult_factor*(0:50:200),...
    'FontSize',fsize, 'FontName',fname);
hold(ax1,'on')

X=randn(size(AntData1,1),1)./15;


plot([mean(X)-0.2 mean(X)+0.2],mean([AntData1.ANT_Alert AntData1.ANT_Alert],'omitnan'),...
    'linewidth',6,'color',[0.75 0.4 0.4])
plot([mean(X)-0.2 mean(X)+0.2]+1,mean([AntData1.ANT_Orient  AntData1.ANT_Orient],'omitnan'),...
    'linewidth',6,'color',[0.75 0.4 0.4])
plot([mean(X)-0.2 mean(X)+0.2]+2,mean([AntData1.ANT_Executive  AntData1.ANT_Executive],'omitnan'),...
    'linewidth',6,'color',[0.75 0.4 0.4])



for j=1:numel(X)

    plot([X(j) X(j)+1 X(j)+2],...
        [AntData1.ANT_Alert(j) AntData1.ANT_Orient(j) AntData1.ANT_Executive(j)],...
        'color',[0.6 0.6 0.6])

end

plot([X X+1 X+2],[AntData1.ANT_Alert AntData1.ANT_Orient AntData1.ANT_Executive],'marker','.',...
    'linestyle','none','markersize',ms,'Color','k')


ylabel('ANT score (ms)')

xlim([-0.5 2.5])
ylim([-1 190])
ax2=axes('parent', fig1, 'TickDir','out', 'position', posits(2,:),...
    'xtick',[0 1 2],'xticklabel',xTL,...
    'FontSize',fsize, 'FontName',fname);
hold(ax2,'on')


X=randn(size(AntData1,1),1)./15;


plot([mean(X)-0.2 mean(X)+0.2],mean([AntData1.Error_Alert_Center_and_Double_Cue AntData1.Error_Alert_Center_and_Double_Cue],'omitnan'),...
    'linewidth',6,'color',[0.75 0.4 0.4])
plot([mean(X)-0.2 mean(X)+0.2]+1,mean([AntData1.Error_Orient_SpatialCue  AntData1.Error_Orient_SpatialCue],'omitnan'),...
    'linewidth',6,'color',[0.75 0.4 0.4])
plot([mean(X)-0.2 mean(X)+0.2]+2,mean([AntData1.Error_Executive_Incongruent AntData1.Error_Executive_Incongruent],'omitnan'),...
    'linewidth',6,'color',[0.75 0.4 0.4])

for j=1:numel(X)
    plot([X(j) X(j)+1 X(j)+2],...
        [AntData1.Error_Alert_Center_and_Double_Cue(j) AntData1.Error_Orient_SpatialCue(j) AntData1.Error_Executive_Incongruent(j)],...
        'color',[0.6 0.6 0.6])
end

plot([X X+1 X+2],[AntData1.Error_Alert_Center_and_Double_Cue AntData1.Error_Orient_SpatialCue AntData1.Error_Executive_Incongruent],'marker','.',...
    'linestyle','none','markersize',ms,'Color','k')

ylabel('ANT error (ms)')

ylim([-1 24])

ax3=axes('parent', fig1, 'TickDir','out', 'position', posits(3,:),...
    'FontSize',fsize, 'FontName',fname);
hold(ax3,'on')

plot(AntData1.ANT_Alert,deltaV_NA_minus_A,'linestyle','none',...
    'Marker','.','markersize',ms,'color',cols(1,:))

plot(AntData1.ANT_Orient,deltaV_NA_minus_A,'linestyle','none',...
    'Marker','o','linewidth',2,'markersize',ms-17,'color',cols(2,:))

plot(AntData1.ANT_Executive,deltaV_NA_minus_A,'linestyle','none',...
    'Marker','+','linewidth',2,'markersize',ms-17,'color',cols(3,:))


text(147,4.2,xTL(1),'fontsize',fsize+1,'color',cols(1,:),'fontweight','bold')
text(147,3.6,xTL(2),'fontsize',fsize+1,'color',cols(2,:),'fontweight','bold')
text(147,3,xTL(3),'fontsize',fsize+1,'color',cols(3,:),'fontweight','bold')
plot(139,4.2,'linestyle','none',...
    'Marker','.','markersize',ms,'color',cols(1,:))
plot(139,3.6,'linestyle','none',...
    'Marker','o','linewidth',2.5,'markersize',ms-17,'color',cols(2,:))
plot(139,3,'linestyle','none',...
    'Marker','+','linewidth',2.5,'markersize',ms-17,'color',cols(3,:))



ylabel({'Attentional Deployment','(Non-Arousing minus Arousing','rating difference)'})
xlabel('ANT score')


title('Valence')
xlim([-1 180])

ax4=axes('parent', fig1, 'TickDir','out', 'position', posits(4,:),...
    'FontSize',fsize, 'FontName',fname);
hold(ax4,'on')

plot(AntData1.ANT_Alert,deltaI_NA_minus_A,'linestyle','none',...
    'Marker','.','markersize',ms,'color',cols(1,:))

plot(AntData1.ANT_Orient,deltaI_NA_minus_A,'linestyle','none',...
    'Marker','o','linewidth',2,'markersize',ms-17,'color',cols(2,:))

plot(AntData1.ANT_Executive,deltaI_NA_minus_A,'linestyle','none',...
    'Marker','+','linewidth',2,'markersize',ms-17,'color',cols(3,:))



xlabel('ANT score (ms)')

title('Intensity')

xlim([-1 180])
ylim([-4.1 4.3])

set(fig1,'PaperUnits','centimeters')
set(fig1, 'PaperPosition', [0 0 25 20])


% print(fig1,[path2saveFigs filesep 'FIG_05_AttDep2023'],'-dpng','-r350')

%% plot V2 (only correlations between ANT and AD measures)

h=0.82;
w=0.37;
left=0.13;
bott1=0.13;


posits=[left bott1 w h;
    left+1.21*w bott1 w h];

fsize=13;
fname='Verdana';

xTL={'Alert','Orientation','Executive'};

ms=25;

mult_factor=1;

% cols=[0.24    0.15    0.66;
%     0.0433    0.7    0.7;
%     0.7    0.54    0];
cols=[0.24    0.15    0.66;
    0.24    0.15    0.66;
    0.24    0.15    0.66];
cols=zeros(3,3);

fig1=figure;

ax1=axes('parent', fig1, 'TickDir','out', 'position', posits(1,:),...
    'FontSize',fsize, 'FontName',fname);
hold(ax1,'on')


plot(AntData1.ANT_Alert,deltaI_NA_minus_A,'linestyle','none',...
    'Marker','.','markersize',ms,'color',cols(1,:))

plot(AntData1.ANT_Orient,deltaI_NA_minus_A,'linestyle','none',...
    'Marker','o','linewidth',2,'markersize',ms-17,'color',cols(2,:))

plot(AntData1.ANT_Executive,deltaI_NA_minus_A,'linestyle','none',...
    'Marker','+','linewidth',2,'markersize',ms-17,'color',cols(3,:))


ylabel({'Attentional Deployment','(Non-Arousing minus Arousing','rating difference)'})
xlabel('ANT score (ms)')

title('Intensity')

xlim([-1 180])
ylim([-4.1 4.3])


ax2=axes('parent', fig1, 'TickDir','out', 'position', posits(2,:),...
    'FontSize',fsize, 'FontName',fname);
hold(ax2,'on')

plot(AntData1.ANT_Alert,deltaV_NA_minus_A,'linestyle','none',...
    'Marker','.','markersize',ms,'color',cols(1,:))

plot(AntData1.ANT_Orient,deltaV_NA_minus_A,'linestyle','none',...
    'Marker','o','linewidth',2,'markersize',ms-17,'color',cols(2,:))

plot(AntData1.ANT_Executive,deltaV_NA_minus_A,'linestyle','none',...
    'Marker','+','linewidth',2,'markersize',ms-17,'color',cols(3,:))

text(144,4.1,xTL(1),'fontsize',fsize+1,'color',cols(1,:),'fontweight','bold')
text(144,3.6,xTL(2),'fontsize',fsize+1,'color',cols(2,:),'fontweight','bold')
text(144,3.1,xTL(3),'fontsize',fsize+1,'color',cols(3,:),'fontweight','bold')
plot(136,4.1,'linestyle','none',...
    'Marker','.','markersize',ms,'color',cols(1,:))
plot(136,3.6,'linestyle','none',...
    'Marker','o','linewidth',2.5,'markersize',ms-17,'color',cols(2,:))
plot(136,3.1,'linestyle','none',...
    'Marker','+','linewidth',2.5,'markersize',ms-17,'color',cols(3,:))


xlabel('ANT score (ms)')

title('Valence')
xlim([-1 180])

set(fig1,'PaperUnits','centimeters')
set(fig1, 'PaperPosition', [0 0 23 13])

print(fig1,[path2saveFigs filesep 'FIG_05'],'-dpng','-r350')

%% stats
disp('%%%%%%%%%%%%%%%%%%')
disp('ANT scores (ms)')
disp(['ALERT: Median = ' num2str(median(AntData1.ANT_Alert,'omitnan')) ', M = ' num2str(mean(AntData1.ANT_Alert,'omitnan')) ', SD = ' num2str(std(AntData1.ANT_Alert,'omitnan')) ])
disp(['ORIENTATION: Median = ' num2str(median(AntData1.ANT_Orient,'omitnan')) ', M = ' num2str(mean(AntData1.ANT_Orient,'omitnan')) ', SD = ' num2str(std(AntData1.ANT_Orient,'omitnan')) ])
disp(['EXECUTIVE: Median = ' num2str(median(AntData1.ANT_Executive,'omitnan')) ', M = ' num2str(mean(AntData1.ANT_Executive,'omitnan')) ', SD = ' num2str(std(AntData1.ANT_Executive,'omitnan')) ])
disp('==')
disp('ANT error rates (%)')
disp(['ALERT: Median = ' num2str(median(AntData1.Error_Alert_Center_and_Double_Cue,'omitnan')) ', M = ' num2str(mean(AntData1.Error_Alert_Center_and_Double_Cue,'omitnan')) ', SD = ' num2str(std(AntData1.Error_Alert_Center_and_Double_Cue,'omitnan')) ])
disp(['ORIENTATION: Median = ' num2str(median(AntData1.Error_Orient_SpatialCue,'omitnan')) ', M = ' num2str(mean(AntData1.Error_Orient_SpatialCue,'omitnan')) ', SD = ' num2str(std(AntData1.Error_Orient_SpatialCue,'omitnan')) ])
disp(['EXECUTIVE: Median = ' num2str(median(AntData1.Error_Executive_Incongruent,'omitnan')) ', M = ' num2str(mean(AntData1.Error_Executive_Incongruent,'omitnan')) ', SD = ' num2str(std(AntData1.Error_Executive_Incongruent,'omitnan')) ])

%%
disp('%%%%%%%%%%%%%%%%%%')
disp('Pearson Correlation between Intensity-based AD and ANT subscales')

corrType='Pearson';

disp('Alert')
[R1,P1]=corr(AntData1.ANT_Alert,deltaI_NA_minus_A,'rows','complete','type',corrType);
disp(['r = ' num2str(R1) ', p = ' num2str(P1)])
disp('Orientation')
[R2,P2]=corr(AntData1.ANT_Orient,deltaI_NA_minus_A,'rows','complete','type',corrType);
disp(['r = ' num2str(R2) ', p = ' num2str(P2)])
disp('Executive')
[R3,P3]=corr(AntData1.ANT_Executive,deltaI_NA_minus_A,'rows','complete','type',corrType);
disp(['r = ' num2str(R3) ', p = ' num2str(P3)])

disp('%%%%%%%%%%%%%')
disp('Pearson Correlation between Valence-based AD and ANT subscales')
disp('Alert')
[R1,P1]=corr(AntData1.ANT_Alert,deltaV_NA_minus_A,'rows','complete','type',corrType);
disp(['r = ' num2str(R1) ', p = ' num2str(P1)])
disp('Orientation')
[R2,P2]=corr(AntData1.ANT_Orient,deltaV_NA_minus_A,'rows','complete','type',corrType);
disp(['r = ' num2str(R2) ', p = ' num2str(P2)])
disp('Executive')
[R3,P3]=corr(AntData1.ANT_Executive,deltaV_NA_minus_A,'rows','complete','type',corrType);
disp(['r = ' num2str(R3) ', p = ' num2str(P3)])
