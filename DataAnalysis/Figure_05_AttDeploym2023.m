%% LOAD ANT DATA
clearvars
AntData1=readtable('ANTData.xlsx');
AntData2=readtable('Nuevos Datos ANT & Tarea DA.xlsx')

% specify directory where the behavioral data is located
behavData_path='';
path2saveFigs='';

AntData1=AntData1(logical(AntData1.Include),:);

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
    
    if AntData1.Include
        
        load([behavData_path filesep AntData1.Filename{thisFile}])
        
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
    
end


% calculate deltas (AD estimate)

%NA: NonArousing
%A: Arousing
%F: Free

deltaI_NA_minus_F_2019=Int_FocusNonArous - Int_FocusFree;
deltaI_NA_minus_A_2019=Int_FocusNonArous - Int_FocusArous;

deltaV_NA_minus_F_2019=Val_FocusNonArous - Val_FocusFree;
deltaV_NA_minus_A_2019=Val_FocusNonArous - Val_FocusArous;





%% plot correlations between ANT and AD measures

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


fig1=figure;

ax1=axes('parent', fig1, 'TickDir','out', 'position', posits(1,:),...
    'FontSize',fsize, 'FontName',fname);
hold(ax1,'on')

plot(AntData1.ANT_Alert,deltaI_NA_minus_A_2019,'linestyle','none',...
    'Marker','.','markersize',ms,'color',cols(1,:))

plot(AntData1.ANT_Orient,deltaI_NA_minus_A_2019,'linestyle','none',...
    'Marker','o','linewidth',2,'markersize',ms-17,'color',cols(2,:))

plot(AntData1.ANT_Executive,deltaI_NA_minus_A_2019,'linestyle','none',...
    'Marker','+','linewidth',2,'markersize',ms-17,'color',cols(3,:))


ylabel({'Attentional Deployment','(Non-Arousing minus Arousing','rating difference)'})
xlabel('ANT score (ms)')

title('Intensity')

xlim([-1 180])
ylim([-4.1 4.3])


ax2=axes('parent', fig1, 'TickDir','out', 'position', posits(2,:),...
    'FontSize',fsize, 'FontName',fname);
hold(ax2,'on')

plot(AntData1.ANT_Alert,deltaV_NA_minus_A_2019,'linestyle','none',...
    'Marker','.','markersize',ms,'color',cols(1,:))

plot(AntData1.ANT_Orient,deltaV_NA_minus_A_2019,'linestyle','none',...
    'Marker','o','linewidth',2,'markersize',ms-17,'color',cols(2,:))

plot(AntData1.ANT_Executive,deltaV_NA_minus_A_2019,'linestyle','none',...
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

% print(fig1,[path2saveFigs filesep 'FIG_05_AttDep2023'],'-dpng','-r350')

%% stats
disp('%%%%%%%%%%%%%%%%%%')
disp('ANT scores (ms)')
disp(['ALERT: Median = ' num2str(median(AntData1.ANT_Alert,'omitnan')) ', M = ' num2str(mean(AntData1.ANTALERTA,'omitnan')) ', SD = ' num2str(std(AntData1.ANTALERTA,'omitnan')) ])
disp(['ORIENTATION: Median = ' num2str(median(AntData1.ANT_Orient,'omitnan')) ', M = ' num2str(mean(AntData1.ANTORIENTACION,'omitnan')) ', SD = ' num2str(std(AntData1.ANTORIENTACION,'omitnan')) ])
disp(['EXECUTIVE: Median = ' num2str(median(AntData1.ANT_Executive,'omitnan')) ', M = ' num2str(mean(AntData1.ANTCEJECUTIVO,'omitnan')) ', SD = ' num2str(std(AntData1.ANTCEJECUTIVO,'omitnan')) ])
disp('==')
disp('ANT error rates (%)')
disp(['ALERT: Median = ' num2str(median(AntData2.ERRORALERTA_Center_DoubleCue_,'omitnan')) ', M = ' num2str(mean(AntData2.ERRORALERTA_Center_DoubleCue_,'omitnan')) ', SD = ' num2str(std(AntData2.ERRORALERTA_Center_DoubleCue_,'omitnan')) ])
disp(['ORIENTATION: Median = ' num2str(median(AntData2.ERRORORIENTACION_SpatialCue_,'omitnan')) ', M = ' num2str(mean(AntData2.ERRORORIENTACION_SpatialCue_,'omitnan')) ', SD = ' num2str(std(AntData2.ERRORORIENTACION_SpatialCue_,'omitnan')) ])
disp(['EXECUTIVE: Median = ' num2str(median(AntData2.ERRORATENCIONEJECUTIVA_Incongruent_,'omitnan')) ', M = ' num2str(mean(AntData2.ERRORATENCIONEJECUTIVA_Incongruent_,'omitnan')) ', SD = ' num2str(std(AntData2.ERRORATENCIONEJECUTIVA_Incongruent_,'omitnan')) ])

%%
disp('%%%%%%%%%%%%%%%%%%')
disp('Pearson Correlation between Intensity-based AD and ANT subscales')

corrType='Pearson';

disp('Alert')
[R1,P1]=corr(AntData1.ANT_Alert,deltaI_NA_minus_A_2019,'rows','complete','type',corrType);
disp(['r = ' num2str(R1) ', p = ' num2str(P1)])
disp('Orientation')
[R2,P2]=corr(AntData1.ANT_Orient,deltaI_NA_minus_A_2019,'rows','complete','type',corrType);
disp(['r = ' num2str(R2) ', p = ' num2str(P2)])
disp('Executive')
[R3,P3]=corr(AntData1.ANT_Executive,deltaI_NA_minus_A_2019,'rows','complete','type',corrType);
disp(['r = ' num2str(R3) ', p = ' num2str(P3)])

disp('%%%%%%%%%%%%%')
disp('Pearson Correlation between Valence-based AD and ANT subscales')
disp('Alert')
[R1,P1]=corr(AntData1.ANT_Alert,deltaV_NA_minus_A_2019,'rows','complete','type',corrType);
disp(['r = ' num2str(R1) ', p = ' num2str(P1)])
disp('Orientation')
[R2,P2]=corr(AntData1.ANT_Orient,deltaV_NA_minus_A_2019,'rows','complete','type',corrType);
disp(['r = ' num2str(R2) ', p = ' num2str(P2)])
disp('Executive')
[R3,P3]=corr(AntData1.ANT_Executive,deltaV_NA_minus_A_2019,'rows','complete','type',corrType);
disp(['r = ' num2str(R3) ', p = ' num2str(P3)])
