%% relevant paths & data files

%run script to get behav and eyetracking data file pairs
map_asc2mat_AttDeployment

%use only participants with both kinds of data
MasterTable_Files=MasterTable_Files(MasterTable_Files.overall,:);


% specify directory where the behavioral data is located
behavData_path='';

% specify directory to save figures
path2saveFigs='';

%% GET DATA (INTENSITY, VALENCE, OTHERS)

nParticip=sum(MasterTable_Files.overall);

nTrials=20;

%preallocate
cell4Excel=cell(nParticip,1);
female=false(nParticip,1);
Int_mean_Mat=nan(5,nParticip); % 5 KINDS OF TRIALS
Val_mean_Mat=nan(5,nParticip);
Int_Mat=nan(nTrials,nParticip);
Val_Mat=nan(nTrials,nParticip);
sessDurat_min=nan(nParticip,1);
RT_Val_Mat=nan(nTrials,nParticip);
RT_Int_Mat=nan(nTrials,nParticip);


partic_ID=cell(1,1);
Stim_Mat=nan(nTrials,nParticip);

for j=1:nParticip
    
    load([behavData_path filesep MasterTable_Files.mat{j}])
    
    if strcmpi(BehavData.info.Subject_Gender,'M')
    else
        female(j)=true;
    end
    
    I=BehavData.vars.ResponseIntensity_seq;
    V=BehavData.vars.ResponseValence_seq;
    Stim=BehavData.vars.Stim_seq;
    
    Stim_Mat(:,j)=Stim;
    Val_Mat(:,j)=V;
    Int_Mat(:,j)=I;
    
    %INTENSITY
    I1=median(I(Stim==1),'omitnan');  % 1 = IMG NEUTRAL, NO FOCUS
    I2=median(I(Stim==2),'omitnan');  % 2 = IMG NEUTRAL, FOCUS NON-AROUSING
    I3=median(I(Stim==3),'omitnan');  % 3 = IMG UNPLEASANT, NO FOCUS
    I4=median(I(Stim==4),'omitnan');  % 4 = IMG UNPLEASANT, FOCUS NON-AROUSING
    I5=median(I(Stim==5),'omitnan');  % 5 = IMG UNPLEASANT, FOCUS AROUSING
    %VALENCE
    V1=median(V(Stim==1),'omitnan');
    V2=median(V(Stim==2),'omitnan');
    V3=median(V(Stim==3),'omitnan');
    V4=median(V(Stim==4),'omitnan');
    V5=median(V(Stim==5),'omitnan');
    
    Int_mean_Mat(:,j)=[I1;I2;I3;I4;I5];
    Val_mean_Mat(:,j)=[V1;V2;V3;V4;V5];
    
    
    partic_ID{j,1}=MasterTable_Files.asc;
    
    
    % RESPONSE TIME
    RT_Val_Mat(:,j)=BehavData.vars.RTvalence_seq;
    RT_Int_Mat(:,j)=BehavData.vars.RTintensity_seq;
    
end

%calculate type of stim sequence (1 out of 3 possibilities)

%the avg of the vector of differences provide a unique ID for the sequence
temp=mean(diff(Stim_Mat));
classes=unique(temp);

stimSeq_label=nan(nParticip,1);

for j=1:nParticip
    for jj=1:numel(classes) %loop through possibilities
        if temp(j)==classes(jj) %if match, label with number 1, 2 or 3
            stimSeq_label(j)=jj;
        end
    end
end
    

%now store intensity and valence measures sorted by stim seq
for i=1:3 %loop through three types of stim sequences
    curr_logi=stimSeq_label==i;
    
    eval(['Int_Mat_stimSeq_' num2str(i) '= Int_Mat(:,curr_logi);'])
    eval(['Val_Mat_stimSeq_' num2str(i) '= Val_Mat(:,curr_logi);'])
    
    S=Stim_Mat(:,curr_logi);
    S=S(:,1);
    
    eval(['Stim_seq_type_' num2str(i) '= S;'])
    
    
end



%delta intensity  (NA:non-arousal; A: Arousal; F: Free)
deltaI_NA_minus_A=Int_mean_Mat(4,:)-Int_mean_Mat(5,:);
deltaV_NA_minus_A=Val_mean_Mat(4,:)-Val_mean_Mat(5,:);

deltaI_NA_minus_F=Int_mean_Mat(4,:)-Int_mean_Mat(3,:);
deltaV_NA_minus_F=Val_mean_Mat(4,:)-Val_mean_Mat(3,:);

disp(['n female = ' num2str(sum(female))])
disp(['n tot = ' num2str(nParticip)])


% 1 = IMG NEUTRAL, NO FOCUS
% 2 = IMG NEUTRAL, FOCUS NON-AROUSING
% 3 = IMG UNPLEASANT, NO FOCUS
% 4 = IMG UNPLEASANT, FOCUS NON-AROUSING
% 5 = IMG UNPLEASANT, FOCUS AROUSING

%calculate mean intensities and valences parsed by stim sequence
% (to check if there was effect of sequence)

%sequence1
Int_mean_Mat_stimSeq_1=nan(5,size(Int_Mat_stimSeq_1,2));
Val_mean_Mat_stimSeq_1=nan(5,size(Val_Mat_stimSeq_1,2));
for thisParticip=1:size(Int_Mat_stimSeq_1,2) % loop through participants
   
   Stim=Stim_seq_type_1; 
    
   I=Int_Mat_stimSeq_1(:,thisParticip);
   
   I1=mean(I(Stim==1),'omitnan');
   I2=mean(I(Stim==2),'omitnan');
   I3=mean(I(Stim==3),'omitnan');
   I4=mean(I(Stim==4),'omitnan');
   I5=mean(I(Stim==5),'omitnan');
   
   Int_mean_Mat_stimSeq_1(:,thisParticip)=[I1;I2;I3;I4;I5];
      
   V=Val_Mat_stimSeq_1(:,thisParticip);
   
   V1=mean(V(Stim==1),'omitnan');
   V2=mean(V(Stim==2),'omitnan');
   V3=mean(V(Stim==3),'omitnan');
   V4=mean(V(Stim==4),'omitnan');
   V5=mean(V(Stim==5),'omitnan');
   
   Val_mean_Mat_stimSeq_1(:,thisParticip)=[V1;V2;V3;V4;V5];
   
end

%sequence2
Int_mean_Mat_stimSeq_2=nan(5,size(Int_Mat_stimSeq_2,2));
Val_mean_Mat_stimSeq_2=nan(5,size(Val_Mat_stimSeq_2,2));
for thisParticip=1:size(Int_Mat_stimSeq_2,2) % loop through participants
   
   Stim=Stim_seq_type_2; 
    
   I=Int_Mat_stimSeq_2(:,thisParticip);
   
   I1=mean(I(Stim==1),'omitnan');
   I2=mean(I(Stim==2),'omitnan');
   I3=mean(I(Stim==3),'omitnan');
   I4=mean(I(Stim==4),'omitnan');
   I5=mean(I(Stim==5),'omitnan');
   
   Int_mean_Mat_stimSeq_2(:,thisParticip)=[I1;I2;I3;I4;I5];
      
   V=Val_Mat_stimSeq_2(:,thisParticip);
   
   V1=mean(V(Stim==1),'omitnan');
   V2=mean(V(Stim==2),'omitnan');
   V3=mean(V(Stim==3),'omitnan');
   V4=mean(V(Stim==4),'omitnan');
   V5=mean(V(Stim==5),'omitnan');
   
   Val_mean_Mat_stimSeq_2(:,thisParticip)=[V1;V2;V3;V4;V5];
   
end


%sequence3
Int_mean_Mat_stimSeq_3=nan(5,size(Int_Mat_stimSeq_3,2));
Val_mean_Mat_stimSeq_3=nan(5,size(Val_Mat_stimSeq_3,2));
for thisParticip=1:size(Int_Mat_stimSeq_3,2) % loop through participants
   
   Stim=Stim_seq_type_3; 
    
   I=Int_Mat_stimSeq_3(:,thisParticip);
   
   I1=mean(I(Stim==1),'omitnan');
   I2=mean(I(Stim==2),'omitnan');
   I3=mean(I(Stim==3),'omitnan');
   I4=mean(I(Stim==4),'omitnan');
   I5=mean(I(Stim==5),'omitnan');
   
   Int_mean_Mat_stimSeq_3(:,thisParticip)=[I1;I2;I3;I4;I5];
   
   
   V=Val_Mat_stimSeq_3(:,thisParticip);
   
   V1=mean(V(Stim==1),'omitnan');
   V2=mean(V(Stim==2),'omitnan');
   V3=mean(V(Stim==3),'omitnan');
   V4=mean(V(Stim==4),'omitnan');
   V5=mean(V(Stim==5),'omitnan');
   
   Val_mean_Mat_stimSeq_3(:,thisParticip)=[V1;V2;V3;V4;V5];
   
end

%% plot


w1=0.38;
w2=0.36;
h=0.29;
left=0.12;
bott0=0.05;
bott1=bott0+2.1*h;

posits=[left bott1 w1 h;
    left+1.3*w1 bott1 w1 h;
    left+0.5*w1 bott0 w2 h*1.1];


ms=14;
lw=7;

tickAngle=18;

% colours=[0.65 0.65 0.1;
%     0.1 0.65 0.65;
%     0.65 0.1 0.65];

colours=[0 0 0;
    0 0 0;
    0 0 0]+0.6;


greyCol=[0.7 0.7 0.7];
orangeCol=[1 0.4 0.15];

conditLabels={'Arousing','Non-Arousing'};

ylimi=[0.5 9.5];
xlimi=[-0.5 1.5];

xlimi_delta=[-0.5 1.5];
ylimi_delta=[-4.1 4.1];

fsize=15;
fname='Verdana';

fig1=figure('color','w');

ax1=axes('parent',fig1,'tickdir','out','position',posits(1,:),...
    'xtick',0:1,'ytick',1:2:9,'xticklabel',conditLabels,...
    'ticklength',[0.02 0.02],'fontsize',fsize,'fontname',fname);
hold(ax1,'all')
xtickangle(tickAngle)

X=randn(nParticip,1).*0.06;

% plot([-0.2 0.2],[nanmean(Int_mean_Mat(3,:)) nanmean(Int_mean_Mat(3,:))],'linewidth',lw,'color',[0 0 0])
plot([-0.2 0.2],[nanmean(Int_mean_Mat(5,:)) nanmean(Int_mean_Mat(5,:))],'linewidth',lw,'color', [0 0 0])
plot([0.8 1.2],[nanmean(Int_mean_Mat(4,:)) nanmean(Int_mean_Mat(4,:))],'linewidth',lw,'color', [0 0 0])

% plot(X, Int_mean_Mat(3,:)', 'marker','.','markersize',ms,'linestyle','none','color',colours(1,:))
plot(X, Int_mean_Mat(5,:)', 'marker','.','markersize',ms,'linestyle','none','color',[0 0 0])
plot(X+1, Int_mean_Mat(4,:)', 'marker','.','markersize',ms,'linestyle','none','color',[0 0 0])
%connector lines
for j=1:size(Int_mean_Mat,2)
   plot([X(j) X(j)+1],[Int_mean_Mat(5,j) Int_mean_Mat(4,j)],'color', colours(1,:)) 
   if Int_mean_Mat(4,j)-Int_mean_Mat(5,j)>2
       plot([X(j) X(j)+1],[Int_mean_Mat(5,j) Int_mean_Mat(4,j)],'color', 'r')
       disp('check on')
       disp(j)
   end
end

annotation(fig1,'textbox', [0.35 0.95 0.38 0.05],'String','Unpleasant Images','fontweight','bold',...
    'HorizontalAlignment','center','FontSize',fsize,'FitBoxToText','off','EdgeColor','none',...
    'fontname',fname);

xlabel('Focus type')
ylabel('Intensity rating')

xlim(xlimi)
ylim(ylimi)

ax2=axes('parent',fig1,'tickdir','out','position',posits(2,:),...
    'xtick',0:2,'ytick',1:2:9,'xticklabel',conditLabels,...
    'ticklength',[0.02 0.02],'fontsize',fsize,'fontname',fname);
hold(ax2,'all')
xtickangle(tickAngle)

X=randn(nParticip,1).*0.05;

% plot([-0.2 0.2],[nanmean(Val_mean_Mat(3,:)) nanmean(Val_mean_Mat(3,:))],'linewidth',lw,'color','k')
plot([-0.2 0.2],[nanmean(Val_mean_Mat(5,:)) nanmean(Val_mean_Mat(5,:))],'linewidth',lw,'color','k')
plot([0.8 1.2],[nanmean(Val_mean_Mat(4,:)) nanmean(Val_mean_Mat(4,:))],'linewidth',lw,'color','k')

% plot(X, Val_mean_Mat(3,:)', 'marker','.','markersize',ms,'linestyle','none','color',colours(1,:))
plot(X, Val_mean_Mat(5,:)', 'marker','.','markersize',ms,'linestyle','none','color',[0 0 0])
plot(X+1, Val_mean_Mat(4,:)', 'marker','.','markersize',ms,'linestyle','none','color',[0 0 0])

%connector lines
for j=1:size(Val_mean_Mat,2)
   plot([X(j) X(j)+1],[Val_mean_Mat(5,j) Val_mean_Mat(4,j)],'color', colours(1,:)) 
end

xlabel('Focus type')
ylabel('Valence rating')

xlim(xlimi)
ylim(ylimi)


ax3=axes('parent',fig1,'position',posits(3,:),'tickdir','out',...
        'xtick',[0 1],'xticklabel',{'Intensity','Valence'},...
    'fontsize',fsize,'fontname',fname);
hold(ax3,'on')
% plot([-0.2 0.2],[mean(deltaI_NA_minus_A) mean(deltaI_NA_minus_A)],...
%     'linewidth',lw+1,'color',orangeCol-min(orangeCol))
% plot([0.8 1.2],[mean(deltaV_NA_minus_A) mean(deltaV_NA_minus_A)],...
%     'linewidth',lw,'color',orangeCol-min(orangeCol))

%plot horizontal line at zero
plot(xlimi_delta,[0 0],'color',[0 0 0])

%produce boxplot
[box,low_whisk,high_whisk]=data4boxplot(deltaI_NA_minus_A);

xlims_leftBoxplot=[-0.45 -0.2];

for j=1:3
    plot(xlims_leftBoxplot,[box(j) box(j)],'linewidth',3,'color',orangeCol)
end

plot([mean(xlims_leftBoxplot) mean(xlims_leftBoxplot)],[box(1) low_whisk],'linewidth',3,'color',orangeCol)
plot([mean(xlims_leftBoxplot) mean(xlims_leftBoxplot)],[box(3) high_whisk],'linewidth',3,'color',orangeCol)

plot([xlims_leftBoxplot(1) xlims_leftBoxplot(1)],[box(1)-0.065 box(3)+0.065],'linewidth',3,'color',orangeCol)
plot([xlims_leftBoxplot(2) xlims_leftBoxplot(2)],[box(1)-0.065 box(3)+0.065],'linewidth',3,'color',orangeCol)


[box,low_whisk,high_whisk,data2plot]=data4boxplot(deltaV_NA_minus_A);

xlims_rightBoxplot=[1.25 1.45];

for j=1:3
    plot(xlims_rightBoxplot,[box(j) box(j)],'linewidth',3,'color',orangeCol)
end

plot([mean(xlims_rightBoxplot) mean(xlims_rightBoxplot)],[box(1) low_whisk],'linewidth',3,'color',orangeCol)
plot([mean(xlims_rightBoxplot) mean(xlims_rightBoxplot)],[box(3) high_whisk],'linewidth',3,'color',orangeCol)

plot([xlims_rightBoxplot(1) xlims_rightBoxplot(1)],[box(1)-0.065 box(3)+0.065],'linewidth',3,'color',orangeCol)
plot([xlims_rightBoxplot(2) xlims_rightBoxplot(2)],[box(1)-0.065 box(3)+0.065],'linewidth',3,'color',orangeCol)




plot([X X+1],[deltaI_NA_minus_A' deltaV_NA_minus_A'],...
    'linestyle','none','marker','.','markersize',ms,'color',orangeCol)








ylabel({'Attentional Deployment','(rating difference)'})
title({'Non-Arousing', 'minus Arousing'})
xlim(xlimi_delta)
ylim(ylimi_delta)


lettW=0.07;
lettH=0.07;

annotation(fig1,'textbox', [0.005 0.93 lettW lettH],'String','A','fontweight','bold',...
    'HorizontalAlignment','center','FontSize',fsize+14,'FitBoxToText','off','EdgeColor','none',...
    'fontname',fname);

annotation(fig1,'textbox', [0.008 0.44 lettW lettH],'String','B','fontweight','bold',...
    'HorizontalAlignment','center','FontSize',fsize+14,'FitBoxToText','off','EdgeColor','none',...
    'fontname',fname);


set(fig1,'PaperUnits','inches')
set(fig1, 'PaperPosition', [0 0 7.5 8])


print(fig1,[path2saveFigs filesep 'FIG_04'],'-dpng','-r450')

%% compare between sequences

% 1 = IMG NEUTRAL, NO FOCUS
% 2 = IMG NEUTRAL, FOCUS NON-AROUSING
% 3 = IMG UNPLEASANT, NO FOCUS
% 4 = IMG UNPLEASANT, FOCUS NON-AROUSING
% 5 = IMG UNPLEASANT, FOCUS AROUSING

for j=1:5

trial_type=j;

y1=Val_mean_Mat_stimSeq_1(trial_type,:);
x1=ones(size(y1));
y2=Val_mean_Mat_stimSeq_2(trial_type,:);
x2=1+ones(size(y2));
y3=Val_mean_Mat_stimSeq_3(trial_type,:);
x3=2+ones(size(y3));

Y=[y1'; y2'; y3';];
X=[x1';x2';x3';];


[p,tbl,stats]=anovan(Y,X,'varnames',{'Sequence type'})


end


%% plot intensity vs valence

fig1=figure;

ax=axes('parent',fig1,'tickdir','out','fontsize',fsize','fontname',fname);
hold(ax,'on')




%% plot delta int vs delta val


fsize=13;
fname='Verdana';

ms=22;

ylimi=[-3 3];
xlimi=[-4 3];

fig=figure;

ax=axes('parent',fig,'tickdir','out','fontsize',fsize,'fontname',fname);
hold(ax,'on')


plot(deltaI_NA_minus_F,deltaV_NA_minus_F,'linestyle','none','marker','.',...
    'markersize',ms,'color',[0 0 0])

plot(xlimi,[0 0],'color',[0 0 0])
plot([0 0],ylimi,'color',[0 0 0])

xlabel('Intensity rating difference')
ylabel('Valence rating difference')
title('Non-Arousing minus Free')
xlim(xlimi)
ylim(ylimi)


print('FIG04_C','-dpng','-r350')
