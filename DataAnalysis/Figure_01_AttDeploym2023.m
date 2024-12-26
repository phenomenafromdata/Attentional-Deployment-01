%% relevant paths & data files

%run script to get behav and eyetracking data file pairs
map_asc2mat_AttDeployment

behavData_path='pwd/data_behavior';

path2saveFigs='pwd';


%use only participants with both kinds of data

%MasterTable_Files=MasterTable_Files(MasterTable_Files.overall,:);




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
    I1=nanmean(I(Stim==1));  % 1 = IMG NEUTRAL, NO FOCUS
    I2=nanmean(I(Stim==2));  % 2 = IMG NEUTRAL, FOCUS NON-AROUSING
    I3=nanmean(I(Stim==3));  % 3 = IMG UNPLEASANT, NO FOCUS
    I4=nanmean(I(Stim==4));  % 4 = IMG UNPLEASANT, FOCUS NON-AROUSING
    I5=nanmean(I(Stim==5));  % 5 = IMG UNPLEASANT, FOCUS AROUSING
    %VALENCE
    V1=nanmean(V(Stim==1));
    V2=nanmean(V(Stim==2));
    V3=nanmean(V(Stim==3));
    V4=nanmean(V(Stim==4));
    V5=nanmean(V(Stim==5));
    
    Int_mean_Mat(:,j)=[I1;I2;I3;I4;I5];
    Val_mean_Mat(:,j)=[V1;V2;V3;V4;V5];
    
    
    partic_ID{j,1}=MasterTable_Files.asc;
    
    %calculate session duration
    ini=BehavData.info.timestart;
    fin=BehavData.info.timefinish;
    tempI=regexp(ini,'-','split');
    tempF=regexp(fin,'-','split');
    
    if numel(tempI)<2 || numel(tempF)<2
    else
        min_ini=str2double(tempI{1})*60 + str2double(tempI{2}); % convert time to mins
        min_fin=str2double(tempF{1})*60 + str2double(tempF{2});
        
        sessDurat_min(j)=min_fin-min_ini;
        
    end
    
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
    
    eval(['Int_Mat_type_' num2str(i) '= Int_Mat(:,curr_logi);'])
    eval(['Val_Mat_type_' num2str(i) '= Val_Mat(:,curr_logi);'])
    
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


% load IAPs data
IAPSdata=readtable('IAPS stimulus data.xlsx');


%% FIGURE

% plot IAPS' intensity and valence distributions

filt_unpl=strcmp(IAPSdata.Type,'Unpleasant');
filt_neut=strcmp(IAPSdata.Type,'Neutral');


w=0.39;
h=0.33;
bott=0.63;
posits=[0.1 bott w h;
    0.57 bott w h];



neutral_col=[0 0.6 0];
unpleas_col=[0.7 0 0];

fname='Verdana';
fsize=17;

fig1=figure;

for j=1:2
    
    if j==1
        unp=IAPSdata.ValenceMean(filt_unpl);
        neu=IAPSdata.ValenceMean(filt_neut);
        xL='IAPS Valence ratings';
    elseif j==2
        unp=IAPSdata.ArousalMean(filt_unpl);
        neu=IAPSdata.ArousalMean(filt_neut);
        xL='IAPS Intensity ratings';
    end
    
    ax=axes('parent',fig1,'tickdir','out','xtick',1:9,'position',posits(j,:),...
        'fontsize',fsize,'fontname',fname);
    hold(ax,'on')
    
    H=histogram(unp);
    H.FaceColor=unpleas_col;
    
    H=histogram(neu);
    H.FaceColor=neutral_col;
    
    xlabel(xL)
    ylabel('Number of images')
    xlim([0.5 9.5])
    if j==2
        text(7,24,'Neutral','FontWeight','bold','color',neutral_col,...
            'FontSize',fsize-1,'fontname',fname)
        text(7,22,'Unpleasant','FontWeight','bold','color',unpleas_col,...
            'FontSize',fsize-1,'fontname',fname)
    end

end


lettW=0.04;
lettH=0.07;

annotation(fig1,'textbox', [0.005 0.93 lettW lettH],'String','C','fontweight','bold',...
    'HorizontalAlignment','center','FontSize',fsize+26,'FitBoxToText','off','EdgeColor','none',...
    'fontname',fname);



% plot: 3 stimuli sequences

ns=[size(Val_Mat_type_1,2) size(Val_Mat_type_2,2) size(Val_Mat_type_3,2)];

w=0.25;
h=0.39;
bott=0.09;
left=0.075;

posits=[left bott w h;
    left+(1.15*w) bott w h;
    left+2*(1.15*w) bott w h];

ms=25;
lw=1.3;

for j=1:3 %loop through 3 types of stim sequences
    
    eval(['Y=Stim_seq_type_' num2str(j) ';'])
    
    if j==1
        yTL={'Free','Non-Arousing','Arousing'};
    else
        yTL='';
    end

    ax=axes('parent',fig1,'Position',posits(j,:),'TickDir','out',...
        'ytick',1:3,'yticklabel',yTL,'ticklength',[0.04 0],...
        'xtick',[1 5 10 15 20],...
        'FontSize',fsize,'FontName',fname);
    hold(ax,"on")
    ytickangle(90)
    
   % 1 = IMG NEUTRAL, NO FOCUS
   % 2 = IMG NEUTRAL, FOCUS NON-AROUSING
    % 3 = IMG UNPLEASANT, NO FOCUS
    % 4 = IMG UNPLEASANT, FOCUS NON-AROUSING
    % 5 = IMG UNPLEASANT, FOCUS AROUSING
    

    Y2=Y;
    Y2(Y2==3)=1; %make 3s into 1s (both are focus free)
    Y2(Y2==4)=2; %make 4s into 2s (both are non-arous)
    Y2(Y2==5)=3; %make 5s into 3s

    logi_neutral=Y<3;
    logi_unpleas=Y>=3;

    X=1:numel(Y);
    plot(X,Y2,'LineWidth',lw,'Color',[0.6 0.6 0.6])
    plot(X(logi_neutral),Y2(logi_neutral),'linestyle','none',...
        'marker','.','markersize',ms,'color',[0 0.65 0])
    plot(X(logi_unpleas),Y2(logi_unpleas),'linestyle','none',...
        'marker','.','markersize',ms,'color',[0.6 0 0])
    
    
    ylim([0.8 3.3])
    xlim([0.5 20.5])

    if j==1
        ylabel('Focus type','FontWeight','bold','FontSize',fsize+3)
        text(1, 3.2,['n=' num2str(ns(j)) ' particip.'],'FontName',fname,'FontSize',fsize-2)
    elseif j==2
        xlabel('Trial','FontWeight','bold','FontSize',fsize+3)
        text(1, 3.2,['n=' num2str(ns(j))],'FontName',fname,'FontSize',fsize-2)
    elseif j==3
        text(18.6,1.2,'Neutral','FontSize',fsize-2,'FontWeight','bold',...
            'color',[0 0.65 0],'FontName',fname)
        text(18.6,1,'Unpleasant','FontSize',fsize-2,'FontWeight','bold',...
            'color',[0.65 0 0],'FontName',fname)    
        text(1, 3.2,['n=' num2str(ns(j))],'FontName',fname,'FontSize',fsize-2)
    end

    
    title(['Sequence ' num2str(j)],'FontWeight','normal')
end



annotation(fig1,'textbox', [0.005 0.48 lettW lettH],'String','D','fontweight','bold',...
    'HorizontalAlignment','center','FontSize',fsize+26,'FitBoxToText','off','EdgeColor','none',...
    'fontname',fname);


set(fig1,'PaperUnits','centimeters')
set(fig1, 'PaperPosition', [0 0 34.9 24.3])


figsPath='/Users/danielrojaslibano/Library/CloudStorage/GoogleDrive-dirl75@gmail.com/My Drive/Matlab_Daniel/Articles/Attentional Deployment/Figures_article';

print(fig1,[figsPath filesep 'FIG_01CD_AttDep2023'],'-dpng','-r330')


