%% relevant paths & data files

%run script to get behav and eyetracking data file pairs
map_asc2mat_AttDeployment

%use only participants with both kinds of data
%MasterTable_Files=MasterTable_Files(MasterTable_Files.overall,:);


% specify data directory
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

Age=nan(nParticip,1);

for j=1:nParticip
    
    load([behavData_path filesep MasterTable_Files.mat{j}])
    
    if strcmpi(BehavData.info.Subject_Gender,'M')
    else
        female(j)=true;
    end
    
    Age(j)=str2double( BehavData.info.Subject_Age);

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
disp(['Age: Mean = ' num2str(mean(Age,'omitnan')) '; SD = ' num2str(std(Age,'omitnan')) ' years' ])

% 1 = IMG NEUTRAL, NO FOCUS
% 2 = IMG NEUTRAL, FOCUS NON-AROUSING
% 3 = IMG UNPLEASANT, NO FOCUS
% 4 = IMG UNPLEASANT, FOCUS NON-AROUSING
% 5 = IMG UNPLEASANT, FOCUS AROUSING
val_unpleas=nan(nParticip,1);
val_neut=nan(nParticip,1);
int_unpleas=nan(nParticip,1);
int_neut=nan(nParticip,1);


for j=1:nParticip
    
    i=Int_Mat(:,j);
    v=Val_Mat(:,j);
    s=Stim_Mat(:,j);
    
    val_unpleas(j)=mean(v(s>=3),'omitnan');
    val_neut(j)=mean(v(s<3),'omitnan');
    int_unpleas(j)=mean(i(s>=3),'omitnan');
    int_neut(j)=mean(i(s<3),'omitnan');
end



countNan_Int=sum(isnan(Int_Mat));
countNan_Val=sum(isnan(Val_Mat));

%% display stats

disp(['N of participants with complete valence responses = ' num2str(sum(countNan_Val==0))])
disp(['N of participants with complete intensity responses = ' num2str(sum(countNan_Int==0))])

disp('N of omitted ratings, for participants with incomplete ratings:')

temp=sum(vertcat(countNan_Int,countNan_Val));
temp=temp(temp~=0);
disp(['1 rating omitted: ' num2str(sum(temp==1)) ' participants'])
disp(['2 ratings omitted: ' num2str(sum(temp==2)) ' participants'])
disp(['3 ratings omitted: ' num2str(sum(temp==3)) ' participants'])
disp(['4 ratings omitted: ' num2str(sum(temp==4)) ' participants'])
disp(['5 ratings omitted: ' num2str(sum(temp==5)) ' participants'])




%% plot


fsize=12;
fname='Verdana';

left=0.15;
bott1=0.66;
bott0=0.11;
w=0.35;
h=0.39;
posits=[left bott1 w*0.8 h*0.75;
    left+1.1*w bott1 w h*0.75;
    left bott0 w*0.95 h;
    left+1.21*w bott0 w*0.95 h];


ms=17;
lw=5;

xlimi_rt=[-0.1 5.3];

ylimi_intval=[0.7 9.1];
xlimi_intval=[-0.5 1.5];

neutral_col=[0 0.6 0];
unpleas_col=[0.7 0 0];

fig1=figure;

ax=axes('parent',fig1,'tickdir','out','position',posits(1,:),...
    'fontsize',fsize,'fontname', fname);
hold(ax,'on')
H=histogram(sessDurat_min,12);
H.FaceColor=[0.5 0.33 0.08];
xlabel('Session duration (min)')
ylabel('Number of particip.')

text(24,15,['n = ' num2str(numel(sessDurat_min)) ' participants'],...
    'fontsize',fsize-4,'FontName',fname)

ylim([1 22.5])

ax2=axes('parent',fig1,'tickdir','out','position',posits(2,:),...
    'fontsize',fsize,'fontname', fname);
hold(ax2,'on')

data4hist=vertcat(RT_Val_Mat(:),RT_Int_Mat(:));

% data4hist=data4hist(data4hist<=5);


H=histogram(data4hist,20);
H.Normalization='probability';
H.FaceColor=[0.9 0.73 0.12];
ylabel('Probability')
xlabel('Response time (s)')

xlim(xlimi_rt);

    text(0.1,0.09,['n = ' num2str(sum(~isnan((data4hist)))) ' ratings'],...
    'fontsize',fsize-4,'FontName',fname)


disp('Response Time:')
disp(['median = ' num2str(median(data4hist, 'omitnan')) ' seconds'])
disp(['mean = ' num2str(mean(data4hist, 'omitnan')) ' seconds'])
disp(['s = ' num2str(std(data4hist, 'omitnan')) ' seconds'])
disp(['n = ' num2str(sum(~isnan((data4hist)))) ' ratings'])

yyaxis right
[F, X] = ecdf(data4hist);
plot(X,F,'linewidth',3)
ylabel('Cumulat. probability')
ylim([0 1.01])



ax3=axes('parent',fig1,'tickdir','out','position',posits(3,:),...
    'xtick',[0 1], 'xticklabel',{'Neutral','Unpleasant'},...
    'ytick',1:2:9,...
    'fontsize',fsize,'fontname', fname);
hold(ax3,'on')


% 
% plot([-0.2 0.2],[mean(int_neut) mean(int_neut)],'linewidth',lw,...
%     'color',[0 0.4 0])
% plot([0.8 1.2],[mean(int_unpleas) mean(int_unpleas)],'linewidth',lw,...
%     'color',[0.4 0 0])


X=randn(nParticip,1)./16;

for j=1:nParticip
    
    plot([X(j) X(j)+1], [int_neut(j) int_unpleas(j)],...
        'color',[0.5 0.5 0.5]) %plot line joining markers
end

for j=1:nParticip
    plot(X(j), int_neut(j),...
        'marker','.','markersize',ms,'color',neutral_col) %plot markers
    
    plot(X(j)+1, int_unpleas(j),...
        'marker','.','markersize',ms,'color',unpleas_col) %plot markers    
end

[box,low_whisk,high_whisk,data2plot]=data4boxplot(int_neut);

plot([-0.4 -0.2],[box(1) box(1)],'linewidth',lw-2,'color',neutral_col)
plot([-0.4 -0.2],[box(2) box(2)],'linewidth',lw-2,'color',neutral_col)
plot([-0.4 -0.2],[box(3) box(3)],'linewidth',lw-2,'color',neutral_col)
plot([-0.3 -0.3],[box(1) low_whisk],'linewidth',lw-2,'color',neutral_col)
plot([-0.3 -0.3],[box(3) high_whisk],'linewidth',lw-2,'color',neutral_col)
plot([-0.39 -0.39],[box(1) box(3)],'linewidth',lw-2,'color',neutral_col)
plot([-0.21 -0.21],[box(1) box(3)],'linewidth',lw-2,'color',neutral_col)



[box,low_whisk,high_whisk,data2plot]=data4boxplot(int_unpleas);

plot([1.2 1.4],[box(1) box(1)],'linewidth',lw-2,'color',unpleas_col)
plot([1.2 1.4],[box(2) box(2)],'linewidth',lw-2,'color',unpleas_col)
plot([1.2 1.4],[box(3) box(3)],'linewidth',lw-2,'color',unpleas_col)
plot([1.3 1.3],[box(1) low_whisk],'linewidth',lw-2,'color',unpleas_col)
plot([1.3 1.3],[box(3) high_whisk],'linewidth',lw-2,'color',unpleas_col)
plot([1.21 1.21],[box(1) box(3)],'linewidth',lw-2,'color',unpleas_col)
plot([1.39 1.39],[box(1) box(3)],'linewidth',lw-2,'color',unpleas_col)

xlim(xlimi_intval)
ylim(ylimi_intval)
xlabel('Stimuli type')
ylabel('Intensity rating')

ax4=axes('parent',fig1,'tickdir','out','position',posits(4,:),...
    'xtick',[0 1], 'xticklabel',{'Neutral','Unpleasant'},...
    'ytick',1:2:9,...
    'fontsize',fsize,'fontname', fname);
hold(ax4,'on')

% plot([-0.18 0.22],[mean(val_neut) mean(val_neut)],'linewidth',lw,...
%     'color',[0 0.4 0])
% plot([0.78 1.22],[mean(val_unpleas) mean(val_unpleas)],'linewidth',lw,...
%     'color',[0.4 0 0])

X=randn(nParticip,1)./17;

for j=1:nParticip
    
    plot([X(j) X(j)+1], [val_neut(j) val_unpleas(j)],...
        'color',[0.5 0.5 0.5]) %plot line joining markers
end

for j=1:nParticip
    plot(X(j), val_neut(j),...
        'marker','.','markersize',ms,'color',neutral_col) %plot markers
    
    plot(X(j)+1, val_unpleas(j),...
        'marker','.','markersize',ms,'color',unpleas_col) %plot markers
    
end

[box,low_whisk,high_whisk,data2plot]=data4boxplot(val_neut);

plot([-0.4 -0.2],[box(1) box(1)],'linewidth',lw-2,'color',neutral_col)
plot([-0.4 -0.2],[box(2) box(2)],'linewidth',lw-2,'color',neutral_col)
plot([-0.4 -0.2],[box(3) box(3)],'linewidth',lw-2,'color',neutral_col)
plot([-0.3 -0.3],[box(1) low_whisk],'linewidth',lw-2,'color',neutral_col)
plot([-0.3 -0.3],[box(3) high_whisk],'linewidth',lw-2,'color',neutral_col)
plot([-0.39 -0.39],[box(1) box(3)],'linewidth',lw-2,'color',neutral_col)
plot([-0.21 -0.21],[box(1) box(3)],'linewidth',lw-2,'color',neutral_col)



[box,low_whisk,high_whisk,data2plot]=data4boxplot(val_unpleas);

plot([1.2 1.4],[box(1) box(1)],'linewidth',lw-2,'color',unpleas_col)
plot([1.2 1.4],[box(2) box(2)],'linewidth',lw-2,'color',unpleas_col)
plot([1.2 1.4],[box(3) box(3)],'linewidth',lw-2,'color',unpleas_col)
plot([1.3 1.3],[box(1) low_whisk],'linewidth',lw-2,'color',unpleas_col)
plot([1.3 1.3],[box(3) high_whisk],'linewidth',lw-2,'color',unpleas_col)
plot([1.21 1.21],[box(1) box(3)],'linewidth',lw-2,'color',unpleas_col)
plot([1.39 1.39],[box(1) box(3)],'linewidth',lw-2,'color',unpleas_col)




xlim(xlimi_intval)
ylim(ylimi_intval)
xlabel('Stimuli type')
ylabel('Valence rating')


lettW=0.07;
lettH=0.07;

annotation(fig1,'textbox', [0.005 0.93 lettW lettH],'String','A','fontweight','bold',...
    'HorizontalAlignment','center','FontSize',fsize+14,'FitBoxToText','off','EdgeColor','none',...
    'fontname',fname);

annotation(fig1,'textbox', [0.008 0.5 lettW lettH],'String','B','fontweight','bold',...
    'HorizontalAlignment','center','FontSize',fsize+14,'FitBoxToText','off','EdgeColor','none',...
    'fontname',fname);




set(fig1,'PaperUnits','inches')
set(fig1, 'PaperPosition', [0 0 7.5 5.5])


figsPath='/Users/danielrojaslibano/Library/CloudStorage/GoogleDrive-dirl75@gmail.com/My Drive/Matlab_Daniel/Articles/Attentional Deployment/Figures_article';

print(fig1,[path2saveFigs filesep 'FIG_02'],'-dpng','-r450')



%% STATS
home
% Intensity: unpleasant vs neutral 
n_perms=1500;
statistic='tstat';
tails='both';

disp('%%%%%% INTENSITY: UNPLEAS VS NEUTRAL')
disp('Intensity: unpleasant images')
disp(['Mdn = ' num2str(median(int_unpleas)) ', M = ' num2str(mean(int_unpleas)) ', SD = ' num2str(std(int_unpleas)) ])
disp('Intensity: neutral images')
disp(['Mdn = ' num2str(median(int_neut)) ', M = ' num2str(mean(int_neut)) ', SD = ' num2str(std(int_neut)) ])
disp('Results of t-stat based permutation test:')
[pval, tstat_no_perm, tstat_dof]=PermutePairedTest(int_unpleas, int_neut, n_perms, statistic, tails);

if pval < 0.01
    pTxt= ', p < 0.01';
else
    pTxt=[', ' num2str(pval)];
end

disp(['t(' num2str(tstat_dof) ') = ' num2str(tstat_no_perm) pTxt])

%[h,p,ci,tstats] =ttest(int_unpleas, int_neut);

gstats=mes(int_unpleas, int_neut,'hedgesg','isDep',1,'nBoot',1000);

disp(['g = ' num2str(gstats.hedgesg)])
disp(['g 95% CI = [' num2str(gstats.hedgesgCi(1)) ', ' num2str(gstats.hedgesgCi(2)) ']'])


disp('  ')
disp('%%%%%% VALENCE: UNPLEAS VS NEUTRAL')
disp('Valence: unpleasant images')
disp(['Mdn = ' num2str(median(val_unpleas)) ', M = ' num2str(mean(val_unpleas)) ', SD = ' num2str(std(val_unpleas)) ])
disp('Valence: neutral images')
disp(['Mdn = ' num2str(median(val_neut)) ', M = ' num2str(mean(val_neut)) ', SD = ' num2str(std(val_neut)) ])
disp('Results of t-stat based permutation test:')

[pval, tstat_no_perm, tstat_dof]=PermutePairedTest(val_unpleas, val_neut, n_perms, statistic, tails);
if pval < 0.01
    pTxt= ', p < 0.01';
else
    pTxt=[', ' num2str(pval)];
end

disp(['t(' num2str(tstat_dof) ') = ' num2str(tstat_no_perm) pTxt])


gstats=mes(val_unpleas, val_neut,'hedgesg','isDep',1,'nBoot',1000);

disp(['g = ' num2str(gstats.hedgesg)])
disp(['g 95% CI = [' num2str(gstats.hedgesgCi(1)) ', ' num2str(gstats.hedgesgCi(2)) ']'])


% Valence: unpleasant vs neutral
