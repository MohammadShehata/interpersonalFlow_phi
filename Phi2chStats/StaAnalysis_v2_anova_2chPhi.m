
%%%%%%%%%%%%%%%%%%%%%%%
%%% Inter-brain Phi %%%
%%%%%%%%%%%%%%%%%%%%%%%
input = ('../phi3/data/');
output = 'output/';

load([input '2ch_diff_perSong_phi3.mat']);

%%% Adjust 2-channel Phi matrix:

% Change zeros to NaN
phis_perSong(:,1,:,:,6)=NaN;
phis_perSong(:,10,:,:,6)=NaN;

% Average across songs
phis_perSub = squeeze(nanmean(phis_perSong,5));

% Mark participant's groups
WantPart1Phi2ch = (1:14)';
WantPart2Phi2ch = (15:28)';
Part1Indx = ismember(networks, WantPart1Phi2ch);
Part2Indx = ismember(networks, WantPart2Phi2ch); 

for pair = 1:10
  
    % Get Inter-brain Phi: 
    InterPartIndx = find(Part1Indx(:,1) .* Part1Indx(:,2) == 0 & Part1Indx(:,1) == 1);
    InterMat1 = networks(InterPartIndx,1);
    InterNetworks(:,:,1) = reshape(InterMat1,14,14);
    InterMat2 = networks(InterPartIndx,2);
    InterNetworks(:,:,2) = reshape(InterMat2,14,14);

    tempInter_raw(:,pair,:,:) = phis_perSub(InterPartIndx,pair,:,:);
end

    Inter_raw = reshape(tempInter_raw,14,14,10,3,4); %groupA x groupB x pair x condition x tau

    % Rearrange groups into 1-7 L; 1-7 R:
    Inter_raw = Inter_raw([3 8 5 12 10 6 1 4 9 11 13 14 7 2],[3 8 5 12 10 6 1 4 9 11 13 14 7 2],:,:,:);
    NewLabel ={'Gp1_L';'Gp2_L';'Gp3_L';'Gp4_L';'Gp5_L';'Gp6_L';'Gp7_L';'Gp1_R';'Gp2_R';'Gp3_R';'Gp4_R';'Gp5_R';'Gp6_R';'Gp7_R'};
    
for pair = 1:10
    P1 = (pair*2)-1;
    P2 = pair*2;
    
    Phi2ch.InterPhi_raw(:,:,P1,:,:) = Inter_raw(:,:,pair,:,:); %groupA x groupB x participant x condition x tau
    I = Inter_raw(:,:,pair,:,:);
    Phi2ch.InterPhi_raw(:,:,P2,:,:) = rot90(fliplr(I),1); %groupA x groupB x participant x condition x tau
    
end

    % Calculate baseline Phi (mean across conditions):
    Phi2ch.allcondavg.InterPhi_base = squeeze(mean(Phi2ch.InterPhi_raw,4)); %groupA x groupB x participant x tau       

    cond={'tm' 'tmrv' 'tmoc'};

for c=1:length(cond)
    
    temp_InterPhi = squeeze(Phi2ch.InterPhi_raw(:,:,:,c,:)); %groupA x groupB x participant x tau
    temp_InterPhi_base = Phi2ch.allcondavg.InterPhi_base; %groupA x groupB x participant x tau
    
    % Average tau 20ms and 40ms
    temp_InterPhi24 = squeeze(mean(temp_InterPhi(:,:,:,1:2),4)); %groupA x groupB x participant
    temp_InterPhi24_base = squeeze(mean(temp_InterPhi_base(:,:,:,1:2),4)); %groupA x groupB x participant
    
    % Average tau 20ms, 40ms, and 60ms
    temp_InterPhi246 = squeeze(mean(temp_InterPhi(:,:,:,1:3),4)); %groupA x groupB x participant
    temp_InterPhi246_base = squeeze(mean(temp_InterPhi_base(:,:,:,1:3),4)); %groupA x groupB x participant

    Phi2ch.raw_BL.(cond{c}).InterPhi24_raw_BL = bsxfun(@minus,temp_InterPhi24,temp_InterPhi24_base); %groupA x groupB x participant
    Phi2ch.raw_BL.(cond{c}).InterPhi246_raw_BL = bsxfun(@minus,temp_InterPhi246,temp_InterPhi246_base); %groupA x groupB x participant
    
    Phi2ch.mean_BL.(cond{c}).InterPhi24_mean_BL = squeeze(mean(Phi2ch.raw_BL.(cond{c}).InterPhi24_raw_BL,3)); %groupA x groupB condition mean
    Phi2ch.mean_BL.(cond{c}).InterPhi246_mean_BL = squeeze(mean(Phi2ch.raw_BL.(cond{c}).InterPhi246_raw_BL,3)); %groupA x groupB condition mean
    
    % Calculate a Grand inter-brain Phi2ch (average across all possible GP-GP connection):      
    Phi2ch.Grandraw_BL.(cond{c}).InterPhi24_Grandraw_BL = squeeze(mean(mean(Phi2ch.raw_BL.(cond{c}).InterPhi24_raw_BL,1),2)); % Grand inter-brain Phi2ch participant
    Phi2ch.Grandraw_BL.(cond{c}).InterPhi246_Grandraw_BL = squeeze(mean(mean(Phi2ch.raw_BL.(cond{c}).InterPhi246_raw_BL,1),2)); % Grand inter-brain Phi2ch participant
        
    Phi2ch.Grandmean_BL.(cond{c}).InterPhi24_Grandmean_BL = squeeze(mean(mean(Phi2ch.mean_BL.(cond{c}).InterPhi24_mean_BL,1),2)); % Grand inter-brain Phi2ch condition mean
    Phi2ch.Grandmean_BL.(cond{c}).InterPhi246_Grandmean_BL = squeeze(mean(mean(Phi2ch.mean_BL.(cond{c}).InterPhi246_mean_BL,1),2)); % Grand inter-brain Phi2ch condition mean

    % Calculate a Grand "L-hemisphere" inter-brain Phi2ch (average across all possible GP-GP connection):      
    tempInterPhi24_raw_BL_Lhemi = Phi2ch.raw_BL.(cond{c}).InterPhi24_raw_BL(1:7,1:7,:);
    Phi2ch.Grandraw_BL_Lhemi.(cond{c}).InterPhi24_Grandraw_BL_Lhemi = squeeze(mean(mean(tempInterPhi24_raw_BL_Lhemi,1),2)); % Grand "L-hemisphere" inter-brain Phi2ch participant
    tempInterPhi246_raw_BL_Lhemi = Phi2ch.raw_BL.(cond{c}).InterPhi246_raw_BL(1:7,1:7,:);
    Phi2ch.Grandraw_BL_Lhemi.(cond{c}).InterPhi246_Grandraw_BL_Lhemi = squeeze(mean(mean(tempInterPhi246_raw_BL_Lhemi,1),2)); % Grand "L-hemisphere" inter-brain Phi2ch participant
        
    % Calculate a Grand "R-hemisphere" inter-brain Phi2ch (average across all possible GP-GP connection):      
    tempInterPhi24_raw_BL_Rhemi = Phi2ch.raw_BL.(cond{c}).InterPhi24_raw_BL(8:14,8:14,:);
    Phi2ch.Grandraw_BL_Rhemi.(cond{c}).InterPhi24_Grandraw_BL_Rhemi = squeeze(mean(mean(tempInterPhi24_raw_BL_Rhemi,1),2)); % Grand "R-hemisphere" inter-brain Phi2ch participant
    tempInterPhi246_raw_BL_Rhemi = Phi2ch.raw_BL.(cond{c}).InterPhi246_raw_BL(8:14,8:14,:);
    Phi2ch.Grandraw_BL_Rhemi.(cond{c}).InterPhi246_Grandraw_BL_Rhemi = squeeze(mean(mean(tempInterPhi246_raw_BL_Rhemi,1),2)); % Grand "R-hemisphere" inter-brain Phi2ch participant
        
    % Calculate a Grand "inter-hemisphere" inter-brain Phi2ch (average across all possible GP-GP connection):      
    tempInterPhi24_raw_BL_Interhemi = Phi2ch.raw_BL.(cond{c}).InterPhi24_raw_BL(1:7,8:14,:);
    Phi2ch.Grandraw_BL_Interhemi.(cond{c}).InterPhi24_Grandraw_BL_Interhemi = squeeze(mean(mean(tempInterPhi24_raw_BL_Interhemi,1),2)); % Grand "inter-hemisphere" inter-brain Phi2ch participant
    tempInterPhi246_raw_BL_Interhemi = Phi2ch.raw_BL.(cond{c}).InterPhi246_raw_BL(1:7,8:14,:);
    Phi2ch.Grandraw_BL_Interhemi.(cond{c}).InterPhi246_Grandraw_BL_Interhemi = squeeze(mean(mean(tempInterPhi246_raw_BL_Interhemi,1),2)); % Grand "inter-hemisphere" inter-brain Phi2ch participant
      
end


%save ([output 'Phi2ch.mat'],'Phi2ch');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Processing 20, 40, and 60 ms Tau average
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotting 

output=('../STAT_Phi2ch/');
NewLabel ={'Gp1_L';'Gp2_L';'Gp3_L';'Gp4_L';'Gp5_L';'Gp6_L';'Gp7_L';'Gp1_R';'Gp2_R';'Gp3_R';'Gp4_R';'Gp5_R';'Gp6_R';'Gp7_R'};
cond={'tm' 'tmrv' 'tmoc'};

subp = 1;
    for c=1:length(cond)

            % plot inter-brain Phi2ch
            figure(1);subplot(1,3,subp);
            lowertri_InterPhi246 = tril(Phi2ch.mean_BL.(cond{c}).InterPhi246_mean_BL);
            imagesc(lowertri_InterPhi246)
            set(gca,'clim',[-0.00009 .00009],'xtick',1:size(NewLabel,1),'ytick',1:size(NewLabel,1),'yticklabel',NewLabel);
            line([7.5 7.5],get(gca,'Ylim'),'Color',[0 1 0],'LineWidth',2)
            line(get(gca,'Xlim'),[7.5 7.5],'Color',[0 1 0],'LineWidth',2)
            axis square
            colorbar
            title([cond{c}  '_InterPhi246_mean_BL_20, 40, and 60ms tau'],'interpreter','none');
            
            subp = subp + 1; 
    end
    
    set(1,'Position', [50, 50, 1600, 400]);
    
    saveas(1,[output 'InterPhi246_mean_BL.png']); saveas(1,[output 'InterPhi246_mean_BL.fig']);

    
%% ANOVA & post-hoc across team/team occulusion/team reverse (inter Phi2ch)
clear;
close all;

input = ('/Users/mohammadshehata/Desktop/MoShehata/1- Projects/F3- Angus-Nao IIT project/Angus/');
output=('/Users/mohammadshehata/Desktop/MoShehata/1- Projects/F3- Angus-Nao IIT project/Angus/STAT_Phi2ch/');

cond={'tm' 'tmrv' 'tmoc'};
posthocType=('bonferroni'); % set posthoc test type: 'hsd' (tukey), 'lsd','bonferroni'

NewLabel ={'Gp1_L';'Gp2_L';'Gp3_L';'Gp4_L';'Gp5_L';'Gp6_L';'Gp7_L';'Gp1_R';'Gp2_R';'Gp3_R';'Gp4_R';'Gp5_R';'Gp6_R';'Gp7_R'};

           
    load([input 'Phi2ch.mat']);
    
    c1=Phi2ch.raw_BL.tm.InterPhi246_raw_BL;
    c2=Phi2ch.raw_BL.tmrv.InterPhi246_raw_BL;
    c3=Phi2ch.raw_BL.tmoc.InterPhi246_raw_BL;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute anova and posthoc test  %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Data arrangment for 3-way anova
            anovamat(:,:,:,1)=c1;
            anovamat(:,:,:,2)=c2;
            anovamat(:,:,:,3)=c3; %groupA x groupB x participant x condition
            
            anovamat1 = permute(anovamat,[3 4 1 2]); %participant x condition x groupA x groupB
    
            anovamat3.Y=reshape(anovamat1,[numel(anovamat1),1]);
            anovamat3.S=repmat([1:20]',588,1); %20subject X (14ROI1 x 14ROI1 x 3conditions)
            F1(1:20,1)={'tm'}; F1(21:40,1)={'tmrv'}; F1(41:60,1)={'tmoc'}; %
            anovamat3.G1=repmat(F1,196,1); %3conditions x (20subject x 14ROI1 X 14ROI2)
            
            for roi1 = 1:14
            F2([1:60],roi1)=roi1; %(20subject x 3conditions)
            end
            F2=reshape(F2,[numel(F2),1]); %14ROI X (20subject x 3conditions)
            anovamat3.G2=repmat(F2,14,1); %14ROI X (20subject x 3conditions) x 14ROI2
            
            for roi1 = 1:14
            F3([1:840],roi1)=roi1; %14ROI X (20subject x 3conditions)
            end
            anovamat3.G3=reshape(F3,[numel(F3),1]);%14ROI X (20subject x 3conditions x 14ROI2)
            
            anovamat3.Gnames={'subject','condition','ROI1','ROI2'};
            
            % Three-way ANOVA
            [p,tbl,stats] = anovan(anovamat3.Y,{anovamat3.S anovamat3.G1 anovamat3.G2 anovamat3.G3},'model','interaction','varnames',anovamat3.Gnames);
            
            INTERPhi2chStat3way.anovamat = anovamat1;
            INTERPhi2chStat3way.stattbl = tbl;
            
            %_________ posthoc
            [c,m,h,gnames]= multcompare(stats,'Dimension',[2 3 4],'CType',posthocType,'Display','off');
            INTERPhi2chStat3way.posthoc123 = c;
            INTERPhi2chStat3way.gnames = gnames;
            
            % Collect the Significance Table across condition (AXCond)
            CCR = 1;
            for roi1 = 1:14
                for roi2 = 1:14
                ROI1num = num2str(roi1);
                ROI2num = num2str(roi2);
                Index1 = find(strcmp(gnames,['condition=tm,ROI1=' ROI1num ',ROI2=' ROI2num]));
                Index2 = find(strcmp(gnames,['condition=tmrv,ROI1=' ROI1num ',ROI2=' ROI2num]));
                IndexF1 = find(c(:,1)==Index1 & c(:,2)==Index2);
                Sign1 = c(IndexF1,6);
                Index3 = find(strcmp(gnames,['condition=tm,ROI1=' ROI1num ',ROI2=' ROI2num]));
                Index4 = find(strcmp(gnames,['condition=tmoc,ROI1=' ROI1num ',ROI2=' ROI2num]));
                IndexF2 = find(c(:,1)==Index3 & c(:,2)==Index4);
                Sign2 = c(IndexF2,6);
                Index5 = find(strcmp(gnames,['condition=tmrv,ROI1=' ROI1num ',ROI2=' ROI2num]));
                Index6 = find(strcmp(gnames,['condition=tmoc,ROI1=' ROI1num ',ROI2=' ROI2num]));
                IndexF3 = find(c(:,1)==Index5 & c(:,2)==Index6);
                Sign3 = c(IndexF3,6);
                INTERPhi2chStat3way.SignTblAXCond(CCR,:) = [roi1,roi2,Sign1,Sign2,Sign3];
                CCR = CCR + 1;
                end 
            end
            
            % Collect the Significance Table across ROI(AXROI)
            CCS = 1;
            for roi1 = 1:14
                for roi2 = 1:14
                    ROI1num = num2str(roi1);
                    ROI2num = num2str(roi2);
                    Indx1 = find(strcmp(gnames,['condition=tm,ROI1=' ROI1num ',ROI2=' ROI2num]));
                    Indx5 = find(strcmp(gnames,['condition=tmrv,ROI1=' ROI1num ',ROI2=' ROI2num]));
                    Indx9 = find(strcmp(gnames,['condition=tmoc,ROI1=' ROI1num ',ROI2=' ROI2num]));
                    
                        for Roi1 = 1:14
                           for Roi2 = 1:14
                                    Roi1num = num2str(Roi1);
                                    Roi2num = num2str(Roi2);
                                    Indx2 = find(strcmp(gnames,['condition=tm,ROI1=' Roi1num ',ROI2=' Roi2num]));
                                    IndxR1 = find(c(:,1)==Indx1 & c(:,2)==Indx2);
                                    SignR1 = c(IndxR1,6);
                                    if isempty(SignR1); SignR1 = 17; end

                                    Indx6 = find(strcmp(gnames,['condition=tmrv,ROI1=' Roi1num ',ROI2=' Roi2num]));
                                    IndxR3 = find(c(:,1)==Indx5 & c(:,2)==Indx6);
                                    SignR3 = c(IndxR3,6);
                                    if isempty(SignR3); SignR3 = 17; end

                                    Indx10 = find(strcmp(gnames,['condition=tmoc,ROI1=' Roi1num ',ROI2=' Roi2num]));
                                    IndxR5 = find(c(:,1)==Indx9 & c(:,2)==Indx10);
                                    SignR5 = c(IndxR5,6);
                                    if isempty(SignR5); SignR5 = 17; end

                                    INTERPhi2chStat3way.SignTblAXROI(CCS,:) = [roi1,roi2,Roi1,Roi2,SignR1,SignR3,SignR5];
                                    CCS = CCS+1;
                           end
                        end
                end
            end
                                    INTERPhi2chStat3way.SignTblAXROI(INTERPhi2chStat3way.SignTblAXROI == 17) = NaN;
            %_______________
            
             save ([output 'INTERPhi2chStat3way.mat'],'INTERPhi2chStat3way');
             
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% One-way ANOVA for the Grand inter-brain Phi 2 channel: 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear;
    close all;

    input = ('/Users/mohammadshehata/Desktop/MoShehata/1- Projects/F3- Angus-Nao IIT project/Angus/');
    output=('/Users/mohammadshehata/Desktop/MoShehata/1- Projects/F3- Angus-Nao IIT project/Angus/STAT_Phi2ch/');

    cond={'tm' 'tmrv' 'tmoc'};
    posthocType=('bonferroni'); % set posthoc test type: 'hsd' (tukey), 'lsd','bonferroni'

    NewLabel ={'Gp1_L';'Gp2_L';'Gp3_L';'Gp4_L';'Gp5_L';'Gp6_L';'Gp7_L';'Gp1_R';'Gp2_R';'Gp3_R';'Gp4_R';'Gp5_R';'Gp6_R';'Gp7_R'};
        
    load([input 'Phi2ch.mat']);
    load([output 'INTERPhi2chStat3way.mat']);
    
    anovamat(:,1) = Phi2ch.Grandraw_BL.tm.InterPhi246_Grandraw_BL;
    anovamat(:,2) = Phi2ch.Grandraw_BL.tmrv.InterPhi246_Grandraw_BL;
    anovamat(:,3) = Phi2ch.Grandraw_BL.tmoc.InterPhi246_Grandraw_BL;
            
            % One-way ANOVA
            [p,tbl,stats] = anova1(anovamat);
            
            INTERPhi2chStat3way.GrandPhi2ch_OneWayANOVA.anovamat = anovamat;
            INTERPhi2chStat3way.GrandPhi2ch_OneWayANOVA.stattbl = tbl;
            
            % posthoc
            c = multcompare(stats,'Dimension',[1 2],'CType',posthocType,'Display','off');

            INTERPhi2chStat3way.GrandPhi2ch_OneWayANOVA.posthoc = c;

            save ([output 'INTERPhi2chStat3way.mat'],'INTERPhi2chStat3way');

%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Two-way ANOVA for the Grand "L-, R-, Inter- hemisphere" inter-brain Phi 2 channel: 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear;
    close all;

    input = ('/Users/mohammadshehata/Desktop/MoShehata/1- Projects/F3- Angus-Nao IIT project/Angus/');
    output=('/Users/mohammadshehata/Desktop/MoShehata/1- Projects/F3- Angus-Nao IIT project/Angus/STAT_Phi2ch/');

    cond={'tm' 'tmrv' 'tmoc'};
    posthocType=('bonferroni'); % set posthoc test type: 'hsd' (tukey), 'lsd','bonferroni'

    NewLabel ={'Gp1_L';'Gp2_L';'Gp3_L';'Gp4_L';'Gp5_L';'Gp6_L';'Gp7_L';'Gp1_R';'Gp2_R';'Gp3_R';'Gp4_R';'Gp5_R';'Gp6_R';'Gp7_R'};
        
    load([input 'Phi2ch.mat']);
    load([output 'INTERPhi2chStat3way.mat']);
    
    % Data arrangment for 2-way anova
            anovamat2(1,1:20)= Phi2ch.Grandraw_BL_Lhemi.tm.InterPhi246_Grandraw_BL_Lhemi;
            anovamat2(1,21:40)=Phi2ch.Grandraw_BL_Rhemi.tm.InterPhi246_Grandraw_BL_Rhemi;
            anovamat2(1,41:60)=Phi2ch.Grandraw_BL_Interhemi.tm.InterPhi246_Grandraw_BL_Interhemi;
            
            anovamat2(2,1:20)= Phi2ch.Grandraw_BL_Lhemi.tmrv.InterPhi246_Grandraw_BL_Lhemi;
            anovamat2(2,21:40)=Phi2ch.Grandraw_BL_Rhemi.tmrv.InterPhi246_Grandraw_BL_Rhemi;
            anovamat2(2,41:60)=Phi2ch.Grandraw_BL_Interhemi.tmrv.InterPhi246_Grandraw_BL_Interhemi;
            
            anovamat2(3,1:20)= Phi2ch.Grandraw_BL_Lhemi.tmoc.InterPhi246_Grandraw_BL_Lhemi;
            anovamat2(3,21:40)=Phi2ch.Grandraw_BL_Rhemi.tmoc.InterPhi246_Grandraw_BL_Rhemi;
            anovamat2(3,41:60)=Phi2ch.Grandraw_BL_Interhemi.tmoc.InterPhi246_Grandraw_BL_Interhemi;
            
            % Two-way ANOVA
            [p,tbl,stats] = anova2(anovamat2',20);
            
            INTERPhi2chStat3way.GrandPhi2ch_WInhemi_TwoWayANOVA.anovamat2 = anovamat2';
            INTERPhi2chStat3way.GrandPhi2ch_WInhemi_TwoWayANOVA.stattbl = tbl;
            
            %_________ posthoc
            c = multcompare(stats,'CType',posthocType,'Display','off');

            INTERPhi2chStat3way.GrandPhi2ch_WInhemi_TwoWayANOVA.posthoc = c;

            save ([output 'INTERPhi2chStat3way.mat'],'INTERPhi2chStat3way');

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Two-way ANOVA for the Grand "L-, R-, Inter- hemisphere" inter-brain Phi 2 channel: 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear;
    close all;

    input = ('/Users/mohammadshehata/Desktop/MoShehata/1- Projects/F3- Angus-Nao IIT project/Angus/');
    output=('/Users/mohammadshehata/Desktop/MoShehata/1- Projects/F3- Angus-Nao IIT project/Angus/STAT_Phi2ch/');

    cond={'tm' 'tmrv' 'tmoc'};
    posthocType=('bonferroni'); % set posthoc test type: 'hsd' (tukey), 'lsd','bonferroni'

    NewLabel ={'Gp1_L';'Gp2_L';'Gp3_L';'Gp4_L';'Gp5_L';'Gp6_L';'Gp7_L';'Gp1_R';'Gp2_R';'Gp3_R';'Gp4_R';'Gp5_R';'Gp6_R';'Gp7_R'};
        
    load([input 'Phi2ch.mat']);
    load([output 'INTERPhi2chStat3way.mat']);
    
    % Data arrangment for 2-way anova
            anovamat2(1,1:20)= Phi2ch.Grandraw_BL_Lhemi.tm.InterPhi246_Grandraw_BL_Lhemi;
            anovamat2(1,21:40)=Phi2ch.Grandraw_BL_Lhemi.tmrv.InterPhi246_Grandraw_BL_Lhemi;
            anovamat2(1,41:60)=Phi2ch.Grandraw_BL_Lhemi.tmoc.InterPhi246_Grandraw_BL_Lhemi;
            
            anovamat2(2,1:20)= Phi2ch.Grandraw_BL_Rhemi.tm.InterPhi246_Grandraw_BL_Rhemi;
            anovamat2(2,21:40)=Phi2ch.Grandraw_BL_Rhemi.tmrv.InterPhi246_Grandraw_BL_Rhemi;
            anovamat2(2,41:60)=Phi2ch.Grandraw_BL_Rhemi.tmoc.InterPhi246_Grandraw_BL_Rhemi;
            
            anovamat2(3,1:20)= Phi2ch.Grandraw_BL_Interhemi.tm.InterPhi246_Grandraw_BL_Interhemi;
            anovamat2(3,21:40)=Phi2ch.Grandraw_BL_Interhemi.tmrv.InterPhi246_Grandraw_BL_Interhemi;
            anovamat2(3,41:60)=Phi2ch.Grandraw_BL_Interhemi.tmoc.InterPhi246_Grandraw_BL_Interhemi;
            
            % Two-way ANOVA
            [p,tbl,stats] = anova2(anovamat2',20);
            
            INTERPhi2chStat3way.GrandPhi2ch_AXhemi_TwoWayANOVA.anovamat2 = anovamat2';
            INTERPhi2chStat3way.GrandPhi2ch_AXhemi_TwoWayANOVA.stattbl = tbl;
            
            %_________ posthoc
            c = multcompare(stats,'CType',posthocType,'Display','off');

            INTERPhi2chStat3way.GrandPhi2ch_AXhemi_TwoWayANOVA.posthoc = c;

            save ([output 'INTERPhi2chStat3way.mat'],'INTERPhi2chStat3way');

%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% N-way ANOVA for the Grand "L-, R-, Inter- hemisphere" inter-brain Phi 2 channel: 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear;
    close all;

    input = ('/Users/mohammadshehata/Desktop/MoShehata/1- Projects/F3- Angus-Nao IIT project/Angus/');
    output=('/Users/mohammadshehata/Desktop/MoShehata/1- Projects/F3- Angus-Nao IIT project/Angus/STAT_Phi2ch/');

    cond={'tm' 'tmrv' 'tmoc'};
    posthocType=('bonferroni'); % set posthoc test type: 'hsd' (tukey), 'lsd','bonferroni'

    NewLabel ={'Gp1_L';'Gp2_L';'Gp3_L';'Gp4_L';'Gp5_L';'Gp6_L';'Gp7_L';'Gp1_R';'Gp2_R';'Gp3_R';'Gp4_R';'Gp5_R';'Gp6_R';'Gp7_R'};
        
    load([input 'Phi2ch.mat']);
    load([output 'INTERPhi2chStat3way_246.mat']);    

            % Data arrangment for 3-way anova (participant x hemisphere x condition)
            anovamat(:,1,1)= Phi2ch.Grandraw_BL_Lhemi.tm.InterPhi246_Grandraw_BL_Lhemi;
            anovamat(:,1,2)=Phi2ch.Grandraw_BL_Lhemi.tmrv.InterPhi246_Grandraw_BL_Lhemi;
            anovamat(:,1,3)=Phi2ch.Grandraw_BL_Lhemi.tmoc.InterPhi246_Grandraw_BL_Lhemi;
            
            anovamat(:,2,1)= Phi2ch.Grandraw_BL_Rhemi.tm.InterPhi246_Grandraw_BL_Rhemi;
            anovamat(:,2,2)=Phi2ch.Grandraw_BL_Rhemi.tmrv.InterPhi246_Grandraw_BL_Rhemi;
            anovamat(:,2,3)=Phi2ch.Grandraw_BL_Rhemi.tmoc.InterPhi246_Grandraw_BL_Rhemi;
            
            anovamat(:,3,1)= Phi2ch.Grandraw_BL_Interhemi.tm.InterPhi246_Grandraw_BL_Interhemi;
            anovamat(:,3,2)=Phi2ch.Grandraw_BL_Interhemi.tmrv.InterPhi246_Grandraw_BL_Interhemi;
            anovamat(:,3,3)=Phi2ch.Grandraw_BL_Interhemi.tmoc.InterPhi246_Grandraw_BL_Interhemi;
            
    
            anovamat3.Y=reshape(anovamat,[numel(anovamat),1]);
            anovamat3.S=repmat([1:20]',9,1); %20subject X (3hemisphere x 3conditions)
            F1(1:60,1)={'tm'}; F1(61:120,1)={'tmrv'}; F1(121:180,1)={'tmoc'};
            anovamat3.G1=F1; %3conditions x (20subject x 3 hemispheres)
            F2(1:20,1)={'L-hemi'}; F2(21:40,1)={'R-hemi'}; F2(41:60,1)={'Inter-hemi'};
            anovamat3.G2=repmat(F2,3,1); %3hemispheres x (20subject x 3 conditions)
            
            anovamat3.Gnames={'subject','condition','hemisphere'};
            
            % Three-way ANOVA
            [p,tbl,stats] = anovan(anovamat3.Y,{anovamat3.S anovamat3.G1 anovamat3.G2},'model','interaction','varnames',anovamat3.Gnames);
            
            INTERPhi2chStat3way.GrandPhi2ch_hemi_NWayANOVA.anovamat = anovamat;
            INTERPhi2chStat3way.GrandPhi2ch_hemi_NWayANOVA.stattbl = tbl;
            
            %_________ posthoc
            [c,m,h,gnames]= multcompare(stats,'Dimension',[2 3],'CType',posthocType,'Display','off');
            INTERPhi2chStat3way.GrandPhi2ch_hemi_NWayANOVA.posthoc123 = c;
            INTERPhi2chStat3way.GrandPhi2ch_hemi_NWayANOVA.gnames = gnames;
            
            save ([output 'INTERPhi2chStat3way_246.mat'],'INTERPhi2chStat3way');
            
            
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Processing 20, and 40 ms Tau average
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plotting 

output=('/Users/mohammadshehata/Desktop/MoShehata/1- Projects/F3- Angus-Nao IIT project/Angus/STAT_Phi2ch/');
NewLabel ={'Gp1_L';'Gp2_L';'Gp3_L';'Gp4_L';'Gp5_L';'Gp6_L';'Gp7_L';'Gp1_R';'Gp2_R';'Gp3_R';'Gp4_R';'Gp5_R';'Gp6_R';'Gp7_R'};
cond={'tm' 'tmrv' 'tmoc'};

subp = 1;
    for c=1:length(cond)

            % plot inter-brain Phi2ch
            figure(1);subplot(1,3,subp);
            lowertri_InterPhi24 = tril(Phi2ch.mean_BL.(cond{c}).InterPhi24_mean_BL);
            imagesc(lowertri_InterPhi24)
            set(gca,'clim',[-0.00009 .00009],'xtick',1:size(NewLabel,1),'ytick',1:size(NewLabel,1),'yticklabel',NewLabel);
            line([7.5 7.5],get(gca,'Ylim'),'Color',[0 1 0],'LineWidth',2)
            line(get(gca,'Xlim'),[7.5 7.5],'Color',[0 1 0],'LineWidth',2)
            axis square
            colorbar
            title([cond{c}  '_InterPhi24_mean_BL_20, and 40 ms tau'],'interpreter','none');
            
            subp = subp + 1; 
    end
    
    set(1,'Position', [50, 50, 1600, 400]);
    
    saveas(1,[output 'InterPhi24_mean_BL.png']); saveas(1,[output 'InterPhi24_mean_BL.fig']);

    
%% ANOVA & post-hoc across team/team occulusion/team reverse (inter Phi2ch)
clear;
close all;

input = ('/Users/mohammadshehata/Desktop/MoShehata/1- Projects/F3- Angus-Nao IIT project/Angus/');
output=('/Users/mohammadshehata/Desktop/MoShehata/1- Projects/F3- Angus-Nao IIT project/Angus/STAT_Phi2ch/');

cond={'tm' 'tmrv' 'tmoc'};
posthocType=('bonferroni'); % set posthoc test type: 'hsd' (tukey), 'lsd','bonferroni'

NewLabel ={'Gp1_L';'Gp2_L';'Gp3_L';'Gp4_L';'Gp5_L';'Gp6_L';'Gp7_L';'Gp1_R';'Gp2_R';'Gp3_R';'Gp4_R';'Gp5_R';'Gp6_R';'Gp7_R'};

           
    load([input 'Phi2ch.mat']);
    
    c1=Phi2ch.raw_BL.tm.InterPhi24_raw_BL;
    c2=Phi2ch.raw_BL.tmrv.InterPhi24_raw_BL;
    c3=Phi2ch.raw_BL.tmoc.InterPhi24_raw_BL;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % compute anova and posthoc test  %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            % Data arrangment for 3-way anova
            anovamat(:,:,:,1)=c1;
            anovamat(:,:,:,2)=c2;
            anovamat(:,:,:,3)=c3; %groupA x groupB x participant x condition
            
            anovamat1 = permute(anovamat,[3 4 1 2]); %participant x condition x groupA x groupB
    
            anovamat3.Y=reshape(anovamat1,[numel(anovamat1),1]);
            anovamat3.S=repmat([1:20]',588,1); %20subject X (14ROI1 x 14ROI1 x 3conditions)
            F1(1:20,1)={'tm'}; F1(21:40,1)={'tmrv'}; F1(41:60,1)={'tmoc'}; %
            anovamat3.G1=repmat(F1,196,1); %3conditions x (20subject x 14ROI1 X 14ROI2)
            
            for roi1 = 1:14
            F2([1:60],roi1)=roi1; %(20subject x 3conditions)
            end
            F2=reshape(F2,[numel(F2),1]); %14ROI X (20subject x 3conditions)
            anovamat3.G2=repmat(F2,14,1); %14ROI X (20subject x 3conditions) x 14ROI2
            
            for roi1 = 1:14
            F3([1:840],roi1)=roi1; %14ROI X (20subject x 3conditions)
            end
            anovamat3.G3=reshape(F3,[numel(F3),1]);%14ROI X (20subject x 3conditions x 14ROI2)
            
            anovamat3.Gnames={'subject','condition','ROI1','ROI2'};
            
            % Three-way ANOVA
            [p,tbl,stats] = anovan(anovamat3.Y,{anovamat3.S anovamat3.G1 anovamat3.G2 anovamat3.G3},'model','interaction','varnames',anovamat3.Gnames);
            
            INTERPhi2chStat3way.anovamat = anovamat1;
            INTERPhi2chStat3way.stattbl = tbl;
            
            %_________ posthoc
            [c,m,h,gnames]= multcompare(stats,'Dimension',[2 3 4],'CType',posthocType,'Display','off');
            INTERPhi2chStat3way.posthoc123 = c;
            INTERPhi2chStat3way.gnames = gnames;
            
            % Collect the Significance Table across condition (AXCond)
            CCR = 1;
            for roi1 = 1:14
                for roi2 = 1:14
                ROI1num = num2str(roi1);
                ROI2num = num2str(roi2);
                Index1 = find(strcmp(gnames,['condition=tm,ROI1=' ROI1num ',ROI2=' ROI2num]));
                Index2 = find(strcmp(gnames,['condition=tmrv,ROI1=' ROI1num ',ROI2=' ROI2num]));
                IndexF1 = find(c(:,1)==Index1 & c(:,2)==Index2);
                Sign1 = c(IndexF1,6);
                Index3 = find(strcmp(gnames,['condition=tm,ROI1=' ROI1num ',ROI2=' ROI2num]));
                Index4 = find(strcmp(gnames,['condition=tmoc,ROI1=' ROI1num ',ROI2=' ROI2num]));
                IndexF2 = find(c(:,1)==Index3 & c(:,2)==Index4);
                Sign2 = c(IndexF2,6);
                Index5 = find(strcmp(gnames,['condition=tmrv,ROI1=' ROI1num ',ROI2=' ROI2num]));
                Index6 = find(strcmp(gnames,['condition=tmoc,ROI1=' ROI1num ',ROI2=' ROI2num]));
                IndexF3 = find(c(:,1)==Index5 & c(:,2)==Index6);
                Sign3 = c(IndexF3,6);
                INTERPhi2chStat3way.SignTblAXCond(CCR,:) = [roi1,roi2,Sign1,Sign2,Sign3];
                CCR = CCR + 1;
                end 
            end
            
            % Collect the Significance Table across ROI(AXROI)
            CCS = 1;
            for roi1 = 1:14
                for roi2 = 1:14
                    ROI1num = num2str(roi1);
                    ROI2num = num2str(roi2);
                    Indx1 = find(strcmp(gnames,['condition=tm,ROI1=' ROI1num ',ROI2=' ROI2num]));
                    Indx5 = find(strcmp(gnames,['condition=tmrv,ROI1=' ROI1num ',ROI2=' ROI2num]));
                    Indx9 = find(strcmp(gnames,['condition=tmoc,ROI1=' ROI1num ',ROI2=' ROI2num]));
                    
                        for Roi1 = 1:14
                           for Roi2 = 1:14
                                    Roi1num = num2str(Roi1);
                                    Roi2num = num2str(Roi2);
                                    Indx2 = find(strcmp(gnames,['condition=tm,ROI1=' Roi1num ',ROI2=' Roi2num]));
                                    IndxR1 = find(c(:,1)==Indx1 & c(:,2)==Indx2);
                                    SignR1 = c(IndxR1,6);
                                    if isempty(SignR1); SignR1 = 17; end

                                    Indx6 = find(strcmp(gnames,['condition=tmrv,ROI1=' Roi1num ',ROI2=' Roi2num]));
                                    IndxR3 = find(c(:,1)==Indx5 & c(:,2)==Indx6);
                                    SignR3 = c(IndxR3,6);
                                    if isempty(SignR3); SignR3 = 17; end

                                    Indx10 = find(strcmp(gnames,['condition=tmoc,ROI1=' Roi1num ',ROI2=' Roi2num]));
                                    IndxR5 = find(c(:,1)==Indx9 & c(:,2)==Indx10);
                                    SignR5 = c(IndxR5,6);
                                    if isempty(SignR5); SignR5 = 17; end

                                    INTERPhi2chStat3way.SignTblAXROI(CCS,:) = [roi1,roi2,Roi1,Roi2,SignR1,SignR3,SignR5];
                                    CCS = CCS+1;
                           end
                        end
                end
            end
                                    INTERPhi2chStat3way.SignTblAXROI(INTERPhi2chStat3way.SignTblAXROI == 17) = NaN;
            %_______________
            
             save ([output 'INTERPhi2chStat3way_24.mat'],'INTERPhi2chStat3way');
             
%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% One-way ANOVA for the Grand inter-brain Phi 2 channel: 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear;
    close all;

    input = ('/Users/mohammadshehata/Desktop/MoShehata/1- Projects/F3- Angus-Nao IIT project/Angus/');
    output=('/Users/mohammadshehata/Desktop/MoShehata/1- Projects/F3- Angus-Nao IIT project/Angus/STAT_Phi2ch/');

    cond={'tm' 'tmrv' 'tmoc'};
    posthocType=('bonferroni'); % set posthoc test type: 'hsd' (tukey), 'lsd','bonferroni'

    NewLabel ={'Gp1_L';'Gp2_L';'Gp3_L';'Gp4_L';'Gp5_L';'Gp6_L';'Gp7_L';'Gp1_R';'Gp2_R';'Gp3_R';'Gp4_R';'Gp5_R';'Gp6_R';'Gp7_R'};
        
    load([input 'Phi2ch.mat']);
    load([output 'INTERPhi2chStat3way_24.mat']);
    
    anovamat(:,1) = Phi2ch.Grandraw_BL.tm.InterPhi24_Grandraw_BL;
    anovamat(:,2) = Phi2ch.Grandraw_BL.tmrv.InterPhi24_Grandraw_BL;
    anovamat(:,3) = Phi2ch.Grandraw_BL.tmoc.InterPhi24_Grandraw_BL;
            
            % One-way ANOVA
            [p,tbl,stats] = anova1(anovamat);
            
            INTERPhi2chStat3way.GrandPhi2ch_OneWayANOVA.anovamat = anovamat;
            INTERPhi2chStat3way.GrandPhi2ch_OneWayANOVA.stattbl = tbl;
            
            % posthoc
            c = multcompare(stats,'Dimension',[1 2],'CType',posthocType,'Display','off');

            INTERPhi2chStat3way.GrandPhi2ch_OneWayANOVA.posthoc = c;

            save ([output 'INTERPhi2chStat3way_24.mat'],'INTERPhi2chStat3way');

%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Two-way ANOVA for the Grand "L-, R-, Inter- hemisphere" inter-brain Phi 2 channel: 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear;
    close all;

    input = ('/Users/mohammadshehata/Desktop/MoShehata/1- Projects/F3- Angus-Nao IIT project/Angus/');
    output=('/Users/mohammadshehata/Desktop/MoShehata/1- Projects/F3- Angus-Nao IIT project/Angus/STAT_Phi2ch/');

    cond={'tm' 'tmrv' 'tmoc'};
    posthocType=('bonferroni'); % set posthoc test type: 'hsd' (tukey), 'lsd','bonferroni'

    NewLabel ={'Gp1_L';'Gp2_L';'Gp3_L';'Gp4_L';'Gp5_L';'Gp6_L';'Gp7_L';'Gp1_R';'Gp2_R';'Gp3_R';'Gp4_R';'Gp5_R';'Gp6_R';'Gp7_R'};
        
    load([input 'Phi2ch.mat']);
    load([output 'INTERPhi2chStat3way_24.mat']);
    
    % Data arrangment for 2-way anova
            anovamat2(1,1:20)= Phi2ch.Grandraw_BL_Lhemi.tm.InterPhi24_Grandraw_BL_Lhemi;
            anovamat2(1,21:40)=Phi2ch.Grandraw_BL_Rhemi.tm.InterPhi24_Grandraw_BL_Rhemi;
            anovamat2(1,41:60)=Phi2ch.Grandraw_BL_Interhemi.tm.InterPhi24_Grandraw_BL_Interhemi;
            
            anovamat2(2,1:20)= Phi2ch.Grandraw_BL_Lhemi.tmrv.InterPhi24_Grandraw_BL_Lhemi;
            anovamat2(2,21:40)=Phi2ch.Grandraw_BL_Rhemi.tmrv.InterPhi24_Grandraw_BL_Rhemi;
            anovamat2(2,41:60)=Phi2ch.Grandraw_BL_Interhemi.tmrv.InterPhi24_Grandraw_BL_Interhemi;
            
            anovamat2(3,1:20)= Phi2ch.Grandraw_BL_Lhemi.tmoc.InterPhi24_Grandraw_BL_Lhemi;
            anovamat2(3,21:40)=Phi2ch.Grandraw_BL_Rhemi.tmoc.InterPhi24_Grandraw_BL_Rhemi;
            anovamat2(3,41:60)=Phi2ch.Grandraw_BL_Interhemi.tmoc.InterPhi24_Grandraw_BL_Interhemi;
            
            % Two-way ANOVA
            [p,tbl,stats] = anova2(anovamat2',20);
            
            INTERPhi2chStat3way.GrandPhi2ch_WInhemi_TwoWayANOVA.anovamat2 = anovamat2';
            INTERPhi2chStat3way.GrandPhi2ch_WInhemi_TwoWayANOVA.stattbl = tbl;
            
            %_________ posthoc
            c = multcompare(stats,'CType',posthocType,'Display','off');

            INTERPhi2chStat3way.GrandPhi2ch_WInhemi_TwoWayANOVA.posthoc = c;

            save ([output 'INTERPhi2chStat3way_24.mat'],'INTERPhi2chStat3way');

    %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Two-way ANOVA for the Grand "L-, R-, Inter- hemisphere" inter-brain Phi 2 channel: 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear;
    close all;

    input = ('/Users/mohammadshehata/Desktop/MoShehata/1- Projects/F3- Angus-Nao IIT project/Angus/');
    output=('/Users/mohammadshehata/Desktop/MoShehata/1- Projects/F3- Angus-Nao IIT project/Angus/STAT_Phi2ch/');

    cond={'tm' 'tmrv' 'tmoc'};
    posthocType=('bonferroni'); % set posthoc test type: 'hsd' (tukey), 'lsd','bonferroni'

    NewLabel ={'Gp1_L';'Gp2_L';'Gp3_L';'Gp4_L';'Gp5_L';'Gp6_L';'Gp7_L';'Gp1_R';'Gp2_R';'Gp3_R';'Gp4_R';'Gp5_R';'Gp6_R';'Gp7_R'};
        
    load([input 'Phi2ch.mat']);
    load([output 'INTERPhi2chStat3way_24.mat']);
    
    % Data arrangment for 2-way anova
            anovamat2(1,1:20)= Phi2ch.Grandraw_BL_Lhemi.tm.InterPhi24_Grandraw_BL_Lhemi;
            anovamat2(1,21:40)=Phi2ch.Grandraw_BL_Lhemi.tmrv.InterPhi24_Grandraw_BL_Lhemi;
            anovamat2(1,41:60)=Phi2ch.Grandraw_BL_Lhemi.tmoc.InterPhi24_Grandraw_BL_Lhemi;
            
            anovamat2(2,1:20)= Phi2ch.Grandraw_BL_Rhemi.tm.InterPhi24_Grandraw_BL_Rhemi;
            anovamat2(2,21:40)=Phi2ch.Grandraw_BL_Rhemi.tmrv.InterPhi24_Grandraw_BL_Rhemi;
            anovamat2(2,41:60)=Phi2ch.Grandraw_BL_Rhemi.tmoc.InterPhi24_Grandraw_BL_Rhemi;
            
            anovamat2(3,1:20)= Phi2ch.Grandraw_BL_Interhemi.tm.InterPhi24_Grandraw_BL_Interhemi;
            anovamat2(3,21:40)=Phi2ch.Grandraw_BL_Interhemi.tmrv.InterPhi24_Grandraw_BL_Interhemi;
            anovamat2(3,41:60)=Phi2ch.Grandraw_BL_Interhemi.tmoc.InterPhi24_Grandraw_BL_Interhemi;
            
            % Two-way ANOVA
            [p,tbl,stats] = anova2(anovamat2',20);
            
            INTERPhi2chStat3way.GrandPhi2ch_AXhemi_TwoWayANOVA.anovamat2 = anovamat2';
            INTERPhi2chStat3way.GrandPhi2ch_AXhemi_TwoWayANOVA.stattbl = tbl;
            
            %_________ posthoc
            c = multcompare(stats,'CType',posthocType,'Display','off');

            INTERPhi2chStat3way.GrandPhi2ch_AXhemi_TwoWayANOVA.posthoc = c;

            save ([output 'INTERPhi2chStat3way_24.mat'],'INTERPhi2chStat3way');

%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% N-way ANOVA for the Grand "L-, R-, Inter- hemisphere" inter-brain Phi 2 channel: 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    clear;
    close all;

    input = ('/Users/mohammadshehata/Desktop/MoShehata/1- Projects/F3- Angus-Nao IIT project/Angus/');
    output=('/Users/mohammadshehata/Desktop/MoShehata/1- Projects/F3- Angus-Nao IIT project/Angus/STAT_Phi2ch/');

    cond={'tm' 'tmrv' 'tmoc'};
    posthocType=('bonferroni'); % set posthoc test type: 'hsd' (tukey), 'lsd','bonferroni'

    NewLabel ={'Gp1_L';'Gp2_L';'Gp3_L';'Gp4_L';'Gp5_L';'Gp6_L';'Gp7_L';'Gp1_R';'Gp2_R';'Gp3_R';'Gp4_R';'Gp5_R';'Gp6_R';'Gp7_R'};
        
    load([input 'Phi2ch.mat']);
    load([output 'INTERPhi2chStat3way_24.mat']);    

            % Data arrangment for 3-way anova (participant x hemisphere x condition)
            anovamat(:,1,1)= Phi2ch.Grandraw_BL_Lhemi.tm.InterPhi24_Grandraw_BL_Lhemi;
            anovamat(:,1,2)=Phi2ch.Grandraw_BL_Lhemi.tmrv.InterPhi24_Grandraw_BL_Lhemi;
            anovamat(:,1,3)=Phi2ch.Grandraw_BL_Lhemi.tmoc.InterPhi24_Grandraw_BL_Lhemi;
            
            anovamat(:,2,1)= Phi2ch.Grandraw_BL_Rhemi.tm.InterPhi24_Grandraw_BL_Rhemi;
            anovamat(:,2,2)=Phi2ch.Grandraw_BL_Rhemi.tmrv.InterPhi24_Grandraw_BL_Rhemi;
            anovamat(:,2,3)=Phi2ch.Grandraw_BL_Rhemi.tmoc.InterPhi24_Grandraw_BL_Rhemi;
            
            anovamat(:,3,1)= Phi2ch.Grandraw_BL_Interhemi.tm.InterPhi24_Grandraw_BL_Interhemi;
            anovamat(:,3,2)=Phi2ch.Grandraw_BL_Interhemi.tmrv.InterPhi24_Grandraw_BL_Interhemi;
            anovamat(:,3,3)=Phi2ch.Grandraw_BL_Interhemi.tmoc.InterPhi24_Grandraw_BL_Interhemi;
            
    
            anovamat3.Y=reshape(anovamat,[numel(anovamat),1]);
            anovamat3.S=repmat([1:20]',9,1); %20subject X (3hemisphere x 3conditions)
            F1(1:60,1)={'tm'}; F1(61:120,1)={'tmrv'}; F1(121:180,1)={'tmoc'};
            anovamat3.G1=F1; %3conditions x (20subject x 3 hemispheres)
            F2(1:20,1)={'L-hemi'}; F2(21:40,1)={'R-hemi'}; F2(41:60,1)={'Inter-hemi'};
            anovamat3.G2=repmat(F2,3,1); %3hemispheres x (20subject x 3 conditions)
            
            anovamat3.Gnames={'subject','condition','hemisphere'};
            
            % Three-way ANOVA
            [p,tbl,stats] = anovan(anovamat3.Y,{anovamat3.S anovamat3.G1 anovamat3.G2},'model','interaction','varnames',anovamat3.Gnames);
            
            INTERPhi2chStat3way.GrandPhi2ch_hemi_NWayANOVA.anovamat = anovamat;
            INTERPhi2chStat3way.GrandPhi2ch_hemi_NWayANOVA.stattbl = tbl;
            
            %_________ posthoc
            [c,m,h,gnames]= multcompare(stats,'Dimension',[2 3],'CType',posthocType,'Display','off');
            INTERPhi2chStat3way.GrandPhi2ch_hemi_NWayANOVA.posthoc123 = c;
            INTERPhi2chStat3way.GrandPhi2ch_hemi_NWayANOVA.gnames = gnames;
                        
            save ([output 'INTERPhi2chStat3way_24.mat'],'INTERPhi2chStat3way');