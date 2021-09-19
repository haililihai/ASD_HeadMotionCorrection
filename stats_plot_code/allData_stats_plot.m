sitename={'caltech','ohsu','olin','pitt','sdsu','yale','301_01','301_02'};
% read xls
T=readtable(fullfile(pwd,'allData_flat.csv'));
newT=[];

for isite=1:length(sitename)
    
    wd=fullfile(pwd,sitename{isite},'plotdir2');

    % select
    rows=(strcmpi(T.Site,sitename{isite}));
    partT=T(rows,:);

    % load data
    load(fullfile(wd,'allData.mat'));


    for denoise=1:4
        iEdge=0;
        for i=2:17
            for j=1:i-1
                iEdge=iEdge+1;

                % add FC
                Name=['Denoise',num2str(denoise),'_FC_',num2str(i),'_',num2str(j)];
                FC=zeros(size(partT,1),1);
                for sub=1:size(partT,1) 
                    FC(sub,1)=allData(denoise).FC(i,j,sub);
                end         
                partT=addvars(partT,FC,'Before',['Denoise',num2str(denoise),'_QCFC_PropSig_corr'],'NewVariableNames',Name);
            end      
        end

        % add GCOR
        Name=['Denoise',num2str(denoise),'_GCOR'];
        partT=addvars(partT,allData(denoise).GCOR,'Before',['Denoise',num2str(denoise),'_QCFC_PropSig_corr'],'NewVariableNames',Name);

        % add QCFC_Vec
        iEdge=0;
        for i=2:17
            for j=1:i-1
                iEdge=iEdge+1;

                % add QCFC for each group
                for group=1:2
                   Name=['Denoise',num2str(denoise),'_QCFC_',num2str(i),'_',num2str(j),'_Group',num2str(group)];
                   QCFC=zeros(size(partT,1),1);
                   for sub=1:size(partT,1) 
                        QCFC(sub,1)=allData(denoise).QCFCVec{1,group}(iEdge,1);
                   end
                   partT=addvars(partT,QCFC,'Before',['Denoise',num2str(denoise),'_QCFC_PropSig_corr'],'NewVariableNames',Name);
                end
            end      
        end

        % split vars
        vars={'QCFC_PropSig_corr','QCFC_PropSig_unc','QCFC_AbsMed','QCFC_DistDep','QCFC_DistDep_Pval','tDOF_mean'};

        for k=1:length(vars)
            orig_name=['Denoise',num2str(denoise),'_',vars{k}];
            partT.(orig_name)=cell2table(cellfun(@(x) str2num(x),partT.(orig_name),'UniformOutput',false),'VariableNames',{orig_name});
            partT.(orig_name)=splitvars(partT.(orig_name),orig_name,'NewVariableNames',{[orig_name,'_Group1'],[orig_name,'_Group2']});
            partT=splitvars(partT,orig_name,'NewVariableNames',{[orig_name,'_Group1'],[orig_name,'_Group2']});

        end

    end
    
    newT=[newT;partT];
end

% QC-FC histogram plot
T=readtable('All_Data_Flat_20210312.csv');
QCFC=cell(4,2);

for i=1:4
    for j=1:2
        for ii=2:17
            for jj=1:ii-1
                QCFC{i,j}=[QCFC{i,j};reshape(unique(T.(['Denoise',num2str(i),'_QCFC_',num2str(ii),'_',num2str(jj),'_Group',num2str(j)])),8,1)];
            end
        end
    end
end

for num=1:4
    subplot(2,2,num);
    histogram(QCFC{num,1},'BinWidth',0.05);
    hold on;
    histogram(QCFC{num,2},'BinWidth',0.05);
    xlim([-1,1]);
    ylim([0,100]);
    legend({'TD','ASD'});
    title(T.(['Denoise',num2str(num),'_Name']){1});
    xlabel('Correlation with head motion');
    ylabel('Count of connections');
end
export_fig(fullfile(pwd,'QCFC_Histogram_highres.tif'),'-r600','-painters','-nocrop');

% Yeo17_3_all correlation (Somatomotor correlation)
% Yeo17_4_all correlation (Somatomotor correlation)
T=readtable('All_Data_Flat_20210312.csv');

Group={'TD','ASD'};
strategy={'6HMP+2Phys','6HMP+2Phys+GSR','ICA-AROMA+2Phys','ICA-AROMA+2Phys+GSR'};
Labels=cell(1,17);
for i=1:length(Labels)
    Labels{i}=num2str(i);
end

% Yeo17_3
roi=16;
mFC=zeros(2,17,4);
for group=1:2
    for i=1:4
        for j=1:roi-1
            mFC(group,j,i)=mean(T.(['Denoise',num2str(i),'_FC_',num2str(roi),'_',num2str(j)])(string(T.Group)==Group{group}));
        end
        %mQCFC(group,roi,i)=1;
        for j=roi+1:17
            mFC(group,j,i)=mean(T.(['Denoise',num2str(i),'_FC_',num2str(j),'_',num2str(roi)])(string(T.Group)==Group{group}));
        end
    end
end

% FC group diff box plot in matrix, two-tailed
T=readtable('All_Data_Flat_20210312.csv');
strategy={'6HMP+2Phys','6HMP+2Phys+GSR','ICA-AROMA+2Phys','ICA-AROMA+2Phys+GSR'};
Group={'TD','ASD'};
addpath(genpath(fullfile(pwd,'export_fig')));
numNetwork=17;
denoise=1;
for denoise=1:3
    p_matrix=readmatrix(fullfile(pwd,['Denoise',num2str(denoise),'_p.txt']));
    tp=tiledlayout(numNetwork,numNetwork,'TileSpacing','none','Position',[10,10,300*numNetwork,300*numNetwork]);
    %fid=fopen(['FC_denoise',num2str(denoise),'_group_diff_p.txt'],'w');
    for m=1:numNetwork
        for n=1:numNetwork
            if m < n
                mm=n;
                nn=m;
            end
            if m > n
                mm=m;
                nn=n;
            end
            if m ~= n
                bb=nexttile((m-1)*numNetwork+n);
                in=T.(['Denoise',num2str(denoise),'_FC_',num2str(mm),'_',num2str(nn)]);
                grp=T.Group;
                in1=in(string(T.Group)==Group{1});
                in2=in(string(T.Group)==Group{2});
                boxplot(in,grp,'GroupOrder',{'TD','ASD'},'Whisker',2);
                set(bb,'XTickLabel',{''},'YTickLabel',{''});
                if n==1
                    ylabel(bb,num2str(m),'FontWeight','bold','Rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right');
                end
                if m==numNetwork
                    xlabel(bb,num2str(n),'FontWeight','bold');
                end
                % sig
                p=p_matrix(m,n);
                if p<=0.05
                   sigstar({[1,2]},[p]);
                   set(bb,'Color','yellow');
                end 
                % left
                %[h,p]=ttest2(in1,in2,'Tail','left');
                %p=p_matrix(m,n);
                %fprintf(fid,'FC_%d_%d denoise %d group diff left p: %f\n',m,n,denoise,p);
                %if h==1
                %   sigstar({[1,2]},[p]);
                %   set(bb,'Color','magenta');
                %end 
                % right
                %[h,p]=ttest2(in1,in2,'Tail','right');
                %p_matrix_right(m,n)=p;
                %fprintf(fid,'FC_%d_%d denoise %d group diff right p: %f\n',m,n,denoise,p);
                %if h==1
                %   sigstar({[1,2]},[p]);
                %   set(bb,'Color','cyan');
                %end
            end
        end
    end
    %fclose(fid);
    tp.XLabel.String=strategy{denoise};
    %save(fullfile(pwd,['Denoise',num2str(denoise),'_p_left.txt']),'p_matrix_left','-ascii');
    %save(fullfile(pwd,['Denoise',num2str(denoise),'_p_right.txt']),'p_matrix_right','-ascii');
    export_fig(fullfile(pwd,['FC_denoise',num2str(denoise),'_group_diff_merged_twotailed_highres.tif']),'-r600','-painters','-nocrop');
end

% FC group diff box plot in matrix
T=readtable('All_Data_Flat_20210312.csv');
strategy={'6HMP+2Phys','6HMP+2Phys+GSR','ICA-AROMA+2Phys','ICA-AROMA+2Phys+GSR'};
Group={'TD','ASD'};
addpath(genpath(fullfile(pwd,'export_fig')));
numNetwork=17;
denoise=1;
for denoise=2:4
    p_matrix_right=zeros(numNetwork,numNetwork);
    p_matrix_left=zeros(numNetwork,numNetwork);
    tp=tiledlayout(numNetwork,numNetwork,'TileSpacing','none','Position',[10,10,300*numNetwork,300*numNetwork]);
    fid=fopen(['FC_denoise',num2str(denoise),'_group_diff_p.txt'],'w');
    for m=1:numNetwork
        for n=1:numNetwork
            if m < n
                mm=n;
                nn=m;
            end
            if m > n
                mm=m;
                nn=n;
            end
            if m ~= n
                bb=nexttile((m-1)*numNetwork+n);
                in=T.(['Denoise',num2str(denoise),'_FC_',num2str(mm),'_',num2str(nn)]);
                grp=T.Group;
                in1=in(string(T.Group)==Group{1});
                in2=in(string(T.Group)==Group{2});
                boxplot(in,grp,'GroupOrder',{'TD','ASD'},'Whisker',2);
                set(bb,'XTickLabel',{''},'YTickLabel',{''});
                if n==1
                    ylabel(bb,num2str(m),'FontWeight','bold','Rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right');
                end
                if m==numNetwork
                    xlabel(bb,num2str(n),'FontWeight','bold');
                end
                % left
                [h,p]=ttest2(in1,in2,'Tail','left');
                p_matrix_left(m,n)=p;
                fprintf(fid,'FC_%d_%d denoise %d group diff left p: %f\n',m,n,denoise,p);
                if h==1
                   sigstar({[1,2]},[p]);
                   set(bb,'Color','magenta');
                end 
                % right
                [h,p]=ttest2(in1,in2,'Tail','right');
                p_matrix_right(m,n)=p;
                fprintf(fid,'FC_%d_%d denoise %d group diff right p: %f\n',m,n,denoise,p);
                if h==1
                   sigstar({[1,2]},[p]);
                   set(bb,'Color','cyan');
                end
            end
        end
    end
    fclose(fid);
    tp.XLabel.String=strategy{denoise};
    save(fullfile(pwd,['Denoise',num2str(denoise),'_p_left.txt']),'p_matrix_left','-ascii');
    save(fullfile(pwd,['Denoise',num2str(denoise),'_p_right.txt']),'p_matrix_right','-ascii');
    export_fig(fullfile(pwd,['FC_denoise',num2str(denoise),'_group_diff_merged_new.png']),'-r300','-painters','-nocrop');
end
% saveas(gca,fullfile(pwd,'FC_Group_diff_merged.png'));


% FDR correction two-tailed
for denoise=1:4
   p_matrix=fullfile(['Denoise',num2str(denoise),'_p.txt']);
   p=readmatrix(p_matrix);
   p_tril=tril(p,-1);
   p_values=reshape(p_tril(p_tril>0),[136 1]);
   q=0.05;
   p_corrected=mafdr(p_values,'BHFDR','true');
   p_fdr=tril(ones(17),-1);
   p_fdr(p_fdr>0)=p_corrected;
   p_fdr_matrix=(p_fdr+p_fdr');
   save(fullfile(pwd,['Denoise',num2str(denoise),'_p_fdr.txt']),'p_fdr_matrix','-ascii');
end

T=readtable('All_Data_Flat_20210312.csv');
strategy={'6HMP+2Phys','6HMP+2Phys+GSR','ICA-AROMA+2Phys','ICA-AROMA+2Phys+GSR'};
Group={'TD','ASD'};
addpath(genpath(fullfile(pwd,'export_fig')));
numNetwork=17;
denoise=1;
for denoise=1:4
    p_matrix=readmatrix(fullfile(pwd,['Denoise',num2str(denoise),'_p_fdr.txt']));
    tp=tiledlayout(numNetwork,numNetwork,'TileSpacing','none','Position',[10,10,300*numNetwork,300*numNetwork]);
    fid=fopen(['FC_denoise',num2str(denoise),'_group_diff_p_fdr.txt'],'w');
    for m=1:numNetwork
        for n=1:numNetwork
            if m < n
                mm=n;
                nn=m;
            end
            if m > n
                mm=m;
                nn=n;
            end
            if m ~= n
                bb=nexttile((m-1)*numNetwork+n);
                in=T.(['Denoise',num2str(denoise),'_FC_',num2str(mm),'_',num2str(nn)]);
                grp=T.Group;
                in1=in(string(T.Group)==Group{1});
                in2=in(string(T.Group)==Group{2});
                boxplot(in,grp,'GroupOrder',{'TD','ASD'},'Whisker',2);
                set(bb,'XTickLabel',{''},'YTickLabel',{''});
                if n==1
                    ylabel(bb,num2str(m),'FontWeight','bold','Rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right');
                end
                if m==numNetwork
                    xlabel(bb,num2str(n),'FontWeight','bold');
                end
                % left
                %[h,p]=ttest2(in1,in2);
                p=p_matrix(m,n);
                fprintf(fid,'FC_%d_%d denoise %d group diff fdr p: %f\n',m,n,denoise,p);
                if p<=0.05
                   sigstar({[1,2]},[p]);
                   set(bb,'Color','yellow');
                end 
              
            end
        end
    end
    fclose(fid);
    tp.XLabel.String=strategy{denoise};
%     save(fullfile(pwd,['Denoise',num2str(denoise),'_p_left.txt']),'p_matrix_left','-ascii');
%     save(fullfile(pwd,['Denoise',num2str(denoise),'_p_right.txt']),'p_matrix_right','-ascii');
    export_fig(fullfile(pwd,['FC_denoise',num2str(denoise),'_group_diff_merged_fdr_twotailed_highres.tif']),'-r600','-painters','-nocrop');
end

% FDR correction one-tailed
for denoise=1:4
   % left
   p_matrix=fullfile(['Denoise',num2str(denoise),'_p_left.txt']);
   p=readmatrix(p_matrix);
   p_tril=tril(p,-1);
   p_values=reshape(p_tril(p_tril>0),[136 1]);
   p_corrected=mafdr(p_values,'BHFDR','true');
   p_fdr=tril(ones(17),-1);
   p_fdr(p_fdr>0)=p_corrected;
   p_fdr_matrix=(p_fdr+p_fdr');
   save(fullfile(pwd,['Denoise',num2str(denoise),'_p_left_fdr.txt']),'p_fdr_matrix','-ascii');
   
   % right
   p_matrix=fullfile(['Denoise',num2str(denoise),'_p_right.txt']);
   p=readmatrix(p_matrix);
   p_tril=tril(p,-1);
   p_values=reshape(p_tril(p_tril>0),[136 1]);
   p_corrected=mafdr(p_values,'BHFDR','true');
   p_fdr=tril(ones(17),-1);
   p_fdr(p_fdr>0)=p_corrected;
   p_fdr_matrix=(p_fdr+p_fdr');
   save(fullfile(pwd,['Denoise',num2str(denoise),'_p_right_fdr.txt']),'p_fdr_matrix','-ascii');
end

denoise=1;
for denoise=1:4
    p_matrix_left=readmatrix(fullfile(pwd,['Denoise',num2str(denoise),'_p_left_fdr.txt']));
    p_matrix_right=readmatrix(fullfile(pwd,['Denoise',num2str(denoise),'_p_right_fdr.txt']));
    tp=tiledlayout(numNetwork,numNetwork,'TileSpacing','none','Position',[10,10,300*numNetwork,300*numNetwork]);
    fid=fopen(['FC_denoise',num2str(denoise),'_group_diff_p_fdr_onetailed.txt'],'w');
    for m=1:numNetwork
        for n=1:numNetwork
            if m < n
                mm=n;
                nn=m;
            end
            if m > n
                mm=m;
                nn=n;
            end
            if m ~= n
                bb=nexttile((m-1)*numNetwork+n);
                in=T.(['Denoise',num2str(denoise),'_FC_',num2str(mm),'_',num2str(nn)]);
                grp=T.Group;
                in1=in(string(T.Group)==Group{1});
                in2=in(string(T.Group)==Group{2});
                boxplot(in,grp,'GroupOrder',{'TD','ASD'},'Whisker',2);
                set(bb,'XTickLabel',{''},'YTickLabel',{''});
                if n==1
                    ylabel(bb,num2str(m),'FontWeight','bold','Rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right');
                end
                if m==numNetwork
                    xlabel(bb,num2str(n),'FontWeight','bold');
                end
                % left
                %[h,p]=ttest2(in1,in2,'Tail','left');
                p=p_matrix_left(m,n);
                fprintf(fid,'FC_%d_%d denoise %d group diff left fdr p: %f\n',m,n,denoise,p);
                if p<=0.05
                   sigstar({[1,2]},[p]);
                   set(bb,'Color','magenta');
                end 
                % right
                %[h,p]=ttest2(in1,in2,'Tail','right');
                p=p_matrix_right(m,n);
                fprintf(fid,'FC_%d_%d denoise %d group diff right fdr p: %f\n',m,n,denoise,p);
                if p<=0.05
                   sigstar({[1,2]},[p]);
                   set(bb,'Color','cyan');
                end
            end
        end
    end
    fclose(fid);
    tp.XLabel.String=strategy{denoise};
    %save(fullfile(pwd,['Denoise',num2str(denoise),'_p_left_fdr.txt']),'p_matrix_left','-ascii');
    %save(fullfile(pwd,['Denoise',num2str(denoise),'_p_right.txt']),'p_matrix_right','-ascii');
    export_fig(fullfile(pwd,['FC_denoise',num2str(denoise),'_group_diff_merged_fdr_onetailed_highres.tif']),'-r600','-painters','-nocrop');
end

% network name, Baker 2014
p_label{1}='1 Visual peripheral';
p_label{2}='2 Visual central';
p_label{3}='3 Somatomotor A';
p_label{4}='4 Somatomotor B';
p_label{5}='5 Dorsal attention A';
p_label{6}='6 Dorsal attention B';
p_label{7}='7 Ventral attention';
p_label{8}='8 Salience';
p_label{9}='9 Limbic A';
p_label{10}='10 Limbic B';
p_label{11}='11 Control C';
p_label{12}='12 Control A';
p_label{13}='13 Control B';
p_label{14}='14 Default D';
p_label{15}='15 Default C';
p_label{16}='16 Default A';
p_label{17}='17 Default B';

% bin p matrix
for denoise=1:4
   % left
   input=fullfile(pwd,['Denoise',num2str(denoise),'_p_left.txt']);
   p=readmatrix(input);
   p(p==0)=10;
   p(p>0.05)=10;
   p(p<=0.001)=3;
   p(p<=0.01)=2;
   p(p<=0.05)=1;
   p(p==10)=0;
   writematrix(p,fullfile(pwd,['Denoise',num2str(denoise),'_p_left_bin.txt']));
   cla
   circularGraph(p,'Colormap',lines(17),'Label',p_label);
   %saveas(gca,fullfile(pwd,['Denoise',num2str(denoise),'_p_left_bin_name.png']));
   export_fig(fullfile(pwd,['Denoise',num2str(denoise),'_p_left_bin_name_highres.tif']),'-r600','-painters','-nocrop');
   
   % right
   input=fullfile(pwd,['Denoise',num2str(denoise),'_p_right.txt']);
   p=readmatrix(input);
   p(p==0)=10;
   p(p>0.05)=10;
   p(p<=0.001)=3;
   p(p<=0.01)=2;
   p(p<=0.05)=1;
   p(p==10)=0;
   writematrix(p,fullfile(pwd,['Denoise',num2str(denoise),'_p_right_bin.txt']));
   cla
   circularGraph(p,'Colormap',lines(17),'Label',p_label);
   %saveas(gca,fullfile(pwd,['Denoise',num2str(denoise),'_p_right_bin_name.png']));
   export_fig(fullfile(pwd,['Denoise',num2str(denoise),'_p_right_bin_name_highres.tif']),'-r600','-painters','-nocrop');

end

% FC denoise diff box plot in matrix
Group={'TD','ASD'};
numNetwork=17;
group=1;
for group=1:2
    tp=tiledlayout(numNetwork,numNetwork,'TileSpacing','none','Position',[10,10,300*numNetwork,300*numNetwork]);
    fid=fopen(['FC_group',num2str(group),'_denoise_diff_p.txt'],'w');
    for m=1:numNetwork
        for n=1:numNetwork
            if m < n
                mm=n;
                nn=m;
            end
            if m > n
                mm=m;
                nn=n;
            end
            if m ~= n
                bb=nexttile((m-1)*numNetwork+n);
                in1=T.(['Denoise1_FC_',num2str(mm),'_',num2str(nn)])(string(T.Group)==Group{group});
                in2=T.(['Denoise2_FC_',num2str(mm),'_',num2str(nn)])(string(T.Group)==Group{group});
                in3=T.(['Denoise3_FC_',num2str(mm),'_',num2str(nn)])(string(T.Group)==Group{group});
                in4=T.(['Denoise4_FC_',num2str(mm),'_',num2str(nn)])(string(T.Group)==Group{group});
                boxplot([in1,in2,in3,in4],'Whisker',2);
                set(bb,'XTickLabel',{''},'YTickLabel',{''});
                if n==1
                    ylabel(bb,num2str(m),'FontWeight','bold','Rotation',0,'VerticalAlignment','middle', 'HorizontalAlignment','right');
                end
                if m==numNetwork
                    xlabel(bb,num2str(n),'FontWeight','bold');
                end
                for i=1:3
                    for j=i+1:4
                        [h,p]=ttest2(T.(['Denoise',num2str(i),'_FC_',num2str(mm),'_',num2str(nn)])(string(T.Group)==Group{group}),...
                            T.(['Denoise',num2str(j),'_FC_',num2str(mm),'_',num2str(nn)])(string(T.Group)==Group{group}));
                        fprintf(fid,'FC_%d_%d group %d denoise %d_%d diff p: %f\n',m,n,group,i,j,p);
                        if h==1
                            sigstar({[i,j]},[p]);
                            %set(bb,'Color','yellow');
                        end
                    end
                end 
            end
        end
    end
    fclose(fid);
    tp.XLabel.String=Group{group};
    export_fig(fullfile(pwd,['FC_group',num2str(group),'_denoise_diff_merged.png']),'-r600','-painters','-nocrop');
end


% radar plot
for i=1:4
    subplot(2,2,i);
    spider_plot_R2019b(mFC(:,:,i),...
                       'AxesLabels',Labels,...
                       'AxesLimits',[-0.55*ones(1,17);1.1*ones(1,17)],...
                       'AxesFontSize',6,...
                       'LabelFontSize',8,...
                       'AxesLabelsEdge','none',...
                       'FillOption','on');
    legend('TD','ASD','Location','southoutside');
    title(strategy{i});
end

% freesurfer plot
hemi={'lh','rh'};
for ihemi=1:length(hemi)
    % annot to mgh command
    mgh_file=fullfile('.',['Yeo17_',hemi{ihemi},'_seg.mgh']);
    command=['mri_annotation2label --subject fsaverage --hemi ',hemi{ihemi},' --annotation Yeo2011_17Networks_N1000 --seg ',mgh_file];
    system(command);
    
    % load mgh
    [mgh_vol,mgh_m]=load_mgh(mgh_file);
    
    % set value
    for group=1:2
        for denoise=1:4
            vol=mgh_vol;
            for i=1:17
                vol(vol==i)=mFC(group,i,denoise);
            end

            % save mgh
            output_mgh_file=fullfile('.',['mFC_network',num2str(roi),'_group',num2str(group),'_denoise',num2str(denoise),'_Yeo17_',hemi{ihemi},'.mgh']);
            save_mgh(vol,output_mgh_file,mgh_m);
            system(['cp -rv ',output_mgh_file,' ~/Documents/test/test_freesurfer/']);
        end
    end   
end

% QCFC_DistDep_Big

% ROI Distance
% REQUIRED
ROIDir = [pwd,'/ROIs/'];
WhichParc='Yeo';
switch WhichParc
    case 'Gordon'
        Parc = 1;
        ROI_Coords = dlmread([ROIDir,'Gordon/Gordon_Centroids.txt']);
        fileName = [ROIDir,'Gordon/CommunityModified.txt'];
    case 'Yeo'
        Parc = 2;
        ROI_Coords = dlmread([ROIDir,'Yeo/Yeo_coordinates.txt']);
        fileName = [ROIDir,'Yeo/Yeo_community.txt'];
end

fileID = fopen(fileName);
ROIStruct = textscan(fileID,'%s'); ROIStruct = ROIStruct{1};
% rearrange by community
[ROIStruct_com,ROI_idx] = sort(ROIStruct);

% Convert text labels to unique integer values
ROILabels = unique(ROIStruct);
ROIStructID = zeros(size(ROIStruct));
numROIComms = length(ROILabels);
for i = 1:numROIComms
    x = find(strcmp(ROIStruct, ROILabels{i}));
    ROIStructID(x) = i;
end

% ------------------------------------------------------------------------------
% Load ROI coordinates
% ------------------------------------------------------------------------------
% Calculate pairwise euclidean distance
ROIDist = pdist2(ROI_Coords,ROI_Coords,'euclidean');

% Flatten distance matrix
ROIDistVec = LP_FlatMat(ROIDist);

% Calculate number of ROIs
numROIs = size(ROIDist,1);

% Calculate number of edges
numConnections = numROIs * (numROIs - 1) / 2;

% Merge

sitename={'caltech','ohsu','olin','pitt','sdsu','yale','301_01','301_02'};
for denoise=1:4
    allData_merged(denoise).metadata=[];
    allData_merged(denoise).FC=[];
    for isite=1:length(sitename)
        wd=fullfile(pwd,sitename{isite},'plotdir2');
        load(fullfile(wd,'allData.mat'));

        allData_merged(denoise).metadata=vertcat(allData_merged(denoise).metadata,allData(denoise).metadata);
        allData_merged(denoise).FC=cat(3,allData_merged(denoise).FC,allData(denoise).FC);
    end
    for group=1:2
       [allData_merged(denoise).QCFCVec{group},allData_merged(denoise).NaNFilter{group},allData_merged(denoise).QCFC_PropSig_corr(group),allData_merged(denoise).QCFC_PropSig_unc(group),allData_merged(denoise).QCFC_AbsMed(group),allData_merged(denoise).QCFC_DistDep(group),allData_merged(denoise).QCFC_DistDep_Pval(group)] = RunQCFC(allData_merged(denoise).metadata.fdJenk_m(allData_merged(denoise).metadata.Diagnosis==group),allData_merged(denoise).FC(:,:,allData_merged(denoise).metadata.Diagnosis==group),ROIDistVec); 
    end
    allData_merged(denoise).noiseOptions=allData(denoise).noiseOptions;
    allData_merged(denoise).noiseOptionsNames=allData(denoise).noiseOptionsNames;
end

% save table as csv
T=table;
T.ROIDistVec=ROIDistVec;
for denoise=1:4
    for group=1:2
        T.(['Denoise',num2str(denoise),'_Group',num2str(group)])=allData_merged(denoise).QCFCVec{group};
    end
end
writetable(T,fullfile(pwd,'allData_QCFCVec_merged.csv'));


% color
tempColors = num2cell([255,105,97;97,168,255;178,223,138;117,112,179;255,179,71]./255,2); 
theColors{1}=tempColors{1};
theColors{2}=tempColors{1};
theColors{3}=tempColors{4};
theColors{4}=tempColors{4};


% ------------------------------------------------------------------------------
    % Figures: QCFC_DistDep
    % ------------------------------------------------------------------------------
    % Initialise figures
    Fig_QCFC_DistDepBig = figure('color','w', 'units', 'centimeters', 'pos', [0 0 25 27], 'name',['Fig_QCFC_DistDepBig']); box('on'); movegui(Fig_QCFC_DistDepBig,'center');

    pipelines2Retain = logical([1 1 1 1]);
    noiseOptions_temp = {allData_merged(pipelines2Retain).noiseOptions};
    noiseOptionsNames_temp = {allData_merged(pipelines2Retain).noiseOptionsNames};
    theColors_temp = theColors(pipelines2Retain);

    numPrePro_temp = length(noiseOptions_temp);
    
    FSize=10;

    for i = 1:numPrePro_temp

        WhichNoise = noiseOptions_temp{i};
        WhichNoiseName = noiseOptionsNames_temp{i};
        idx = strmatch(WhichNoise,{allData_merged(:).noiseOptions},'exact');

        % ------------------------------------------------------------------------------
        % Plot: distance dependence
        % ------------------------------------------------------------------------------
        figure(Fig_QCFC_DistDepBig)
        subplot(4,ceil(numPrePro_temp/4),i)
        set(gca,'FontSize',FSize)

        % Bin QCFC data by distance and generate means and stds for each
        numThresholds = numConnections; % bins
        % BF_PlotQuantiles(ROIDistVec(allData(idx).NaNFilter{k}),allData(idx).QCFCVec{k},numThresholds,0,0,theColors_temp{i})
        BF_PlotQuantiles(ROIDistVec(allData_merged(idx).NaNFilter{1}),allData_merged(idx).QCFCVec{1},numThresholds,0,0,theColors_temp{1},'^')
        BF_PlotQuantiles(ROIDistVec(allData_merged(idx).NaNFilter{2}),allData_merged(idx).QCFCVec{2},numThresholds,0,0,theColors_temp{4},'v')
        hold on
        % plot([0:200],zeros(1,201),'--','Color','k')
        plot([0:160],zeros(1,161),'--','Color','k')
        % text(140,1,['\delta','-TD'],'Color',theColors_temp{1})
        % text(140,0.8,['\nabla','-ASD'],'Color',theColors_temp{4})
        
        xlabel('Distance (mm)')
        ylabel('QC-FC')

        title(WhichNoiseName,'Interpreter', 'none','FontSize',FSize,'FontWeight','normal')
    end

    saveas(Fig_QCFC_DistDepBig,fullfile(pwd,'QCFC_DistDepBig_merged.bmp'));
    export_fig(fullfile(pwd,'QCFC_DistDepBig_merged_highres.tif'),'-r600','-painters','-nocrop');

    % ------------------------------------------------------------------------------
    % Figures: QCFC bins
    % ------------------------------------------------------------------------------
    % Initialise figures
    Fig_QCFC_DistDepBig_bins = figure('color','w', 'units', 'centimeters', 'pos', [0 0 25 27], 'name',['Fig_QCFC_DistDepBig_bins']); box('on'); movegui(Fig_QCFC_DistDepBig_bins,'center');

    pipelines2Retain = logical([1 1 1 1]);
    noiseOptions_temp = {allData(pipelines2Retain).noiseOptions};
    noiseOptionsNames_temp = {allData(pipelines2Retain).noiseOptionsNames};
    theColors_temp = theColors(pipelines2Retain);

    numPrePro_temp = length(noiseOptions_temp);

    for i = 1:numPrePro_temp

        WhichNoise = noiseOptions_temp{i};
        WhichNoiseName = noiseOptionsNames_temp{i};
        idx = strmatch(WhichNoise,{allData(:).noiseOptions},'exact');

        % ------------------------------------------------------------------------------
        % Plot: distance dependence
        % ------------------------------------------------------------------------------
        figure(Fig_QCFC_DistDepBig_bins)
        subplot(4,ceil(numPrePro_temp/4),i)
        set(gca,'FontSize',FSize)

        % Bin QCFC data by distance and generate means and stds for each
        numThresholds = 11; % bins
        % BF_PlotQuantiles(ROIDistVec(allData(idx).NaNFilter{k}),allData(idx).QCFCVec{k},numThresholds,0,0,theColors_temp{i})
        BF_PlotQuantiles(ROIDistVec(allData(idx).NaNFilter{1}),allData(idx).QCFCVec{1},numThresholds,0,0,theColors_temp{1},'^')
        BF_PlotQuantiles(ROIDistVec(allData(idx).NaNFilter{2}),allData(idx).QCFCVec{2},numThresholds,0,0,theColors_temp{4},'v')
        hold on
        % plot([0:200],zeros(1,201),'--','Color','k')
        plot([0:160],zeros(1,161),'--','Color','k')
        
        xlabel('Distance (mm)')
        ylabel('QC-FC')

        title(WhichNoiseName,'Interpreter', 'none','FontSize',FSize,'FontWeight','normal')
    end

    saveas(Fig_QCFC_DistDepBig_bins,fullfile(pwd,'QCFC_DistDepBig_bins_merged.bmp'));
    export_fig(fullfile(pwd,'QCFC_DistDepBig_bins_merged_highres.tif'),'-r600','-painters','-nocrop');
 
 % Group diff stats
    
 T=readtable('allData_QCFCVec_merged.csv');
 for denoise=1:4
     subplot(2,2,denoise)
     boxplot([T.(['Denoise',num2str(denoise),'_Group1']),T.(['Denoise',num2str(denoise),'_Group2'])],'Labels',{'TD','ASD'},'Whisker',2);
     [h,p]=ttest2(T.(['Denoise',num2str(denoise),'_Group1']),T.(['Denoise',num2str(denoise),'_Group2']));
     fprintf('Denoise %d Group diff p: %f\n',denoise,p)
     if h==1
         sigstar({[1,2]},[p]);
     end
     
     xlabel(allData_merged(denoise).noiseOptionsNames);
 end
 saveas(gca,fullfile(pwd,'QCFC_Group_diff.bmp'));

 % Denoise diff stats
 Group={'TD','ASD'};
for group=1:2
    subplot(1,2,group);
    boxplot([T.(['Denoise1_Group',num2str(group)]),T.(['Denoise2_Group',num2str(group)]),T.(['Denoise3_Group',num2str(group)]),T.(['Denoise4_Group',num2str(group)])],'Labels',{'1','2','3','4'});
    for i=1:3
        for j=i+1:4
            [h,p]=ttest2(T.(['Denoise',num2str(i),'_Group',num2str(group)]),T.(['Denoise',num2str(j),'_Group',num2str(group)]));
            fprintf('Group %d denoise %d_%d p: %f\n',group,i,j,p)
            if h==1
                sigstar({[i,j]},[p]);
            end
        end
    end
    xlabel(Group{group});
end
 saveas(gca,fullfile(pwd,'QCFC_Denoise_diff.bmp'));
 
% ------------------------------------------------------------------------------
    % Figures: QCFC PropSig
    x = {allData_merged(:).noiseOptionsNames};
    % xy = cell(x);
    x_new = repelem(x,1,2);
    x_new2(1:2:7) =  cellfun(@(c) [c,' (TD)'],x_new(1:2:7),'uni',false);
    x_new2(2:2:8) =  cellfun(@(c) [c,' (ASD)'],x_new(2:2:8),'uni',false);
    xy=cell(x_new2);
    for i = 1:length(x)
        xy{i} = x_new2{i};
    end
    
    for denoise=1:4
        strs = strsplit(allData_merged(denoise).noiseOptions,'+');
        if any(strmatch('GSR',strs,'exact')) == 1 | any(strmatch('4GSR',strs,'exact')) == 1 | any(strmatch('2GSR',strs,'exact')) == 1
            theLines{denoise} = ':';
            % theLines{i} = '--';
        else
            theLines{denoise} = '-';
        end
    end
    
    
    Fig_QCFC_Dist = figure('color','w', 'units', 'centimeters', 'pos', [0 0 16 9], 'name',['Fig_QCFC_Dist']); box('on'); movegui(Fig_QCFC_Dist,'center');
    sp = subplot(1,3,1);
    pos = get(sp,'Position');
    % set(gca,'Position',[pos(1)*2.75, pos(2)*1.2, pos(3)*1.4, pos(4)*1]); % [left bottom width height]
    set(gca,'Position',[pos(1)*3.25, pos(2)*1.2, pos(3)*1.2, pos(4)*1]); % [left bottom width height]

    % Create data
    % data = {[allData_merged(:).QCFC_PropSig_corr]'};
    numGroups=2;
    numPrePro=4;
    if numGroups > 1
        QCFC_PropSig_unc=[allData_merged(:).QCFC_PropSig_unc];
        data = {[[QCFC_PropSig_unc(1:4),QCFC_PropSig_unc(5:8)]]'};
        data_std = cell(1,length(data)); [data_std{:}] = deal(zeros(size(data{1})));
        % data{1}=QCFC_PropSig_unc(1:2:7);
        % data{2}=QCFC_PropSig_unc(2:2:8);
        % data_std = cell(1,length(data{1})); [data_std{:}] = deal(zeros(size(data{1})));
        % data = {[QCFC_PropSig_unc(1:2:7)';QCFC_PropSig_unc(2:2:8)']'};
        % data_std=cell(1,length(data)); [data_std{:}] = deal(zeros(size(data{1})));
        % data1 = {QCFC_PropSig_unc(1:2:7)};
        % data2 = {QCFC_PropSig_unc(2:2:8)};
        % data_std1 = cell(1,length(data1)); [data_std1{:}] = deal(zeros(size(data1{1})));
        % data_std2 = cell(1,length(data2)); [data_std2{:}] = deal(zeros(size(data2{1})));
    else
        data = {[allData_merged(:).QCFC_PropSig_unc]'};
        data_std = cell(1,length(data)); [data_std{:}] = deal(zeros(size(data{1})));
    end
    

    % Create table
    % T = table(data{1},'RowNames',{allData_merged(:).noiseOptionsNames}','VariableNames',{'QCFC_PropSig'})

    % Create bar chart
    clear extraParams
    extraParams.xTickLabels = xy;
    % extraParams.xTickLabels = '';
    extraParams.xLabel = ''; % 'Pipeline'
    % extraParams.yLabel = 'QC-FC (%)';
    extraParams.yLabel = 'QC-FC uncorrected (%)';
    extraParams.theColors = theColors;
    extraParams.theLines = theLines;
    extraParams.yLimits = [0 110];
    % extraParams.xLimits = [0 2*size(data{1},1)+1];
    % extraParams.yLimits = [-110 110];
    extraParams.makeABS = true;

    TheBarChart(data,data_std,false,extraParams)
    % % TheBarChart(data,data_std,true,extraParams)
    % TheBarChart(data1,data_std1,false,extraParams)
    % TheBarChart(data2,data_std2,false,extraParams)

    % ------------------------------------------------------------------------------
    % QCFC distributions
    % ------------------------------------------------------------------------------
    sp = subplot(1,3,3);
    pos = get(sp,'Position');
    set(gca,'Position',[pos(1)*1, pos(2)*1.2, pos(3)*1.4, pos(4)*1]); % [left bottom width height]
    if numGroups > 1
        for jj=1:4
            data2{2*jj-1}=allData_merged(jj).QCFCVec{1};
        end
        for jj=1:4
            data2{2*jj}=allData_merged(jj).QCFCVec{2};
        end
    else
        data2 = {allData_merged(:).QCFCVec};
    end
    clear extraParams
    % extraParams.theLabels = {allData_merged(:).noiseOptionsNames};
    extraParams.customSpot = '';
    extraParams.add0Line = true;
    extraParams.theColors = theColors;
    BF_JitteredParallelScatter(data2,1,1,0,extraParams);
    ax = gca;

    % Set axis stuff
    ax.FontSize = FSize;
    % ax.XTick = [1:size(data{1},1)];
    ax.XTick = [];
    ax.XTickLabel = [];
    % ax.XLim = ([0 numPrePro+1]);
    ax.XLim = ([0 2*numPrePro+1]);
%     if ismember('NYU_2',WhichProject,'rows') | ismember('OCDPG',WhichProject,'rows')
%         ax.YLim = ([-0.8 1.25]);
%     else
        % ax.YLim = ([-0.6 1]);
        ax.YLim = ([-0.8 0.8]);
%     end
    ylabel('QC-FC (Pearson''s r)')

    % add text
    TextRotation = 0;
    strprec = '%0.2f';
    if numGroups > 1
        for jj=1:4
            data3{2*jj-1}=allData_merged(jj).QCFC_AbsMed(1);
        end
        for jj=1:4
            data3{2*jj}=allData_merged(jj).QCFC_AbsMed(2);
        end
    else
        data3 = {allData_merged(:).QCFC_AbsMed};
    end
    fi

    text(1:size(data3,2),repmat(ax.YLim(2) - ax.YLim(2)*.05,1,size(data3,2)),num2str([data3{1,:}]',strprec),... 
    'HorizontalAlignment','right',... 
    'VerticalAlignment','middle',...
    'Color','black',...
    'FontSize', FSize,...
    'Rotation',TextRotation)

    view(90,90)

    saveas(Fig_QCFC_Dist,fullfile(pwd,'QCFC_Dist_merged.png'));
    export_fig(fullfile(pwd,'QCFC_Dist_merged_highres.tif'),'-r600','-painters','-nocrop');
    
    
   % ------------------------------------------------------------------------------
    % QC-FC corrected
    % ------------------------------------------------------------------------------
    % Create data
    % data = {[allData(:).QCFC_PropSig_corr]'};
    % data_std = cell(1,length(data)); [data_std{:}] = deal(zeros(size(data{1})));
    if numGroups > 1
        QCFC_PropSig_corr=[allData_merged(:).QCFC_PropSig_corr];
        data = {[[QCFC_PropSig_corr(1:4)],[QCFC_PropSig_corr(5:8)]]'};
        data_std = cell(1,length(data)); [data_std{:}] = deal(zeros(size(data{1})));
    else
        data = {[allData_merged(:).QCFC_PropSig_corr]'};
        data_std = cell(1,length(data)); [data_std{:}] = deal(zeros(size(data{1})));
    end

    % Create table
    % T = table(data{1},'RowNames',{allData_merged(:).noiseOptionsNames}','VariableNames',{'QCFC_PropSig_corr'})

    % Create bar chart
    clear extraParams
    % extraParams.xTickLabels = {allData_merged(:).noiseOptionsNames};
    extraParams.xTickLabels = xy;
    extraParams.xLabel = ''; % 'Pipeline'
    extraParams.yLabel = 'QC-FC FDR-corr. (%)';
    extraParams.theColors = theColors;
    extraParams.theLines = theLines;
    extraParams.yLimits = [0 110];

    Fig_QCFC_FDR = figure('color','w', 'units', 'centimeters', 'pos', [0 0 10.5 9], 'name',['Fig_QCFC_FDR']); box('on'); movegui(Fig_QCFC_FDR,'center');
    sp = subplot(1,1,1);
    pos = get(sp,'Position');
    set(gca,'Position',[pos(1)*4.9, pos(2)*1.2, pos(3)*0.425, pos(4)*1]); % [left bottom width height]

    TheBarChart(data,data_std,false,extraParams)

    saveas(Fig_QCFC_FDR,fullfile(pwd,'QCFC_FDR_merged.bmp'));
    export_fig(fullfile(pwd,'QCFC_FDR_merged_highres.tif'),'-r600','-painters','-nocrop');
    
    %--------------------------
    % QCFC DistDep
    % data_std = cell(1,length(data)); [data_std{:}] = deal(zeros(size(data{1})));
    if numGroups > 1
        QCFC_DistDep=[allData_merged(:).QCFC_DistDep];
        data = {[[QCFC_DistDep(1:4)],[QCFC_DistDep(5:8)]]'};
        data_std = cell(1,length(data)); [data_std{:}] = deal(zeros(size(data{1})));
    else
        data = {[allData_merged(:).QCFC_DistDep]'};
        data_std = cell(1,length(data)); [data_std{:}] = deal(zeros(size(data{1})));
    end
    
    % Create table
    % T = table(data{1},'RowNames',{allData_merged(:).noiseOptionsNames}','VariableNames',{'QCFC_DistDep'})

    % Create bar chart
    extraParams.xTickLabels = xy;
    extraParams.xLabel = '';
    % extraParams.yLabel = 'QC-FC distance dependence (Spearman''s rho)';
    extraParams.yLabel = 'QC-FC distance dependence';
    extraParams.theColors = theColors;
    extraParams.theLines = theLines;
    extraParams.yLimits = [-0.5 0.5];

    Fig_QCFC_DistDep = figure('color','w', 'units', 'centimeters', 'pos', [0 0 10.5 9], 'name',['Fig_QCFC_DistDep']); box('on'); movegui(Fig_QCFC_DistDep,'center');
    sp = subplot(1,1,1);
    pos = get(sp,'Position');
    set(gca,'Position',[pos(1)*4.9, pos(2)*1.2, pos(3)*0.425, pos(4)*1]); % [left bottom width height]

    TheBarChart(data,data_std,false,extraParams)

    saveas(Fig_QCFC_DistDep,fullfile(pwd,'QCFC_DistDep_merged.bmp'));
    export_fig(fullfile(pwd,'QCFC_DistDep_merged_highres.tif'),'-r600','-painters','-nocrop');