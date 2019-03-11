d0= '/nifti_toolboxdir/';
d1='/control_scans_dir/'; %choose path to where control lesions are.
d2= '/patient_scans_dir/'; %choose path to where patient lesions are.
addpath(d0); addpath(d1); addpath(d2);
filelistC = dir([d1 '/*.nii']);
filelistP = dir([d2 '/*.nii']);

img_binsumC=zeros(182,218,182);
img_binsumP=zeros(182,218,182);


min_overlap=16; %min overlap for binarized image
thr=.0005; %Normalized tract threshold below which to 0 values
[hdr,T1_1mm_mask]=read_nifti('YourMask.nii');



for a=1:length(filelistC)
    [hdr,imgNC]= read_nifti(filelistC(a).name);
    [hdr2,imgNP]= read_nifti(filelistP(a).name);
    %binarize images
    img_binC=imgNC;
    %threshold img
    img_binC(img_binC<=thr)=0;
    img_binC(img_binC>0)=1;
    img_binsumC=img_binC+img_binsumC;
    img_binP=imgNP;
    %threshold img
    img_binP(img_binP<=thr)=0;
    img_binP(img_binP>0)=1;
    img_binsumP=img_binP+img_binsumP;
end
bin_ind=find((img_binsumC>=min_overlap)&(T1_1mm_mask>0));

%now generate the ttest matrix
for a=1:length(filelistC)
    [hdr,imgNC]= read_nifti(filelistC(a).name);
    [hdr2,imgNP]= read_nifti(filelistP(a).name);
    valsC(a,:)=imgNC(bin_ind);
    valsP(a,:)=imgNP(bin_ind);
end

%calculate wilcoxan at each voxel
for a=1:size(valsC,2)
    [p,hR,statsR]=ranksum(valsC(:,a),valsP(:,a));
    Z_true(a)=statsR.zval;
    RS_true(a)=statsR.ranksum;
    pR_true(a)=p;
end


%PERMUTATION BASED, non-parametric CORRECTION FOR MULT COMPARISONS
%FIRST GENERATE PERMAUTIONS FOR PERM TEST
nperms=1000;  %set number of permuattions to generate, seeing as 32Choose16 is prohibitively large
bins=50;%set number of batches of t val maxes to calculate (so as to avoid generating a memory-killing giant matrix)
allvals=[valsC;valsP];
n_C=size(valsC,1);
n_P=size(valsP,1);
valsC_permd=zeros(n_C,size(allvals,2));
valsP_permd=zeros(n_P,size(allvals,2));
count=0;%count variable
Zval=zeros(1,size(valsC,2));
RS=zeros(1,size(valsP,2));
RS_max=zeros(1,nperms);
Z_max=zeros(1,nperms);
%permute the indices of patient/control
for a=1:nperms
    label_ind=randperm(32);
    vals_permd=allvals(label_ind,:);
    valsC_permd=vals_permd(1:16,:);
    valsP_permd=vals_permd(17:end,:);
    for b=1:size(valsC,2)
    [p1,h1,stats1]=ranksum(valsC_permd(:,b),valsP_permd(:,b));
    Zval(b)=stats1.zval;
    RS(b)=stats1.ranksum;
    end
    RS_max(a)=max(RS);
    Z_max(a)=max(abs(Zval));
    count=count+1
end
RS_max=[RS_max max(RS_true)];
[f,x]=ecdf(RS_max);
crit_val=min(x(f>=.95));
RS_img_perm=zeros(91,109,91);
%plot on img RS/Z values greater than threshold
supra_thresh_RS=RS_true((RS_true>crit_val)&(Z_true>0));
% t_img_perm(bin_ind(abs(tvals)>crit_val))=supra_thresh_t;
% write_nifti(hdr,t_img_perm,'t_img_perm.nii');

%correct the p values, as below:
%what is probability that a random t value is at least as extreme as the one from
%our image
for a=1:length(supra_thresh_RS)
    inv_p(a)=min(f(x>=supra_thresh_RS(a)));
end
p_inv_img_perm=zeros(182,218,182);
p_inv_img_perm(bin_ind(RS_true>crit_val))=inv_p;
%write to nifti
write_nifti(hdr,p_inv_img_perm,'p_inv_img_permsig.nii.gz');


%CORRECTION FOR MULTIPLE COMPARISONS, PARAMETRIC
%benjamini hochberg FDR correction
[p_sort,temp_ind]=sort(pR_true);
p_sort_adj=p_sort.*(repmat(length(bin_ind),1,length(bin_ind))./[1:length(bin_ind)]);
for a=2:max(sort_adj_less05_ind)
    if p_sort_adj(a)<p_sort_adj(a-1)
        p_sort_adj(p_sort_adj(1:a-1)>p_sort_adj(a))=p_sort_adj(a);
    end
end
sort_adj_less05_ind=find(p_sort_adj<.05);
%write to nifti
p_img_HochFDR=zeros(182,218,182);
p_img_HochFDR(bin_ind(temp_ind(sort_adj_less05_ind)))=p_sort_adj((sort_adj_less05_ind));
