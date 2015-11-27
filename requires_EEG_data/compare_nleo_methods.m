%-------------------------------------------------------------------------------
% compare_nleo_methods: Compare frequency-weighted energy measures (nonlinear energy operators)
%
% Syntax: [rocs_compare_st,corr_cc]=compare_nleo_methods(data_type,method1,method2,fout)
%
% Inputs: 
%     data_type  - EEG data type either 'normterms' or 'preterms'
%     method1,method2  - methods to compare from the list:
%                        a) Teager-Kaiser            'teager'  
%                        b) Agarwal-Gotman           'agarwal' 
%                        c) envelope-derivative      'envelope_diff'
%                        b) positive Agarwal-Gotman  'palmu'   
%                        e) envelope                 'env_only'
%     fout - handle to write the output (default is stdout)
%
% Outputs: 
%     [rocs_compare_st,corr_cc] - 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 04-04-2014
%
% last update: Time-stamp: <2014-06-05 17:36:01 (otoolej)>
%-------------------------------------------------------------------------------
function [roc_summary_stats,corr_cc]=compare_nleo_methods(data_type,method1,method2,fout,MA_FILT)
if(nargin<1 || isempty(data_type)), data_type='normterms'; end
if(nargin<2 || isempty(method1)), method1='teager'; end
if(nargin<3 || isempty(method2)), method2='envelope_diff'; end
if(nargin<4 || isempty(fout)), fout=1; end
if(nargin<5 || isempty(MA_FILT)), MA_FILT=[]; end



corr_cc=[];
DBplot=0; 
CROSS_VALIDATION=0;
TRAIN_ALL_DATA=0;

%---------------------------------------------------------------------
% 0. method's parameters are in here:
%---------------------------------------------------------------------
nleo_parameters;


%---------------------------------------------------------------------
% 1. load data from .mat file
%---------------------------------------------------------------------
b=load([D_NLEO_SIM 'eeg_segments_' data_type '.mat']);
x_data=b.eeg_data; x_mask=b.anno; Fs=b.Fs{1};
L=length(x_data);


roc_mod_st(L)=struct('auc',[],'ss_eq',[],'thres_eq',[]);
roc_kirsi_st(L)=struct('auc',[],'ss_eq',[],'thres_eq',[]);


if(LOW_PASS_DATA>0)
    if(HIGH_PASS_DATA==0), HIGH_PASS_DATA=[]; end
    for n=1:L
        x_data{n}=do_bandpass_filtering(x_data{n},Fs,LOW_PASS_DATA,HIGH_PASS_DATA);
    end
end

% override whats in the 'nleo_parameters' file:
if(~isempty(MA_FILT)),  DO_MA_FILTER=MA_FILT; end 


ma_len=[];
if(DO_MA_FILTER), ma_len=WIN_LENGTH*Fs; end


%---------------------------------------------------------------------
% IF doing cross-validation (optimise parameters)
%---------------------------------------------------------------------
if(CROSS_VALIDATION)
    fw_fprintf(fout,'| %s | %s | difference |\n',method1,method2);
    for n=1:L
        if(iscell(x_data))
            x_train=x_data; x_train(n)=[];
            mask_train=x_mask; mask_train(n)=[];
            
            x_test=x_data{n}; mask_test=x_mask{n};
        else
            % not yet implemented:
            x=x_data(n,:); mask=x_mask(n,:);
        end
        
        % a. find best parameters on training set:
        [mod_params,kirsi_params]=cross_validation_training(x_train,mask_train,Fs, ...
                                                          method1,method2,data_type,fout);
        
        % b. test on left-one-out:
        x_nleo_mod  =cal_freqweighted_energy(x_test,Fs,method1,ma_len,mod_params);
        x_nleo_kirsi=cal_freqweighted_energy(x_test,Fs,method2,ma_len,kirsi_params);        
        
% $$$         [x_nleo_mod,Fs_down]=hfunc_alt_method(x_test,Fs,DO_MA_FILTER, ...
% $$$                                  mod_params,TEAGER,0);
% $$$         [x_nleo_kirsi,Fs_down_k]=plain_nleo(x_test,Fs,DO_MA_FILTER,kirsi_params);
% $$$         mask_downsampled=decimate_mask(mask_test,Fs_down,Fs);
% $$$         mask_k=decimate_mask(mask_test,Fs_down_k,Fs);        
        
        % c. testing results:
        x_nleo_kirsi(isnan(mask_test))=NaN;
        x_nleo_mod(isnan(mask_test))=NaN;    
        
        roc_mod_st(n)=assess_performance(x_nleo_mod,mask_test);
        roc_kirsi_st(n)=assess_performance(x_nleo_kirsi,mask_test);
        
        fw_fprintf(fout,'TESTING: | %g | %g | %g |\n',roc_mod_st(n).auc, ...
                roc_kirsi_st(n).auc,roc_mod_st(n).auc-roc_kirsi_st(n).auc);
    end
    
else
    
    %---------------------------------------------------------------------
    % .... or training on all data (to find best parameters):
    %---------------------------------------------------------------------
    if(TRAIN_ALL_DATA)
        [mod_params,kirsi_params]=cross_validation_training(x_data,x_mask,Fs, ...
                                                          method1,method2,data_type,fout);
        
        fw_fprintf(fout,'OPT params, mod: [%g,%g]\n',mod_params(1),mod_params(2));
        fw_fprintf(fout,'OPT params, mod: [%g,%g]\n',kirsi_params(1),kirsi_params(2));
    else
        mod_params=NLEO_params_opt.(data_type).(method1);
        kirsi_params=NLEO_params_opt.(data_type).(method2);        
    end
    

    %---------------------------------------------------------------------
    % .... or just test using given parameters:
    %---------------------------------------------------------------------
    fw_fprintf(fout,'| %s | %s | difference |\n',method1,method2);
    corr_cc=zeros(1,L);
    for n=1:L
        if(iscell(x_data))
            x=x_data{n}; mask=x_mask{n};
        else
            x=x_data(n,:); mask=x_mask(n,:);
        end

        x_nleo_mod  =cal_freqweighted_energy(x,Fs,method1,ma_len,mod_params);
        x_nleo_kirsi=cal_freqweighted_energy(x,Fs,method2,ma_len,kirsi_params);        
        
% $$$         [x_nleo_mod,Fs_down]=hfunc_alt_method(x,Fs,DO_MA_FILTER,mod_params,TEAGER,DBplot);
% $$$         x_nleo_kirsi=plain_nleo(x,Fs,DO_MA_FILTER,kirsi_params);
% $$$         mask_downsampled=decimate_mask(mask,[],Fs);        
        
        x_nleo_kirsi(isnan(mask))=NaN;
        x_nleo_mod(isnan(mask))=NaN;    

        
        corr_cc(n)=do_corr(x_nleo_kirsi,x_nleo_mod);
        
        
        if(DBplot)
            figure(1); clf; hold all;

% $$$         plotyy(1:N,x_nleo_kirsi,1:N,x_nleo_mod);
% $$$         plot(x_nleo_kirsi); plot(x_nleo_mod);        
            N=length(x_nleo_kirsi);
            t1=(0:N-1)./Fs; 
% $$$             if(~isempty(Fs_down))
% $$$                 t2=(0:length(x_nleo_mod)-1)./Fs_down;
% $$$             else
            t2=(0:length(x_nleo_mod)-1)./Fs;
% $$$             end
            h1=subplot(211); hold all;
            plot(t1,x_nleo_kirsi./max(x_nleo_kirsi));
            if(isreal(x_nleo_mod))
                plot(t2,x_nleo_mod./max(x_nleo_mod));
            else
                plot(t2,real(x_nleo_mod)./max(abs(x_nleo_mod)));
                plot(t2,imag(x_nleo_mod)./max(abs(x_nleo_mod)));            
            end
            
            m=max(max(x_nleo_mod),max(x_nleo_kirsi));                
            plot(t1,mask);
            legend(method1,method2,'mask');        
            
            h2=subplot(212);
            plot(t1,x);
            linkaxes([h1 h2],'x');
            disp('--- paused; hit key to continue ---'); pause;
        end

        
        roc_mod_st(n)=assess_performance(x_nleo_mod,mask);
        roc_kirsi_st(n)=assess_performance(x_nleo_kirsi,mask);
        
        fw_fprintf(fout,'| %g | %g | %g |\n',roc_mod_st(n).auc, ...
                roc_kirsi_st(n).auc,roc_mod_st(n).auc-roc_kirsi_st(n).auc);
    end
end



median_mod_auc=median([roc_mod_st.auc]);
median_kirsi_auc=median([roc_kirsi_st.auc]);
mean_mod_auc=mean([roc_mod_st.auc]);
mean_kirsi_auc=mean([roc_kirsi_st.auc]);

iqr_mod_auc=iqr([roc_mod_st.auc]);
iqr_kirsi_auc=iqr([roc_kirsi_st.auc]);


roc_summary_stats.median_methodOne_auc=median_mod_auc;
roc_summary_stats.mean_methodOne_auc=mean_mod_auc;
roc_summary_stats.median_methodTwo_auc=median_kirsi_auc;
roc_summary_stats.mean_methodTwo_auc=mean_kirsi_auc;
roc_summary_stats.iqr_methodOne_auc=prctile([roc_mod_st.auc],[25,75]);
roc_summary_stats.iqr_methodTwo_auc=prctile([roc_kirsi_st.auc],[25,75]);


fw_fprintf(fout,'-------------------------------------\n');
fw_fprintf(fout,'median (IQR): | %g (%g) | %g (%g)|\n', ...
        median_mod_auc,iqr_mod_auc,median_kirsi_auc,iqr_kirsi_auc);
fw_fprintf(fout,'mean (IQR): | %g (%g) | %g (%g)|\n', ...
        mean_mod_auc,iqr_mod_auc,mean_kirsi_auc,iqr_kirsi_auc);


% paired t-test:
[h,p]=ttest([roc_mod_st.auc]-[roc_kirsi_st.auc]);
fw_fprintf(fout,'|t-test: [%g,%g]\n',h,p);

effect_size_st=mes([[roc_mod_st.auc]-[roc_kirsi_st.auc]]',0,'g1');
fw_fprintf(fout,'effect size (CI): %g (%g -- %g)\n',effect_size_st.g1,...
        effect_size_st.g1Ci(1),effect_size_st.g1Ci(2));


if(~isempty(corr_cc))
    fw_fprintf(fout,'Correlation between %s and %s is [range]: %g [%g--%g]\n',...
               method1,method2,mean(corr_cc),min(corr_cc),max(corr_cc));
end


    
function roc_st=assess_performance(nleo_stat,mask)
%---------------------------------------------------------------------
% generate ROC and metrics from that
%---------------------------------------------------------------------
[R,auc,ss_eq,thres_eq]=do_time_roc(nleo_stat,mask);
sens=R(2,:); spec=1-R(1,:);

roc_st.auc=auc;
roc_st.ss_eq=ss_eq;
roc_st.thres_eq=thres_eq;    
    


function x_upsample=int_upsample(x,Fs,Fs_new)
%---------------------------------------------------------------------
% upsample signal by increase sampling frequency by factor a (a integer)
%---------------------------------------------------------------------
a=Fs_new/Fs;
if( ~rem(Fs,Fs_new) ), error('not integer sampling relation.'); end
    

N=length(x); M=a*N; Nh=ceil(N/2);
X=fft(x);

X_up=zeros(1,M);
X_up(1:Nh)=X(1:Nh);
X_up(M:-1:(M-Nh+2))=conj(X(2:Nh));
x_upsample=a.*ifft(X_up);

% $$$ x_upsample=interp1(1:a:M,x,1:M);
% $$$ figure(80); clf; hold all;
% $$$ plot((0:N-1)./Fs,x,'-+'); plot((0:M-1)./Fs_new,x_upsample,'-o');
% $$$ disp('--- paused; hit key to continue ---'); pause;


function mask_downsampled=decimate_mask(mask_test,Fs_down,Fs)
%---------------------------------------------------------------------
% downsample mask
%---------------------------------------------------------------------
if(~isempty(Fs_down) & Fs_down~=Fs)
    mask_downsampled=round( resample(mask_test,Fs_down,Fs) );
else
    mask_downsampled=mask_test;
end

mask_downsampled(mask_downsampled>1)=1;
mask_downsampled(mask_downsampled<0)=0;    

if( any(mask_downsampled>1) | any(mask_downsampled<0) )
    '>>>?'
    keyboard;
end

    
function cc=do_corr(x,y)
%---------------------------------------------------------------------
% remove NaNs before doing correlation
%---------------------------------------------------------------------
a=find(isnan(x)); b=find(isnan(y)); inans=union(a,b);

if(~isempty(inans))
    x(inans)=[];  y(inans)=[]; 
end

if( size(x,1)==1 ), x=x.'; end
if( size(y,1)==1 ), y=y.'; end

cc=corr(x,y);


function [mod_best_params,kirsi_best_params]=...
    cross_validation_training(x,mask,Fs,method1,method2,data_type,fout)
%---------------------------------------------------------------------
% find best parameters on training data
%---------------------------------------------------------------------

nleo_parameters;
L=length(x);


ma_len=[];
if(DO_MA_FILTER), ma_len=WIN_LENGTH*Fs; end


mod_params_all=NLEO_params.(data_type).(method1);
kirsi_params_all=NLEO_params.(data_type).(method2);

P=length(mod_params_all{1});
Q=length(mod_params_all{2});
if(length(mod_params_all)>2)
    R=length(mod_params_all{3});
else
    R=1;
end

N=length(x);

if( P~=length(kirsi_params_all{1}) | ...
    Q~=length(kirsi_params_all{2}) )
    error('parameter sets need to be same size');
end

if(R>1)
    mod_auc=zeros(P,Q,R,N); kirsi_auc=zeros(P,Q,R,N);
else
    mod_auc=zeros(P,Q,N); kirsi_auc=zeros(P,Q,N);
end


for n=1:N
    fprintf('-%d',n);
    x_train=x{n}; mask_train=mask{n};

    for r=1:R
        for p=1:P
            for q=1:Q
                mod_params=[mod_params_all{1}(p) mod_params_all{2}(q)];
                kirsi_params=[kirsi_params_all{1}(p) kirsi_params_all{2}(q)];
                if(R>1)
                    mod_params=[mod_params_all mod_params_all{3}(r)];
                    kirsi_params=[kirsi_params_all kirsi_params_all{3}(r)];
                end

                x_nleo_mod  =cal_freqweighted_energy(x_train,Fs,method1,ma_len,mod_params);
                x_nleo_kirsi=cal_freqweighted_energy(x_train,Fs,method2,ma_len,kirsi_params); 
                    
% $$$                 [x_nleo_mod,Fs_down]=hfunc_alt_method(x_train,Fs,DO_MA_FILTER,mod_params,TEAGER);
% $$$                 [x_nleo_kirsi,Fs_down_k]=plain_nleo(x_train,Fs,DO_MA_FILTER,kirsi_params);
% $$$                 
% $$$                 mask_downsampled=decimate_mask(mask_train,Fs_down,Fs);        
% $$$                 mask_k=decimate_mask(mask_train,Fs_down_k,Fs);                            

                x_nleo_kirsi(isnan(mask_train))=NaN;
                x_nleo_mod(isnan(mask_train))=NaN;    

                roc_mod_st=assess_performance(x_nleo_mod,mask_train);
                roc_kirsi_st=assess_performance(x_nleo_kirsi,mask_train);

                if(R>1)
                    mod_auc(p,q,r,n)=roc_mod_st.auc;
                    kirsi_auc(p,q,r,n)=roc_kirsi_st.auc;        
                else
                    mod_auc(p,q,n)=roc_mod_st.auc;
                    kirsi_auc(p,q,n)=roc_kirsi_st.auc;        
                end
            end
        end
    end
end
fprintf('\n');

if(R>1), dim=4; else dim=3; end

mod_auc_mean=mean(mod_auc,dim);
kirsi_auc_mean=mean(kirsi_auc,dim);


figure(33); clf; hold all;
plot(mod_params_all{2},mod_auc_mean.','-o'); 
plot(mod_params_all{2},kirsi_auc_mean.','-+');
drawnow;

figure(34); clf; hold all;
imagesc(mod_params_all{1},mod_params_all{2}, mod_auc_mean.');
colorbar; axis('tight');
figure(35); clf; hold all;
imagesc(kirsi_params_all{1},kirsi_params_all{2}, kirsi_auc_mean.');
colorbar; axis('tight');
drawnow;


[m,im]=max(mod_auc_mean(:)); 

if(R>1)
    [a,b,c]=ind2sub(size(mod_auc_mean),im);
    mod_best_params=[mod_params_all{1}(a) mod_params_all{2}(b) mod_params_all{3}(c)];
else
    [a,b]=ind2sub(size(mod_auc_mean),im);
    mod_best_params=[mod_params_all{1}(a) mod_params_all{2}(b)];
end

[m,im]=max(kirsi_auc_mean(:)); 
if(R>1)
    [a,b,c]=ind2sub(size(kirsi_auc_mean),im);
    kirsi_best_params=[kirsi_params_all{1}(a) kirsi_params_all{2}(b) kirsi_params_all{3}(c)];
else
    [a,b]=ind2sub(size(kirsi_auc_mean),im);
    kirsi_best_params=[kirsi_params_all{1}(a) kirsi_params_all{2}(b)];
end


if(R>1)
    p_mod_str=sprintf('[%g,%g,%g]\n',mod_best_params(1),mod_best_params(2),...
                  mod_best_params(3));
    p_kirsi_str=sprintf('[%g,%g,%g]\n',kirsi_best_params(1),kirsi_best_params(2),...
                  kirsi_best_params(3));
else
    p_mod_str=sprintf('[%g,%g]\n',mod_best_params(1),mod_best_params(2));
    p_kirsi_str=sprintf('[%g,%g]\n',kirsi_best_params(1),kirsi_best_params(2));    
end


fw_fprintf(fout,'OPT: mod-NLEO best parameters=%s',p_mod_str);
fw_fprintf(fout,'OPT: Kirsi-NLEO best parameters=%s',p_kirsi_str);

