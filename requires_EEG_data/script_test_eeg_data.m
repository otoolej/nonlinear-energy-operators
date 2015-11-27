%-------------------------------------------------------------------------------
% script_test_eeg_data: Test the different methods with the EEG data and write results
% to file. 
%
% **NB: EEG data set is REQUIRED for this to run.**
%
% Syntax: []=script_test_eeg_data()
%
% Inputs: 
%      PRINT_PLOTS - print or not, either 0 or 1; 
%
% Example:
%     
%

% John M. O' Toole, University College Cork
% Started: 04-04-2014
%
% last update: Time-stamp: <2014-06-05 17:50:09 (otoolej)>
%-------------------------------------------------------------------------------
function []=script_test_eeg_data(PRINT_PLOTS)
if(nargin<1 || isempty(PRINT_PLOTS)), PRINT_PLOTS=0; end

nleo_parameters;

if(~exist(DATA_DIR,'dir'))
    error('EEG data is needed for this function.');
end
    

methods_all={'teager','agarwal','abs_teager','envelope_diff','palmu','env_only'};
fout=fopen([PIC_DIR 'eeg_results_out.org'],'w');

figure(1); clf; hold all;
fprintf(fout,'* NORMTERMS\n');
fprintf('\n\n *** NORMTERMS *** \n\n');
hax(1)=do_one_dataset('normterms',methods_all,fout,1);

fprintf(fout,'\n\n* PRETERMS\n');
fprintf('\n\n *** PRETERMS *** \n\n');
hax(2)=do_one_dataset('preterms',methods_all,fout,2);

fclose(fout);


if(PRINT_PLOTS)
    set_gca_fonts('Arial',12,hax(1));
    set_gca_fonts('Arial',12,hax(2));    
    print2eps([PIC_DIR 'burst_detection_results.eps']);
end


function hax=do_one_dataset(data_type,methods_all,fout,figno)
%---------------------------------------------------------------------
% 
%---------------------------------------------------------------------
L=length(methods_all);

for n=1:L
    fprintf(fout,'| %s\t',methods_all{n});
end
fprintf(fout,'|\n');
fprintf(fout,'|------------------------------------------------------------|\n');

LWIDTH=8; lshift=0.3;
lcolor={[34 85 51],[68 187 204],[136 221 221],[187 238 255],[0 85 187]};
lcolor={[0 47 47],[105 0 17],[40 138 207],[242 97 1],[84 84 84],[84 84 84]};
ll=lines;
ll=ll.*256;
lcolor={ll(1,:),ll(2,:),ll(3,:),ll(4,:),ll(5,:),ll(6,:)};

hax=subtightplot(2,1,figno,[0.07,0.2],[0.15,0.1],[0.1,0.1]);
hold all;

for p=1:2
    for n=1:2:(L-rem(L,2))
        method1=methods_all{n};   method2=methods_all{n+1};

        roc_st=compare_nleo_methods(data_type,method1,method2,1,p-1);

        fprintf(fout,'| %0.2f (%0.2f--%0.2f) | %0.2f (%0.2f--%0.2f) ', ...
                roc_st.median_methodOne_auc, roc_st.iqr_methodOne_auc(1), ...
                roc_st.iqr_methodOne_auc(2),roc_st.median_methodTwo_auc, ...
                roc_st.iqr_methodTwo_auc(1),roc_st.iqr_methodTwo_auc(2));
        
        pshift=(p-1).*8;
        hl(n)=plot(ones(1,20).*(pshift+n),linspace(roc_st.iqr_methodOne_auc(1),...
                                                roc_st.iqr_methodOne_auc(2),20));
        hp(n)=line([pshift+n-lshift pshift+n+lshift],[roc_st.median_methodOne_auc ...
                            roc_st.median_methodOne_auc]);        
        
        hl(n+1)=plot(ones(1,20).*(pshift+n+1),linspace(roc_st.iqr_methodTwo_auc(1),...
                                                roc_st.iqr_methodTwo_auc(2),20));
        hp(n+1)=line([pshift+n+1-lshift pshift+n+1+lshift],[roc_st.median_methodTwo_auc ...
                            roc_st.median_methodTwo_auc]);        
    end
    if(rem(L,2)~=0)
        method1=methods_all{L};   method2=methods_all{L-1};

        roc_st=compare_nleo_methods(data_type,method1,method2,1,p-1);

        fprintf(fout,'| %0.2f (%0.2f--%0.2f) ',roc_st.median_methodOne_auc, ...
                roc_st.iqr_methodOne_auc(1),roc_st.iqr_methodOne_auc(2));
        
        hl(L)=plot(ones(1,20).*(pshift+L),linspace(roc_st.iqr_methodOne_auc(1),...
                                                roc_st.iqr_methodOne_auc(2),20));
        hp(L)=line([pshift+L-lshift pshift+L+lshift],[roc_st.median_methodOne_auc ...
                            roc_st.median_methodOne_auc]);        
    end

    set(hp,'color','k','linewidth',2);    
    for q=1:L
        set(hl(q),'linewidth',LWIDTH,'color',lcolor{q}./256);        
    end
    
    fprintf(fout,'|\n');
end

if(figno==1)
    ylim([0.45 1]);
else
    ylim([0.5 1]);        
end
set(gca,'ytick',[0.5:0.1:1]);
set(gca,'xtick',[]);
ylabel('AUC');

set(hax,'xcolor',[1 1 1]);

if(figno==1)
    hleg=legend(hl,'a: Teager-Kaiser','b: Agarwal-Gotman','c: abs. Teager-Kaiser', ...
                'd: abs. Agarwal-Gotman','e: envelope-derivative', ...
                'f: instantaneous energy','location','southeast');
    legend('boxoff');

    yloc=0.50;  yloc2=0.97;
    text(1,yloc,'a.','horizontalAlignment','center'); 
    text(9, yloc2,'a.','horizontalAlignment','center'); 
    text(2,yloc,'b.','horizontalAlignment','center'); 
    text(10,yloc2,'b.','horizontalAlignment','center');    
    text(3,yloc,'c.','horizontalAlignment','center'); 
    text(11,yloc2,'c.','horizontalAlignment','center');        
    text(4,yloc,'d.','horizontalAlignment','center'); 
    text(12,yloc2,'d.','horizontalAlignment','center');            
    text(5,yloc,'e.','horizontalAlignment','center'); 
    text(13,yloc2,'e.','horizontalAlignment','center');                
    text(6,yloc,'f.','horizontalAlignment','center'); 
    text(14,yloc2,'f.','horizontalAlignment','center');                
else
    yloc=0.54;   yloc2=0.98;
    text(1,yloc,'a.','horizontalAlignment','center'); 
    text(9, yloc2,'a.','horizontalAlignment','center'); 
    text(2,yloc,'b.','horizontalAlignment','center'); 
    text(10,yloc2,'b.','horizontalAlignment','center');    
    text(3,yloc,'c.','horizontalAlignment','center'); 
    text(11,yloc2,'c.','horizontalAlignment','center');        
    text(4,yloc,'d.','horizontalAlignment','center'); 
    text(12,yloc2,'d.','horizontalAlignment','center');            
    text(5,yloc,'e.','horizontalAlignment','center'); 
    text(13,yloc2,'e.','horizontalAlignment','center');                
    text(6,yloc,'f.','horizontalAlignment','center'); 
    text(14,yloc2,'f.','horizontalAlignment','center');                
    
    text(1.23,1.06,'no post-processing');
    text(8.55,1.06,'with moving-average filter');    
    text(1.23,0.47,'no post-processing');
    text(8.55,0.47,'with moving-average filter');    
end





