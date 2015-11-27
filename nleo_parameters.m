DATA_DIR='~/ucc/software/preterm_IBI_detector/data/';
PIC_DIR='~/ucc/software/NLEOs/pics/';
D_NLEO_SIM=[DATA_DIR 'nleo_testing/'];

DO_MA_FILTER=0;
WIN_LENGTH=1.5; % MA window length (in seconds)

% low-pass filter the data (both preterms and normterms data); 
% in Hz, zero to turn off
LOW_PASS_DATA=20; 
HIGH_PASS_DATA=0;


% $$$ lcolor={[0 47 47]./256,[105 0 17]./256,[40 138 207]./256, ...
% $$$         [242 97 1]./256,[84 84 84]./256};
ll=lines;
lcolor={ll(1,:),ll(2,:),ll(3,:),ll(4,:),ll(5,:),ll(6,:)};


%---------------------------------------------------------------------
% the REST is not used (for EMBC-2014 paper)
%---------------------------------------------------------------------
% $$$ return;

%---------------------------------------------------------------------
% parameters for the methods
%---------------------------------------------------------------------
% now testing four different methods:
% a) Teager-Kaiser            'teager'
% b) Agarwal-Gotman           'agarwal'
% c) envelope-derivative      'envelope_diff'
% b) positive Agarwal-Gotman  'palmu'
NLEO_params_opt.normterms.teager          =[];  %[]; %[0,3.5];
NLEO_params_opt.normterms.agarwal         =[];  %[]; %[0,3.5];
NLEO_params_opt.normterms.envelope_diff   =[];  %[]; %[0,3.5];
NLEO_params_opt.normterms.palmu           =[];  %[];
NLEO_params_opt.normterms.abs_teager      =[];
NLEO_params_opt.normterms.env_only        =[];

NLEO_params_opt.preterms.teager          =[];
NLEO_params_opt.preterms.agarwal         =[];
NLEO_params_opt.preterms.envelope_diff   =[];
NLEO_params_opt.preterms.palmu           =[];  
NLEO_params_opt.preterms.abs_teager      =[];
NLEO_params_opt.preterms.env_only        =[];



%---------------------------------------------------------------------
% range of parameters (if using cross-validation)
%---------------------------------------------------------------------
% $$$ CROSS_VALIDATION=0;
% $$$ TRAIN_ALL_DATA=0;
% $$$ 
% $$$ if(CROSS_VALIDATION), TRAIN_ALL_DATA=0; end
% $$$ 
% $$$ % for term data:
% $$$ NLEO_params.normterms.teager       ={ [0:0.2:1],  [1:2:12] };  
% $$$ NLEO_params.normterms.agarwal      ={ [0:0.2:1],  [1:2:12] };   
% $$$ NLEO_params.normterms.envelope_diff={ [0:0.2:1],  [1:2:12] };   
% $$$ NLEO_params.normterms.palmu        ={ [0:0.2:1],  [1:1:4] };   
% $$$ NLEO_params.normterms.env_only     ={ [0:0.2:1],  [1:1:4] };   
% $$$ 
% $$$ % $$$ NLEO_mod_params_range.normterms  ={ [0:0.2:1],  [1:1:4] };
% $$$ % $$$ NLEO_kirsi_params_range.normterms={ [0:0.2:1],  [1:1:4] };
% $$$ 
% $$$ % for preterm data:
% $$$ % $$$ NLEO_params.preterms.down_samp      ={ [0:0.1:0.5],  [1:0.5:14] };
% $$$ % $$$ NLEO_params.preterms.plain_nleo     ={ [0:0.1:0.5],  [1:1:14] };
% $$$ % $$$ NLEO_params.preterms.hilbert_diff_op={ [0:0.1:0.5],  [1:1:14] };
% $$$ % $$$ NLEO_params.preterms.env_only       ={ [0:0.1:0.5],  [1:0.5:14] };
% $$$ 
% $$$ NLEO_params.preterms.teager       ={ [0:0.2:2],  [2:2:8] };  
% $$$ NLEO_params.preterms.agarwal      ={ [0:0.2:2],  [2:2:8] };
% $$$ NLEO_params.preterms.envelope_diff={ [0:0.2:2],  [2:2:8] };   
% $$$ NLEO_params.preterms.palmu        ={ [0:0.2:2],  [2:2:8] };   
% $$$ NLEO_params.preterms.env_only     ={ [0:0.2:2],  [2:2:8] };   
% $$$ 
% $$$ 
% $$$ % $$$ NLEO_mod_params_range.preterms  ={ [0:0.1:0.5],  [1:0.5:14] };
% $$$ % $$$ NLEO_kirsi_params_range.preterms={ [0:0.1:0.5],  [1:0.5:14] };
% $$$ 
% $$$ % $$$ NLEO_mod_params_range.preterms  ={ [0:0.2:0.5],  [1:1:14] };
% $$$ % $$$ NLEO_kirsi_params_range.preterms={ [0:0.2:0.5],  [1:1:14] };
% $$$ 
% $$$ 
% $$$ 
% $$$ 
