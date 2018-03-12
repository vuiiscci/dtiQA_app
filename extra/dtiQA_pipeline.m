try
    % Set job directory path
    job_dir_path = '/OUTPUTS/';
    
    % Set FSL path
    fsl_path = '/extra/fsl_5_0_10_eddy_5_0_11';
    
    % Set camino path
    camino_path = '/extra/camino';

    % Read config file
    dtiqa_config = read_config('/INPUTS/dtiQA.conf');
    
    % Set dwmri_info - this will set base path to nifti/bvec/bval, phase encoding direction, and readout times
    if isfield(dtiqa_config,'dwmri_info_base_path')
        if iscell(dtiqa_config.dwmri_info_base_path)
            for i = 1:length(dtiqa_config.dwmri_info_base_path)
                dwmri_info(i).base_path = dtiqa_config.dwmri_info_base_path{i}; %#ok<SAGROW>
            end
        else
            dwmri_info.base_path = dtiqa_config.dwmri_info_base_path;
        end
    end
    if isfield(dtiqa_config,'dwmri_info_pe_dir')
        if iscell(dtiqa_config.dwmri_info_pe_dir)
            for i = 1:length(dtiqa_config.dwmri_info_pe_dir)
                dwmri_info(i).pe_dir = dtiqa_config.dwmri_info_pe_dir{i};
            end
        else
            dwmri_info.pe_dir = dtiqa_config.dwmri_info_pe_dir;
        end
    end
    if isfield(dtiqa_config,'dwmri_info_scan_descrip')
        if iscell(dtiqa_config.dwmri_info_scan_descrip)
            for i = 1:length(dtiqa_config.dwmri_info_scan_descrip)
                dwmri_info(i).scan_descrip = dtiqa_config.dwmri_info_scan_descrip{i};
            end
        else
            dwmri_info.scan_descrip = dtiqa_config.dwmri_info_scan_descrip;
        end
    end
    if isfield(dtiqa_config,'dwmri_info_readout_time')
        if iscell(dtiqa_config.dwmri_info_readout_time)
            for i = 1:length(dtiqa_config.dwmri_info_readout_time)
                dwmri_info(i).readout_time = dtiqa_config.dwmri_info_readout_time{i};
            end
        else
            dwmri_info.readout_time = dtiqa_config.dwmri_info_readout_time;
        end
    end
    
    % Make sure at least one "scan" exists
    if ~any(strcmp('scan',{dwmri_info.scan_descrip}))
        error('At least one dwmri must be a "scan"');
    end

    % Test if first scan is a "b0"; if it is, set the first dwmri to be the
    % first "scan". This is because first diffusion image should correspond to
    % the first b0 for eddy.
    if strcmp(dwmri_info(1).scan_descrip,'b0')
        first_scan_idx = find(strcmp({dwmri_info.scan_descrip},'scan'),1);
        disp(['Setting "scan": ' dwmri_info(first_scan_idx).base_path ' ' ...
              'in front because first dwmri is "b0": '  ...
              dwmri_info(1).base_path '. First dwmri should be a scan for eddy.']);
        % Rearrange
        dwmri_info = dwmri_info([first_scan_idx 1:first_scan_idx-1 first_scan_idx+1:length(dwmri_info)]);
    end
    
    % BET params
    bet_params = dtiqa_config.bet_params;
        
    % ADC fix - apply it for Philips scanner
    ADC_fix = dtiqa_config.ADC_fix;
    
    % zero_bval_thresh - will set small bvals to zero
    zero_bval_thresh = dtiqa_config.zero_bval_thresh;
    
    % prenormalize - will prenormalize data prior to eddy
    prenormalize = dtiqa_config.prenormalize;
    
    % use all b0s for topup
    use_all_b0s_topup = dtiqa_config.use_all_b0s_topup;
    
    % topup params
    topup_params = dtiqa_config.topup_params;
    
    % Sometimes name of eddy is 'eddy', 'eddy_openmp', or 'eddy_cuda'
    eddy_name = dtiqa_config.eddy_name;
    
    % use b0s in eddy
    use_b0s_eddy = dtiqa_config.use_b0s_eddy;
    
    % eddy params
    eddy_params = dtiqa_config.eddy_params;
    
    % normalize - will normalize data and output a single B0
    normalize = dtiqa_config.normalize;
    
    % sort scans - will sort scans by b-value
    sort_scans = dtiqa_config.sort_scans;
    
    % Set number of threads (only works if eddy is openmp version)
    setenv('OMP_NUM_THREADS',num2str(dtiqa_config.OMP_NUM_THREADS));
    
    % Perform preprocessing
    [dwmri_path, bvec_path, bval_path, mask_path, movement_params_path, topup_eddy_pdf_path] = ...
        topup_eddy_preprocess(job_dir_path, ...
                              dwmri_info, ...
                              fsl_path, ...
                              ADC_fix, ...
                              zero_bval_thresh, ...
                              prenormalize, ...
                              use_all_b0s_topup, ...
                              topup_params, ...
                              eddy_name, ...
                              use_b0s_eddy, ...
                              eddy_params, ...
                              normalize, ...
                              sort_scans, ...
                              bet_params);
                              
    % CSF info
    csf_info.label_path = '/extra/dti_stats/csf_label/mni_icbm152_csf_tal_nlin_sym_09a_trunc.nii.gz';
    csf_info.template_path = '/extra/dti_stats/csf_label/mni_icbm152_t2_tal_nlin_sym_09a.nii.gz';
    csf_info.template_masked_path = '/extra/dti_stats/csf_label/mni_icbm152_t2_tal_nlin_sym_09a_mask.nii.gz';
    
    % Perform dti stats pipeline
    dti_stats_pdf_path = dti_stats(job_dir_path, ...
                                   dwmri_path, ...
                                   bvec_path, ...
                                   bval_path, ...
                                   mask_path, ...
                                   fsl_path, ...
                                   camino_path, ...
                                   csf_info, ...
                                   bet_params, ...
                                   movement_params_path);
    
    % Merge PDFs
    pdf_path = fullfile(job_dir_path,'PDF','dtiQA.pdf');
    system_utils.system_with_errorcheck(['gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -dPDFSETTINGS=/prepress -sOutputFile=' pdf_path ' ' topup_eddy_pdf_path ' ' dti_stats_pdf_path ],'Failed to merge output PDFs');
    
    % Remove single-paged pdfs
    delete(topup_eddy_pdf_path);
    delete(dti_stats_pdf_path);
catch e
    disp(['Matlab script failed. Reason: ' getReport(e)]);
    exit(1);
end
exit(0);
