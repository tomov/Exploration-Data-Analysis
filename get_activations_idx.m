function act_idx = get_activations_idx(EXPT, onset, session, nTRs)

    % find peak of HRF
    hrf = spm_hrf(0.001);
    [~,hrf_offset] = max(hrf);
    hrf_offset = hrf_offset / 1000;

    TR = EXPT.TR;
    trs = TR/2 : TR : nTRs * TR;

    act_idx = [];
    for i = 1:length(onset)
        [~, idx] = min(abs(trs - (onset(i) + hrf_offset)));
        act_idx = [act_idx; idx + nTRs * (session(i) - 1)];
    end
end
