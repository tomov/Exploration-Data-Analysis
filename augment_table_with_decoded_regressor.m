function tbl = augment_table_with_decoded_regressor(tbl, regressor, dec, standardize, exclude, V_all)

    % given a table generated with data2table
    % and a decoded regressor dec, append it to the table
    %
    assert(~any(isnan(dec(~exclude))), 'There are NaNs in the decoded regressor');

    switch regressor
        case 'RU'
            decRU = dec;
            if standardize == 1
                decRU(~exclude) = zscore(decRU(~exclude));
            elseif standardize == 2
                decRU(~exclude) = decRU(~exclude) / norm(decRU(~exclude));
            end
            tbl = [tbl table(decRU)];

            % orthogonalized version
            tmp = spm_orth([tbl.RU(~exclude), decRU(~exclude)]);
            decRU_orth = decRU;
            decRU_orth(~exclude) = tmp(:,end);
            if standardize == 1
                decRU_orth(~exclude) = zscore(decRU_orth(~exclude));
            elseif standardize == 2
                decRU_orth(~exclude) = decRU_orth(~exclude) / norm(decRU_orth(~exclude));
            end
            tbl = [tbl table(decRU_orth)];

        case 'TU'
            VdecTU = V_all ./ dec;
            if standardize == 1
                VdecTU(~exclude) = zscore(VdecTU(~exclude));
            elseif standardize == 2
                VdecTU(~exclude) = VdecTU(~exclude) / norm(VdecTU(~exclude));
            end
            tbl = [tbl table(VdecTU)];

            % orthogonalized version
            tmp = spm_orth([tbl.VTU(~exclude), VdecTU(~exclude)]);
            VdecTU_orth = VdecTU;
            VdecTU_orth(~exclude) = tmp(:,end); 
            if standardize == 1
                VdecTU_orth(~exclude) = zscore(VdecTU_orth(~exclude));
            elseif standardize == 2
                VdecTU_orth(~exclude) = VdecTU_orth(~exclude) / norm(VdecTU_orth(~exclude));
            end
            tbl = [tbl table(VdecTU_orth)];

        case 'V'
            decV = dec;
            if standardize == 1
                decV(~exclude) = zscore(decV(~exclude));
            elseif standardize == 2
                decV(~exclude) = decV(~exclude) / norm(decV(~exclude));
            end
            tbl = [tbl table(decV)];

            % orthogonalized version
            tmp = spm_orth([tbl.V(~exclude), decV(~exclude)]);
            decV_orth = decV;
            decV_orth(~exclude) = tmp(:,end);
            if standardize == 1
                decV_orth(~exclude) = zscore(decV_orth(~exclude));
            elseif standardize == 2
                decV_orth(~exclude) = decV_orth(~exclude) / norm(decV_orth(~exclude));
            end
            tbl = [tbl table(decV_orth)];


        case 'DV'
            decDV = dec;
            if standardize == 1
                decDV(~exclude) = zscore(decDV(~exclude));
            elseif standardize == 2
                decDV(~exclude) = decDV(~exclude) / norm(decDV(~exclude));
            end
            tbl = [tbl table(decDV)];

            % orthogonalized version
            tmp = spm_orth([tbl.V(~exclude), tbl.RU(~exclude), tbl.VTU(~exclude), decDV(~exclude)]);
            decDV_orth = decDV;
            decDV_orth(~exclude) = tmp(:,end);
            if standardize == 1
                decDV_orth(~exclude) = zscore(decDV_orth(~exclude));
            elseif standardize == 2
                decDV_orth(~exclude) = decDV_orth(~exclude) / norm(decDV_orth(~exclude));
            end
            tbl = [tbl table(decDV_orth)];



        otherwise
            assert(false);
    end
