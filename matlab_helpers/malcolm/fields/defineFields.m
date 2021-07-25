function [inField, outField] = defineFields(fr, fr_dist, cent, idx_spike, posx,params)
% fixed bug in requiring fields to be at least 3 bins 
% 5/31/17 MGC

    nbins = length(fr);
    
    % calculate firing rate percentiles
    % UniversalParams.xlsx has thresholds for defining fields
    % (percentiles of a shuffled distribution)
    if ~exist('params','var')
        params = readtable('UniversalParams.xlsx');
    end
    frp = nan(3,nbins);
    frp(1,:) = prctile(fr_dist,params.InFieldThreshold);
    frp(2,:) = prctile(fr_dist,params.FieldExtensionThreshold);
    frp(3,:) = prctile(fr_dist,params.OutFieldThreshold);
    
    inField = fr > frp(1,:);
    outField = fr <= frp(3,:); % <= in case this percentile is 0
    
    % outField: keep only groups of 2 or more adjacent bins
    if outField(1)
        outField(1) = outField(2);
    end
    if outField(end)
        outField(end) = outField(end-1);
    end
    for k = 2:nbins-1
        if outField(k)
            if outField(k+1)
            elseif outField(k-1)
            else
                outField(k) = 0;
            end
        end
    end
    
    % inField: keep only groups of 3 or more adjacent bins
    if inField(1)
        inField(1) = inField(2) && inField(3);
    end
    if inField(2)
        inField(2) = (inField(1) && inField(3)) || (inField(3) && inField(4));
    end
    if inField(end)
        inField(end) = inField(end-1) & inField(end-2);
    end
    if inField(end-1)
        inField(end-1) = (inField(end) && inField(end-2)) || (inField(end-2) && inField(end-3));
    end
    for j = 3:nbins-2
        if inField(j)
            if inField(j+1)
                if inField(j-1)
                elseif inField(j+2)
                else
                    inField(j) = 0;
                end
            elseif inField(j-1)
                if inField(j-2)
                else
                    inField(j) = 0;
                end
            elseif ~inField(j-1) && ~inField(j+1)
                inField(j) = 0;
            end
        end
    end
    
    % extend fields by 1 bin to left and/or right if those bins > 70th
    % percentile
    newFields = zeros(1,nbins);
    for m = 2:nbins-1
        if inField(m)
            if ~inField(m-1)
                if fr(m-1) > frp(2,m-1)
                    newFields(m-1) = 1;
                end
            elseif ~inField(m+1)
                if fr(m+1) > frp(2,m+1)
                    newFields(m+1) = 1;
                end
            end
        end
    end
    inField = inField | newFields;
    
    % must have spikes on > 20% of traversals (eliminates fields with just
    % one burst of spikes in periods lasting 5 trials)
    
    % find indices of traversal bounderies
    idx_traversal_boundaries = [1; find(diff(posx)<-100); length(posx)+1];
    numTraversals = length(idx_traversal_boundaries)-1;
    numFields = inField(1) + sum(diff(inField)==1);
    inFieldBins = find(inField);
    traversalsWithSpikes = zeros(numFields,1);
    
    % calculate number of traversals through each field
    % depends on whether field was entered in final traversal or not
    % (if desired, change to require full field traversal)
    numTraversalsPerField = repmat(numTraversals,numFields,1);
    fieldStartingPoints = find(diff(inField)==1)+1;
    if inField(1)
        fieldStartingPoints = [1,fieldStartingPoints];
    end
    fieldEndingPoints = find(diff(inField)==-1);
    if inField(end)
        fieldEndingPoints = [fieldEndingPoints,nbins];
    end
    numTraversalsPerField = numTraversalsPerField - [cent(fieldStartingPoints)>posx(end)]';
    
    % iterate through traversals, counting those where spikes occur
    for n = 1:numTraversals
        idx_spike_traversal = idx_spike(idx_spike>=idx_traversal_boundaries(n) & idx_spike<idx_traversal_boundaries(n+1));
        spike_cts_traversal = hist(posx(idx_spike_traversal),cent)>0;
        counter = 1;
        spikesInThisField = 0;
        for p = 1:sum(inField)
            if spike_cts_traversal(inFieldBins(p))
                spikesInThisField = 1;
            end
            if p<sum(inField)
                if inFieldBins(p+1) - inFieldBins(p) > 1
                    traversalsWithSpikes(counter) = traversalsWithSpikes(counter) + spikesInThisField;
                    spikesInThisField = 0;
                    counter = counter+1;
                end
            elseif p == sum(inField)
                traversalsWithSpikes(counter) = traversalsWithSpikes(counter) + spikesInThisField;
            end
        end
    end
    % remove fields with spikes on fewer than 20% of runs
    keepField = traversalsWithSpikes./numTraversalsPerField > 0.2;
    counter = 1;
    numInField = sum(inField);
    for i = 1:numInField
        inField(inFieldBins(i)) = keepField(counter);
        if i<numInField
            if inFieldBins(i+1) - inFieldBins(i) > 1
                counter = counter+1;
            end
        end
    end
    
    % combine fields if there is only one out of field bin between them
    for i = 2:nbins-1
        if inField(i-1)==1 && inField(i+1)==1
            inField(i)=1;
        end
    end

end