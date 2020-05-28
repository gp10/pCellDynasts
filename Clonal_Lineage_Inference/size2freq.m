function [Freq, Freq_rel] = size2freq(clonesizes,timepoints,vartype,clonesizes_ref,cutoff1,clonesizes_ref2,cutoff2)
%% CLONE-SIZES-TO-FREQUENCIES CONVERTER:
% It computes the histogram counts or clone size frequencies from the
% individual no. of basal cells per clone

% from Piedrafita et al, 2020

%% Input:
% clonesizes: cell array {:,timepoints} or matrix [:,timepoints] containing clone sizes
% timepoints: row-vector containing time points
% vartype: type of format of 'clonesizes' data: 1=cell array | 2=matrix
% clonesizes_ref: cell array or matrix of clone sizes {:,timepoints} used for cutoff1 (usually the TOTAL clone sizes)
% cutoff1: minimum clone size considered for the histogram (by default=0: all sizes considered)
% clonesizes_ref2: cell array or matrix of clone sizes {:,timepoints} used for cutoff2 (usually the BASAL clone sizes) (optional cutoff)
% cutoff2: minimum clone size considered for the histogram (by default=0: all sizes considered)

%% Output:
% Freq: cell array {:,timepoints} or matrix [:,timepoints] containing clone size frequencies. Each row r contains absolute number of clones with r-1 basal cells (subject to excluding rules -cutoffs- above)
% Freq_rel: cell array {:,timepoints} or matrix [:,timepoints] containing relative clone size frequencies. Each row r contains relative number of clones with r-1 basal cells (subject to excluding rules -cutoffs- above)

%%
if (nargin < 4)
    clonesizes_ref = clonesizes;
    cutoff1 = 0;
end

switch vartype
    
    case 1 % clonesizes given in a cell_array format
        Freq = cell(1,size(timepoints,2));
        Freq_rel = cell(1,size(timepoints,2));
        for aa = 1:size(timepoints,2)
            % Apply cutoffs:
            loc_prolif = find(clonesizes_ref{:,aa}>=cutoff1);
            if (nargin > 5)
                loc_prolif = loc_prolif(find(clonesizes_ref2{:,aa}(loc_prolif,1)>=cutoff2));
            end
            % Retrieve histograms of bin counts:
            Freq{:,aa} = histc(clonesizes{:,aa}(loc_prolif,1),[0:1:max(clonesizes{:,aa})]);
            Freq_rel{:,aa} = Freq{:,aa}./sum(Freq{:,aa});
        end
        
    case 2 % clonesizes given in a matrix format
        Freq = zeros(size([0:max(max(clonesizes))],2),size(timepoints,2));
        Freq_rel = zeros(size([0:max(max(clonesizes))],2),size(timepoints,2));
        for at = 1:size(timepoints,2)
            % Apply cutoffs:
            loc_prolif = find(clonesizes_ref(:,at)>=cutoff1);
            if (nargin > 5)
                loc_prolif = loc_prolif(find(clonesizes_ref2(loc_prolif,at)>=cutoff2));
            end
            % Retrieve histograms of bin counts:
            Freq(:,at) = histc(clonesizes(loc_prolif,at),[0:1:max(max(clonesizes))]);
            Freq_rel(:,at) = Freq(:,at)./sum(Freq(:,at));
        end

end
