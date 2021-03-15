function [ d ] = hcompare_KL( h1,h2 )
%This routine evaluates the Kullback-Leibler (KL) distance between histograms. 
%             Input:      h1, h2 - histograms
%             Output:    d – the distance between the histograms.
%             Method:    KL is defined as: 
%             Note, KL is not symmetric, so compute both sides.
%             Take care not to divide by zero or log zero: disregard entries of the sum      for which with H2(i) == 0.

% temp = sum(h1 .* log(h1 ./ h2));
% temp( isinf(temp) ) = 0; % this resloves where h1(i) == 0 
% d1 = sum(temp);
% 
% temp = sum(h2 .* log(h2 ./ h1)); % other direction of compare since it's not symetric
% temp( isinf(temp) ) = 0;
% d2 = sum(temp);
% 
% d = d1 + d2;


%%
% %# you may want to do some input testing, such as whether h1 and h2 are
% %# of the same size
% 
% %# preassign the output
% % d = zeros(size(h1));
% 
% %# create an index of the "good" data points
% goodIdx = h1>0 & h2>0; %# bin counts <0 are not good, either
% 
% d1 = sum(h1(goodIdx) .* log(h1(goodIdx) ./h2(goodIdx)));
% d2 = sum(h2(goodIdx) .* log(h2(goodIdx) ./h1(goodIdx)));
% 
% %# overwrite d only where we have actual data
% %# the rest remains zero
% % d(goodIdx) = d1 + d2;
% d = d1 + d2

%%

d1 = sum((h1+eps) .* log((h1+eps) ./(h2+eps)));
d2 = sum((h2+eps) .* log((h2+eps) ./(h1+eps)));

d = d1 + d2;

%%

% d=sum(h1.*log2(h1+eps)-h1.*log2(h2+eps)); %KL(h1,h2)
% % d=sum(h2.*log2(h2+eps)-h2.*log2(h1+eps)) %KL(h2,h1)


%% ASIF
% P_vec = h1(:);
% Q_vec = h2(:);
% Dkl = P_vec.*(log2(P_vec+eps)-log2(Q_vec+eps));
% Dkl_mean = mean(Dkl);



end