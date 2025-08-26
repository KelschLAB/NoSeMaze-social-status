function DS_data = compute_DS_from_match_matrix(match_matrix)
%% David's score
%%% Pij = proportion of wins. number of time that i defeats j (aij) divided by the total number
%%% of interaction btween i and j (nij). Pji is the proportion of losses

%%% For each mouse one can deduce then a David's score DS=w+w2-l-l2

%%% w is the sum of i's Pij values and w2 is the sum of w values (weighted for the corresponding Pij) of those
%%% animals that with which i interacted (l and l2 are calculated in similar way)

for i=1:size(match_matrix,1)
    for j=1:size(match_matrix,2)
        int(i,j)=match_matrix(i,j)+match_matrix(j,i);
    end
end


P=match_matrix./int;
W=sum(P,2,'omitnan')'; %%% w of each animal
L=sum(P,1,'omitnan');  %%% l for each animal

for i=1:length(W)
    W2(i,:)=P(i,:).*W;
end

for i=1:length(L)
    L2(i,:)=P(:,i)'.*L;
end

w2=sum(W2,2,'omitnan')';
l2=sum(L2,2,'omitnan')';

DS=W+w2-L-l2;
%%%%% David's score plot
[DSv,DSi]=sort(DS,'descend');


%% parse output
DS_data.DS = DS;
DS_data.DS_sorted = DSv;
DS_data.DS_sortedIndex = DSi;
DS_data.match_matrix = match_matrix;


