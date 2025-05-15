%% p values from figure 3K,L
yaxis_height = 26.777;

% K 0-25ms PINK
K25p = [...
11.851 / yaxis_height;
6.364  / yaxis_height;
4.263  / yaxis_height];
% K 0-25ms BLUE
K25b = [...
5.05   / yaxis_height;
3.232  / yaxis_height;
0.013  / yaxis_height];
% K 25-50ms PINK
K50p = [...
21.24  / yaxis_height;
11.96  / yaxis_height;
11.588 / yaxis_height];
% K 25-50ms BLUE
K50b = [...
4.356  / yaxis_height;
4.071  / yaxis_height;
2.864  / yaxis_height];
% K 50-100ms PINK
K100p = [...
18.395 / yaxis_height;
11.856 / yaxis_height;
10.327 / yaxis_height];
% K 50-100ms BLUE
K100b = [...
7.466 / yaxis_height;
5.165 / yaxis_height;
2.236 / yaxis_height];


% L 0-25ms PINK
L25p = [...
11.029 / yaxis_height;
8.138  / yaxis_height;
4.696  / yaxis_height];
% L 0-25ms BLUE
L25b = [...
7.392  / yaxis_height;
3.953  / yaxis_height;
0.921  / yaxis_height];
% L 0-50ms PINK
L50p = [...
22.482 / yaxis_height;
13.128 / yaxis_height;
11.35  / yaxis_height];
% L 0-50ms BLUE
L50b = [...
8.795  / yaxis_height;
3.657  / yaxis_height;
2.976  / yaxis_height];
% L 0-100ms PINK
L100p = [...
18.724 / yaxis_height;
11.592 / yaxis_height;
10.994 / yaxis_height];
% L 0-100ms BLUE
L100b = [...
7.923  / yaxis_height;
5.415  / yaxis_height;
4.700  / yaxis_height];

%% t tests
[h_K25,p_K25] = ttest2(K25p,K25b);
[h_K50,p_K50] = ttest2(K50p,K50b);
[h_K100,p_K100] = ttest2(K100p,K100b);

[h_L25,p_L25] = ttest2(L25p,L25b);
[h_L50,p_L50] = ttest2(L50p,L50b);
[h_L100,p_L100] = ttest2(L100p,L100b);

RejectNull = [h_K25 h_K50 h_K100 h_L25 h_L50 h_L100]';
Pvalue = [p_K25 p_K50 p_K100 p_L25 p_L50 p_L100]';
RN = {'K 0-25ms','K 0-50ms','K 0-100ms','L 0-25ms','L 0-50ms','L 0-100ms'}';
% make a table
T = table(RejectNull,Pvalue,RowNames=RN)

%% now do a one tailed ttest assuming RFA inact has larger effect
% (p)ink is RFA, (b)lue is CFA
[h_K25,p_K25] = ttest2(K25p,K25b,Tail='right');
[h_K50,p_K50] = ttest2(K50p,K50b,Tail='right');
[h_K100,p_K100] = ttest2(K100p,K100b,Tail='right');

[h_L25,p_L25] = ttest2(L25p,L25b,Tail='right');
[h_L50,p_L50] = ttest2(L50p,L50b,Tail='right');
[h_L100,p_L100] = ttest2(L100p,L100b,Tail='right');

RejectNull = [h_K25 h_K50 h_K100 h_L25 h_L50 h_L100]';
Pvalue = [p_K25 p_K50 p_K100 p_L25 p_L50 p_L100]';
RN = {'K 0-25ms','K 0-50ms','K 0-100ms','L 0-25ms','L 0-50ms','L 0-100ms'}';
% make a table
T = table(RejectNull,Pvalue,RowNames=RN)


%%