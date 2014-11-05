function EIs = importEIs()

%Read in the info in the excel file
[~,genes] = xlsread('genefxn_EI.xlsx','McN work','A7:A1785');
sums = xlsread('genefxn_EI.xlsx','McN work','K7:K1785');

%The sums is the number of trials where it's essential. Separate into
%buckets
EIs.four = lower(genes(sums>=4));
EIs.three = lower(genes(sums>=3));
EIs.two = lower(genes(sums>=2));
EIs.one = lower(genes(sums>=1));
EIs.none = lower(genes(sums==0));
EIs.all = lower(genes);
