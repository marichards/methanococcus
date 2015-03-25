function idx = findIndex(array,phrase)

%Given an array that contains a phrase, return the index of that phrase
%Phrase can be a list
[~,idx] = intersect(array,phrase);