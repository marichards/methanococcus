function drawTheirVitalMap(theirs,sol)

%Take in the M.maripaludis model and its solution, print out the Paint4Net
%biographs for both the Biosynthetic and Catabolic pathways

%Catabolic
%Make a list of catabolic reactions
catabolic = {'R76'...
            ,'R41'...
            ,'R489'...
            ,'R57'...
            ,'R304'...
            ,'R305'...
            ,'R306'...
            ,'R45'...
            ,'R55'...
            ,'R46'...
            ,'R53'...
            ,'R56'...
            ,'R535'...
            ,'R442'...
            ,'R44'...
            ,'R42'...
            ,'R67'...
            ,'R65'...
            ,'R54'...
            }';

%Draw the catabolic biograph
draw_by_rxn(theirs,catabolic,'true','struc',{''},{''},sol.x);
%Biosynthetic
biosynthetic = {'R54'...
                ,'R41'...
                ,'R75'...
                ,'R122'...
                ,'R124'...
                ,'R125'...
                ,'R24'...
                ,'R495'...
                ,'R23'...
                ,'R473'...
                ,'R493'...
                ,'R494'...
                ,'R71'...
                ,'R26'...
                ,'R25'...
                ,'R15'...
                ,'R16'...
                ,'R27'...
                ,'R4'...
                ,'R8'...
                ,'R20'...
                ,'R466'...
                ,'R29'...
                ,'R17'...
                }';

%Draw the biosynthetic biograph
draw_by_rxn(theirs,biosynthetic,'true','struc',{''},{''},sol.x);