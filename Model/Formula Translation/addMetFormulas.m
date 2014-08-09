function model = addMetFormulas(model)
%
%Take in our model, use metabolite names to draw out formulas
%

%First use the two lists. Load them
load ('name_to_formula_files.mat')

%Create the cell array for the model
model.metFormulas = cell(length(model.mets),1);

%Loop through all the metabolites
for i=1:length(model.mets)
   %If it's in the name->ID translation file
   if ismember(model.mets{i},seed_name_to_id(:,1))
        %Find its index so we can find the ID
        [~,idx]=intersect(seed_name_to_id(:,1),model.mets{i});
        %Assign the ID and remove any "_e0" or "_c0"
        ID = regexprep(seed_name_to_id{idx,2},'_(c|e)0','');
        %If it's in the dictionary
        if ismember(ID,seed_id_to_formula(:,1))
            %Find the index in the other dictionary
            [~,idx]=intersect(seed_id_to_formula(:,1),ID);
            %Finally, assign to the metFormulas field for the entry
            model.metFormulas{i} = seed_id_to_formula{idx,2};
        else
            model.metFormulas{i}='';
        end
   else
       model.metFormulas{i}='';
   end
end
   
%Remainder must be done manually...there are 251 on 08/6/2014

%3 are 'None': Oxidizedferredoxin, Reducedferredoxin, apo-ACP
[~,idx]=intersect(model.mets,'Oxidizedferredoxin_c0');
model.metFormulas{idx}='Fe2R4S6';
[~,idx]=intersect(model.mets,'Reducedferredoxin_c0');
model.metFormulas{idx}='Fe2R4S6';
[~,idx]=intersect(model.mets,'apo-ACP_c0');
model.metFormulas{idx}='HOR';

%248 are missing...boo!
%Put in the ferredoxins I inserted all at once, they're the same
[~,idx]=intersect(model.mets,{'Fdox*1_c0','Fdred*1_c0','Fdox*2_c0','Fdred*2_c0'});
model.metFormulas(idx,:)={'Fe2R4S6'};

%List of things...find them
idx = find(ismember(model.metFormulas,''));
%List off the formulas
model.metFormulas(idx,:)={'C8H15N2O3S','C7H14N2O4','C51H72O2','C51H74O2',...
    'C6H9N2O5','C36H68NO10P','C42H75N3O15P2','C35H68O8P','C44H79N3O15P2',...
    'C48H87N3O15P2','C39H76O8P','C5H10N2O3S','C38H72NO10P','C44H79N3O15P2',...
    'C8H13N2O5','C115H224O135P24','C307H536N24O255P24','C15H29O2','C36H61N7O17P3S',...
    'C95H154N9O27P2','','C313H548N24O255P24','C121H236O135P24','',...%Up to DNA Replication
    'C6H11N3O4','C46H83N3O15P2','C37H72O8P','C49H76O4','C40H67N3O15P2',...
    'C31H56O8P','C34H60NO10P','C8H16N2O3','C35H68O8P','C48H83N3O15P2',...
    'C39H72O8P','C42H80NO10P','C5H10N2O3','C119H232O135P24','C311H544N24O255P24',...
    'C123H240O135P24','C315H552N24O255P24','C7H11N2O5','C16H31O2','C37H63N7O17P3S',...
    'C36H68NO10P','C42H75N3O15P2','C7H13N3O4','C40H71N3O15P2','C31H60O8P',...
    '','C38H65N7O17P3S','C17H33O2','C34H64NO10P','C40H71N3O15P2',...
    'C40H76NO10P','C11H14N2O3','C31H60O8P','C9H14N4O3','C14H27O2',...
    'C35H59N7O17P3S','C7H11N2O5','C8H15N3O4','C46H83N3O15P2','C40H76NO10P',...
    'C15H29O2','C36H61N7O17P3S','C33H64O8P','C119H232O135P24','C311H544N24O255P24',...
    'CH4','C115H224O135P24','C307H536N24O255P24','C38H65N7O17P3S','C17H33O2',...
    'C36H63N3O15P2','C30H56NO10P','C11H14N2O4','C313H548N24O255P24','C121H236O135P24',...
    'C94H153N8O26P2','C94H155N9O25P2','C24H42O21','C7H12N2O3','C9H18N2O3',...
    'C27H52O8P','C34H64NO10P','C18H33O2','C39H65N7O17P3S','C7H14N2O3S',...
    'C117H228O135P24','C309H540N24O255P24','C50H70O2','C50H72O2','C42H76NO10P',...
    'C38H68NO10P','C44H75N3O15P2','C33H64O8P','C38H72NO10P','C35H64O8P',...
    'C37H61N7O17P3S','C16H29O2','C309H540N24O255P24','C117H228O135P24','C37H72O8P',...
    'C80H125N16O42R','C77H148O17P2','C81H156O17P2','C39H78NO8P','C40H78O10P',...
    '','C40H63N8O21R','C39H78NO8P','C40H78O10P','C176H307N2O100P4',...
    'C42H82O10P','C77H148O17P2','C41H82NO8P','C20H35N2O8PRS','C20H37N2O9PRS',...
    'C131H234N2O60P3','C124H222N2O54P3','C22H39N2O8PRS','C22H41N2O9PRS','C18H33N2O8PRS',...
    'C18H31N2O8PRS','C11H8O5','C11H10O6','C20H35N2O9PRS','C18H31N2O9PRS',...
    'C16H29N2O8PRS','C36H40N4O8','C8H13O11P','C5H10O8P','C26H49N2O8PRS',...
    'C26H47N2O8PRS','C40H78O13P2','C12H23O2','C12H23O2','C15H25N2O9PRS',...
    'C29H55N2O8PRS','C17H22N4O9P','C2H6NO','C6H12O6','C28H53N2O8PRS',...
    'C22H41N2O8PRS','C22H39N2O8PRS','C24H45N2O9PRS','C24H43N2O8PRS','C24H43N2O9PRS',...
    'C20H40O7P','C11H15N2O8P','C124H222N2O51P2','C7H5NO4','C28H51N2O9PRS',...
    'C28H51N2O8PRS','C22H41N2O8PRS','C24H43N2O9PRS','C22H39N2O9PRS','C18H33N2O9PRS',...
    'C26H49N2O9PRS','C26H49N2O8PRS','C26H47N2O8PRS','C18H31N2O8PRS','C18H33N2O8PRS',...
    'C26H49N2O9PRS','C55H77CoN15O11','C3H10NO','C10H30N4','C22H39N2O9PRS',...
    'C22H41N2O9PRS','C10H30N4','C28H53N2O9PRS','C46H64O2','C42H82O13P2',...
    'C151H265N2O79P4','C157H275N2O84P4','C18H33N2O9PRS','C28H53N2O8PRS','C9H15N2O5',...
    'C55H69CoN11O15','C45H57CoN6O12','C169H295N2O94P4','C163H285N2O89P4','C40H78O13P2',...
    'C20H37N2O8PRS','C20H35N2O8PRS','C24H45N2O8PRS','C24H43N2O8PRS','C9H11NO6',...
    'C117H210N2O45P2','C24H45N2O9PRS','C20H40O7P','C145H255N2O74P4','C139H245N2O70P4',...
    'C131H234N2O63P4','C18H31N2O9PRS','C16H29N2O8PRS','C32H40N7O20P3S','C24H45N2O8PRS',...
    'C34H38N4O4','C10H13N5O10P2','C34H32N4O4','C26H47N2O9PRS','C18H33O2',...
    'C28H51N2O9PRS','C20H37N2O8PRS','C46H64O2','C21H42O7P','C26H47N2O9PRS',...
    'C20H30N6O12S2','C28H51N2O8PRS','C48H73CoN11O8','C10H8O6','C28H53N2O9PRS',...
    'C20H30N6O12S2','C20H37N2O9PRS','C20H35N2O9PRS','C5H7O3S','C6H10O6PS',...
    'C32H39N7O19P3S','C3H7NO2','C8H13N2O5S','C2H4O5P','NH3',...
    'C25H42O7P2','CH33N','CH33NO','C14H13O9','C10H12N5O10PS',...
    'C2H3O4S','C5H8NO5S2','C5H10NO5S2','C6H5O7','C6H3O6',...
    'C6H5O7','C5H4O5','C3H2O6S','C3H4O6S','C5H6O5',...
    'C5H4O4','C5H6O5','C2H3O2','C11H8NO3','C31H40N8O17P3S'};

end