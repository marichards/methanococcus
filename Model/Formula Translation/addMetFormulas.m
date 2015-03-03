function model = addMetFormulas(model)
%
%Take in our model, use metabolite names to draw out formulas
%
%Do it all manually

model.metFormulas={...
'C6H5O4';...
'H';...
'CO2';...
'C7H4O6';...
'C10H17N4O6';...
'C4H2O4';...
'C6H15N4O2';...
'C6H14NO8P';...
'C8H15NO9P';...
'C21H33N7O16P3S';...
'C23H35N7O17P3S';...
'C6H7NO2R2S2';...
'C6H11NO4';...
'C4H9NO2S';...
'HO3S';...
'C2H3O2';...
'S2O3';...
'C6H9NO2R2S2';...
'C5H11NO2S';...
'H2O';...
'C3H7NO2';...
'C8H15N2O3S';...
'C4H9NO3';...
'C7H14N2O4';...
'C44H47N4O16';...
'C21H26N7O17P3';...
'C44H49N4O16';...
'C21H27N7O17P3';...
'Mg';...
'Mg';...
'C4H4O4';...
'C51H72O2';...
'C51H74O2';...
'C7H11O10P';...
'C7H9O6';...
'HO4P';...
'C9H15N5O14P3';...
'C10H15N5O15P3';...
'CHO2';...
'C29H33N9O12';...
'C30H35N9O12';...
'C4H4O5';...
'C4H2O5';...
'C21H27N7O14P2';...
'C21H26N7O14P2';...
'C4H7NO3';...
'C4H9NO3';...
'C13H18N4O6';...
'C9H16N4O6';...
'C17H20N4O6';...
'C9H13N3O9P';...
'C10H13N5O13P3';...
'C10H13N5O10P2';...
'C8H14N3O7P';...
'CHO3';...
'C6H9N2O5';...
'C4H6NO4';...
'C2H5NO2';...
'C6H15N2O2';...
'C6H15N2O2';...
'C43H44CoN4O16';...
'C45H47CoN4O16';...
'C2H4O';...
'C5H7O3';...
'C7H10O5';...
'H2S';...
'C36H68NO10P';...
'C3H7NO3';...
'C9H13N3O8P';...
'C42H75N3O15P2';...
'C42H36FeN4O16';...
'C42H38N4O16';...
'Fe';...
'C3H3O6P';...
'C3H5O7P';...
'C35H68O8P';...
'C44H79N3O15P2';...
'H2O7P2';...
'C9H13N3O14P3';...
'C48H87N3O15P2';...
'C39H76O8P';...
'C3H3O7P';...
'C3H5O7P';...
'C10H13N5O14P3';...
'HO4P';...
'C7H3NO4';...
'C4H3NO4';...
'C3H6O6P';...
'C3H4O9PS';...
'C8H16N3O8P';...
'C6H12O9P';...
'C6H12O9P';...
'C4H14N2';...
'C11H15N5O3S';...
'C7H22N3';...
'C14H24N6O3S';...
'NH4';...
'C5H8NO4';...
'C5H10N2O3';...
'C7H14O10P';...
'C3H6O6P';...
'C5H10O8P';...
'C5H10O8P';...
'C10H13N5O12P3';...
'C10H12N4O5';...
'C10H12N4O8P';...
'C3H7NO2S';...
'C5H10N2O3S';...
'C9H12N2O14P3';...
'C9H12N2O11P2';...
'C38H72NO10P';...
'C44H79N3O15P2';...
'C10H14N2O14P3';...
'C10H14N2O11P2';...
'C21H33N7O13P2S';...
'C5H6N2O5';...
'CH3NO5P';...
'C6H9O3';...
'C5H4O5';...
'C6H13NO2';...
'C15H21N5O14P2';...
'C15H21N5O20P4';...
'C11H22N2O7PS';...
'C12H21N2O9PS';...
'C15H23N5O15P2';...
'C15H23N5O15P2';...
'C6H12O12P2';...
'C7H10NO4';...
'C7H14N2O3';...
'C8H13N2O5';...
'H';...
'K';...
'K';...
'C10H13N5O7P';...
'HO3PSe';...
'H2Se';...
'C9H12N2O9P';...
'C10H11N2O11P';...
'C14H20N6O5S';...
'C42H40N4O16';...
'C15H23N6O5S';...
'C43H44N4O16';...
'C45H50CoN4O16';...
'C45H53CoN4O14';...
'C4H4N2O2';...
'C4H5N3O';...
'C115H224O135P24';...
'C9H12N2O12P2';...
'C307H536N24O255P24';...
'C17H25N3O17P2';...
'C15H29O2';...
'C36H61N7O17P3S';...
'C9H10N2';...
'C6H4NO2';...
'C11H13NO9P';...
'C14H18N2O7P';...
'C95H154N9O27P2';...
'C95H152N8O28P2';...
'Fe2R4S6';...
'Fe2R4S6';...
'C43H42CoN4O16';...
'C42H40CoN4O16';...
'';...
'C94H155N9O25P2';...
'C94H153N8O26P2';...
'C11H13NO6P';...
'C8H7N';...
'C9H16N4O9P';...
'C9H14N4O9P';...
'C313H548N24O255P24';...
'C121H236O135P24';...
'C5H10O7P2';...
'C20H34O7P2';...
'C15H26O7P2';...
'C5H9O4';...
'C4H9NO6P';...
'C4H8O2';...
'C5H7O4';...
'C9H12N2O15P3';...
'C6H14NO8P';...
'C6H9N3O4P';...
'C6H9N3O';...
'C9H7O4';...
'C9H11NO3';...
'NO3';...
'NO3';...
'C45H53CoN4O14';...
'';...
'Co';...
'C7H14N2O4';...
'C7H14N2O4';...
'C4H8N2O3';...
'C6H11N3O4';...
'C9H13N5O13P3';...
'C10H13N5O9P2';...
'C31H44N6O16P';...
'C3H7O3S2';...
'C30H42N6O16P';...
'C2H5O3S2';...
'C3H7NO2';...
'C5H5N5';...
'C6H12O7PS';...
'C46H83N3O15P2';...
'C37H72O8P';...
'C9H16NO8P';...
'C6H9O3';...
'C6H11O4';...
'C6H8N2O5P';...
'C6H12N3O4P';...
'C27H31N9O15P2';...
'C49H76O4';...
'C49H74O4';...
'C27H33N9O15P2';...
'C44H46CoN4O16';...
'C12H14NO9P';...
'C55H90O7P2';...
'C55H90O4P';...
'C24H37N7O17P3S';...
'C3H5O2';...
'C4H8O7P';...
'C18H22N4O23P4';...
'C6H13N3O3';...
'C5H10O8P';...
'C40H67N3O15P2';...
'C31H56O8P';...
'C9H13N3O11P2';...
'C7H7O5';...
'C20H34O7P2';...
'C7H8O5';...
'C7H10O5';...
'C4H9NO4';...
'C4H9NO7P';...
'C9H12N2O8P';...
'C5H9NO3';...
'C10H13N2O4';...
'C45H53N4O14';...
'C11H12N2O2';...
'C9H16NO5';...
'C34H60NO10P';...
'C6H12O9P';...
'C6H9O4';...
'C9H13N3O10P2';...
'C9H13N3O13P3';...
'C8H16N2O3';...
'C10H14N4O9P';...
'C17H25N3O17P2';...
'C5H7O4';...
'C9H13N3O7P';...
'C44H45N4O17';...
'C45H47N4O17';...
'C10H8O6';...
'C10H12NO5';...
'C35H68O8P';...
'C48H83N3O15P2';...
'C39H72O8P';...
'C8H9NO6P';...
'Na';...
'Na';...
'C5H8NO2';...
'C5H8NO2';...
'C42H80NO10P';...
'C3H6O10P2';...
'C6H11O7P';...
'C6H11O4';...
'C9H7O3';...
'C13H17N5O8P';...
'C12H17N4O4PS';...
'C6H9N3O7P2';...
'C6H9NO4PS';...
'C28H42N7O19P3S';...
'C7H10O4';...
'C5H10N2O3';...
'C119H232O135P24';...
'C311H544N24O255P24';...
'C5H10O14P3';...
'C123H240O135P24';...
'C315H552N24O255P24';...
'C7H11N2O5';...
'C16H31O2';...
'C37H63N7O17P3S';...
'C9H15N5O14P3';...
'C3H6O3';...
'C3H8O3';...
'C21H24N6O15P2';...
'C9H14N4O8P';...
'C13H16N4O12P';...
'C4H5O3';...
'C25H37N7O18P3S';...
'C3H3O3';...
'C7H6NO2';...
'C10H8O6';...
'C5H10N2O3S';...
'C8H15NOS2';...
'C8H17NOS2';...
'C36H68NO10P';...
'C42H75N3O15P2';...
'C9H12N2O9P';...
'C7H13N3O4';...
'H2O2';...
'C10H10O10P';...
'C14H21N4O8P2S';...
'C6H9O4';...
'C12H17N4O7P2S';...
'C4H5O3';...
'C5H5N2O4';...
'C5H3N2O4';...
'O2';...
'C40H71N3O15P2';...
'C31H60O8P';...
'C37H50N2O';...
'C11H20NO7PS';...
'C13H23NO10PS3';...
'C37H52N2O';...
'C35H39N4O16';...
'C34H40N4O15';...
'C40H36N4O16';...
'C9H15N5O8P';...
'C2H7NO3S';...
'C2H7NO3S';...
'C6H12O7PS';...
'';...
'C4H4O6';...
'C5H9NO3';...
'N2';...
'H2';...
'C29H34N5O18P';...
'C29H32N5O18P';...
'C7H5O3';...
'C8H14N2O9P';...
'C8H15NO6';...
'C38H65N7O17P3S';...
'C17H33O2';...
'C5H10N2O3';...
'C10H14N2O8P';...
'C7H9O8P';...
'C6H13NO2';...
'C34H64NO10P';...
'C40H71N3O15P2';...
'C5H11NO2';...
'C5H11NO2';...
'C15H21N5O14P2';...
'C10H12N4O14P3';...
'C10H12N4O11P2';...
'C5H6N2O2';...
'C5H7N3O';...
'C40H76NO10P';...
'C9H11NO2';...
'C11H14N2O3';...
'C5H4N4O';...
'C5H4N4O';...
'C7H8O4';...
'C5H10O7P2';...
'C10H13N5O3';...
'C10H17N2O3';...
'C10H15N2O3S';...
'S';...
'C14H15N5O11P';...
'C19H19N7O6';...
'C20H21N7O6';...
'C10H13N5O13P3';...
'C6H10N3O';...
'C6H9N3O2';...
'C18H35O2';...
'C39H67N7O17P3S';...
'C31H60O8P';...
'Fe';...
'C6H12N3O';...
'C4H9NO6P';...
'C3H10NO4P';...
'C9H15N5O3';...
'C9H13N5O3';...
'C45H57CoN6O12';...
'C9H14N4O3';...
'C14H27O2';...
'C35H59N7O17P3S';...
'C44H48CoN4O16';...
'Cl';...
'Cl';...
'CHO2';...
'C7H11N2O5';...
'C44H44CoN4O16';...
'C5H12NO7P';...
'C40H38N4O17';...
'C6H10N2O6P';...
'C10H12N4O9P';...
'C8H15N3O4';...
'C4H8O6P';...
'C31H42N6O16P';...
'C9H17NO3';...
'C11H15N2O8P';...
'C46H83N3O15P2';...
'C40H76NO10P';...
'C15H29O2';...
'C36H61N7O17P3S';...
'C33H64O8P';...
'C25H36N7O19P3S';...
'C119H232O135P24';...
'C311H544N24O255P24';...
'C7H9NO5';...
'C15H22N2O17P2';...
'C15H19N2O18P2';...
'C7H15N2O8P';...
'C7H9O5';...
'C2H6O';...
'C10H12N5O3R';...
'C15H19N6O6R';...
'CH4';...
'C115H224O135P24';...
'C307H536N24O255P24';...
'C12H14NO9P';...
'C38H65N7O17P3S';...
'C17H33O2';...
'C35H59N7O17P3S';...
'C14H27O2';...
'C4H5N3O';...
'C7H10NO8P';...
'C36H63N3O15P2';...
'C30H56NO10P';...
'C8H12NO6';...
'C43H44N4O17';...
'C10H10NO5';...
'C11H14N2O4';...
'C3H6O2';...
'C6H12O8P';...
'C5H11NO2';...
'C37H63N7O17P3S';...
'C16H31O2';...
'C313H548N24O255P24';...
'C121H236O135P24';...
'C10H13N5O11P2';...
'C17H25N3O17P2';...
'C10H11N5O6P';...
'C10H13N5O4';...
'C27H40N7O20P3S';...
'O4S';...
'O4S';...
'C4H7NO7P';...
'C10H18O7P2';...
'C94H153N8O26P2';...
'C94H155N9O25P2';...
'C6H12O9P';...
'C24H42O21';...
'C30H52O26';...
'C45H58N6O12';...
'C45H54N4O14';...
'C9H21N2O2';...
'C7H12N2O3';...
'C31H41N6O16P';...
'C31H42N6O17P';...
'C9H18N2O3';...
'C27H52O8P';...
'C15H22N2O17P2';...
'C4H4N2O2';...
'C34H64NO10P';...
'C18H33O2';...
'C39H65N7O17P3S';...
'C7H14N2O3S';...
'C117H228O135P24';...
'C309H540N24O255P24';...
'C50H70O2';...
'C50H72O2';...
'C5H11NO3S';...
'C42H76NO10P';...
'C15H19N5O6S';...
'Co';...
'C38H68NO10P';...
'C44H75N3O15P2';...
'C16H22N2O15P2';...
'C16H26N3O14P2';...
'CHRS2';...
'CH3RS2';...
'C33H64O8P';...
'C38H72NO10P';...
'C5H13N2O2';...
'C35H64O8P';...
'C63H103NO12P2';...
'C10H13N5O8P';...
'C10H13N5O5';...
'C6H12O9P';...
'C6H12O9P';...
'C12H17N4OS';...
'C12H17N4OS';...
'C10H13N5O10P2';...
'CH4N2O';...
'CH4N2O';...
'C37H61N7O17P3S';...
'C16H29O2';...
'C5H4O6';...
'C5H8NO5';...
'C6H9O4';...
'C6H11O4';...
'C6H9NOS';...
'C9H13N3O5';...
'C3H7NO6P';...
'C309H540N24O255P24';...
'C117H228O135P24';...
'C37H72O8P';...
'C5H16N4';...
'Fe';...
'Fe';...
'C5H4N4O4';...
'C4H6N4O3';...
'C4H6N4O3';...
'H2O';...
'CO2';...
'O2';...
'C5H9NO3';...
'C5H6NO2';...
'C4H7NO7P';...
'C3H8NO5P';...
'C7H10NO5';...
'C7H7NO4';...
'CH2NO2';...
'C4H7NO3';...
'C3H8NO';...
'C5H5N4O5';...
'C6H11NO3';...
'C6H8NO2';...
'C6H8NO2';...
'C6H11NO3';...
'Mn';...
'C80H125N16O42R';...
'C77H148O17P2';...
'C81H156O17P2';...
'Cu';...
'C39H78NO8P';...
'C40H78O10P';...
'C10H16N3O6S';...
'HOR';...
'';...
'C40H63N8O21R';...
'C39H78NO8P';...
'Ca';...
'C48H73CoN11O8';...
'C40H78O10P';...
'C20H23N7O6';...
'C19H21N7O6';...
'C176H307N2O100P4';...
'C42H82O10P';...
'C72H101CoN18O17P';...
'C77H148O17P2';...
'C11H21N2O7PRS';...
'C41H82NO8P';...
'C20H21N7O7';...
'C34H30FeN4O4';...
'Zn';...
'C20H35N2O8PRS';...
'C20H37N2O9PRS';...
'C47H69O3';...
'C46H70O';...
'C131H234N2O60P3';...
'C17H25N5O16P2';...
'C124H222N2O54P3';...
'C22H39N2O8PRS';...
'C22H41N2O9PRS';...
'C18H33N2O8PRS';...
'C18H31N2O8PRS';...
'C11H8O5';...
'C11H10O6';...
'HO10P3';...
'C14H22N2O10PRS';...
'C20H35N2O9PRS';...
'Zn';...
'C18H31N2O9PRS';...
'C16H29N2O8PRS';...
'C34H65NO12P';...
'C43H75N3O20P2';...
'C36H40N4O8';...
'C8H13O11P';...
'C5H10O8P';...
'C26H49N2O8PRS';...
'C26H47N2O8PRS';...
'C40H78O13P2';...
'C12H23O2';...
'C12H23O2';...
'C15H25N2O9PRS';...
'C29H55N2O8PRS';...
'C17H22N4O9P';...
'C2H6NO';...
'C17H20N4O9P';...
'C6H12O6';...
'C28H53N2O8PRS';...
'C22H41N2O8PRS';...
'C22H39N2O8PRS';...
'C24H45N2O9PRS';...
'C24H43N2O8PRS';...
'C3H8O6P';...
'C24H43N2O9PRS';...
'C58H85CoN16O11';...
'C20H40O7P';...
'C11H15N2O8P';...
'C16H24N2O16P2';...
'C16H22N2O15P2';...
'C16H24N2O15P2';...
'C124H222N2O51P2';...
'C7H5NO4';...
'C48H72O3';...
'C48H71O4';...
'C19H17N7O6';...
'C19H17N7O6';...
'C28H51N2O9PRS';...
'C28H51N2O8PRS';...
'C22H41N2O8PRS';...
'C24H43N2O9PRS';...
'C22H39N2O9PRS';...
'C47H70O3';...
'C18H33N2O9PRS';...
'Cu';...
'C30H50O7P2';...
'C26H49N2O9PRS';...
'C76H139N2O30P2';...
'C68H128N2O23P2';...
'C17H24N3O15P';...
'C110H198N2O39P2';...
'C96H172N2O38P2';...
'C25H47N2O8PRS';...
'C26H49N2O8PRS';...
'C26H47N2O8PRS';...
'C26H41N7O17P3S';...
'C18H31N2O8PRS';...
'C18H33N2O8PRS';...
'C26H49N2O9PRS';...
'C58H84CoN16O14P';...
'C28H39N5O23P2';...
'C35H51N7O26P2';...
'C55H77CoN15O11';...
'C3H10NO';...
'C10H30N4';...
'C35H58O7P2';...
'C40H66O7P2';...
'C25H47N2O9PRS';...
'C29H50N3O18P2';...
'C22H39N2O9PRS';...
'C22H41N2O9PRS';...
'C10H30N4';...
'C26H41N7O17P3S';...
'C68H128N2O20P';...
'C28H53N2O9PRS';...
'C45H62O2';...
'C46H64O2';...
'C42H82O13P2';...
'C23H33N4O20P2';...
'C5H8NO4';...
'C151H265N2O79P4';...
'C157H275N2O84P4';...
'C18H33N2O9PRS';...
'C28H53N2O8PRS';...
'Mn';...
'C9H15N2O5';...
'C55H69CoN11O15';...
'C45H57CoN6O12';...
'C169H295N2O94P4';...
'C163H285N2O89P4';...
'Ca';...
'C68H96CoN21O21P2';...
'C40H78O13P2';...
'C20H26N3O19P2';...
'C20H37N2O8PRS';...
'C20H35N2O8PRS';...
'C24H45N2O8PRS';...
'C24H43N2O8PRS';...
'C23H43N2O8PRS';...
'C84H150N2O37P2';...
'C47H72O2';...
'C21H33N7O13P2S';...
'C9H11NO6';...
'C117H210N2O45P2';...
'C24H45N2O9PRS';...
'C20H40O7P';...
'C7H14O13P2';...
'C7H14O10P';...
'C145H255N2O74P4';...
'C139H245N2O70P4';...
'C131H234N2O63P4';...
'C7H14O10P';...
'C46H70O2';...
'C18H31N2O9PRS';...
'C16H29N2O8PRS';...
'C32H40N7O20P3S';...
'C24H45N2O8PRS';...
'C34H38N4O4';...
'C10H13N5O10P2';...
'C6H12N2O3';...
'C31H51N3O19P2';...
'C14H18N2O4';...
'C34H32N4O4';...
'C26H47N2O9PRS';...
'C18H33O2';...
'C20H28N3O19P2';...
'C8H13O8';...
'C28H51N2O9PRS';...
'C20H37N2O8PRS';...
'C87H139N7O23P2';...
'C46H64O2';...
'C21H42O7P';...
'C11H7O4';...
'C26H47N2O9PRS';...
'C20H30N6O12S2';...
'C17H25N5O16P2';...
'C28H51N2O8PRS';...
'C48H73CoN11O8';...
'C41H61N9O28P2';...
'C10H8O6';...
'C18H35O2';...
'C28H53N2O9PRS';...
'C20H30N6O12S2';...
'C20H37N2O9PRS';...
'C20H35N2O9PRS';...
'C14H27O2';...
'C5H7O3S';...
'C6H10O6PS';...
'C32H39N7O19P3S';...
'C24H34N7O19P3S';...
'C3H7NO2';...
'C8H13N2O5S';...
'C2H4O5P';...
'NH4';...
'C6H11O10P2';...
'C25H42O7P2';...
'C3H10N';...
'C3H9NO';...
'C14H13O9';...
'C10H12N5O10PS';...
'C2H3O4S';...
'C5H8NO5S2';...
'C5H10NO5S2';...
'C6H5O7';...
'C6H3O6';...
'C6H5O7';...
'C3H2O6S';...
'C3H4O6S';...
'C5H6O5';...
'C5H4O4';...
'C5H6O5';...
'C2H3O2';...
'C11H8NO3';...
'C31H40N8O17P3S';...
};

end