AAString    					= "ACDEFGHIKLMNPQRSTVWY";

ExecuteAFile 					(HYPHY_BASE_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "Utility" + DIRECTORY_SEPARATOR + "GrabBag.bf");
ExecuteAFile 					(HYPHY_BASE_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "Utility" + DIRECTORY_SEPARATOR + "AncestralMapper.bf");
ExecuteAFile 					(HYPHY_BASE_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "Utility" + DIRECTORY_SEPARATOR + "ReadDelimitedFiles.bf");

SetDialogPrompt 				("Select the protein alignment to process");

DataSet							baseAlignment = ReadDataFile (PROMPT_FOR_FILE);
basePath						= LAST_FILE_PATH;
trunkPath						= splitFilePath (LAST_FILE_PATH);
extensionStack					= {};
extensionStack[Abs(extensionStack)] = trunkPath["EXTENSION"];
trunkPath						= trunkPath["DIRECTORY"] + trunkPath["FILENAME"];

reducedTrunk					= splitFilePath (trunkPath);

ACCEPT_ROOTED_TREES				= 1;

while (Abs(reducedTrunk["EXTENSION"]))
{
	extensionStack[Abs(extensionStack)] = reducedTrunk["EXTENSION"];
	reducedTrunk					= reducedTrunk["DIRECTORY"] + reducedTrunk["FILENAME"];
	reducedTrunk					= splitFilePath (reducedTrunk);
}

trunkPath									  = reducedTrunk["DIRECTORY"] + reducedTrunk["FILENAME"];
DataSetFilter					baseFilter	  = CreateFilter (baseAlignment,1);

fprintf							(stdout, "[SET TRUNK PATH TO ", trunkPath, "]\n");

bySiteFreqs						= {};

COUNT_GAPS_IN_FREQUENCIES = 0;

adder = {1,20}["1"];

for	(s = 0; s < baseFilter.sites; s=s+1)
{
	DataSetFilter 		siteFilter = CreateFilter (baseAlignment,1,siteIndex == s);
	siteF	 = {20,1};
	currentM = siteFilter.species;
	for (z = 0; z < currentM; z = z+1)
	{
		GetDataInfo (charInfo, siteFilter, z, 0);
		k = (adder*charInfo)[0];
		if (k == 1)
		{
			siteF = siteF + charInfo;
		}
	}
	bySiteFreqs [s] = siteF;
}


fprintf							(stdout, "[COMPUTED SITE-BY-SITE BASE COMPOSITIONS]\n");
FEL_matrix			=			(ReadCSVTable (trunkPath + ".fel", 1))[1];
fprintf							(stdout, "[READ FEL RESULTS]\n");
ExecuteAFile					(basePath + ".base");
fprintf							(stdout, "[RELOADED BASE AA FIT]\n");

rootStateBySite					= {};

felSummary						 = {2,1};
for	(s = 0; s < Rows(FEL_matrix); s=s+1)
{
	if (FEL_matrix[s][5] < 0.05)
	{
		s2 = FEL_matrix[s][0] > FEL_matrix[s][1];
		felSummary[s2] = felSummary[s2] + 1;
	}
}
fprintf							(stdout, "[RECONSTRUCTED FEL SUMMARY]\n", felSummary, "\n");
psOutput	= trunkPath + ".ps";

fprintf							(psOutput, CLEAR_FILE, KEEP_OPEN);
SetDialogPrompt					("Save the resulting HTML file to");
fprintf							(PROMPT_FOR_FILE, CLEAR_FILE, KEEP_OPEN,
								"<!DOCTYPE HTML PUBLIC '-//W3C//DTD HTML 4.01 Transitional//EN'>\n<html>\n\t<head>\n\t\t<META http-equiv='Content-Style-Type' content='text/css'>\n\t\t<meta http-equiv='Content-Type' content='text/html; charset=ISO-8859-1'>\n\t\t<META NAME='description' CONTENT='Explore your sequence alignments for evidence of adaptive or purifying evolution, build NJ trees and select appropriate evolutionary models using fast likelihood based methods on our server..'>\t\n\t\t<title>\n\t\t\tAdaptive Evolution Server @ Datamonkey.org\n\t\t</title>\n\t\t<LINK REL=STYLESHEET TYPE='text/css' HREF='http://www.datamonkey.org/2007/styles'>\n\t\t\n</head><body background = 'http://www.hyphy.org/images/aquastripe.gif' bgcolor = '#FFFFFF'>",
								"\n<H1 class = 'SuccessCap'>DEPS report for <span style='font-size: 14px'>",basePath,"</span></H1>\n"
								);
								
htmlOutPath			=			LAST_FILE_PATH;

DEPS_Results		=			basePath + "_bysite.csv";
fscanf							(DEPS_Results, "Lines", DIR_matrix);
pvString			= 			splitOnRegExp (DIR_matrix[1], "\\,");

fprintf							(stdout, "[READ DEPS RESULTS]\n");

test_p_values					 = {20,2}["_MATRIX_ELEMENT_ROW_"];
for	(s = 0; s < 20; s=s+1)
{
	test_p_values[s][0] = 0+pvString[s+1];
}
test_p_values      				 = test_p_values % 0;
rejectedHypotheses   			 = {};
parameterEstimates	 			 = {20,3}["1"];
sitesByResidue					 = {};

ancCachesByResidue  = {};
rootStatesByResidue = {};
rootIndexByResidue	= {};

for (k=0; k<20; k=k+1)
{
	if (test_p_values[k][0] < (0.05/(20-k)))
	{
		residueIndex								 = test_p_values[k][1];
		rejectedHypotheses  [residueIndex]           = 1;
		rejectedHypotheses  [AAString[residueIndex]] = 1;
		parameterEstimates  [residueIndex][0]		 = test_p_values[k][0];
		fprintf (stdout, "[WORKING ON RESIDUE ",AAString[residueIndex],"]\n"); 		
		ExecuteAFile		(basePath + "." + AAString[residueIndex]);
		ancCacheID 						= _buildAncestralCache ("lfb", 0);
		ancCachesByResidue[k] 			= ancCacheID;
		parameterEstimates  [residueIndex][1]		 = rateBiasTo;
		parameterEstimates  [residueIndex][2]		 = P_bias;
		rootStatesByResidue[residueIndex]			 = {};
		rootIndexByResidue[residueIndex]			 = {};
		for	(s = 0; s < baseFilter.sites; s=s+1)
		{
			_rss = _rootState (ancCacheID,s);
			(rootStatesByResidue[residueIndex])[s] = _rss["CHAR"];
			(rootIndexByResidue[residueIndex])[s]  = _rss["INDEX"];
		}		
	}
	else
	{
		break;
	}
}

for (s=2; s<Columns(DIR_matrix); s=s+1)
{
	siteInfo 	= splitOnRegExp(DIR_matrix[s],"\\,");
	siteIndex 	= 0+siteInfo[0];
	for (r = 1; r < 21; r=r+1)
	{
		thisBF = 0 + siteInfo[r];
		if (thisBF >= 100.)
		{
			sitesByResidue[r-1] = sitesByResidue[r-1] + 1;
		}
	}
}

pCount = 0;

baselineBL						= BranchLength (givenTree,-1);
referenceL						= (baselineBL * (Transpose(baselineBL)["1"]))[0];

bySiteAVL					= {};

for (s=2; s<Columns(DIR_matrix); s=s+1)
{
	siteInfo 	= splitOnRegExp(DIR_matrix[s],"\\,");
	siteIndex 	= 0+siteInfo[0];
	myBFs		= {};
	tRI			= 0;
	for (r = 1; r < 21; r=r+1)
	{
		thisBF = 0 + siteInfo[r];
		if (thisBF >= 100.)
		{
			myBFs[r-1] = thisBF;
			tRI = r-1;
		}
	}
	if (Abs(myBFs))
	{
		thisResidue = {};
		
		thisResidue [0] = returnSiteSignature (bySiteFreqs[siteIndex-1]);
		thisResidue [1] = (rootStatesByResidue[tRI]) [siteIndex-1];
		
		thisSite 		= _substitutionsBySite (ancCachesByResidue[tRI],siteIndex-1);
		thisSitePS 		= _mapSubstitutionsBySite (ancCachesByResidue[tRI],siteIndex-1,1);
		thisSitePS      = thisSitePS ^ {{"showpage"}{"/Times-Roman findfont 12 scalefont setfont 500 770 moveto (Site " + siteIndex +") show stroke\nshowpage"}};
		
		fprintf			(psOutput, "\n" , thisSitePS);
		subBufferer = ""; subBufferer * 128;
		
		doComma = 0;
		for (char1 = 0; char1 < 20; char1 = char1+1)
		{
			for (char2 = char1+1; char2 < 20; char2 = char2+1)
			{
				if (char1 != char2 && ((thisSite["COUNTS"])[char1][char2]+(thisSite["COUNTS"])[char2][char1]))
				{	
					ccountf = (thisSite["COUNTS"])[char1][char2];
					ccountb = (thisSite["COUNTS"])[char2][char1];
					if (doComma)
					{
						subBufferer *("<br>");
					}
					doComma = 1;
					subBufferer * (AAString[char1] + "<sub>" + ccountb + "</sub>&harr;<sub>" + ccountf + "</sub>" + AAString[char2]);
				}
			}
		}

		
		subBufferer*0;
		thisResidue [2] = subBufferer;
		subBufferer * 128; 
		subBufferer2 = ""; subBufferer2 * 128;
		char1 = 0;
		
		for (r = 0; r < 20; r=r+1)
		{
			if (myBFs[r])
			{
				if (char1)
				{
					subBufferer * "<br>";
					subBufferer2 *"<br>";
				}
				subBufferer  * (AAString[r] + ":");
				subBufferer2 * (AAString[r] + ":");
				
				if (myBFs[r] > 10000)
				{
					subBufferer *  ("&gt;10<sup>5</sup>");
				}
				else
				{
					subBufferer *  (Format(myBFs[r],6,1));				
				}
				subBufferer2 * guessSelectionKind ((rootIndexByResidue[tRI]) [siteIndex-1], r,bySiteFreqs[siteIndex-1],thisSite["COUNTS"]);
				char1 = char1+1;
			}
		}
		subBufferer*0; subBufferer2 * 0;
		thisResidue [3] = subBufferer;
		thisResidue [4] = subBufferer2;
		thisResidue [5] = _standardizeRatio(FEL_matrix[siteIndex-1][0],FEL_matrix[siteIndex-1][1]);
		thisResidue [6] = Format(FEL_matrix[siteIndex-1][5],8,4);
		bySiteAVL[siteIndex] = thisResidue;
	}
}

byResidueAVL					= {};
labelKeysR						= {};
labelKeysR[0] = "Residue";
labelKeysR[1] = "p-value";
labelKeysR[2] = "Bias term";
labelKeysR[3] = "Proportion of affected sites";
labelKeysR[4] = "Directionally evolving sites";

for (k=0; k<20; k=k+1)
{
	if (parameterEstimates[k][0] != 1)
	{
		residueReport		= {};
		residueReport		[0] = Format(parameterEstimates[k][0],8,4);
		residueReport		[1] = Format(parameterEstimates[k][1],5,2);
		residueReport		[2] = Format(100*parameterEstimates[k][2],5,2) + "%";
		residueReport		[3] = sitesByResidue[k];
	
		byResidueAVL[AAString[k]] = residueReport;
	}
}


summaryAVL												= {};
summaryAVL	["Sequences"]								= baseFilter.species;
summaryAVL	["Sites"]									= baseFilter.sites;
summaryAVL	["Tree Length (subs/site)"]					= Format (referenceL, 5,2);
summaryAVL	["Positively selected sites (FEL)"]		    = felSummary[1];
summaryAVL	["Negatively selected sites (FEL)"]		    = felSummary[0];
summaryAVL	["Directionally selected residues (DEPS)"]  = Abs(byResidueAVL);
summaryAVL	["Directionally selected sites (DEPS)"]		= Abs(bySiteAVL);

fprintf							(htmlOutPath, "<DIV CLASS = 'RepClass' style = 'font-weight: bolder;'>Analysis statistics</DIV><DIV CLASS = 'RepClassSM'>",echoAVLAsTable (summaryAVL, 0, 0, 0, 12), "</DIV>\n");


fprintf							(htmlOutPath, "<DIV CLASS = 'RepClass' style = 'font-weight: bolder;'>Selected residue report </DIV><DIV CLASS = 'RepClassSM'>",echoAVLAsTable (byResidueAVL, labelKeysR, 1, 0, 12), "</DIV>\n");

labelKeys						= {};
labelKeys[0] = "Site";
labelKeys[1] = "Composition";
labelKeys[2] = "MRCA Residue";
labelKeys[3] = "Inferred Substitutions";
labelKeys[4] = "DEPS EBF";
labelKeys[5] = "Selection kind";
labelKeys[6] = "FEL dN/dS";
labelKeys[7] = "FEL p";



fprintf							(htmlOutPath, "<DIV CLASS = 'RepClass' style = 'font-weight: bolder;'>Selected sites report </DIV><DIV CLASS = 'RepClassSM'>",echoAVLAsTable (bySiteAVL, labelKeys, 1, 0, 12), "</DIV>\n");


GetString (timeStamp, TIME_STAMP, -1);
fprintf							(htmlOutPath, "<DIV CLASS = 'RepClassSM' style = 'text-align:right;'>Generated on ", timeStamp, " GMT</DIV>\n</body></html>", CLOSE_FILE);
fprintf							(psOutput, CLOSE_FILE);

/* ------------------------------------------------------------------------------- */

function						guessSelectionKind (root, target, freqs, subCounts)
{
	toTarget = subCounts[-1][target]; toTarget     = ((Transpose(toTarget["_MATRIX_ELEMENT_ROW_!=target"]))*toTarget)[0];
	fromTarget = subCounts[target][-1]; fromTarget = (fromTarget*(Transpose(fromTarget["_MATRIX_ELEMENT_COLUMN_!=target"])))[0];

	
	if (toTarget > fromTarget) /* more to than from */
	{
		if (freqs[target] >= toTarget*2)
		{
			return "(Partial) selective sweeps";			
		}
		else
		{
			return "Convergent evolution/Repeated Substitutions";
		}
	}

	if (freqs[target] >= toTarget*2)
	{
		return "Frequency dependent selection";			
	}
	
	return "Variable site/toggling";
}

/* ------------------------------------------------------------------------------- */

function 						returnSiteSignature (countVector)
{	
	countVector2 = {20,2}["_MATRIX_ELEMENT_ROW_"];
	for (_cc = 0; _cc < 20; _cc = _cc+1)
	{
		countVector2 [_cc][0] = countVector[_cc];
	}
	countVector2 = countVector2 % 0;
	
	sigString = "";
	sigString * 128;
	for (_cc = 19; _cc >=0 && countVector2[_cc][0]; _cc = _cc-1)
	{
		sigString * (AAString[countVector2[_cc][1]] + "<sub>" + countVector2[_cc][0] + "</sub>");
	}
	sigString * 0;
	return sigString;
}

/*--------------------------------------------------------------------------------------------------------------------*/

function		echoAVLAsTable (theData, theKeys, avlKind, avlMap, tFontSize)
{
	_resStr = ""; _resStr * 128;
	if (doPlainText == 0)
	{
		_resStr * "<TABLE BORDER = '0' cellspacing = '1px' cellpadding = '5px' style = 'padding:10px'>";
	}
	if (Abs(theKeys))
	{
		_rc = Abs (theKeys)-1;
		if (doPlainText == 0)
		{
			_resStr * "\n<TR CLASS = 'HeaderClassSM'>";
			for (_k=0; _k <= _rc; _k=_k+1)
			{
				_resStr * ("<TH>" + theKeys[_k] + "</TH>");
			}
			_resStr * ("</TR>\n");
		}
		else
		{
			fprintf (stdout, theKeys[0]);
			for (_k=1; _k <= _rc; _k=_k+1)
			{
				_resStr * ("," + theKeys[_k]);
			}		
		}
	}

	_avlSize = Abs  (theData);
	_avlKeys = Rows (theData);
	
	if (avlKind == 0)
	{
		for (_k=0; _k<_avlSize;_k=_k+1)
		{
			_resStr * ("<TR CLASS = 'ModelClass" + (1+_k%2) + "' style = 'padding:3px; text-align: left; font-size: " + tFontSize + "px;'><TH>"+
							 _avlKeys[_k]+
							 "</TH><TD>"+
							 theData[_avlKeys[_k]]+
							 "</TR>");
		}	
	}
	else
	{
	
		for (_k=0; _k<_avlSize;_k=_k+1)
		{
			if (doPlainText == 0)
			{
				_resStr * ("<TR CLASS = 'ModelClass" + (1+_k%2) + "' style = 'padding:3px; text-align: left; font-size: " + tFontSize+"px; text-indent: 0px;'>");	
			}
				
			if (avlKind == 1)
			{	
				if (doPlainText == 0)
				{
					_resStr * ("<TH>"+_avlKeys[_k]+"</TH>");
					for (_k2 = 0; _k2 < _rc; _k2 = _k2+1)
				    {
				    	_resStr *  ("<TD>"+
									 (theData[_avlKeys[_k]])[_k2]+
									 "</TD>");
									 	
					}
				}
				else
				{
					fprintf (stdout,"\n",_avlKeys[_k]);
					for (_k2 = 0; _k2 < _rc; _k2 = _k2+1)
				    {
				    	_resStr * ("," + (theData[_avlKeys[_k]])[_k2]);
									 	
					}	
				}
			}
			else
			{
		    	_resStr * ("<TH>" + (_k+1) + "</TH>");
		    	
				for (_k2 = 0; _k2 < _rc; _k2 = _k2+1)
			    {
			    	if (_k2 == 0 && Abs(avlMap["LINK"]))
			    	{
			    		theLink = avlMap["LINK"] ^ {{"THELINK",(theData[_avlKeys[_k]])[avlMap[0]]}};
			    		_resStr * ("<TD><a href = '", theLink, "'>" + (theData[_avlKeys[_k]])[avlMap[_k2]] + "</a></TD>");
			    	}
			    	else
			    	{			
			    		_resStr * ("<TD>" + (theData[_avlKeys[_k]])[avlMap[_k2]] + "</TD>");
					}	
				}	
					
			}
							 
			if (doPlainText == 0)
			{
				_resStr * ("</TR>\n");
			}
		}	
		
	}
	if (doPlainText == 0)
	{
		_resStr * ("</TABLE>\n");
	}
	_resStr * 0;
	return _resStr;
}