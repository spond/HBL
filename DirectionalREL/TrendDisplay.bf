AAString    					= "ACDEFGHIKLMNPQRSTVWY";

ExecuteAFile 					(HYPHY_BASE_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "Utility" + DIRECTORY_SEPARATOR + "GrabBag.bf");
ExecuteAFile 					(HYPHY_BASE_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "Utility" + DIRECTORY_SEPARATOR + "AncestralMapper.bf");
ExecuteAFile 					(HYPHY_BASE_DIRECTORY + "TemplateBatchFiles" + DIRECTORY_SEPARATOR + "Utility" + DIRECTORY_SEPARATOR + "ReadDelimitedFiles.bf");

ChoiceList 						(plotType, "Plot Type", 1, SKIP_NONE, 
									"Stacked bars","Stacked bars for all residue frequencies",
									"Frequency trend","Track the frequency of a given residue over");

if (plotType < 0)
{
	return ;
}
if (plotType == 0)
{
	ExecuteAFile					("StackedBarPlot.bf");
}
else
{
	ExecuteAFile					("FrequencyTrendPlot.bf");
}

SetDialogPrompt 				("Select the protein alignment to process");


DataSet							baseAlignment = ReadDataFile (PROMPT_FOR_FILE);
basePath						= LAST_FILE_PATH;
trunkPath						= splitFilePath (LAST_FILE_PATH);
initialFPath					= LAST_FILE_PATH;
extensionStack					= {};
extensionStack[Abs(extensionStack)] = trunkPath["EXTENSION"];
trunkPath						= trunkPath["DIRECTORY"] + trunkPath["FILENAME"];

reducedTrunk					= splitFilePath (trunkPath);

while (Abs(reducedTrunk["EXTENSION"]))
{
	extensionStack[Abs(extensionStack)] = reducedTrunk["EXTENSION"];
	reducedTrunk					= reducedTrunk["DIRECTORY"] + reducedTrunk["FILENAME"];
	reducedTrunk					= splitFilePath (reducedTrunk);
}

trunkPath	= reducedTrunk["DIRECTORY"] + reducedTrunk["FILENAME"];

DataSetFilter					baseFilter	  = CreateFilter (baseAlignment,1);

fprintf							(stdout, "[SET TRUNK PATH TO ", trunkPath, "]\n");

bySiteFreqs						= {};

COUNT_GAPS_IN_FREQUENCIES = 0;
GetString (seqIDs, baseFilter, -1);

for	(s = 0; s < baseFilter.sites; s=s+1)
{
	DataSetFilter 		siteFilter = CreateFilter (baseAlignment,1,siteIndex == s);
	HarvestFrequencies  (siteF, siteFilter, 1, 1, 1);
	
	currentM = siteFilter.species;
	z = 0;
	while (z < 20)
	{
		bySiteFreqs[s] 	= siteF*currentM;
		for (z = 0; z < 20; z=z+1)
		{
			if (((bySiteFreqs[s])[z]+0.5)$1 != (bySiteFreqs[s])[z])
			{
				currentM = currentM - 1;
				break;
			}
		}	
	}
}

fprintf							(stdout, "[COMPUTED SITE-BY-SITE BASE COMPOSITIONS]\n");


ChoiceList 						(yrMap, "Year Mapping", 1, SKIP_NONE, "Direct","Extract from sequence IDs","Premapped","Reload from a previously computed ID->year file");
if (yrMap < 0)
{
	return;
}

if (yrMap > 0)
{
	ExecuteAFile					(basePath + ".map");
	yearMap							= _hyphyAssociativeArray;
	fprintf							(stdout, "[RELOADED SEQID-YEAR MAP]\n");
}
else
{
	yearMap = {};
	theRegExp  = ".*\\|([0-9]+)$";
		
	for (specIndex = 0; specIndex < baseFilter.species; specIndex = specIndex + 1)
	{
		specName  = seqIDs[specIndex];
		specMatch = specName $ theRegExp;
		
		if (specMatch[0]>=0)
		{
			if (Rows (specMatch) > 2)
			{
				key = 0+specName[specMatch[2]][specMatch[3]];
			}
			else
			{
				key = 0+specName[specMatch[0]][specMatch[1]];
			}
		}
		else
		{
			key = 0;
		}
		yearMap[specName] = key;
	}
}

minYear   = 3000;
maxYear	  = 0;

for (seq = 0; seq < siteFilter.species; seq = seq+1)
{
	curYear = yearMap[seqIDs[seq]];
	if (curYear < minYear)
	{
		minYear = curYear;	
	}
	if (curYear > maxYear)
	{
		maxYear = curYear;
	}
}

fprintf (stdout, minYear, "--", maxYear, "\n");
yearSpan = maxYear-minYear+1;


summer	  = {1,20}["1"];
id_taker  = {1,20}["_MATRIX_ELEMENT_COLUMN_"];


haveBaseFit = initialFPath + ".base";

if (!haveBaseFit)
{
	ExecuteAFile (haveBaseFit);
	ancCacheID 						= _buildAncestralCache ("lf", 0);
	rootState						= {};
	
	for	(s = 0; s < baseFilter.sites; s=s+1)
	{
		rootState[s] = (_rootState (ancCacheID,s))["CHAR"];
	}
	
	fprintf (stdout, "[RELOADED BASELINE MODEL FIT]\n");
	haveBaseFit = 1;
}
else
{
	fprintf (stdout, "[BASELINE MODEL FIT IS MISSING]\n");
	haveBaseFit = 0;
}


while (1)
{
	fprintf							(stdout, "Which site (1-based, <=0 to finish)?");
	fscanf							(stdin,"Number", whichSite);
	whichSite = whichSite $ 1;
	if (whichSite < 1)
	{
		break;
	}
	byResidue = {};
	DataSetFilter 		siteFilter = CreateFilter (baseAlignment,1,siteIndex == (whichSite-1));
	for (seq = 0; seq < siteFilter.species; seq = seq+1)
	{
		GetDataInfo (dInfo, siteFilter, seq, 0);
		id = (summer*dInfo)[0];
		if (id == 1)
		{
			id 		  = AAString[(id_taker*dInfo)[0]];
			if (Abs (byResidue[id]) == 0)
			{
				byResidue[id] = {};
			}
			yearID 						= yearMap[seqIDs[seq]];
			(byResidue[id])[yearID] 	= (byResidue[id])[yearID]+1;
		}
	}
	residueCount 	= Abs (byResidue);
	resID		 	= Rows (byResidue);
	rawCounts	 	= {yearSpan,residueCount+1};

	sumByResidue	= {residueCount,1};

	for (r = 0; r < residueCount; r=r+1)
	{
		for (y = minYear; y <= maxYear; y=y+1)
		{
			yi = y-minYear;
			rawCounts [yi][r] 			 = (byResidue[resID[r]])[y];
			rawCounts [yi][residueCount] = rawCounts [yi][residueCount] + rawCounts [yi][r];
			sumByResidue[r] = sumByResidue[r] + rawCounts [yi][r]; 
		}
	}

	sum1 	   = {1,residueCount};
	totalCount = (sum1["1"]*sumByResidue)[0];
	
	if (haveBaseFit)
	{
		fprintf (stdout, "Root state: ", rootState[whichSite-1], "\n");
	}
	
	residuesFound = Rows (byResidue);
	choices		  = {residueCount,2};
	for (r = 0; r < residueCount; r=r+1)
	{
		choices[r][0] = residuesFound[r];
		choices[r][1] = "Residue " + choices[r][0] + "(" + Format (sumByResidue [r]/totalCount * 100, 6,2) + "%)";
	}
	ChoiceList (whichResidue, "Target residue", 1, SKIP_NONE, choices);
	if (whichResidue < 0)
	{
		return 0;
	}
	
	nonZeroYears = 0;
	for (yi = 0; yi < yearSpan; yi=yi+1)
	{
		if (rawCounts [yi][residueCount] > 0)
		{
			nonZeroYears = nonZeroYears + 1;
		}
	}
	
	yearLabels = {nonZeroYears,1};
	resLabels  = {residueCount,1};
	compressed = {nonZeroYears,residueCount+1+2*plotType};
	segmentColors = {nonZeroYears,3};
	
	resLabels[0] = residuesFound[whichResidue];
	for (r = 0; r < residueCount; r=r+1)
	{
		if (r<whichResidue)
		{
			resLabels[r+1] = residuesFound[r];
		}
		else
		{
			if (r>whichResidue)
			{
				resLabels[r] = residuesFound[r];
			}
		}
	}

	step = 1/(residueCount);
	for (r = 1; r < residueCount; r=r+1)
	{
		segmentColors[r][0] = step*r;
		segmentColors[r][1] = step*r;
		segmentColors[r][2] = step*r;
	}

	if (plotType)
	{
		alpha = 0.95;
		FindRoot (chiCut, CChi2(x,residueCount-1)-alpha, x, 0, 1000);
	}

	nonZeroYears = 0;
	for (yi = 0; yi < yearSpan; yi=yi+1)
	{
		if (rawCounts [yi][residueCount] > 0)
		{
			yearLabels[nonZeroYears] = "" + (minYear+yi);
			compressed[nonZeroYears][0] = rawCounts[yi][whichResidue];
			for (r = 0; r < residueCount; r=r+1)
			{
				if (r<whichResidue)
				{
					compressed[nonZeroYears][r+1] = rawCounts[yi][r];
				}
				else
				{
					if (r>whichResidue)
					{
						compressed[nonZeroYears][r] = rawCounts[yi][r];
					}
				}
			}
			compressed[nonZeroYears][residueCount] = rawCounts[yi][residueCount];
			if (plotType)
			{
				ci = estimateMNCI (compressed[nonZeroYears][0],compressed[nonZeroYears][residueCount],chiCut); 
				compressed[nonZeroYears][residueCount+1] = ci[0];
				compressed[nonZeroYears][residueCount+2] = ci[1];
			}
			nonZeroYears = nonZeroYears + 1;
		}
	}
	
	/*
	fprintf (stdout, yearLabels, resLabels, compressed);
	*/
	
	if (haveBaseFit)
	{
		thisSite = _substitutionsBySite (ancCacheID,whichSite-1);
	
		for (char1 = 0; char1 < 20; char1 = char1+1)
		{
			for (char2 = 0; char2 < 20; char2 = char2+1)
			{
				if (char1 != char2 && (thisSite["COUNTS"])[char1][char2])
				{	
					ccount = (thisSite["COUNTS"])[char1][char2];
					fprintf (stdout,  AAString[char1], "->", AAString[char2], " (", ccount, ")\n");
				}
			}
		}
	}

	fileOut = trunkPath + "_" + whichSite + "_" + resLabels[0] + ".ps";
	if (plotType == 0)
	{
		fprintf 		(fileOut,CLEAR_FILE,ProportionBars("compressed",{{400,400,16}},segmentColors,{{"Total samples/year","Sampling year","Observed residue proportions"}},resLabels,yearLabels,1));
	}
	else
	{
		segmentColors = {{0,0,0}{0.8,0.8,0.8}};
		fprintf 		(fileOut,CLEAR_FILE,FrequencyTrend("compressed",{{400,400,9}},segmentColors,{{"Total samples/year","Sampling year","Observed residue proportions"}},resLabels,yearLabels,-1));
	}
}

function estimateMNCI (n,N,A)
{
	p = n/N;
	upper = (A+2*n+Sqrt(A^2+4*N*A*(1-p)*p))/(2N+2A);
	lower = (A+2*n-Sqrt(A^2+4*N*A*(1-p)*p))/(2N+2A);
	return {{lower__,upper__}};
}

