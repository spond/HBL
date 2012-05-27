/*--------------------------------------------------------*/

function determineCoordinateTicks (x1,x2) 
{
	_range 	   = x2-x1;
	_log10Tick = Log(_range)/Log(10) /* round to the next smallest integer */
				 $1;
				 
	_log10Tick = Exp(_log10Tick*Log(10));
	if (_range/_log10Tick < 4)
	{
		_log10Tick = _log10Tick / (((4*_log10Tick/_range)+0.5)$1);
	}
	return _log10Tick;
}


/*--------------------------------------------------------*/

ExecuteAFile  ( HYPHY_BASE_DIRECTORY+"TemplateBatchFiles"+DIRECTORY_SEPARATOR+"Utility"+DIRECTORY_SEPARATOR+"PostScript.bf");


function FrequencyTrend		 (xy&, 			/* 
												NxK matrix with values to plot
											   	each row reprsents a bar
											   	columns are as follows
											   		0 - (K-2)	: contributing values for the label
											   		K-1			: total # of points
											   		
											   		C = K-1 is the number of components 
											*/
											
							  plotDim, 		/* 1x3 matrix {{width, height,font_size}} of the plot in points */
							  colors, 		/* Cx3 matrix of RGB colors to plot each class with */
							  labels, 		/* 1x3 matrix of strings: plot-label, x-axis label, y-axis label*/
							  clabels,  	/* Cx1 string matrix of labels for each class */
							  xlabels,		/* Nx1 string matrix of x-axis labels */
							  trackResidue 
							  )	
{
	
	psDensityPlot = ""; psDensityPlot*1024;
	
	plotHeight = Max (100, plotDim[1]);
	plotWidth  = Max (100, plotDim[0]);
	
	plotOriginX = 4.5*plotDim[2];
	plotOriginY = 3.5*plotDim[2];
	
	numberOfComponents = Columns (xy)-3;
	_x 				   = Rows    (xy);
	
	yMin			   = 0;
	yMax			   = 1;
	
	legendWidth		   = 0;
	
	for (_dataPoint = 0; _dataPoint < 1; _dataPoint = _dataPoint + 1)
	{
		px = Abs (clabels[_dataPoint]) * plotDim[2];
		if (px > legendWitdh)
		{
			legendWitdh = px;
		}
	}
	
	if (legendWitdh)
	{
		legendWitdh = legendWitdh + 2*plotDim[2];
	}
	
	
	psDensityPlot * _HYPSPageHeader (plotWidth + 5*plotDim[2] + legendWitdh, plotHeight + 7*plotDim[2], "Density Plot");
	psDensityPlot * "\n";
	psDensityPlot * _HYPSSetFont ("Times-Roman", plotDim[2]);
	psDensityPlot * "\n";
	psDensityPlot * _HYPSTextCommands(0);
	
	psDensityPlot * "\n 1 setlinewidth 1 setlinecap 0 setlinejoin 0 0 0 setrgbcolor";
	psDensityPlot * ("\n " + plotOriginX + " " + plotOriginY + " " + plotWidth + " " + plotHeight + " rectstroke\n");
	
	
	px 		   		= plotWidth / _x; /* raw width per bar; the bars will occupy bar_width proportion of this space */
	py				= plotHeight;
	bar_width  		= 0.75*px;
	padding_width	= (px-bar_width)/2;
	
	barHeight  		= plotHeight;
	
	totalSum		= 0;
	trs				= 0;
	
	for (_dataPoint = 0; _dataPoint < _x; _dataPoint = _dataPoint + 1)
	{
		xLoc		  = plotOriginX+_dataPoint*px+padding_width;
		yLoc		  = plotOriginY;
		
		for (_segment = 0; _segment < 1; _segment = _segment + 1)
		{
			blah = xy[_dataPoint][_segment];
		
			_seg_height = (xy[_dataPoint][numberOfComponents+2]-xy[_dataPoint][numberOfComponents+1]) * py;

			psDensityPlot * ("" + colors[1][0] + " " + colors[1][1] + " " + colors[1][2] + " setrgbcolor\n");
			
			psDensityPlot * ("newpath " + xLoc + " " 
								+ (yLoc + xy[_dataPoint][numberOfComponents+1] * py) + " " 
								+ bar_width + " " 
								+ _seg_height + " " 
								+ " rectfill\n");
																
								

			psDensityPlot * ("" + colors[0][0] + " " + colors[0][1] + " " + colors[0][2] + " setrgbcolor\n");
			_estimate 	  = plotOriginY + (xy[_dataPoint][0]/xy[_dataPoint][numberOfComponents]) * py;
			
			psDensityPlot * ("newpath " + xLoc + " " 
								+ (_estimate - 1) + " " 
								+ (bar_width) + " " 
								+ 2 + " " 
								+ " rectfill\n");

			psDensityPlot * ("newpath " + (xLoc+bar_width/2) + " " 
								+ (_estimate - 2) + " moveto " 
								+ 0 + " " 
								+ 4 + " " 
								+ " rlineto stroke\n");

			yLoc = yLoc + _seg_height;
			
			if (trackResidue)
			{
				totalSum = totalSum + xy[_dataPoint][numberOfComponents];
				trs		 = trs 		+ xy[_dataPoint][0];
			}
		}
	
		psDensityPlot * ("" + (plotOriginX+(_dataPoint+0.5)*px) + " " + (plotOriginY + barHeight + 0.5*plotDim[2]) +" (" + xy[_dataPoint][numberOfComponents] + ") vrighttext\n");
		psDensityPlot * ("" + (plotOriginX+(_dataPoint+0.5)*px) + " " + (2.5*plotDim[2]) +" (" + (xlabels[_dataPoint])[2][3] + ") centertext\n");
	}
	
	if (trackResidue >= 0)
	{
		psDensityPlot * ("newpath [2] 1 setdash " + plotOriginX + " " + (plotOriginY + trs/totalSum * py) + " moveto " + _x*px + " 0 rlineto stroke [] 0 setdash\n");
	}
	
	
	yscaler = determineCoordinateTicks (yMin,yMax);
	_y	= ((yMin/yscaler)$1)*yscaler;
	while (_y <= yMax)
	{
		yStep = (plotOriginY + py*(_y-yMin));
		psDensityPlot * ("" +  (4*plotDim[2]) + " " + yStep + " (" + Format(_y,4,2) + ") righttext\n");  
		psDensityPlot * ("" +  plotOriginX    + " " + yStep + " moveto "+(0.25*plotDim[2]) +" 0 rlineto stroke\n");  
		_y = _y + yscaler;
	}

	
	psDensityPlot * ("" + (plotOriginX+plotWidth/2) + " " + (0.5*plotDim[2]) +" (" + labels[1] + ") centertext\n");
	psDensityPlot * ("" + (plotOriginX+plotWidth/2) + " " + (plotOriginY + py + 2*plotDim[2]) +" (" + labels[0] + ") centertext\n");
	psDensityPlot * ("" + (plotOriginY+plotHeight/2) + " " + (-1.5*plotDim[2]) +" ("+ labels[2] + ") vcentertext\n");
	
	if (legendWitdh)
	{
		yLoc = plotOriginY + py - 1.5*plotDim[2];
		xLoc = plotOriginX +  _x * px + 0.5*plotDim[2];
		
		for (_segment = 0; _segment < 1; _segment = _segment + 1)
		{
			psDensityPlot * ("" + colors[_segment][0] + " " + colors[_segment][1] + " " + colors[_segment][2] + " setrgbcolor\n");
			psDensityPlot * ("newpath " + xLoc + " " 
									+ yLoc + " " 
									+ plotDim[2] + " " 
									+ plotDim[2] + " " 
									+ " rectfill\n");
			psDensityPlot * ("0 0 0 setrgbcolor\n");
			psDensityPlot * ("newpath " + xLoc + " " 
									+ yLoc + " " 
									+ plotDim[2] + " " 
									+ plotDim[2] + " " 
									+ " rectstroke\n");
									
			psDensityPlot * ("newpath " + (xLoc + plotDim[2]*1.5) + " " + (yLoc+plotDim[2]/6) + " moveto (" + clabels[_segment] + ") show stroke\n");

			yLoc = yLoc - plotDim[2] * 1.5;
			
		}
	}
	
	psDensityPlot * "\nshowpage\n";
	psDensityPlot * 0;
	
	return psDensityPlot;
}

