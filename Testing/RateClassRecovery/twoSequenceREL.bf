/* simulation settings */

sim_settings = {"Globals" :
					{"AC" : 1,
				     "AG" : 1,
				     "AT" : 1,
				     "CG" : 1,
				     "CT" : 1,
				     "GT" : 1},
				"Frequencies" : 
						{{0.25,0.25,0.25}
						{0.25,0.25,0.25}
						{0.25,0.25,0.25}
						{0.25,0.25,0.25}},
				"Branch Length" : 1,
				"Codons" : 10000,
				"omegas" : {{0.01,   0.5},
						     {0.1, 0.25}
						     {1, 0.15}
						     {10, 0.10}
							}
				}; 

/* end simulation settings */

LoadFunctionLibrary ("chooseGeneticCode", {"0" : "Universal"});
LoadFunctionLibrary ("LocalMGREV");
LoadFunctionLibrary ("CF3x4");
LoadFunctionLibrary ("GrabBag");

// define MG94 local w/ CF3x4

corrected_nucleotide_frequencies = CF3x4 (sim_settings["Frequencies"], GeneticCodeExclusions);
PopulateModelMatrix ("MGLocalQ", corrected_nucleotide_frequencies);
codon_frequencies = BuildCodonFrequencies (corrected_nucleotide_frequencies);
Model MGLocal = (MGLocalQ, codon_frequencies, 0);

Tree simulation_tree = (ancestral,derived);

(sim_settings["Globals"])["setParameterValues"][""];

omega_rate_count = Rows (sim_settings["omegas"]);
generate_gdd_freqs (omega_rate_count, "omega_weights", "omega_mixing", "omega_pi", 0);
setFrequenciesToValues ("omega_pi", (sim_settings["omegas"])[-1][1]);


Tree inference_tree  = (ancestral,derived);


defineGDD ("omega_cat", omega_weights, (sim_settings["omegas"])[-1][0]);
ExecuteCommands ("inference_tree.derived.nonSynRate := omega_cat*inference_tree.derived.synRate");

FindRoot (t0, forBLSolver (z) - sim_settings ["Branch Length"], z, 0, 1000);
simulation_tree.derived.synRate = t0;
inference_tree.derived.synRate = t0;

omegas = sim_settings["omegas"];

characters_universal_code = {{"A","C","G","T"}
			  			   {"3",GeneticCodeExclusions,"",""}};

normalizer = sim_settings["Codons"]/(+((sim_settings["omegas"])[-1][1]));

for (segment = 0; segment < Rows (omegas); segment += 1) {
	simulation_tree.derived.nonSynRate = omegas[segment][0] *simulation_tree.derived.synRate;
	
	if (segment>0) {
		DataSet sim   = Simulate 	(simulation_tree,codon_frequencies,characters_universal_code,(omegas[segment][1]*normalizer)$1,0);
		DataSet bigDS = Concatenate (bigDS, sim);
	}
	else
	{
		DataSet bigDS =  Simulate 	(simulation_tree,codon_frequencies,characters_universal_code,(omegas[segment][1]*normalizer)$1,0);
	}
}			  

DataSetFilter codonDF = CreateFilter (bigDS,3,"","",GeneticCodeExclusions);

USE_LAST_RESULTS 			= 1;
VERBOSITY_LEVEL				= 1;

LikelihoodFunction verifier = (codonDF, inference_tree);

LFCompute(verifier,LF_START_COMPUTE);
LFCompute(verifier,logL);
LFCompute(verifier,LF_DONE_COMPUTE);

fprintf (stdout, "N = ", sim_settings ["Codons"], " codons \n");
fprintf (stdout, "Log (L) under the true parameter values = ", logL, "\n");
fprintf (stdout, "Simulated branch length = ", sim_settings["Branch Length"], "\n");

Optimize (results, verifier);
fprintf (stdout, "\nLog (L) under the MLEs = ", results[1][0], "\n");
fprintf (stdout, "Inferred branch length = ", BranchLength (inference_tree, 0), "\n\n");

GetInformation (inferred_rates, omega_cat);

for (k = 0; k < Rows (omegas); k+=1) {
	fprintf (stdout, "Class\t", k+1, " omega ", Format(inferred_rates[0][k], 10,8), " (simulated ", Format(omegas[k][0], 10,8), "), weight ", Format (inferred_rates[1][k], 8, 7), " (simulated ", Format(omegas[k][1], 8,7), ")\n");
}


//----------------------------------------------------------------------------------------


function forBLSolver (value) {
	inference_tree.derived.synRate = value;
	return BranchLength (inference_tree, 0);
}


function setParameterValues (key, value) {
	ExecuteCommands ("`key` = " + value);
}

function setFrequenciesToValues (prefix, values) {
	values   = values * (1/(+values));
	count = Rows (values);
	left_over = 1;
	for (k = 0; k < count-1; k += 1) {
		ExecuteCommands ("`prefix`_" + (k+1) + " = " + values[k] / left_over);
		left_over += (-values[k]);
	}	
}


function defineGDD (prefix, frequencies, omegas) {
	count = Rows (frequencies);
	fm = "`prefix`_weights";
	fr = "`prefix`_rates";
	rm = {count,1};
	for (k = 0; k < count; k+=1) {
		rm[k] = "`prefix`.omega_" + k;
		ExecuteCommands ("global " + rm[k] + " = " + omegas[k]);
		
	}
	
	defineMatrixFromStrings (fm, frequencies);
	defineMatrixFromStrings (fr, rm);
	
	ExecuteCommands ("category `prefix`  = (count, `fm` , MEAN, ,`fr`, 0, 1e25)");
}

function defineMatrixFromStrings (target, source) {
	count = Rows (source);
	ExecuteCommands ("`target` = {count, 1}");
	for (k = 0; k < count; k+=1) {
		ExecuteCommands ("`target`[k] := " + source[k]);
	}
}