A self-contained HBL to simulate two sequences under the MG94 REV model with multiple rate categories on &omega; 
and then try two recover them using a REL model.

Edit the ```sim_settings``` global to change simulation settings.

Example output:

```
N = 10000 codons
	      1232 unique pairs
	      7057 (  70.570%) are invariable
	      1840 (  18.400%) have 1 nucleotide differences
		      1144 (   0.622%) are synonymous
		       696 (   0.378%) are non-synonymous
	       691 (   6.910%) have 2 nucleotide differences
		        25 (   0.036%) are synonymous
		       666 (   0.964%) are non-synonymous
	       412 (   4.120%) have 3 nucleotide differences
		         1 (   0.002%) are synonymous
		       411 (   0.998%) are non-synonymous

Log (L) under the true parameter values = -55915.42597225038
Simulated branch length = 0.25

Log (L) under the MLEs = -55911.21018871951
Inferred branch length = 0.2407338058239624

Class	1 omega 0.09619144 (simulated 0.01000000), weight 0.5103205 (simulated 0.5000000)
Class	2 omega 0.16319436 (simulated 0.10000000), weight 0.2469205 (simulated 0.2500000)
Class	3 omega 0.23286364 (simulated 1.00000000), weight 0.1250200 (simulated 0.1500000)
Class	4 omega 8.00657195 (simulated 10.00000000), weight 0.1177390 (simulated 0.1000000)
```
