# F8_standards_map
Generating a variant to barcode map for the Factor 8 standards library with Pacybara intermediate files
This requires that Pacybara was already attempted on mapping the variants to barcodes. The above intermediate files were pulled from the Pacybara directory after the run failed at the clustering step.

Running instructions:
The barcode, genotype, and parameter files were placed in the directory with pacybara_barcode_genotypes.py and run the script was run to generate the map. 

Output columns: 
"barcode" is the 16bp barcode sequence
"reads" is the number of times the barcode was seen
"geno" is the most common ORF-indexed genotype for that barcode
"frac_exact_genotype" is the fraction of times that a particular barcode is paired with the exact match for the listed variant. The non-exact reads are almost always sequencing errors for this set
"min_variant_support" is the fraction of reads of that barcode that contain the rarest of the nucleic acid changes for a codon change
"aa_construct" is the amino acid change in the trucated ORF lacking the B domain
"aa_fl_precursor" is the amino acid change for full length F8
"aa_legacy" is the amino acid change for the F8 orf after celavage of the signal peptide
"kind" is the type of mutation


