# codon-project-zaitlen-lab
The following files are used in the Zaitlen Lab for the Codon Project. Sorry in advance if these files are messy.

GWAS data was gathered from the Neale Lab:
https://docs.google.com/spreadsheets/d/1kvPoupSzsSFBNSztMzl04xMoSC3Kcx3CrjVf4yBmESU/edit#gid=178908679

The scripts were run in the order they appear.

codon_project_script.R - script run to combine the variant file with the phenotype information 
				for the GWAS data. Filters the data to only include synonymous variants. 
				Data is named according to their phenotype. Ex. "combined_height_data.tsv" for height phenotype

rest_api.R - script used to query the Ensembl database. Takes all the synonymous SNPs from GWAS data and gets their annotation,
		in particular, information about the ancestral allele, codons, and SNP location. Removed SNPs with more than two alleles. 
		Data is written out as "ancestral_and_derived_from_rsid.tsv".

codon_frequency.R - script written to calculate the relative codon frequency of the each codon relative to their amino acid. 
			  Takes the cds of each gene from Ensembl build GRCh37 and counts each codon to get the total occurance in the genome. 
			  Resulting codon usage table is named "codon_frequency_table.tsv" 

delta.R - splits the codons from rest_api into the ancestral and derived codon along with the frequency of the ancestral and derived codon. 
	    calculates the change in codon frequency and appends it to the data. Data is written out to "ancestral_and_derived_with_delta.tsv"

lin_regress_scripts folder - folder containing the scripts for the linear regression models for (effect size ~ delta codon freq) run on the individual phenotypes. 
				     Ran linear regressions while Partitioning with smaller p-values of the phenotypes to look at more significantly associated SNPs. 
				     These were the original scripts before automating them into a table. Changed the sign of beta if the 
				     ref allele != ancestral allele before running the regressions. This is because the beta is based on the minor allele and 
				     thus the if the minor allele == ancestral allele, then switch the sign to get the right direction of effect. 

plotting.R - script used to output graphs for our data and models used for our capstone presentation. plots codon usage table, linear regressions for 
		 significant relationships and histogram of change in codon frequency for throughout the entire dataset. 

table.R - the same as the lin_regress scripts folder, but the coefficients of the linear model are inputted into a table. 
	    These values are saved as "lin_regress_data.tsv"

derived_AF_graph.R - script used to make plot for change in codon frequency based on the derived allele frequency. The allele frequency 
			   was split into quantiles before graphing. Allele freq was based off of the variant file. 

percentage_table.R - linear regressions run on effect size against the percentage change of the codon frequency. The coefficients of the linear regression were 
			   inputted into a table. These values are saved as "lin_regress_change_percent.tsv"

abs_change_table.R - linear regressions run on absolute value of effect size against the absolute value of percentage change. The coefficients of the linear 
			   regression were inputted into a table. These values are saved as "lin_regress_abs_change_percent.tsv" 
