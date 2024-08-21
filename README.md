## MEMHC 1.0 (minimal-epitope-for-maximum-MHC-coverage version 1.0)

#Creator; Mohammad Arabpour, mohammad.arabpour.sinior@gmail.com

## Objective(minimal-epitope-for-maximum-MHC-coverage version 1.0)

This software aims to generate a minimal number of peptides that provide maximum MHC and HLA coverage, specifically for use in vaccine and immunotherapy applications targeting T cell epitopes. It takes an affinity file for a  list of HLA alleles against an arraye of peptides that is already generated in affinity prediction models such as mhcflurry, NetMHC, and others as inputs. The script predicts and selects the ones that cover the widest range of HLA alleles with high affinity. The result is a CSV file listing the selected peptides with their corresponding coverage scores an example is showed in the following. This tool is adaptable to various MHC affinity prediction methods, such as mhcflurry, NetMHC, and others.

# Pipeline & run  
The code  provided is a Python script that processes a CSV file containing MHC (Major Histocompatibility Complex) ranking data for a set of peptides. It performs various operations such as filtering, binary conversion of Kd scores, calculating coverage, and generating output files.
Using this code you would be able to use a csv file including columnes with 'id' = prptide ID that could be peptide index and or any attribution , 'peptide' = peptide sequence in AA,'peptide_length'= lenght of peptide, 'start'= starting amino acid position in protein,  the rest f columns are the MHC alleles headings for correspondent peptides and kd(nM) score for each peptide.

The output is a ranking file that give the minimal epitopes(peptides) needed to cover maximum HLA alleles provided.
the applications encompases but not limited to  any immunoassay and vaccine developments need the minimal epitope with maximal HLA coverage.

Here is an example of input csv file :

A part of  protein inflrunza A/Sendai/TU66/2008(H1N1)):   mkvkllvllctft
peptide sequences 8-11 mers
HLA type I supertype binding affinity Kd(nM): HLA-A01:01,HLA-B15:01,HLA-A02:01,HLA-A03:01,HLA-A24:02,HLA-A26:01,HLA-B07:02,HLA-B08:01,HLA-B27:05,HLA-B39:01,HLA-B40:01,HLA-B58:01

id	peptide	peptide_length	start	HL+D1:BM18A-A01:01	HLA-A02:01	HLA-A03:01	HLA-A24:02	HLA-A26:01	HLA-B07:02	HLA-B08:01	HLA-B15:01	HLA-B27:05	HLA-B39:01	HLA-B40:01	HLA-B58:01
Sequence1	MKVKLLVLLC	10	0	37290.2266	24836.6055	29509.4492	38739.0938	40020.6367	36766.1914	17256.5215	24450.7969	12728.3242	15750.3096	31257.5508	25245.1699
.
.
Sequence14	LVLLCTFT	8	5	36421.332	19980.8965	33236.9688	44220.8711	45617.0742	38217.8594	26305.9082	33364.5117	38968.2148	42882.7188	44307.5664	34264.1367

the output csv  file for above would  like :
peptide sequence	# HLA allel hits	accumulativcoverage(%)	Peptide lenght
MKVKLLVLL	1	33.33333333	9
MKVKLLVL	1	66.66666667	8
KLLVLLCTF	1	100	9

## License

This software is licensed under the [CC BY-NC 4.0 License](LICENSE) for non-commercial use.

For commercial use, please contact Mohammad Arabpour at mohammad.arabpour.sinior@gmail.com to obtain a commercial license.
