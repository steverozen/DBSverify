# DBSverify

This an alpha version.

DBSverify stands for "Doublet Base Substitution verify". 


DBSverify stands for "Doublet Base Substitution verify". The
main function is Read_DBS_VCF_and_BAMs_to_verify_DBSs. This
reads a VCF (variant call format) file containing somatic 
DBS mutations and examines the supporting reads in the 
corresponding tumor and normal BAM files to assess wither 
each DBS is likely real, i.e. arose at one time on a single
chromosome. This is as opposed to being two adjacent single 
base mutations on homologous chromosomes, or adjacent 
mutations that occurred on the same chromosome but at 
detectably different times. Many somatic mutation callers do
not call DBSs at all, and those that do sometime make 
obvious errors, such as calling DBSs that consist of a 
germline SNP next to a somatic single base mutation.

The main function is *Read_DBS_VCF_and_BAMs_to_verify_DBSs*
if you have a VCF with already-called DBSs. If you have a variant
caller that does not call DBSs you can use
*Read_SBS_VCF_and_BAMs_to_verify_DBSs* to combine adjacent
SBSs to be analyzed as DBSs.

To install

remotes::install_github("steverozen/DBSverify")

