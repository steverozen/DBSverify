# DBSverify

This an alpha version.

DBSverify stands for "Doublet Base Substitution verify". 
The main function is *ReadVCFAndBAMsAndVerifyDBS*.

DBSverify stands for "Doublet Base Substitution verify". 
The main function is ReadVCFAndBAMsAndVerifyDBS. This reads 
a VCF (variant call format) file for somatic DBS mutations
and examines the supporting reads in the corresponding tumor and
normal BAM files to assess wither each DBS is likely real,
i.e. arose at one time on a single chromosome. This is
as opposed to being two adjacent single base mutations on
homologous chromosomes, or adjacent mutations that occurred
on the same chromosome but a detectably different times. Many
somatic mutation callers do call DBSs at all, and those that
do sometime make obvious errors, such as calling DBSs that
consist of a germline SNP next to a somatic single base 
mutation.

To install

remotes::install_github("steverozen/DBSverify", ref = "master")

