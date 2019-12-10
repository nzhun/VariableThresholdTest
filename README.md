# VariableThresholdTest

This script required the following package:
Math::GSL::CDF
Math::GSL::Randist
usage:
## hard threshold command line:
perl VT_binomial.pl $pwd/test/test.vcf.gz $pwd/test/test.ped $pwd/test/test.transcript.txt  test.OUT
## soft threshold command line:
perl VT_weight.pl $pwd/test/test.vcf.gz $pwd/test/test.ped $pwd/test/test.transcript.txt  test.OUT [K]

notes: 
vcf file has to be bgzipped and tabixed, it is not a regular vcf format, please follow format shown in test.vcf.gz.
ped file has 6 columns.
transcript format: transcriptName[tab]chr:start-end[space]chr:start2-end2
   it starts with transcript name, tab-delimited, coding regions delimited by a space.
 