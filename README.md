# VariableThresholdTest

usage:
perl VT_binomial.pl $pwd/test/test.vcf.gz $pwd/test/test.ped $pwd/test/test.transcript.txt  test.OUT

notes: 
vcf file has to be bgzipped and tabixed
ped must only contains ID,gender and status
transcript format: transcriptName[tab]chr:start-end[space]chr:start2-end2
   it starts with transcript name, tab-delimited, coding regions delimited by a space.
 