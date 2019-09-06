#!/usr/bin/perl;
package binomtest;
use warnings;
use strict;
use Math::GSL::CDF qw/:binomial/;
use Math::GSL::Randist qw/:binomial/;

sub binomtest {
	my ($num,$total,$mp)=@_;
    	my $gsl_p=gsl_ran_binomial_pdf($num,$mp,$total);
    	my $gsl_cdf=1-gsl_cdf_binomial_P($num,$mp,$total);
    	my $r=$gsl_p+$gsl_cdf;
    	return ($r);
}

1
