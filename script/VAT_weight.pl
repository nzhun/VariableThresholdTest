#!/usr/bin/perl;
use warnings;
use strict;
our $keycol=6;
our $start=0.2;
our $step=0.05;
our $k=5;
our $NR=(1-$start)/$step+1;


sub generate_arr {
  my $cnt=$_[0];
  my $mp=$_[1];
  my @arr=();
  for(my $i=0;$i<$cnt;$i++){
     my $r=rand(1);
	 if($r>$mp){push(@arr,0);}else{push(@arr,1)}
  }

  return (\@arr);
}


sub clean {
	 my ($mat,$phe)=@_;
	 my @matrix2=@$mat;
	 my @phe2=@$phe;
	 my @newmat=();
	 my @newphe=();
	 for (my $i=0;$i<@matrix2;$i++){
		 if($matrix2[$i]==0){next}
		 push(@newmat,$matrix2[$i]);
		 push(@newphe,$phe2[$i]);
	 }
	 return (\@newmat,\@newphe);
}

sub get_best {
 	my $cutoff=0;
	my $p=1;
	my ($new_mat,$new_phe,$mp)=@_;
	my @mat=@$new_mat;
	my @phe=@$new_phe;
	my @cs=();
	my @ct=();
	push @cs, 0 foreach (0..($NR-1));
	push @ct, 0 foreach (0..($NR-1));
#	print "selected K: ".$k."\n";
	my $xbest=0;
	my $tbest=0;
	for(my $x0=0;$x0<1;$x0=$x0+0.1){
		my $sum=0;
		my $dum=0;
	#	print @mat."\n";
		for(my $i=0;$i<@mat;$i++){
			my $weight=1/(1+exp (-$k*($mat[$i]-$x0)));
			$sum+=$weight*($phe[$i]-$mp);
			$dum+=$weight*$weight;
		}
		my $T=$sum/sqrt($dum);
	#	print $x0."\t".$T."\t".$sum."\t".$sum2."\n";
		if(abs($T) > abs($tbest)){
			$xbest=$x0;
			$tbest=$T;
		}
	}
	my @r=($tbest,$xbest);

     return (\@r);
}

sub deal_gene {
	 my ($matd,$phed,$mp)=@_; 
	 my @matrix=@$matd;
	 my @phe=@$phed;
	 my ($new_mat,$new_phe)=clean(\@matrix,\@phe);
	 my $LEFT=@$new_phe;
	 my @rs=(1,1,1,1);
	 if($LEFT==0){return \@rs;}
	 my ($p,$cutoff)=@{get_best($new_mat,$new_phe,$mp)};
	 $rs[0]=$cutoff;
	 $rs[1]=$p;
	 #@rs=($cutoff,$p);
	 my $T1=1000;
	 my $Te=$T1;
	 my $Ts=0;
	 my $pcount=0;
	 my $MAXT=1000000;
	 my $permup=1;
	 LOOP:
	 while($Te < $MAXT+1){
		 for(my $ccct=$Ts;$ccct<$Te;$ccct++){
			  my $new_phe2=generate_arr($LEFT,$mp);
			  my ($p2,$cutoff2)=@{get_best($new_mat,$new_phe2,$mp)};
			  if(abs($p2) <abs($p)){
				  next;
			  }
			  $pcount+=1;
		 }
		 $permup=($pcount+1)/($Te+1);
	#	print $pcount."\t".$Te."\t$p\n";
		 if($pcount >4||$Te==$MAXT){
			 last LOOP;
	      }else{
			  $Ts=$Te;
			  $Te=$Te*10;
	      }
    }
		 my $permup=($pcount+1)/($Te+1);
		 $rs[2]=$permup;
		 $rs[3]=$Te;
		 #push(@rs,$permup);
		 #push(@rs,$Te);
	
	
	return \@rs;
}
sub gen_best {
        my $fin=$_[0];
        my $region=$_[1];
        my @phe=@{$_[2]};
        my $mp=$_[3];
        my $keycol=$_[4];
		  my @includes=@{$_[5]};
        my $lastgene="";
		  my @matrix=();
		  push @matrix, 0 foreach (0..(@phe-1));
		  my @init_matrix=();
		  push @init_matrix, 0 foreach (0..(@phe-1));
        my %result=();
        my $bestT=0;
        my $bestZ=1;
        my $N=10; ##Gt column
		  my $call=0;
        open my $IN, "tabix $fin $region|";
        while(my $line=<$IN>){
				chomp($line);
				if ($line =~ /^#/){next;}
				my @sets=split(/\t+/,$line);
				my @c=@sets[@includes];
				my $site_num=0;
				my $site_denum=0;
				my $revel=$sets[$keycol];
				my $gene=$sets[5];
				if($revel eq "." || $revel<$start){next;}
				if($lastgene eq ""){$lastgene=$gene}
				for(my $i=0;$i<@c;$i++){
					if($c[$i]==0){next;}
					if($matrix[$i] <$revel) {$matrix[$i]=$revel}
					$call=1;
	  		 }		 
 		}
		close $IN;
		my ($cutoff,$p,$permutp,$Te)= @{deal_gene(\@matrix,\@phe,$mp)};
		return($lastgene."\t".$cutoff."\t".$p."\t".$permutp."\t".$Te);
}



sub main {
    #print "input format: inputFile pedfile(case+control)  region  OutPut_prefix\n ";
#	print "input format pedfile: sampleName	Affected_state(1 or 2)\n ";
	if(@_<4){print "check the input parameters please!\n";exit;}
    my ($fin,$fped,$ftrans,$prefix)=@_;
	 print join("\t",@_)."\n";
	my %ped=();
	my @phe=();
	open my $PED, "$fped";
	while(my $line=<$PED>){
		chomp($line);
		if($line =~ /^#/){next;}
		my @sets=split(/\t+/,$line);
		if(@sets <6){
			print  $line." is incomplete pedigree\n";
			exit;
		}
	   $ped{$sets[1]}=$sets[5]-1;
	 }
	close $PED;
	open my $IN,"tabix $fin  -H |" or die "$fin cannot find!\n";
	my $line=<$IN>;
	chomp($line);
	my @sets=split(/\t+/,$line);
	my $outkey=$sets[$keycol];
	my $fout="$prefix.$outkey.txt";
	close $IN;
	my $N=10;
	my $mp=0;
	my @includes=();
	for(my $i=$N;$i<@sets;$i++){
	    if(exists($ped{$sets[$i]})){
	      push(@phe,$ped{$sets[$i]});
		  push(@includes,$i);
           
	    }else{
	      next;
	    }
	    $mp+=$ped{$sets[$i]};
	  }
	$mp=$mp/(@phe);
	my $permut=1000;#1000; #00;
	open my $OUT,">$fout";
	print "the result is outputted to $fout\n";
	print $OUT "#Group\tcutoff\t$p.obv\tp.empirical\tTimes\n";
	open my $TIN, "$ftrans" or die "cannot find $ftrans";
	while(my $tline=<$TIN>){
		chomp($tline);
		my @tsets=split(/\s+/,$tline);
		my $region=join(" ",@tsets[1..(@tsets-1)]);
		print $OUT gen_best($fin,$region,\@phe,$mp,$keycol,\@includes)."\n";
	}
	close $OUT;
	close $TIN; 

}


print "Totally bins ".$NR."\n";
print "Input format: vcf.gz pedgree_file transcript_file output_prefix [K]\n";
print "tabixed vcf.gz.\n"
print "Pedigree format: http://zzz.bwh.harvard.edu/plink/data.shtml#ped\n";
print "transcript format: group[tab]1:start1-end1 1:start2-end2...\n";
print "output_prefix: prefix of the output result\n";
print "K: optional, default:5, used for the weight:1/(1+exp (-k*(x-x0)));\n";
print "run example: perl VAT_weight.pl ../test/test.vcf.gz ../test/test.ped ../test/test.transcript.txt test 10\n;"

if(@ARGV<4){
	print "Input format: fin fped region outpredix\n";
	
	exit;
}

print "Input: ".join(" ",@ARGV)."\n";
my $fin=$ARGV[0]; 
my $fped=$ARGV[1];
my $ftrans=$ARGV[2];
my $prefix=$ARGV[3];
if(@ARGV>4){$k=$ARGV[4];}
main($fin,$fped,$ftrans,$prefix);