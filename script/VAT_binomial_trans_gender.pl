#!/usr/bin/perl;
use warnings;
use strict;
use binomtest;

our %self_dist=();
our $keycol=6;
our $start=0.1;
our $step=0.02;
our $NR=(1-$start)/$step+1;
print "Totally bins ".$NR."\n";
if(@ARGV<4){
	print "required more parameters \n
	Input format: fin fped transcript out_predix\n";
	exit;
}

print "Input: ".join(" ",@ARGV)."\n";
my $fin=$ARGV[0]; 
my $fped=$ARGV[1]; 
my $ftrans=$ARGV[2];
my $prefix=$ARGV[3];
main($fin,$fped,$ftrans,$prefix);
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

sub selectN {
	my ($N, $total)=@_;
	my %sm=();
	my @arr= (0) x $total;
	while(scalar(keys %sm) <$N ){
		my $k=rand($total);
		 $arr[$k]=1;
		 $sm{$k}=1;
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


sub Btest {
	my $num=$_[0];
	my $denum=$_[1];
	my $mp=$_[2];
	my $r=-1;
	if(@_ <3){print "$num:$denum:$mp\nERROR\n";return ($r);}
	if(exists($self_dist{$num.":".$denum})){
		$r=$self_dist{$num.":".$denum}
	}else{
		$r=binomtest::binomtest($num,$num+$denum,$mp);
		$self_dist{$num.":".$denum}=$r;
	}
	return ($r);
}

sub best_chose {
  my @denum=@{$_[0]};
  my @num=@{$_[1]};
  my $mp=$_[2];
  my $bestZ=1;
  my $bestT=$start;
  my $bestk=0;
  my @r=();
  for(my $k=0;$k<$NR;$k++){
    if($denum[$k]==0 && $num[$k]==0){last;}
	my $z0=Btest($num[$k],$denum[$k],$mp);
     if($z0>$bestZ){next;}
     $bestZ=$z0;
     $bestT=$start+$k*$step;
	 $bestk=$k;

   }
if($num[$bestk]==0){$num[$bestk]+=0.5;}
@r=($bestZ,$bestT,$denum[$bestk],$num[$bestk]);
   return \@r;
}

sub get_best {
 	my $cutoff=0;
	my $p=1;
	my ($new_mat,$new_phe,$mp)=@_;
	my @mat=@$new_mat;
	my @nphe=@$new_phe;
	my @cs=();
	my @ct=();
	push @cs, 0 foreach (0..($NR-1));
	push @ct, 0 foreach (0..($NR-1));
	for(my $i=0;$i<@mat;$i++){
		my $hit=int(($mat[$i]-$start)/$step);
		for(my $h=0;$h<$hit+1;$h++){
				 if($nphe[$i]==0){
				 	$ct[$h]=$ct[$h]+1;
				 }else{
				 	$cs[$h]=$cs[$h]+1;
				 }
		 	}

		}
	 my @r=@{best_chose(\@ct,\@cs,$mp)};
     return (\@r);
}

sub deal_gene {
	 my ($matd,$phed,$mp)=@_;
	 my @matrix=@$matd;
	 my @sphe=@$phed;
	 my ($new_mat,$new_phe)=clean(\@matrix,\@sphe);
	 my $LEFT=@$new_phe;
	 if($LEFT==0){return 1;}
	 ## compute the best cutoff and p-value
	 my ($p,$cutoff,$N1,$N2)=@{get_best($new_mat,$new_phe,$mp)};
	 if($N1==0){$N1=0.05;}
	 my $OR=($N2/$N1)*(1/$mp-1);
	 my @rs=($cutoff,$p);
	 my $T1=100;
	 my $Te=$T1;
	 my $Ts=0;
	 my $pcount=0;
	 my $MAXT=10000000;
	 LOOP:
	 while($Te < $MAXT+1){
		 for(my $ccct=$Ts;$ccct<$Te;$ccct++){
			 ##permut new_phe; ## resign phenotype based on the freq, ignore the acutal number of case-control in new_phe
			  my $new_phe2=generate_arr($LEFT,$mp);
			  my ($p2,$cutoff2,$n1,$n2)=@{get_best($new_mat,$new_phe2,$mp)};
			  if($p2 >$p){
				  next;
			  }
			  $pcount+=1;
		 }
		 if($pcount >9||$Te==$MAXT){
			 last LOOP;
	      }else{
			  $Ts=$Te;
			  $Te=$Te*10;
			  print "increase to ".$Te."\n";
	      }
    }
		 my $permup=($pcount+1)/($Te+1);
		 push(@rs,$permup);
		 push(@rs,$Te);
		 push(@rs,$OR);
		 push(@rs,$N2);
		 push(@rs,$N1);


	return \@rs;
}
sub gen_best {
	my ($fin,$gene,$chr,$adphe,$mp,$keycol,$adinclude)=@_;
	my @includes=@$adinclude;
	my @phe=();
	foreach my $ph(@$adphe){push(@phe,$ph); };
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
	my $ostr="";
        print "size: ".@matrix."\t".@phe."\n";
	print join(":",@phe)."\n";
	open IN, "tabix $fin $chr |";
	#print $fin."\t".$chr."\n";
        while(my $line=<IN>){
		chomp($line);
	  	if($line =~ /^#/){next;}
          	my @sets=split(/\t+/,$line);
          	my @c=@sets[@includes];
          	my $site_num=0;
          	my $site_denum=0;
          	my $revel=$sets[$keycol];
          	my $gene=$sets[5];
	#	  if($gene =~ /;/){next}
         	 if($revel eq "." || $revel<$start){next;}
	  	for(my $i=0;$i<@c;$i++){
			  if($c[$i]==0){next;}
			  if($matrix[$i] <$revel) {$matrix[$i]=$revel}
			  $call=1;
		  }
	}
         close IN;
	if($call ==1){
		      my @ret=@{deal_gene(\@matrix,\@phe,$mp)};
		      if(@ret <5){print "$gene went wrong!\n";return($ostr);}
		      my ($cutoff,$p,$permutp,$Te,$OR,$Nc,$Ns)= @ret; #{deal_gene(\@matrix,\@phe,$mp)};
			 return ($gene."\t".$cutoff."\t".$p."\t".$permutp."\t".$Te."\t".$OR."\t$Nc\t$Ns");
		}
		  return ($ostr);
}
#close OUT;



sub main {
    print "input format: inputFile pedfile(case+control)  region  OutPut_prefix\n ";
	print "input format pedfile: sampleName	Affected_state(1 or 2)\n ";
	if(@_<4){print "check the input parameters please!\n";exit;}
    my ($fin,$fped,$trans,$prefix)=@_;
	my %gender=();
	my %ped=();
	my @aphe=();
	open PED, "$fped";
	while(my $line=<PED>){
	   chomp($line);
	   if($line =~ /^#/){next;}
	   my @sets=split(/\t+/,$line);
	   $ped{$sets[0]}=$sets[@sets-1]-1;
	   if(@sets>=3){
		$gender{$sets[0]}=$sets[@sets-2];
	   }
	 }
	close PED;



	open IN,"tabix $fin  -H |" or die "$fin cannot find!\n";
	print "tabix $fin -H\n";#open OUT,">$fout";
	my $line=<IN>;
	chomp($line);
	my @sets=split(/\t+/,$line);
	my $outkey=$sets[$keycol];
	my $fout="$prefix.$outkey.txt";
	close IN;
	my $N=10;
	my $mp=0;
	my @includes=();
	my $total=0;
	my $NX=0;
	for(my $i=$N;$i<@sets;$i++){
	    if(exists($ped{$sets[$i]})){
	    	  push(@aphe,$ped{$sets[$i]});
		  push(@includes,$i);
		  if(exists($gender{$sets[$i]})){
			$total=$total+$gender{$sets[$i]};
			$NX=$NX+$gender{$sets[$i]}*$ped{$sets[$i]};	
		  }
	    }else{
	      next;
	    }
	   #print $sets[$i]."\t".$ped{$sets[$i]}."\n";
	    $mp+=$ped{$sets[$i]};
		
	  }

	print "Total samples:".@aphe.", cases: $mp\n";
	my $Ncase=$mp;
	$mp=$mp/(@aphe);
        my $Amp=$mp;
	my $permut=1000;#1000; #00;
	#print $mp."\t".join(":",@phe)."\n";
	open OUT, ">$fout" ;
	open TRIN,$trans or die "$trans cannot find!\n";
	while(my $tl=<TRIN>){
		chomp($tl);
		my @tsets=split(/\s+/,$tl);
		my $gene=$tsets[0];
		my $region=join("\t",@tsets[1..(scalar(@tsets)-1)]);
		if($region =~/X/){$mp=$NX/$total;print "MP in X:$mp\t$gene\n";}else{$mp=$Amp; print "MP in autosomal:$mp\t$gene\n";}
		my $st=gen_best($fin,$gene,$region,\@aphe,$mp,$keycol,\@includes);
		if($st ne "") {print OUT $st."\n";}
	}
	close	TRIN;
	close OUT;
}

## /home/local/ARCS/nz2274/PAH/PAH_10032017/Result/Type1_Error/data/refGene_mRNA_protein_coding_hg37_Ensemble95_4VT.cannoical.bed
