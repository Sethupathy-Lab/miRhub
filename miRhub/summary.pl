#!/usr/bin/perl -w
use Math::Complex;
use List::Util qw[min max];
#Calculate avg/ sd and pValue for each miR
#usage: summary OutName geneList


$file =shift;
$input = shift;

open(FILE,$file) or die;
$header = <FILE>; # get miR list
chomp($header);
@hparts = split(/\t/,$header);

$test=<FILE>; # get test scores
chomp($test);
@tparts = split(/\t/,$test);
for ($i=1; $i < scalar(@hparts) ; $i++){ # calculate stats for test 
   $testH{$hparts[$i]} = $tparts[$i]; #hash score by miRNA name
   $randH{$hparts[$i]} = 1; #ititialze count to 1
   $s1H{$hparts[$i]} = 0; #init sum
   $s2H{$hparts[$i]} = 0; #init sum^2
}
$N=0;
while(<FILE>) {
   chomp();
   $N++; #count iterations
   @parts=split(/\t/,$_);
   if( $_ !~ m/,/) {
      for ($i=1; $i < scalar(@hparts) ; $i++){
	 $randH{$hparts[$i]} ++ if ($parts[$i]>= $testH{$hparts[$i]}); # count # time random score >= test score
	 $s1H{$hparts[$i]} += $parts[$i]; #sum scores
	 $s2H{$hparts[$i]} += ($parts[$i] * $parts[$i]); # sum square of scores
	 $maxS{$hparts[$i]} = $parts[$i] if (!exists($maxS{$hparts[$i]}) || ($parts[$i] > $maxS{$hparts[$i]})); # record max score
      }
   }
   else {
      for ($i=1; $i < scalar(@parts) ; $i++){
	 $gL{$hparts[$i]} = $parts[$i];
      }
   }
}
foreach $k (keys(%randH)) {
   push(@{$rLU{$randH{$k}}},$k);
}
@pVarr = sort {$a <=> $b} keys(%rLU);
@sigMirs=();
foreach $k (@pVarr) {
   @tmp =sort {$testH{$b} <=> $testH{$a}} @{$rLU{$k}};
   push (@sigMirs,@tmp);
}

close(FILE);

open (OUT,">$file.out2.txt");
print OUT "miR\tScore\tAvgRandScore\tsd(randScore)\tmax(randScore)\tRank\tcount>=\tpValue\tFDR(0.05) thresh\tBH-sig\tpv_corrected\tpv_adjusted\n"; # print header
#@sigMirs = sort {$randH{$a} <=> $randH{$b}} keys(%testH); # sort miRNAs by pValue
@revSigMirs = reverse(@sigMirs);

$V=scalar(@sigMirs);
$alpha=0.05;

$rank=0;
$fdr_k=0;
foreach $s (@sigMirs) {
   $rank++;
   $pval{$s}= $randH{$s} / ($N+1); # calculate pValue 
   $miRrank{$s} = $rank;

   $fdr_thresh{$s}= ($rank)*$alpha/($V);
   $pv_cor{$s} = $pval{$s} * $V / $rank;
   if ($pval{$s} <= $fdr_thresh{$s}) {
      $fdr_k = $rank;
   }
}
$ref = 1;
foreach $r (@revSigMirs) {
   # Yekutieli & Benjamini (1999), equ #3: Adjusted p-value (monatonic)
   $pv_adj{$r} = min($ref,$pval{$r}*$V/$miRrank{$r});
   $ref = $pv_adj{$r};
}


foreach $m (@sigMirs){
   $testSc=$testH{$m}; 
   $avSc = $s1H{$m} / $N; # calculate avg score = sum / count
   $sdSc = sqrt( $N*$s2H{$m} - ($s1H{$m}*$s1H{$m})) / $N; #calculate sd(score) 
   $sig=0;
   $sig=1 if ($miRrank{$m}<=$fdr_k);
   $line = join ("\t",$m,sprintf("%.4f",$testSc),sprintf("%.4f",$avSc),sprintf("%.4f",$sdSc),sprintf("%.4f",$maxS{$m}),$miRrank{$m},sprintf("%d",$randH{$m}),sprintf("%f",$pval{$m}),sprintf("%f",$fdr_thresh{$m}),sprintf("%d",$sig),sprintf("%f",$pv_cor{$m}),sprintf("%f",$pv_adj{$m}));
   print OUT "$line";

   if (exists($gL{$m})) {
      print OUT "\t$gL{$m}" ;
   }

   print OUT "\n";
}

print OUT "\n\n";

close (OUT);
