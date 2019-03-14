#!/usr/bin/perl -w
use POSIX;
use List::Util qw(min);

# file: miRhub.pl
# Arguments: 
# 
# $mirlist: list of miRNAs to scan for target enrichment
# $list: list of genes to examine
# $Name0: Output file name
# $ITER: # of iterations to run
# $SPECIES: 9606(human) or 10090(mouse)
# $CONS: # of species (N<=5) for which the targets site must be conserved to count
# $NET: Mode of operation (list=0, net = 1)
# $ALPHA: Value of Alpha in scoring scheme (weight of protein interacions in scoring).


if ($#ARGV < 4) {
   print "Usage: miRhub.pl miRlist.txt GeneList.txt OutName ITER CONS (NET) (ALPHA)\n";
   print "miRList.txt: text file containing miR list (1 perl line).\n";
   print "GeneList.txt: text file containing gene list (1 perl line).\n";
   print "OutName: Base name of output files\n";
   print "ITER: number of iterations you want to run\n";
   print "SpeciesID: 9606(human) or 10090(mouse) or 10116(rat)\n";
   print "CONS: {0,1 .. 4} is the required number of species the target site must be present in.\n";
   print "Net: {0,1} if you want a network score use 1, a genelist score is 0\n";
   print "Alpha: coefficient on log10(degree) of gene\n";
   exit();
}
$nArgs = scalar(@ARGV);
$mirlist=shift;
$list=shift;
$PPI=shift;
$Name0=shift;
$ITER=shift;
$species=shift;
$CONS=shift;
$NET =0; 
$ALPHA=3;

print "miRNA list: $mirlist\n";
print "geneList: $list\n";
print "Protein interactions: $PPI\n";
print "Output: $Name0\n";
print "# Random simulations: $ITER\n";
print "Required conservation: Human + $CONS others\n";
print "Mode: ";
print "gene list\n";


#$miRConLev='betamiR_cons_leve2.txt';
$miRConLev='/home/pr46_0001/shared/lib/python3.6/miRhub/miR_family_conservation.txt';
open (MIR,$miRConLev)  or die;
while (<MIR>) {
   chomp();
   my($a,$b) =split(/\t/,$_);
   #record the conservation level of each miRNA
   $mirCL{($a)} = $b-1;
}
close(MIR);
print $mirlist;
open (MIR,$mirlist) or die;
while (<MIR>) {
   chomp();
   #determine if the miRNA is conserved in the minimum number of species
   if ((exists($mirCL{$_}))) {
	 if ($mirCL{$_}>=$CONS) {
	    $miRNAset{$_} =1;
	 }
   }
}
close(MIR);
if ($species==9606) {
   $ScoreCard='/home/pr46_0001/shared/lib/python3.6/miRhub/scorecard_HSA.txt';
}
elsif ($species==10116) {
   $ScoreCard='/home/pr46_0001/shared/lib/python3.6/miRhub/scorecard_RNO.txt';
}
else {
   $ScoreCard='/home/pr46_0001/shared/lib/python3.6/miRhub/scorecard_MMU.txt';
}

$info=`date`;
print "$info\t$ScoreCard\n";

open (OUT,">$Name0");

print $ScoreCard;
open (SC,$ScoreCard) or die;
# Read in predicted targets for all genes/miRNAs

$six_c=0;
while (<SC>) {
   chomp();
   
   my($gene,$transcript,$utrLen,$miRgroup,$type,$pos,$cons) = split(/\t/,$_);
   if ($type eq '6mer') {
       $six_c ++;
       next;
      }
   if (exists ($miRNAset{$miRgroup})) {
      if ($type eq '7mer-m8') {
	 $score =1.25;
      }
      elsif ($type eq '8mer-1a') {
	 $score = 1.5
      }
      elsif ($type eq '7mer-1a') { 
	 $score = 1;
      }
      else {
          print "Not recognize binding site type: $type\n"
      }
      @conList = split(":",$cons);
      $nCon = scalar(@conList) -1; # -1 to remove human (9606);

      if ($nCon >= $CONS) { # only record targeting info in genTgt if at conservation threshold

	 if (!exists($miRtgts{$miRgroup}{$gene})) {
	    $miRtgts{$miRgroup}{$gene}=0;
	 }
         #score without clustering for debug only (overwritten later)
	 $miRtgts{$miRgroup}{$gene}+=($score);

	 my ($st,$ed) = split(":",$pos);

         #individual target score for each gene/miRNA/base position
	 $genTgt{$gene}{$miRgroup}{$st} = $score;

      }
   }
   #record the longest UTR length of each gene
   if (exists($gutrLen{$gene})) {
      $gutrLen{$gene} = $utrLen  if ($utrLen>$gutrLen{$gene});
   }
   else {
      $gutrLen{$gene} = $utrLen ;
   }
}
close(SC);

#Cluster miRNA target scores by position within the UTR
foreach $g (keys(%genTgt)) {
   foreach $m (keys(%{$genTgt{$g}})) {
      # sort start posiont in numerical order
      @starts = sort { $a <=> $b } keys(%{$genTgt{$g}{$m}});

      # start with the first position
      $prevSt=$starts[0];
      $score = $genTgt{$g}{$m}{$starts[0]};

      $k=1; # position within the cluster
      for ($ip =1 ; $ip<scalar(@starts); $ip++) {
	 $dist = $starts[$ip] - $prevSt;
	 if (($dist > 8) && ($dist < 60)) {
	    $mod = 1.5;
	    $mod = 1 if ($k==0);
	    $score += ($genTgt{$g}{$m}{$starts[$ip]}*$mod);
	    $k++;
	    $prevSt=$starts[$ip];
	 }
	 elsif ($dist>=60) {
	    $prevSt=$starts[$ip];
	    $score +=$genTgt{$g}{$m}{$starts[$ip]};
	    $k=1;
	 }
	 else { # closer than 8
            # re-score the close targets
	    if ($genTgt{$g}{$m}{$starts[$ip]} > $genTgt{$g}{$m}{$prevSt}) {
	       $k--; # get previous position
	       $mod = 1.5; # set modifier
	       $mod = 1 if ($k==0); #modifier is 0 for first score in cluster
	       $score -= ($genTgt{$g}{$m}{$prevSt}*$mod); # remove old score for first target site
	       $score += 0.5*($genTgt{$g}{$m}{$prevSt}); # add a reduced score for the first target site

	       $prevSt = $starts[$ip]; # move to next target site
	       $score += ($genTgt{$g}{$m}{$prevSt}*$mod); # Add score to previous cluster

	       $k++; # move to next position in cluster
	    }
	    else { # newest site has a smaller score
	       $score += 0.5*($genTgt{$g}{$m}{$starts[$ip]}); # add reduced score for second target site
               #do not advance cluster position
	    }
	 } # end closer than 8
      } # end cluster loop

      $miRtgts{$m}{$g} = $score ; # store accumulated score for miRNA/Gene pair

   } #end miRNA loop
} #end gene loop


%bComp=(); # results hash

# Check to see if gene in list is in the score card
open (GLIST,$list) or die "Cant open genelist!"; # open input gene list
while (<GLIST>) {
   chomp();
   if (exists($gutrLen{$_})) { # only add if the gene is in our list of gene (has a UTR)
      $g0{$_}=1; # Orig gene list
      $glist0{$_}=1; # full node list
      $ginter0{$_}=0; # interaction count : initialize 
      $beta0{$_}=1000 / $gutrLen{$_}; # normalizer
      $GINTER{$_}{NA} =1; #Full list of interactions : initialize
   }
   elsif (exists($gutrLen{uc($_)})) {
      $geneSymbol=uc($_);
      print "Warning $_ converted to $geneSymbol\n";
      $g0{$geneSymbol}=1; # Orig gene list
      $glist0{$geneSymbol}=1; # full node list
      $ginter0{$geneSymbol}=0; # interaction count : initialize 
      $beta0{$geneSymbol}=1000 / $gutrLen{$geneSymbol}; # normalizer
      $GINTER{$geneSymbol}{NA} =1; #Full list of interactions : initialize
   }
   else {
      print "Warning discarding gene: $_ (no UTR found)\n";
   }
}
close(GLIST);
if (scalar(keys(%g0)) ==0) {
   print "ERROR no UTR found for any input gene: Are your gene symbols correct ?\n";
   exit();
}

@Agenes=();
# open protein protein interaction file
open (PPIFILE,$PPI) or die "Cant open protein network file!";
while(<PPIFILE>) {
   chomp();
   my ($A,$D,$B) = split("\t",$_); # read gene, Degree(not used) , list of genes
   #$D = degree in origional PPI network without restriction to high confidence nodes
   if (exists ($gutrLen{$A}) ){ # only add if the gene is in our list of gene (has a UTR)
      $GINTER{$A}{NA} =1; #initialize
      push (@Agenes,$A);
      if ($B eq "") { # if this gene (A) doesn't interact with any others
	 $Ggroup{0}{$A} =1; # add gene to the "0" group (Has 0 interactions)
      }
      else {

	 @genePartners = split(":",$B); # split genelist into individual genes

	 $Ncat=0;

	 foreach $g(@genePartners) {

	    if (exists($gutrLen{$g})) { # if partner has UTR
	       $GINTER{$A}{$g} =1; # Add interaction A->g into full network

	       if (exists($g0{$A})){# if gene A in input gene list
		  $glist0{$g} = 2 if (!exists($g0{$g})); # add node g to node list
		  $ginter0{$A} = 0 if (!exists($ginter0{$A}));
		  $ginter0{$A} += 1; # add to # of interactions for node A
	       } # end if A is source
	       $Ncat++; #increment interaction counter for node A
	    }
	 } # end for each partner
	 $Ggroup{$Ncat}{$A} = 1; # Add Node A to the (count = Ncat ) group
      } # end else
   } #end if has a UTR
}
print "PPI file opening complete!\n";
close(PPIFILE);



@degCats = sort { $a <=>$b } keys (%Ggroup); # sort degree categories (full list of genes with degree N)


# score graph 0 -> This is the input gene network/list
scoreGraph (\%miRtgts, \%glist0, \%ginter0, \%bComp ,1); #\%beta0);

#Output results
print "Working on writing the output results!\n";
@miRarry = keys(%miRtgts);
foreach $m (@miRarry){
   print OUT "\t$m";
}
print OUT "\n";
foreach $m (@miRarry){
   $val = $bComp{$m}{ref};
   print OUT "\t$val";
}
print OUT "\n";


#Iterate generation and scoring of random graphs
print "Entering Monte-Carlo simulation\n";
for ($i=0; $i<$ITER; $i++) {
   print "Iteration #: $i\n";
   %gI=(); #Random gene list
   %glistI=(); #Random node list (orig + neighbors)
   %ginterI=(); #Interaction count
   %ginfoI=(); #Info hash

   $utr0=0;
   $utrI=0;

   foreach $g(keys(%g0)){ # loop through input gene list

      $Ncat=scalar(keys(%{$GINTER{$g}}))-1;  # Get the number of partners of gene g (-1 to remove NA)


      @gNames=();
      if ($NET==1) {
	 push(@gNames,keys(%{$Ggroup{$Ncat}})); # get all gene names of genes with Ncat partners
      }
      else{
	 push (@gNames,@Agenes);

      }
      $nRandHubs=scalar(@gNames); #count genes of degree N (with Ncat partners)
      $nHubs = $nRandHubs;
      my $index = 0;
      my $max_dist =0;
      # get index (in degree array) of target degree
      ++$index until $degCats[$index] == $Ncat or $index > $#degCats;
      $ii=1;

      while ($nHubs < 20) { #while choices < 20 add genes from bins on both sides
	 if (($index+$ii) <= $#degCats){
	    $Npos=$degCats[$index+$ii] ; #pos side: index + ii
	    push (@gNames,keys(%{$Ggroup{$Npos}})); # add gene names
	    $max_dist = abs($Npos-$Ncat) if (abs($Npos-$Ncat) > $max_dist); # record distance
	 }
	 if (($index-$ii) >= 0 ){ #Neg side: index - ii
	    $Nneg=$degCats[$index-$ii] ;
	    push (@gNames,keys(%{$Ggroup{$Nneg}})); #add genes
	    $max_dist = abs($Nneg-$Ncat) if (abs($Nneg-$Ncat) > $max_dist); # record distance
	 }
	 $nHubs=scalar(@gNames); # count genes in list
	 $ii++; #increment bin width
      }
      $idx=int(rand($nHubs)); # choose a random gene from the list
      while (exists($gI{$gNames[$idx]})) {
	 $idx=int(rand($nHubs)); # choose a random gene from the list
      }

      if ($nRandHubs<20) {
	 print "Warning: only $nRandHubs genes of degree $Ncat\t$g: expanded to $nHubs \t$max_dist\n" if ($i ==0);
	 print "$idx / $nHubs: $gNames[$idx] \n" if($i==0); 
         # warning message for low count bins
      }


      unless($nHubs<1) { 
	 $gI{$gNames[$idx]}=1; # add random gene to list 
	 $glistI{$gNames[$idx]} = 1; # add random gene to full node list

         #$betaI{$gNames[$idx]} = 1000 / $gutrLen{$gNames[$idx]}; # normalizer
	 $utr0+=$gutrLen{$g}; # counter for input utr gene length
	 $utrI+=$gutrLen{$gNames[$idx]}; # counter for random gene length

#
	 $ginfoI{$gNames[$idx]}{source} = $g; # record info
	 $ginfoI{$gNames[$idx]}{degree} = $Ncat;
	 $ginfoI{$gNames[$idx]}{groupSize} = $nHubs;
      }
   } # end random selection

   # add partners
   foreach $gi(keys(%gI)) {
      @genePartners = keys(%{$GINTER{$gi}}); # partners of gene gi
      foreach $pi(@genePartners) {
	 unless ($pi eq "NA") {
	    $glistI{$pi} = 2 if (!exists($glistI{$pi})); # add all partners to node list
	 }
      }
   }
   foreach $gi(keys(%glistI)) {
      $ginterI{$gi} = scalar(keys(%{$GINTER{$gi}}))-1; # add iteraction count
   }

   # score graph I
scoreGraph (\%miRtgts, \%glistI, \%ginterI, \%bComp, $utr0/$utrI) ; #\%betaI); 

   # output results
   foreach $m (@miRarry){
      $val = $bComp{$m}{rand};
      print OUT "\t$val";
   }
   print OUT "\n";


   foreach $g(keys (%ginterI))
   {
      delete $ginterI{$g};
      delete $glistI{$g};
      delete $ginfoI{$g};
      delete $gI{$g} if (exists($gI{$g}));
   }
} #end iteration
foreach $m (@miRarry){
   $val ='';
   $val = $bComp{$m}{tc} if (exists($bComp{$m}{tc}));
   print OUT "\t$val";
   print "$m\t$val\n";
}
print OUT "\n";

$info=`date`;
print "$info\n";


close(OUT);

#scoring function
sub scoreGraph {
   my $miRtgt = $_[0]; # first imput target hash
   my $paramList = $_[1]; # geneList
   my $paramGraph = $_[2]; # degree count
   my $paramBeta = $_[4]; # normalizer
   my %miRtgtH = %$miRtgt; # conversion to hash
   my %LiH = %$paramList; # conversion to hash
   my %GrH = %$paramGraph; # conversion to hash


   foreach $mirG (keys(%miRtgtH)) { # loop through miRNAs

      $tgtCount=0; # 
      $count=0;
      $geneCount=0;
      $tgtList='';
      foreach $g(keys(%LiH)){ # loop thourgh genes
	 if (exists($miRtgtH{$mirG}{$g})) { # if miRNA m targets gene g
	    if ($LiH{$g} ==1) { # only score nodes in orig List
	       $nPPI = $GrH{$g}; # get degree of gene g
	       $ppi_part = 0;
	       if ($NET==1) {
		  $ppi_part = $ALPHA*(log10($nPPI)) if ($nPPI>0); # if network informed add log 10 of degree
	       }

	       $isc = ($miRtgtH{$mirG}{$g} *(1+ $ppi_part));
               $isc*=$paramBeta; #$utrLH{$g}; #normalize score
	       $count += $isc;
	       $tgtCount++; # count targets of miRNA m
               $tgtList=join(',',$tgtList,$g); # record list of targets 

	    }
	 }
	 $geneCount++ if ($LiH{$g} ==1) ; # count genes

      } # end gene loop
      $geneScore=($count/$geneCount); 
      if (!exists (${$_[3]}{$mirG})) { # output score
	 ${$_[3]}{$mirG}{ref}=$geneScore;
	 print "$mirG\t$tgtList\n";
	 ${$_[3]}{$mirG}{tc}=$tgtList; # output target count of miRNA
      }
      else{
	 ${$_[3]}{$mirG}{rand}=$geneScore;
      }
   }
}
