#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Pod::Usage;
use POSIX;
my $man = '';
my $help= '';
my $query_seq = '';
my $hit_seq = '';
my $skip_find_best_hit='F';
my $gene_pair='';
my $gff = '';
my $spe = '';
my $blat_dir = 'blat';
my $blast2seq = 'bl2seq';
my $is_genome='F';
my $qSeq_sep=' ';
my $hSeq_sep=' ';
my $chrom = 'all';
#my $codeml = 'codeml';#updated,it can also be yn00 or KaKs_Calculator
my $program = 'codeml';#updated,it can also be yn00 or KaKs_Calculator
my $icode="0"; #updated, if the $codeml is codeml or yn00, then it can be 0,1~10; if it is KaKs_Calculator, it should be 1~23;
my @KaKsCms;#updated, it can be all the options that in the KaKs_Calculator
my $blat_option="";
my $bl2seq_option="";
my $fastacmd = '';
my $problem_loc = 'problem_locs.log';
my $detail ='DNA2AA_alignment.log';
my $blat_out = 'blat_rel.query';
my $blat_done = 'F';
my $kaks_file = 'kaks_rel.txt';
my $min_match_length = '100';
my $min_identity = '0.9';
my $version = '';

GetOptions('help|?' => \$help,
	      man => \$man,
	      'query_seq=s' => \$query_seq,
	      'hit_seq=s' => \$hit_seq,
	      'gff=s' =>\$gff,
	      'spe=i' => \$spe,
	      'min_length=i' => \$min_match_length,
	      'min_identity=f' => \$min_identity,
	      'is_genome=s' => \$is_genome,
	      'sfbh=s' => \$skip_find_best_hit,
	      'gene_pair=s' => \$gene_pair,
	      'qSeq_sep=s' => \$qSeq_sep,
	      'hSeq_sep=s' => \$hSeq_sep,
	      'blat=s' =>\$blat_dir,
	      'fastacmd=s' =>\$fastacmd,
	      'blast2seq=s'=> \$blast2seq,
	      'bl2seq_option=s'=>\$bl2seq_option,
   	      'blat_option=s'=>\$blat_option,
	      'chrom=s' => \$chrom,
	      'program=s' => \$program,
	      'icode=s' =>\$icode,
	      'KaKsCms=s'=>\@KaKsCms,
	      'problem_loc=s' =>  \$problem_loc,
	      'detail=s' => \$detail,
	      'blat_out=s' => \$blat_out,
	      'blat_done=s' => \$blat_done,
	      'kaks_file=s' => \$kaks_file,
	      'version'=> \$version
    )
    or pod2usage(2);
pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

if ($version){
    print "GKas Version 1.0. This version is the version when submit paper.\n";
    print "GKas Version 1.1. Add yn00 and KaKs_Calculator for -codeml option, add two more options=>icode and KaKsCms. Not released.\n";
    print "GKas Version 1.2. Add -blat_option and -bl2seq_option, and change option -codeml to -program. Current Version.\n";
	print "gKaKs Version 1.2.1. We renamed the pipe line, since the GKas sounds like 'jackass'; and \$kaks_file format are updated with title.\n";
    exit;
}

pod2usage(1) unless ($query_seq && $hit_seq && $gff && $spe && $chrom);

if ($bl2seq_option=~/-[i|j|p|o|a|T|I|J|D|A]/){
    print "options -i|-j|-p|-o|-a|-T|-I|-J|-D|-A can't be used in bl2seq in this pipe line\n";
    exit;
}

if ($blat_option=~/-[t|q|prot|noHead|fastMap|out]/){
    print "options -t|-q|-prot|-noHead|-fastMap|-out can't be used in blat in this pipe line\n";
    exit;
}


my $KaKsC_methods;

if (substr($program,-6,6) eq "codeml" or substr($program,-4,4) eq "yn00"){
    if ($icode!~/\d/ or $icode>10){
	print "the option icode must be in the range 0~10 for option program = codeml/yn00\n";
	print qq(
* Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,
* 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt.,
* 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt.,
* 10: blepharisma nu.
* These codes correspond to transl_table 1 to 11 of GENEBANK.
);
	exit;
    }
}else{
    $icode="1" if $icode eq "0";
    if ($icode!~/\d/ or $icode>23){
	print "the option icode must be in the range 1~23 for option program = KaKs_Calculator\n";
	print qq[
        -c      Genetic code table (Default = 1-Standard Code).
                  1-Standard Code                         2-Vertebrate Mitochondrial Code
                  3-Yeast Mitochondrial Code                      4-Mold Mitochondrial Code
                  5-Invertebrate Mitochondrial Code                       6-Ciliate, Dasycladacean and Hexamita Code
                  9-Echinoderm and Flatworm Mitochondrial Code                    10-Euplotid Nuclear Code
                  11-Bacterial and Plant Plastid Code                     12-Alternative Yeast Nuclear Code
                  13-Ascidian Mitochondrial Code                          14-Alternative Flatworm Mitochondrial Code
                  15-Blepharisma Nuclear Code                     16-Chlorophycean Mitochondrial Code
                  21-Trematode Mitochondrial Code                         22-Scenedesmus obliquus mitochondrial Code
                  23-Thraustochytrium Mitochondrial Code
                  (More information about the Genetic Codes: http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi?mode=c)

];
	exit;
    }
    $KaKsC_methods="";
    foreach (my $kaksi=0;$kaksi<=$#KaKsCms;$kaksi++){
	$KaKsC_methods.="-m $KaKsCms[$kaksi] ";
    }
    $KaKsC_methods="-m MA" if $KaKsC_methods eq "";
    print $KaKsC_methods,"\n\n";
}

if ($is_genome eq "T"){
    if (! $fastacmd){
	print "please make sure you have format the genome data by formatdb with option -oT and the index name is same with $hit_seq\n";
	pod2usage(1);
    }
}

if ($skip_find_best_hit eq "T"){
    if (! $gene_pair){
	print "you must give out the gene name pairs' file\nEXAMPLE:\nLOC_Os01g01090.1 glaberrima_3s-4902445-4905356\n Make sure the name is unique and same with the $query_seq and $hit_seq\n";
	pod2usage(1);
    }
}

my %multiple;
my %chimeric;
my %locs_seq1 = &read_fasta_file("$query_seq",$qSeq_sep);
my %locs_seq2;
if ($is_genome eq "F"){
    %locs_seq2 = &read_fasta_file("$hit_seq",$hSeq_sep);
}
my %codon=&codon_table;
my @both_stop;
my @hit_stop;
my @query_stop;
my %predict_cds_seq;
my %query_cds_seq;
my $jap_cds_id_info = get_gff($gff,$spe,$chrom);
my %jap_cds_id_info = %{$jap_cds_id_info};

my @reverse;
my @frame_shift;
open (ERRORS,">$problem_loc") or die "cat not open $problem_loc, please make sure you have the write permissions";
open (DET,">$detail") or die "cat not open $detail, please make sure you have the write permissions";
if (-e $kaks_file){
    my $old_kaks_file=$kaks_file."_back_".time();
    my $cmd="mv $kaks_file $old_kaks_file";
    system $cmd;
    print "your previous result was backup as $old_kaks_file\n";
    sleep (5);
}

if ($skip_find_best_hit eq "F"){
    if ($blat_done eq "F"){
	if (-e $blat_out){
	    my $old_blat_out=$blat_out."_back_".time();
	    my $cmd="mv $blat_out $old_blat_out";
	    system $cmd;
	    print "your previous blat results was backup as $old_blat_out\n";
	    sleep (5);
	}
	print "blat may cost a lot of time, please wait!\n";
	`$blat_dir $blat_option $hit_seq $query_seq $blat_out`;
    }
    my $new_blat_output = best_hit($blat_out);
    open (NEWKS,"$new_blat_output") or die "cat not open file $new_blat_output";
    <NEWKS>;
    <NEWKS>;
    <NEWKS>;
    <NEWKS>;
    <NEWKS>;
}else{
    open (NEWKS,"$gene_pair") or die "cat not open file $gene_pair";
}
if (substr($program,-4,4) eq "yn00"){
	open(OUTPUTt,">>$kaks_file");
	print OUTPUTt  "compare_id\tseq.\tseq.\tS\tN\tt\tkappa\tomega\tdN +- SE\tdS +- SE\n";
	close OUTPUTt;
}elsif(substr($program,-6,6) eq "codeml"){
}else{
	open(OUTPUTt,">>$kaks_file");
	print OUTPUTt  "Sequence	Method	Ka	Ks	Ka/Ks	P-Value(Fisher) Length  S-Sites N-Sites Fold-Sites(0:2:4)	Substitutions   S-Substitutions N-Substitutions Fold-S-Substitutions(0:2:4)     Fold-N-Substitutions(0:2:4)     Divergence-Time Substitution-Rate-Ratio(rTC:rAG:rTA:rCG:rTG:rCA/rCA)   GC(1:2:3)       ML-Score        AICc    Akaike-Weight   Model\n";
	close OUTPUTt;
}
while (<NEWKS>){
    my $model_id;
    my $model_id2;
    if ($skip_find_best_hit eq "F"){
	my @lines=split/\t/,$_;
	next if $lines[0]<$min_match_length;
	next if $lines[0]/($lines[0]+$lines[1])<$min_identity;
	next if not exists $jap_cds_id_info{$lines[9]};
	$model_id=$lines[9];
	$model_id2=$lines[13];
    }else{
	my $genepair=$_;
	chomp $genepair;
	while($genepair=~s/\t/ /){};
	while($genepair=~s/  / /){};
	my @gene_pair=split/ /,$genepair;
	$model_id=$gene_pair[0];
	$model_id2=$gene_pair[1];
    }
    my $flag="";
    my $two_ids=$model_id."_".$model_id2;
    foreach my $cds (sort {$a<=>$b} keys %{$jap_cds_id_info{$model_id}}){
	my $start_point=0 if $cds eq 1;
	my $tmp=$cds-1;
        if(($cds > 1) && (exists $jap_cds_id_info{$model_id}{$tmp})){
	    $start_point+=$jap_cds_id_info{$model_id}{$tmp} ;
        }elsif(($cds > 1) && (!(exists $jap_cds_id_info{$model_id}{$tmp}))){
	    $start_point = 0;
	    print "$two_ids\twarning at line 159\n";
	    exit;
        }
	my $tmp_cds_seq=substr($locs_seq1{$model_id},$start_point,$jap_cds_id_info{$model_id}{$cds});
	$flag=&get_aligment($tmp_cds_seq,$locs_seq2{$model_id2},$model_id,$model_id2,$cds,$blast2seq);
	last if $flag eq "both_rev";
	last if $flag eq "bl2seq_wrong"
    }
    goto ERR if $flag eq "both_rev";
    goto ERR if $flag eq "bl2seq_wrong";
    my $predict_gene_seq="";
    my $aligned_gene_seq="";
    foreach my $cds (sort {$a<=>$b} keys %{$predict_cds_seq{$two_ids}}){
	$predict_gene_seq.=$predict_cds_seq{$two_ids}{$cds};
	$aligned_gene_seq.=$query_cds_seq{$two_ids}{$cds};
    }
    my ($all_tmp1,$all_tmp2,$all_aa1,$all_aa2,$yn_tmp1,$yn_tmp2)=&calculate_yn00($two_ids,$predict_gene_seq,$aligned_gene_seq);
    my $alignment_len=length($locs_seq1{$model_id});
    my $alignment_len1=length($locs_seq2{$model_id2});
    my $alignment_len2=length($yn_tmp1);
    print DET ">$two_ids\n$all_tmp1\n$all_tmp2\n$all_aa1\n$all_aa2\n";
    if (substr($program,-6,6) eq "codeml"){
	open (YN,">codeml3") or die "wrong";
	print  YN "\t2\t$alignment_len2\nquery\n$yn_tmp1\nhit\n$yn_tmp2\n";
	close YN;
	open(CTL,">codeml.ctl");
	my $ctl_text = qq{
    seqfile = codeml3
    treefile = codeml.tree
    outfile = codeml_rel
    noisy = 0  * 0,1,2,3,9: how much rubbish on the screen
    verbose = 0 * 0: concise; 1: detailed, 2: too much
    runmode = -2
    seqtype = 1
    CodonFreq = 2 
    clock = 0
    model = 0
    NSsites = 0
    icode = $icode
    fix_kappa = 0
    kappa = 2
    fix_omega = 0
    omega = 0.1
    fix_alpha = 1
    alpha = .0
    Malpha = 0
    ncatG = 10
    getSE = 0
    RateAncestor = 0
    method = 0
    Small_Diff = .5e-6
};
	print CTL "$ctl_text";
	system "$program codeml.ctl";
	unlink("2ML.dN");
	unlink("2ML.dS");
	unlink("2ML.t");
	unlink("2NG.dN");
	unlink("2NG.dS");
	unlink("2NG.t");
	system qq( tail -n 1 codeml_rel > temp_codeml_output); #bug report by Jose Guillen, previous version is:  system qq( tail codeml_rel -n 1  > temp_codeml_output);

#############################################################################
#Otherwise the bash gives:
#tail: -n: No such file or directory
#tail: 1: No such file or directory
#CODONML in paml version 4.7, January 2013

#and kaks_rel.txt looks like this:
#A0weG_0401_contig108_orf00001    ==> codeml_rel <==
#############################################################################	

	open(INPUTt,"temp_codeml_output");
	my $line_codeml = <INPUTt>;
	close INPUTt;
	chomp $line_codeml;
	open(OUTPUTt,">>$kaks_file");
	print OUTPUTt  "$two_ids\t$line_codeml\n"; 
	close OUTPUTt;
    }elsif(substr($program,-4,4) eq "yn00") {
	open (YN,">yn00") or die "wrong";
	print  YN "\t2\t$alignment_len2\nquery\n$yn_tmp1\nhit\n$yn_tmp2\n";
	close YN;    
	open(CTL,">yn00.ctl");
	my $ctl_text = qq{
seqfile      = yn00
outfile      = yn00_rel
verbose      = 1  * 1: detailed output (list sequences), 0: concise output
icode        = $icode  * 0:universal code; 1:mammalian mt; 2-10:see below
weighting    = 0  * weighting pathways between codons (0/1)
commonf3x4   = 0  * use one set of codon freqs for all pairs (0/1)
};
	print CTL "$ctl_text";

	system "$program";
	system qq( grep "kappa" yn00_rel -A 2 |tail -n1 > temp_codeml_output);
	open(INPUTt,"temp_codeml_output");
	my $line_codeml = <INPUTt>;
	close INPUTt;
	chomp $line_codeml;
	open(OUTPUTt,">>$kaks_file");
	print OUTPUTt  "$two_ids\t$line_codeml\n";
	close OUTPUTt;

    }else{
	open (YN,">kaks.axt") or die "wrong";
	print  YN "$two_ids\n$yn_tmp1\n$yn_tmp2\n";
	close YN;    
	system qq($program -i kaks.axt -o kaks.axt.kaks -c $icode $KaKsC_methods);
	system qq(cat kaks.axt.kaks |grep -v Sequence >> $kaks_file);
    }
    ERR:
}
print ERRORS "\n===reverse===\n".join("\n",@reverse);
print ERRORS "\n===shift===\n".join("\n",@frame_shift);
print ERRORS "\n===both-stop===\n".join("\n",@both_stop);
print ERRORS "\n===hit-stop===\n".join("\n",@hit_stop);
print ERRORS "\n===query-stop===\n".join("\n",@query_stop);
close ERRORS;
close DET;
close NEWKS;



sub codon_table
{
    my %codon;
    $codon{"TTT"}="F";
    $codon{"TTC"}="F";
    $codon{"TTA"}="L";
    $codon{"TTG"}="L";
    $codon{"CTA"}="L";
    $codon{"CTG"}="L";
    $codon{"CTT"}="L";
    $codon{"CTC"}="L";
    $codon{"ATA"}="I";
    $codon{"ATC"}="I";
    $codon{"ATT"}="I";
    $codon{"ATG"}="M";
    $codon{"GTA"}="V";
    $codon{"GTC"}="V";
    $codon{"GTT"}="V";
    $codon{"GTG"}="V";
    $codon{"TCA"}="S";
    $codon{"TCT"}="S";
    $codon{"TCG"}="S";
    $codon{"TCC"}="S";
    $codon{"CCT"}="P";
    $codon{"CCA"}="P";
    $codon{"CCC"}="P";
    $codon{"CCG"}="P";
    $codon{"ACA"}="T";
    $codon{"ACT"}="T";
    $codon{"ACG"}="T";
    $codon{"ACC"}="T";
    $codon{"GCT"}="A";
    $codon{"GCA"}="A";
    $codon{"GCC"}="A";
    $codon{"GCG"}="A";
    $codon{"TAT"}="Y";
    $codon{"TAC"}="Y";
    $codon{"TAA"}="*";
    $codon{"TAG"}="*";
    $codon{"CAT"}="H";
    $codon{"CAC"}="H";
    $codon{"CAA"}="Q";
    $codon{"CAG"}="Q";
    $codon{"AAT"}="N";
    $codon{"AAC"}="N";
    $codon{"AAG"}="K";
    $codon{"AAA"}="K";
    $codon{"GAT"}="D";
    $codon{"GAC"}="D";
    $codon{"GAG"}="E";
    $codon{"GAA"}="E";
    $codon{"TGG"}="W";
    $codon{"TGA"}="*";
    $codon{"TGC"}="C";
    $codon{"TGT"}="C";
    $codon{"CGT"}="R";
    $codon{"CGC"}="R";
    $codon{"CGG"}="R";
    $codon{"CGA"}="R";
    $codon{"AGA"}="R";
    $codon{"AGG"}="R";
    $codon{"AGT"}="S";
    $codon{"AGC"}="S";
    $codon{"GGG"}="G";
    $codon{"GGC"}="G";
    $codon{"GGA"}="G";
    $codon{"GGT"}="G";
    return %codon;
}

sub read_fasta_file()
{
    my %hash;
    my ($file_name,$sep)=@_;
    my $aa;
    my $id;
    my @lines;
    my $pep;
    open (PEP,"$file_name") or die "can not open $file_name, please check the path.";
    my $first_line=<PEP>;
    chomp $first_line;
    if($first_line =~ s/>//){
	while($first_line=~s/\t/ /){}
	@lines = split(/$sep/,$first_line);
    }
    $id = $lines[0];
    while ($pep=<PEP>){
	chomp $pep;
	if ($pep=~s/>//){
	    $aa=~s/\*//;
	    $hash{$id}=$aa;
	    $aa="";    
	    while($pep=~s/\t/ /){}
	    @lines = split(/$sep/,$pep);
	    $id=$lines[0];
	}else{
	    $aa.=$pep;
	}
    }
    $aa=~s/\*//;
    $hash{$id}=$aa;
    return %hash;
}


sub calculate_yn00(){
    my ($two_ids,$predict_seq,$aligned_seq)=@_;
    $predict_seq=~ tr/[a-z]/[A-Z]/;
    $aligned_seq=~ tr/[a-z]/[A-Z]/;
    my $yn_length=length($predict_seq);
    my $yn2_length=length($aligned_seq);
    if ($yn_length ne $yn2_length){
	print ERRORS "$two_ids\t===length-not-equal===\n";
	print ERRORS "$predict_seq\n$aligned_seq\n";
	print ERRORS "the alignment length is not equal to each other, please check\n";
    }
    my $new_predict_seq;
    my $new_aligned_seq;
    for (my $i=0;$i<$yn2_length;$i++){
	my $a=substr($aligned_seq,$i,1);
	my $b=substr($predict_seq,$i,1);
	if ($a eq "+"){
	    my $new_a=substr($aligned_seq,$i,3);
	    $a="";
	    $b="";
	    if ($new_a eq "+++"){
		$i=$i+2;
	    }else{
		push (@frame_shift,$two_ids);
	    }
	}
	$new_predict_seq.=$b;
	$new_aligned_seq.=$a;
    }
    if (substr($new_aligned_seq,-3)!~/-/){
	if ($codon{substr($new_aligned_seq,-3)} eq "*"){
	    $new_predict_seq=substr($new_predict_seq,0,-3);
	    $new_aligned_seq=substr($new_aligned_seq,0,-3);
	}
    }
    $yn_length=length($new_predict_seq);
    print ERRORS "$two_ids can't devided by 3, wrong length\n" if ($yn_length/3 ne int ($yn_length/3));
    my $yn_tmp1;
    my $yn_tmp2;
    my $all_tmp1;
    my $all_tmp2;
    my $all_aa1;
    my $all_aa2;
    for (my $i=0;$i< $yn_length;$i=$i+3){
	my $tmp1=substr($new_predict_seq,$i,3);
	my $tmp2=substr($new_aligned_seq,$i,3);
	last if length($tmp1) ne "3" and length($tmp2) ne "3";
	if (($tmp1 =~/\-|\~|n|N/ and $tmp1!~/\+/) or $tmp2 =~/\-|\~|N/ ){
    next;
	}elsif($tmp1 =~/\+/){
	    my $shift_a=0;
	    for (my $j=0;$j<3;$j++){
		my $shift=substr($tmp1,$j,1);
		$shift_a++ if $shift eq "+";
	    }
	    if ($shift_a eq 3){
		next;
	    }else{
		push (@frame_shift,$two_ids);
		next;
	    }
	}
	if ($codon{$tmp1} eq "*" or $codon{$tmp2} eq "*"){
	    if ($codon{$tmp1} eq "*" and $codon{$tmp2} eq "*"){
		push (@both_stop,$two_ids);
	    }elsif($codon{$tmp1} eq "*" and $codon{$tmp2} ne "*"){
		push (@hit_stop,$two_ids);
	    }elsif($codon{$tmp1} ne "*" and $codon{$tmp2} eq "*"){
		push (@query_stop,$two_ids);
	    }
	}
	if ($codon{$tmp1} ne "*" and $codon{$tmp2} ne "*"){
	    $yn_tmp1.=$tmp1;
	    $yn_tmp2.=$tmp2;
	}
	$all_tmp1.=$tmp1;
	$all_tmp2.=$tmp2;
	$all_aa1.="-".$codon{$tmp1}."-";
	$all_aa2.="-".$codon{$tmp2}."-";
    }
    return ($all_tmp1,$all_tmp2,$all_aa1,$all_aa2,$yn_tmp1,$yn_tmp2);
}


sub get_aligment()
{
    my ($dna_seq1,$dna_seq2,$mid1,$mid2,$cds,$blast2seq)=@_;
    open (TMP,">tmpdna1.fa") or die "wrong";
    print TMP ">$mid1"."_query\n$dna_seq1\n";
    close TMP;
    open (TMP,">tmpdna2.fa") or die "wrong";
    print TMP ">$mid2"."_hit\n$dna_seq2\n";
    close TMP;
    system "rm tmp.fas.txt" if (-e "tmp.fas.txt");
    system "$blast2seq $bl2seq_option -p blastn -i tmpdna1.fa -j tmpdna2.fa -D 1 -o tmp.fas.txt";
    my $split_hit=0;
    my $two_ids=$mid1."_".$mid2;
    my @multiple_identity;
    undef @multiple_identity;    
    my $max_identity=0;
    my $max_hit=0;
    my $max_match_length=0;
    my $max_len_hit=0;
    my $max_match_bit=0;
    my $max_len_bit=0;
    my $max_bit=0;
    my $newrev_flag=0;
    NEWREV:
    if(-e "tmp.fas.txt"){
	open (TMP,"tmp.fas.txt") or die "wrong";
    }else{
	print ERRORS "$mid1\t$mid2\tcan't open file tmp.fas.txt\n";
	return("bl2seq_wrong");
    }
    while(<TMP>){
	if ($_=~/_hit/){
	        my @lines=split/\t/,$_;
		$split_hit++ if $lines[10] < 4e-5;
		next if $lines[10] > 4e-5;
		$multiple_identity[$split_hit-1]=$lines[2] if $split_hit>0;
		if ($split_hit eq 1){
		    $max_identity =$lines[2];
		    $max_match_length=$lines[3];
		    $max_hit=$split_hit;
		    $max_len_hit=$split_hit;
		    $max_match_bit=$lines[11];
		    $max_len_bit=$lines[11];
		}else{
		    if  ($max_identity <$lines[2]){
			$max_identity =$lines[2];
			$max_hit=$split_hit;
			$max_match_bit=$lines[11];
		    }
		    if ($max_match_length<$lines[3]){
			$max_match_length=$lines[3];
			$max_len_hit=$split_hit;
			$max_len_bit=$lines[11];
		    }
		}
	}
    }
    my $follow="F";
    if ($split_hit>1){
	print ERRORS "more than one hit\t$mid1\t$mid2\t$cds\t$split_hit\n";
	foreach my $tmp_line2 (@multiple_identity){
	    if(abs($max_identity - $tmp_line2)>5){
		$chimeric{$two_ids}{$cds}++;
	    }
	}
	if (abs($max_match_bit-$max_len_bit)<50){
	    if (exists $chimeric{$two_ids}{$cds}){
		$max_bit=$max_match_bit;
	    }else{
		$max_bit= $max_len_bit;
	    }
	}else{
	    $max_bit=$max_match_bit;
	    if ($max_match_bit<$max_len_bit){
		$max_bit=$max_len_bit;
	    }
	    $follow="T";
	}
    }

    close TMP; 
    if(-e "tmp.fas.txt"){
	open (TMP,"tmp.fas.txt") or die "wrong";
    }else{
	print ERRORS "$mid1\t$mid2\tcan't open file tmp.fas.txt\n";
	return("bl2seq_wrong");
    }

    my $new_tmp_split=0;
    while(<TMP>){
	if ($_=~/_hit/){
	    my @lines=split/\t/,$_;
	    $new_tmp_split++ if $lines[10] < 4e-5;
	    next if $lines[10] > 4e-5;
	    if ($follow eq "T"){
		next if $lines[11] ne "$max_bit";
	    }else{
		next if (exists $chimeric{$two_ids}{$cds} and $new_tmp_split ne $max_hit);
		next if (not exists $chimeric{$two_ids}{$cds} and $new_tmp_split ne $max_len_hit);
	    }
	    if ($lines[6]<$lines[7] and $lines[9]<$lines[8]){
		my $rcvalue = reverse $locs_seq2{$mid2};
		$rcvalue =~ tr/ACGTacgt/TGCAtgca/;  
		open (TMP,">tmpdna2.fa") or die "wrong";
		print TMP ">$mid2"."_hit\n$rcvalue\n";
		close TMP;
		system "$blast2seq $bl2seq_option -p blastn -i tmpdna1.fa -j tmpdna2.fa -D 1 -o tmp.fas.txt";
		$newrev_flag++;
		if ($newrev_flag eq "1"){
		    goto NEWREV;
		}else{
		    print ERRORS "both direction have hit\t$mid1\t$mid2\t$cds\t$split_hit\n";
		    return "both_rev";
		}
	    }
	}
    }
    close TMP;

    if ($split_hit eq "0"){
	print ERRORS "no good hit\t$mid1\t$mid2\t$cds\t$split_hit\n";
	for (my $i=1;$i<=length($dna_seq1);$i++){  
	    $predict_cds_seq{$two_ids}{$cds}.="-";
	    $query_cds_seq{$two_ids}{$cds}.="-";
	}
	return "not_hit";
    }
    my $rev_flag="F";
    REV:
    my $new_split=0;
    if(-e "tmp.fas.txt"){
	open (TMP,"tmp.fas.txt") or die "wrong";
    }else{
	print ERRORS "$mid1\t$mid2\tcan't open file tmp.fas.txt\n";
	return("bl2seq_wrong");
    }

    while(<TMP>){
	if ($_=~/_hit/){
	    my @lines=split/\t/,$_;
	    $new_split++ if $lines[10] < 4e-5;
	    next if $lines[10] > 4e-5;
	    next if $new_split>1;
	    $max_bit=$lines[11];
	    if ($lines[6]<$lines[7] and $lines[9]<$lines[8]){
		my $rcvalue = reverse $locs_seq2{$mid2};
		$rcvalue =~ tr/ACGTacgt/TGCAtgca/;  
		open (TMP,">tmpdna2.fa") or die "wrong";
		print TMP ">$mid2"."_hit\n$rcvalue\n";
		close TMP;
		system "$blast2seq $bl2seq_option -p blastn -i tmpdna1.fa -j tmpdna2.fa -D 1 -o tmp.fas.txt";
		$locs_seq2{$mid2}=$rcvalue;
		push (@reverse,$two_ids);
		if ($rev_flag eq "F"){
		    $rev_flag = "T";
		    goto REV;
		}else{
		    print ERRORS "both direction have hit\t$mid1\t$mid2\t$cds\t$split_hit\n";
		    return "both_rev";
		}
	    }
	    if ($lines[5] >0){
		my $read_line=int(5*((($lines[7]-$lines[6])/60)+1)+11)+1;
		my $read_line2=int(5*((($lines[9]-$lines[8])/60)+1)+11)+1;
		$read_line=$read_line2 if $read_line<$read_line2;
		&have_gap($dna_seq1,$dna_seq2,$read_line,$cds,$two_ids,$lines[6],$max_bit,$blast2seq);
		return "gap";
	    }
	    open(INPUTtemp,"tmpdna2.fa");
	    <INPUTtemp>;
	    my $temp_seq2 = <INPUTtemp>;
	    my $value=substr($temp_seq2,($lines[8]-1),($lines[9]-$lines[8]+1));
	    my $temp_start = $lines[8]-1;
	    my $temp_end = $lines[9] - $lines[8] +1;
	    for (my $i=1;$i<$lines[6];$i++){
		$value="-".$value;
	    }
	    while (length($value)<length($dna_seq1)){
		$value.="-";
	    }
	    $predict_cds_seq{$two_ids}{$cds}=$value;
	    $query_cds_seq{$two_ids}{$cds}=$dna_seq1;
	}
    }
    return "ok";
}


sub have_gap()
{
    my ($dna_seq1,$dna_seq2,$read_line,$cds,$two_ids,$dash,$max_bit,$blast2seq)=@_;
    chomp $max_bit;
    while($max_bit=~s/ //){}
    system "$blast2seq $bl2seq_option -p blastn -i tmpdna1.fa -j tmpdna2.fa -D 0 -o tmp.fas.gap.txt";
    open (TMP,"tmp.fas.gap.txt") or die "wrong";
    my $count_line=0;
    my $query_seq;
    my $hit_seq;
    my $al_flag="F";
    while(<TMP>){
	if ($max_bit ne "0"){
	    if ($_=~/Score/){
		my @arr=split/=/,$_;
		while ($arr[1]=~s/  / /){}
		my @arrb=split/ /,$arr[1];
		if ($arrb[1] eq $max_bit){
		    $read_line-=8;
		    $al_flag="T";
		}else{
		    next;
		}
	    }
	    next if $al_flag eq "F";
	}
	$count_line++;
	last if $count_line>$read_line;
	while($_=~s/  / /){}
	my @lines=split/ /;
	if ($_=~/Sbjct/){
	    $hit_seq.=$lines[2];
	}
	if ($_=~/Query:/){
	    $query_seq.=$lines[2];
	}
    }
    my $tmp_seq=$query_seq;
    my $tmp_dash=0;
    while ($tmp_seq=~s/-//){
	$tmp_dash++;
    }
    $hit_seq=~s/-/+/g;
    $query_seq=~s/-/+/g;
    for (my $i=1;$i<$dash;$i++){
	$hit_seq="-".$hit_seq;
	$query_seq="-".$query_seq;
    }
    while(length($query_seq)<(length($dna_seq1)+$tmp_dash)){
	$hit_seq.="-";
	$query_seq.="-";
    }
    $predict_cds_seq{$two_ids}{$cds}=$hit_seq;
    $query_cds_seq{$two_ids}{$cds}=$query_seq;
}

sub get_gff{
    my ($gff_file,$spe,$chrom)=(@_);
    my %jap_cds_id_info;
    open(GFF,$gff_file) or die "can't find the file $gff_file";
    if (($spe ==1) || ($spe==2) || ($spe ==3)){
	my %tmp_cds_hash;
	while (<GFF>){
	    next if $_!~/CDS/;
my @lines=split/\s+/,$_;
	    if($chrom eq "all"){
		goto line;
	    }else{
		if($lines[0] eq $chrom){
		  line:my $tmp_id=substr($lines[11],1,(length($lines[11])-3));
		    $tmp_cds_hash{$tmp_id}++;
		    my $cds_id = $tmp_cds_hash{$tmp_id};
		    $jap_cds_id_info{$tmp_id}{$cds_id}=abs($lines[3]-$lines[4])+1;
		}
	    }
	}
    }elsif($spe==4){
	while (<GFF>){
	    next if $_!~/CDS/;
	    if($chrom eq "all"){
		goto line1;
	    }else{
		next if $_!~/$chrom/;
	    }
	  line1:my @lines=split/\t/,$_;
	    my @tmp_id=split/\:/,$lines[8];
	    my @model_id=split/=/,$tmp_id[0];
	    my @cds_id=split/;/,$tmp_id[1];
	    $cds_id[0]=~s/cds_//;
	    $jap_cds_id_info{$model_id[1]}{$cds_id[0]}=abs($lines[3]-$lines[4])+1;
	}
    }elsif($spe==5){
	my %tmp_cds_hash;
	while (<GFF>){
	    next if $_!~/CDS/;
	    if($chrom eq "all"){
		goto line2;
	    }else{
		next if $_!~/$chrom/;
	    }
	  line2:my @lines = split(/\s+/,$_);
	    my @tmp_id = split(/[=|\,|;]/,$lines[8]);
	    my $model_id = $tmp_id[1];
	    $tmp_cds_hash{$model_id}++;
	    my $cds_id= $tmp_cds_hash{$model_id};
	    $jap_cds_id_info{$model_id}{$cds_id}=abs($lines[3]-$lines[4])+1;
	}
    }elsif($spe==6){
	my %tmp_cds_hash;
	while (<GFF>){
	    next if $_!~/CDS/i;
	    if($chrom eq "all"){
		goto line3;
	    }else{
		next if $_!~/$chrom/;
}
#Here is the Danaus line that we need to process
#DPSCF300466	OGS2	CDS	40587	41969	.	+	0	ID=DPOGS200002-TA.cds-8;Parent=DPOGS200002-TA;
	  line3:my @lines = split(/\s+/,$_);
	    my @tmp_id = split(/[=|;]/,$lines[8]);  # so split the 9th column on = and ; 
	    my $model_id = $tmp_id[3];
	    
	    #print "my ID: $model_id\tGFF: $_"; 
#	    my $model_id = $lines[8];
	    $tmp_cds_hash{$model_id}++;
	    my $cds_id = $tmp_cds_hash{$model_id};
	    $jap_cds_id_info{$model_id}{$cds_id}=abs($lines[3]-$lines[4])+1;
	}
    }

    close GFF;
    return(\%jap_cds_id_info);
}
sub best_hit{
    my $blat_out = shift(@_);
    my $new_blat_out = "new_blatout";
    open(OUTPUTb,">$new_blat_out") or die "cat not open $new_blat_out, please make sure you have the write permissions\n";
    open(INPUTb,$blat_out) or die "cat not open $blat_out, please make sure the blat result is there\n";
    my %hash;
    for(my $i=0; $i<=4; $i++){
	my $temp_line = <INPUTb>;
	chomp $temp_line;
	print  OUTPUTb "$temp_line\n";
    }
    my %count;
    while(my $line = <INPUTb>){
	chomp $line;
	my @info = split(/\s+/,$line);
	my ($qcover,$tcover,$iden) = pslCal($line);
	my $iden_score = 100.0 - $iden * 0.1;
	my $blat_score = blat_score($line);
	if(!(exists $count{$info[9]})){
	    $count{$info[9]}=0;
	}
	$hash{$info[9]}->{$count{$info[9]}}->{iden_score} = $iden_score;
	$hash{$info[9]}->{$count{$info[9]}}->{score} = $blat_score;
	$hash{$info[9]}->{$count{$info[9]}}->{qcover} = $qcover;
	$hash{$info[9]}->{$count{$info[9]}}->{tcover} = $tcover;
	$hash{$info[9]}->{$count{$info[9]}}->{line} = $line;
	$count{$info[9]}++;
    }
    my %count_key;
    for my $key (sort keys %hash){
	my %sub_hash;
	undef %sub_hash;
	%sub_hash = %{$hash{$key}};
	for my $sub_key (sort {$sub_hash{$b}->{score} <=> $sub_hash{$a}->{score} || $sub_hash{$b}->{iden_score} <=> $sub_hash{$a}->{iden_score} } keys %sub_hash){
	    if(!(exists $count_key{$key})){
		if ($is_genome eq "F"){
		    print  OUTPUTb "$sub_hash{$sub_key}->{line}\n";
		}else{
		    my $tmp_line=$sub_hash{$sub_key}->{line};
		    my @tmp_lines=split/\t/,$tmp_line;
		    my $cmd="$fastacmd -d $hit_seq -s $tmp_lines[13] -L $tmp_lines[15],$tmp_lines[16] >tmp_fastacmd.fasta";
		    system $cmd;
		    $tmp_lines[13]=$tmp_lines[13]."-".$tmp_lines[15]."-".$tmp_lines[16];
		    open (FASTACMD,"tmp_fastacmd.fasta") or die "can not open file tmp_fastacmd.fasta, wrong\n";
		    <FASTACMD>;
		    my $tmp_fastacmd_seq="";
		    while (my $tmp_seq=<FASTACMD>){
			chomp $tmp_seq;
			$tmp_fastacmd_seq.=$tmp_seq;
		    }
		    close FASTACMD;
		    $locs_seq2{$tmp_lines[13]}=$tmp_fastacmd_seq;
		    print OUTPUTb join("\t",@tmp_lines),"\n";
		}
		$count_key{$key} =1;
	    }
	}
    }
    close OUTPUTb;
    return("new_blatout");
}

sub blat_score{
    my $psl = shift(@_);
    my @info = split(/\s+/,$psl);
#    my $score = ($info[0] + $info[2]/2) - $info[1] - $info[5] - $info[7];
    my $score = ($info[0] + $info[2]/2) - $info[1] - $info[4] - $info[6];

    return($score);
}

sub pslCal{
    my $psl = shift(@_);
    my $sizeMul = 1;
    my ($qAliSize,$tAliSize,$aliSize);
    my $milliBad=0;
    my $sizeDif;
    my $insertFactor; 
    my $total;
    my @info = split(/\s+/,$psl);
    $qAliSize = $sizeMul * ($info[12] - $info[11]);
    $tAliSize = $info[16] - $info[15];
    my $qcover = $qAliSize/$info[10];
    my $tcover = $tAliSize/$info[14];
    if($qAliSize <= $tAliSize){
	$aliSize = $qAliSize;
    }else{
	$aliSize = $tAliSize;
    }

    if($aliSize <= 0){
	goto line1;
    }
    $sizeDif = $qAliSize - $tAliSize;
    if($sizeDif < 0){
#	$sizeDif = -$sizeDif;
	$sizeDif = 0;

    }
#    $insertFactor = $info[5];
#    $insertFactor += $info[7];
    $insertFactor = $info[4];
        $total = ($sizeMul * ($info[0] + $info[1]+$info[2]));
    if($total !=0){
	$milliBad = (1000*($info[1]*$sizeMul+$insertFactor + ceil(3*log(1+$sizeDif))))/$total;
      line1: return ($qcover,$tcover,$milliBad);
    }
}

__END__

=head1 NAME

gKaKs.pl - align DNA seq to compute ka/ks at genome level

=head1 SYNOPSIS

gKaKs.pl [options]

 Options:
    -help Help command

    -man Explanation command

    -query_seq *CDS sequences of reference genome in fasta format

    -gff *gff3 file of reference genome

    -hit_seq *Target genome sequences in fasta format (can be chromosomes, scaffolds, contigs and short genomic sequences)

    -spe *Reference species, the value ranges from 1 to 6. 1:human 2:mouse 3:fruitfly 4:rice 5:arabidopsis 6:others 

    -chrom *The chromosome or contig ID of the reference CDS sequences (-query_seq value) should be provided, which should be corresponding to identifier in gff file. If the reference CDS sequences come from several or all the chromosomes, the option should be set to "all". The default value is "all".

    -qSeq_sep The separation symbol of gene id in fasta file in query CDS file. The default value is gap. 

    -hSeq_sepThe separation symbol of sequence id in fasta file in target sequence file. The default value is gap.

    -sfbhSkip finding the best hits using blat. The default value is F. If the value is T, it indicates that the best match between ref CDS and target sequence has been generated and don’t need blat to search for the best match. And if it is set to T, the parameter for -gene_pair option should be provided.

    -gene_pairWhen -sfbh has been set as T, it is necessary to give this parameter. And this parameter is a file listing the paired CDS ID and target sequence ID. 

    -blat_doneThe default is F. When the value is T, it will skip blat step and it should guarantee the existence of blat output file (the parameter of -blat_out). 

    -is_genomeIndicate whether the target sequences are chromosome-wise whole genome sequences. The default value is F (the target sequence file includes multiple sequences, which can be small contigs or just homolog genomic sequences). If the value is T, the whole genome sequences should be formatted by "formatdb" and use option -o T. when the value is T, the directory of fastacmd should be given. 

    -fastacmdThe directory of fastacmd. It is needed, only when -is_genome is set as T.

    -blatThe directory of blat

    -blast2seqThe directory of bl2seq

    -program The command for Ka/Ks calculate program, it can be codeml/yn00/KaKs_Calculator, default is codeml
    
    -icode Genetic code type, for codeml and yn00, it can be 0~10; for KaKs_Calculator, it can be 1~23, the default is -Standard Code

    -KaKsCms the option for KaKs_Calculator, it can be used for more than onec, for example: -KaKsCms=NG -KaKsCms=LWL, their're 11 options: NG, LWL, LPB, MLWL, MLPB, GY, YN, MYN, MS, MA, ALL, the default is MA, but this will cost a lot of time to calculate

    -blat_option It can add options to the blat, use "" in the command, for example: -blat_option="-minMatch=1 -minIdentity=30", some of options are forbid in the pipe line, like "-t|-q|-prot|-noHead|-fastMap|-out".

    -bl2seq_option It can add options to the bl2seq, use "" in the command, for example: -bl2seq_option="-e=1 -mT", some of options are forbid in the pipe line, like "-i|-j|-p|-o|-a|-T|-I|-J|-D|-A".
    
    -problem_loc The file records the CDS ID with problem

    -detailThe file records the DNA and protein sequence alignment used to compute Ka/Ks 

    -blat_outThe file contains the output of blat program

    -kaks_fileThe file records the final kaks ratio for each homolog sequence pair

    -min_lengthThe minium match length, the default value is 100bp. 

    -min_identity The minium identity value, should be <=1, and is suggested to >= 0.7. The default is 0.9. 

    -version Version changes


=head1 DESCRIPTION

B<gKaKs> 

The main purpose of this program is to align the CDSs from a well-annotated genome to a target genome, which hasn’t been annotated yet. It uses blat to find the best match between the CDSs and target genome sequences, and then uses bl2seq to align every CDS to the blat-identified target genome region. After merging the aligned sequences and removing gaps according to reference CDS codon, it uses codeml of PAML to compute the ka/ks ratio between CDSs and their homolog sequences in the target genome. Also this program can compute Ka/Ks ratio for two lists of homologous DNA sequences, with one list as CDS sequences and the other list as genomic sequences.

This program can be used to compute ka/ks ratio between the genes in one well-annotated genome and their ortholog sequences in another closely related genome, which hasn’t been annotated. The result a) can be used to compute the diverge time between two species through estimating average Ks and mutation rate; b) can be used to estimate how many ortholog sequence pairs are under functional constraints; c) can be used as evidence to annotate genes.

Note

1. This program works for closely related species. If species splits from each other too long ago, the sequences diverge strikingly, bl2seq cannot generate good results.

2. Because this program deletes all the gaps that can’t be aligned in codon level, the Ks generated by this approach tend to be smaller.

3. In –problem_loc records, the Ks/Ka cannot be computed for the CDS id under ===reverse=== category. 

=head1 AUTHOR 

Chengjun Zhang @ University of Chicago
JunWang@ Wayne state university

=head1 BUGS

none.

=head1 COPYRIGHT 

This program is free software. You may copy or redistribute it under the same terms as Perl itself.

=cut

################################################

