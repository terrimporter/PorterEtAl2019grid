#!/usr/bin/perl
# Teresita M. Porter, Oct. 9, 2019
# Script to fix badly formatted sample names from Lisa's 2016 Grid Experiment sequence files.  Use mapping file based on 2016 metadata spreadsheet.
# USAGE perl reviseGRDInames.plx LV2016_1.txt GRDIname.map

use strict;
use warnings;

#declare var
my $i=0;
my $line;
my $key;#old name
my $value;#new correctly formatted name
my $sampleField;
my $sampleName;
my $new;
my $layer;
my $updated;
my $site;
my $date;
my $exptPart;
my $cores;
my $extractions;
my $gridSampleLayer;
my $sample;
my $gridCoord;
my $scalar;
my $marker_OTU;
my $marker;
my $OTU;

#declare array
my @report;
my @map;
my @line;
my @sampleField;
my @new;
my @marker_OTU;

#declare hash
my %map;#key=oldname, value=new GOOD name

open (IN1, "<", $ARGV[0]) || die "Error cannot open first infile: $!\n";
@report = <IN1>;
close IN1;

open (IN2, "<", $ARGV[1]) || die "Error canot open second infile: $!\n";
@map = <IN2>;
close IN2;

open (OUT, ">>", "LV2016_2.csv") || die "Error cannot open outfile: $!\n";
print OUT "marker,OTU,marker_OTU,GRDIname,site,date,cores,extractions,gridcoord,sample,layer,reads,strand,cellularOrganisms,colabel,coBP,superkingdom,sklabel,skBP,kingdom,klabel,kBP,phylum,plabel,pBP,class,clabel,cBP,order,olabel,oBP,family,flabel,fBP,genus,glabel,gBP,species,slabel,sBP\n";

#put mapping file in a hash
while($map[$i]) {
	$line = $map[$i];
	chomp $line;

	@line = split(/\t/,$line);
	$key = $line[0];
	$value = $line[1];
	$map{$key} = $value;
	$i++;
}
$i=0;

#remap GRDI sample names
while ($report[$i]) {
	$line = $report[$i];
	chomp $line;

	@line = split(/\t/,$line);
	$marker_OTU = $line[0];
	@marker_OTU = split(/_/,$marker_OTU);
	$marker = shift @marker_OTU;
	$OTU = shift @marker_OTU;
	$sampleField = $line[1];
	@sampleField = split(/_/,$sampleField);
	pop @sampleField;
	pop @sampleField;
	$sampleName = join ("_",@sampleField);

	if (exists $map{$sampleName}) {
		$new = $map{$sampleName};
		$line[1] = $new;
		@new = split("_",$new);
		$site = $new[0];
		$date = $new[1];
		$exptPart = $new[2];
		if ($exptPart =~ /\d{1,2}C(1|3)E/) {
			$exptPart =~ /(\d{1,2})C(1|3)E/;
			$cores = $1;
			$extractions = $2;
		}
		else {
			print "Error parsing exptPart\n";
			print $exptPart."\n";
		}
		$gridSampleLayer = $new[3];
		if ($gridSampleLayer =~ /S\d{1}(B|M|O)/) {
			$gridSampleLayer =~ /S(\d{1})(B|M|O)/;
			$sample = $1;
			$layer = $2;
			$gridCoord = "-";
		}
		elsif ($gridSampleLayer =~ /\d{2}(B|M|O)/) {
			$gridSampleLayer =~ /(\d{2})(B|M|O)/;
			$gridCoord = $1;
			$layer = $2;
			$sample = "-";
		}
		else {
			print "Error parsing gridSampleLayer\n";
			print $gridSampleLayer."\n";
		}
			
		splice @line, 2, 0, $site, $date, $cores, $extractions, $gridCoord, $sample, $layer;
#		$scalar = scalar @line;
#		print $cores."\n"; #test
		$updated = join(",",@line);
		print OUT $marker.",".$OTU.",".$updated."\n";
	}
	else {
		print "Problem finding $sampleName in hash\n";
	}
	$i++;
}
$site=();
$date=();
$cores=();
$extractions=();
$sample=();
$gridCoord=();
$layer=();
$i=0;
close OUT;
