#! /usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use Bio::SeqIO;
use File::Path qw(make_path);
use File::Util;
use File::Basename qw(basename);

##Program Info
# Name: convert_gbk.pl
# Based on scripts extractUpStreamDNA.pl, gbk2faaAV.pl, gbk2ffn.pl. and GBKtoPTT.pl
# Function: Takes a Genbank file as input. Parses for every CDS.Optionally outputs a .ptt format conversion, a .ffn (fasta format nucleotide),
#           a .faa (fasta format amino acid), or a .xls feature table. Can also produe an .ffn of 
#           a given upstream length + 3 to include the start codon. This only works for prokaryotic sequences
#           because it does not handle Splits or Joins found in eukaryotes. This only works for linearized genomes.
# Author: Max Hargreaves
# Copyright (c) University of Guelph, 2017
#
# Licence: This script may be used freely as long as no fee is charged
#    for use, and as long as the author/copyright attributions
#    are not removed. It remains the property of the copyright holder.
# History:
# 26 March 2017     - Finished initial version

my ( $gbkPath, $outDir);
my ( $ptt, $ffn, $faa, $xls, $l ) = ( 0, 0, 0, 0, 0);

GetOptions(
	'i|input=s' => \$gbkPath,
    'o|out-dir=s' => \$outDir,
    'ptt' => \$ptt,
    'ffn' => \$ffn,
    'faa' => \$faa,
    'xls' => \$xls,
    'l|length=i' => \$l
);

sub usage {
	print <<USAGE;
convert_gbk.pl Takes a Genbank file or a directory containing Genbank files. Parses through each file for 
	       every CDS. Optionally outputs a .ptt format conversion, a .ffn (fasta format nucleotide),
               a .faa (fasta format amino acid), or a .xls feature table. Can also produe an .ffn of 
               a given upstream length + 3 to include the start codon. This only works for prokaryotic sequences
               because it does not handle Splits or Joins found in eukaryotes. This only works for linearized genomes.

Usage:	[option]
    -i | --input    Path to input Genbank file or directory
    -o | --out-dir  Path to output directory. Creates directory if it doesn't exist
    --ptt           Convert CDS to .ptt format
    --ffn           Convert CDS to .ffn (fasta format nucleotide) 
    --xls           Create and associated feature table (.xls) for .ffn file
    --faa           Convert CDS to .faa (fasta format amino acid)
    -l | --length   length of upstream region to extract (outputs to separate .ffn)

Examples:
perl gbk_to_faa.pl -i input --out-dir dir --ptt --ffn --ffa --xls -l 35

USAGE
	exit;
}

usage() unless $gbkPath && $outDir && ( $ptt || $ffn || $faa || $xls || $l );

if (-d $gbkPath) {
	opendir(DIR0, $gbkPath) or die "Can't open $gbkPath: $!\n";
	while (defined(my $gbkName = readdir(DIR0))) {
		next unless $gbkName =~ /\.gbk$/;
        produce_output(parse_gbk(join(File::Util->SL, $gbkPath, $gbkName), $l), $gbkName, $outDir, $ptt, $ffn, $xls, $faa, $l);
	}
    closedir(DIR0);
} elsif (-e $gbkPath && $gbkPath =~ /\.gbk$/) {
    produce_output(parse_gbk($gbkPath, $l), $gbkPath, $outDir, $ptt, $ffn, $xls, $faa, $l,);
} else {
	print "Input file(s) do not meet requirements.\n\n";
	usage();
}

sub parse_gbk {
    my ($gbkName, $l ) = @_;

	my $gbk_object = Bio::SeqIO->new(-file => , "<$gbkName", -format => 'genbank');

    my %seqsData;
    while (my $seq_object = $gbk_object->next_seq) {
        my $accession = $seq_object->accession_number;
        my @organism = "";
        my @strain = "";
        my @type = "";
        my $range;

        my %cdsData;
        for my $feat_object ($seq_object->get_SeqFeatures) {
            my $feat_type = $feat_object->primary_tag;
            if ($feat_type eq "source") {
                if ($feat_object->has_tag("organism")) {
                    @organism = $feat_object->get_tag_values("organism");
                } 
                if ($feat_object->has_tag("strain")) {
                    @strain = $feat_object->get_tag_values("strain");
                }
                if ($feat_object->has_tag("mol_type")) {
                    @type = $feat_object->get_tag_values("mol_type");
                }
                $range = $feat_object->start."..".$feat_object->end;
            } elsif ($feat_type eq 'CDS') {

                my ($nucSeq, $loc, $upStreamSeq, $strand, $length, $geneInfo,
                $locusTag, $productInfo, $noteInfo, $dbXref, $protSeq) = get_data($feat_object, $l);

                $cdsData{$locusTag} = [$nucSeq, $loc, $upStreamSeq, $strand, $length,
                                       $geneInfo, $productInfo, $noteInfo, $dbXref, $protSeq]
            }
        }
        $seqsData{join(":", $accession, $organism[0], $strain[0], $type[0], $range)} = \%cdsData;
    }
    return(\%seqsData);
}

sub get_data {
    my ($feat_object, $l) = @_;
    
    ## sequence
    my @nucSeq = $feat_object->spliced_seq->seq;

    ## location range
    my $location = $feat_object->location;
    my $loc = $location->start."..".$location->end."\t";

    ## strand
    my $strand = $feat_object->strand >= 0 ? '+' : '-';

    ## upstram seq
    my $upStreamSeq;
    if ($strand eq "+" && defined $l) {
        $upStreamSeq = $feat_object->entire_seq()->subseq($location->start - $l, $location->start + 2)
    } elsif (defined $l) {
        $upStreamSeq = $feat_object->entire_seq()->subseq($location->end - 2, $location->end + $l)
    }
    
    ## length
    my $length = $feat_object->length;

    ## gene info
    my @geneInfo;
    if ($feat_object->has_tag('gene')) {
        @geneInfo = $feat_object->get_tag_values('gene');
    } else {
        @geneInfo = "-";
    }

    ## locus tag
    my @locusTag;
    if ($feat_object->has_tag('locus_tag')) {
        @locusTag = $feat_object->get_tag_values('locus_tag');
    } else {
        @locusTag = "-";
    }

    ## product tag
    my @productInfo;
    if ($feat_object->has_tag('product')) {
        @productInfo = $feat_object->get_tag_values('product');
    } else {
        @productInfo = "-";
    }

    ## notes
    my @noteInfo;
    if ($feat_object->has_tag('note')) {
        @noteInfo = $feat_object->get_tag_values('note');
    } else {
        @noteInfo = "-";
    }

    ## db xref
    my @dbXref;
    if ($feat_object->has_tag('db_xref')) {
        @dbXref = $feat_object->get_tag_values('db_xref')
    } else {
        @dbXref = "-";
    }

    ## protein seq
    my @protSeq;
    if ($feat_object->has_tag('translation')) {
        @protSeq = $feat_object->get_tag_values('translation')
    } else {
        @protSeq = "-";
    }

    return ($nucSeq[0], $loc, $upStreamSeq, $strand, $length, $geneInfo[0],
            $locusTag[0], $productInfo[0], $noteInfo[0], $dbXref[0], $protSeq[0]);
}

sub produce_output {
    my ($seqsData, $gbkName, $outDir, $ptt, $ffn, $xls, $faa, $l) = @_;

    $outDir =~ s/[\/\\]+$//;
    if (!-e $outDir) { make_path($outDir); }
    if (!-e $outDir) { die "Can't find output directory \"$outDir\"\n"; }
    $outDir .= File::Util->SL;

    my $base = basename($gbkName);

    my @key = keys %$seqsData;
    my ($accession, $organism, $strain, $type, $range ) = split(":", $key[0]);

    my ($pttOut, $ffnOut, $xlsOut, $faaOut, $upstreamOut);
    if ($ptt) {
        open ($pttOut, ">", join("/", $outDir, substr($base,0,-4).".ptt"));
        print $pttOut join(" ", $organism, $strain, $type, $range)."\n";
        print $pttOut (scalar keys %{ $$seqsData{$key[0]} })." proteins\n";
        print $pttOut "Location\tStrand\tLength\tPID\tGene\tSynonym Code\tCOG\tProduct\n";
    }
    if ($ffn) {
        open ($ffnOut, ">", join("/", $outDir, substr($base,0,-4).".ffn"));
    }
    if ($xls) {
        open ($xlsOut, ">", join("/", $outDir, substr($base,0,-4).".xls"));
        print $xlsOut "Coordinates\tStrand\tLength\tGene\tLocus\tProduct\tMiscellaneous\n";
    }
    if ($faa) {
        open ($faaOut, ">", join("/", $outDir, substr($base,0,-4).".faa"));
    }
    if ($l) {
        open ($upstreamOut, ">", join("/", $outDir, substr($base,0,-4)."-upstream-$l"."bp.ffn"));
    }

    foreach my $locusTag (sort keys %{ $$seqsData{$key[0]} }) {
        my ($nucSeq, $loc, $upstreamSeq, $strand, $length, $geneInfo,
            $productInfo, $noteInfo, $dbXref, $protSeq) = @{ ${ $$seqsData{$key[0]} }{$locusTag} };

        if ($ptt) {
            print $pttOut "$loc\t$strand\t$length\t$dbXref\t$geneInfo\t$locusTag\t-\t$productInfo\n";
        }
        if ($ffn | $faa) {
            my $fastaHeader = fasta_header($accession, $organism, $strain, $type, $loc, $productInfo, $locusTag, $geneInfo, $noteInfo);
            if ($ffn) {
                print $ffnOut "$fastaHeader$nucSeq\n";
            }
            if ($faa) {
                print $faaOut "$fastaHeader$protSeq\n";
            } 
        }
        if ($l) {
            my $upstreamHeader = upstream_header($accession, $organism, $strain, $type, $l, $loc, $strand, $productInfo, $locusTag, $geneInfo, $noteInfo);
            print $upstreamOut "$upstreamHeader$upstreamSeq\n";
        }
        if ($xls) {
            print $xlsOut "$loc\t$strand\t$length\t$geneInfo\t$locusTag\t$productInfo\t$noteInfo\n";
        }
    }

    if ($ptt) {
        close $pttOut;
    }
    if ($ffn) {
        close $ffnOut;
    }
    if ($xls) {
        close $xlsOut;
    }
    if ($faa) {
        close $faaOut;
    } 
    if ($l) {
        close $upstreamOut;
    }   
}

sub fasta_header {
    my ($accession, $organism, $strain, $type, $loc, $productInfo, 
        $locusTag, $geneInfo, $noteInfo) = @_;

    my $fastaHeader;
    if ( !($locusTag eq "-") and ($productInfo eq "-") ) {
        $fastaHeader = ">$accession|$locusTag|$organism $strain|$type|From $loc\n";
    } elsif ( ($locusTag eq "-") and !($productInfo eq "-") ) {
        $fastaHeader = ">$accession|$productInfo|$organism $strain|$type|From $loc\n";
    } elsif ( !($locusTag eq "-") and !($productInfo eq "-") ) {
        $fastaHeader = ">$accession|$locusTag|$productInfo|$organism $strain|$type|From $loc\n";
    } else {
        $fastaHeader = ">$accession|$organism $strain|$type|From $loc\n";
    }
    return($fastaHeader);
}

sub upstream_header {
    my ($accession, $organism, $strain, $type, $l, $loc, $strand,
        $productInfo, $locusTag, $geneInfo, $noteInfo) = @_;

    my @locs = split("[..]", $loc);
    my ($start, $end);
    if ($strand eq "+") {
        $start = $locs[0] - $l;
        $end = $locs[0] + 2;
    } else {
        $start = $locs[2] - 2;
        $end = $locs[2] + $l;
    }
    
    my $fastaHeader;
    if ( !($locusTag eq "-") and ($productInfo eq "-") ) {
        $fastaHeader = ">$accession|$locusTag|$organism $strain|$type|From $start..$end\n";
    } elsif ( ($locusTag eq "-") and !($productInfo eq "-") ) {
        $fastaHeader = ">$accession|$productInfo|$organism $strain|$type|From $start..$end\n";
    } elsif ( !($locusTag eq "-") and !($productInfo eq "-") ) {
        $fastaHeader = ">$accession|$locusTag|$productInfo|$organism $strain|$type|From $start..$end\n";
    } else {
        $fastaHeader = ">$accession|$organism $strain|$type|From $start..$end\n";
    }
    return($fastaHeader);
}