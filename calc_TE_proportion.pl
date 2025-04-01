#! /usr/bin/perl

use strict;
use warnings;

die "file?\n$0 GFF FA" if (@ARGV!=2);

my $gff=shift;
my $fa=shift;

my (%seq_len,%TE_len,%gene_trans_id);
open GFF, "<$gff";
while (<GFF>){
	chomp;
	next if (/^\#\#gff-version/);
	if (/^\#\#sequence-region/) {
		my @a=split/\s+/,$_;
		$seq_len{$a[1]}=$a[-1];
	} elsif ($_!~m/^\#\#/){
		my @a=split/\t/,$_;
		$TE_len{$a[0]}+=$a[4]-$a[3];
	} else {
		die "wrong format: $_\n";
	}
}
close GFF;

open FA, "<$fa";
#>anno1.g21253.t1 gene=g_23503 seq_id=h1tg000190l type=cds
while (<FA>){
	chomp;
	if ($_=~s/^>//){
		my @a=split/\s+/,$_;
		$a[1]=~s/gene=//;
		$gene_trans_id{$a[0]}=$a[1];
	} else {
		next;
	}
}
close FA;


for my $i (keys %seq_len){
	my $gene_length=$seq_len{$i};
	my $te_length=$TE_len{$i};
	my $prop=$te_length/$gene_length;
	print "$i\t$gene_trans_id{$i}\t$te_length\t$gene_length\t$prop\n";
}
