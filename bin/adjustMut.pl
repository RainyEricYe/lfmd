#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Cwd qw/abs_path/;
use File::Basename;

my ( $Help, $freq,)
=  ( 0,     1e-7, );
GetOptions(
    'help|?' => \$Help,
    'freq=f' => \$freq,
);
die `pod2text $0` if $Help or @ARGV < 1;

my ($inf, ) = @ARGV;

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#&showLog("START");

my @ps;
open my $IN, "<$inf" or die $!;
LINE:
  while ( <$IN> ) {
    chomp;
    my @t = split;

    for my $x ( @t[12..$#t] ) {
        if ( $x =~ /^[ID]/ ) {
            next LINE;
        }
        else {
            my $p = (split /:/, $x)[1];
            next if $p == 0;
            push @ps, $p;
        }
    }
}
close $IN;

my $num = scalar @ps;
my %h;
my $rank = 1;
for my $k (sort {$a<=>$b} @ps) {
    push @{ $h{$k} }, $k * $num / $rank;
    $rank++;
}

$" = "\t";

my @nt = qw/T C G A/;

open my $RE, "<$inf" or die $!;
my $title = <$RE>;
print $title;
my %nearIndel;
open my $TMP, ">$inf.tmp" or die $!;

while ( <$RE> ) {
    chomp;
    my @t = split;

    if ( /[ID]/ && $t[4] > 0 ) {
        for my $i (0..3) {
            if ( $t[$i+5] > 0 && $nt[$i] ne $t[1] ) {
                $t[4] -= $t[$i+5];
                $t[$i+5] = 0;
            }
        }
    }

    my %snv = ();
    my $out = "";
    my %indel;
    for my $x ( @t[12..$#t] ) {
        if ( $x =~ /^[ID]/ ) {
            my ($t, $num) = split /:/, $x;
            $indel{$t} = $num;
        }
        else {
            my ($base, $p) = split /:/, $x;
            next if $base eq "*";

            if ( !defined $p ) {
                next;
            }
            elsif ( $p == 0 ) {
                $out .= "\t$x";
                $snv{$base}++;
            }
            elsif ( !/[ID]/ ) {
                my $adjust_p = pop @{$h{$p}};
                $out .= sprintf "\t$base:%.2e", $adjust_p;

                if ( $adjust_p > $freq ) {
                    $out .= "F";
                }
                else {
                    $snv{$base}++;
                }
            }
            else {
                1;
            }
        }
    }


    # skip snv which don't have pvalue
    for my $i (0..3) {
        if ( $t[$i+5] > 0 && !defined $snv{ $nt[$i] } ) {
            $t[4] -= $t[$i+5];
            $t[$i+5] = 0;
        }
    }


    my @ks = sort { $indel{$b} <=> $indel{$a} } keys %indel;
    if ( scalar @ks > 1 ) {
        for my $k (@ks[1..$#ks]) {
            if ( $k =~ /^I/ ) {
                $t[9] -= $indel{ $k };
            }
            else {
                $t[10] -= $indel{ $k };
            }
        }
    }

    print $TMP "@t[0..11]$out";
    if ( defined $ks[0] ) {
        print $TMP "\t$ks[0]:$indel{ $ks[0] }";

        if ( $ks[0] =~ /I/ ) {
            $nearIndel{ $t[0] }{ $t[2] } = $t[9]/$t[3];
        }
        else {
            $nearIndel{ $t[0] }{ $t[2] } = $t[10]/$t[3];
        }
    }
    print $TMP "\n";
}
close $RE;
close $TMP;

open my $TH, "<$inf.tmp" or die $!;
while (<$TH>) {
    chomp;
    my @t = split;
    my $nearby = 0;


    if ( defined $nearIndel{ $t[0] } ) {
        for my $k (keys %{ $nearIndel{ $t[0] } }) {

            if ( (( $k - $t[2] ) ** 2) <= 25 ) {  # distance 5 bp
                if ( $nearIndel{ $t[0] }{$k} > $nearby ) {
                    $nearby = $nearIndel{ $t[0] }{$k};
                }
            }
        }
    }

    if ( $nearby > 0 ) {
        if ( /[ID]/ ) {
            print "$_\n";
        }
        else {
            for my $i (0..3) {
                if ( $t[4] > 0 && $t[$i+5] > 0 && ($t[$i+5]/$t[3]) <= $nearby ) {
                    $t[4] -= $t[$i+5];
                    $t[$i+5] = 0;
                }
            }

            print "@t\n";
        }
    }
    else {
        print "$_\n";
    }
}
close $TH;

`rm -f $inf.tmp`;

&showLog("END");

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

sub showLog {
    my @t = localtime();
    printf STDERR "[%04d-%02d-%02d %02d:%02d:%02d]\t%s\n", $t[5] + 1900, $t[4] + 1, @t[3,2,1,0], $_[0];
}

__END__

=head1 Function
    adjust the pvalues of mutations detected by lhmut

=head1 Usage
    perl $0 in.lhmut > out.mut.adj

=head1 Options
    -f [f]  lowest frequency [1e-7]
    -h      help

=head1 Author
    yerui;    yerui@connect.hku.hk

=head1 Version
    v4.0;    2022-7-8

=cut
