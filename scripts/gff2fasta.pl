#!/usr/bin/env perl
###############################################################################
#
#    gff2fasta.pl
#     
#    Convert a GFF file to fasta
#
#    Copyright (C) Michael Imelfort
#    Edited: Kranti Konganti
#    Adapted by Sishuo Wang from The University of British Columbia (sishuowang@hotmail.ca)
#
#    Updates
#    2015-09-16
#    New features:
#      bugs related to --attribute fixed
#    2014-10-15
#    New features:
#      --with_ambiguous_posi:  include information of positions (start + end) with '<' or '>
#      --with_ambiguous_start: include information of start positions with '<'
#      --with_ambiguous_end:   include information of end positions with '>
#    2014-09-02
#      Problems of extracting sequences from the Crick strand have been fixed.
#    2014-08-18
#    New features:
#      --upstream, --downstream, --no_middle
#      --combine_exons
#      --first_feature
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################

#pragmas
use 5.010;
use strict;
#use warnings;

#core Perl modules
use Getopt::Long;
use Carp;
use List::Util qw/max min/;
no autovivification;

#CPAN modules
use Bio::SeqIO;
use Bio::Perl;
use Bio::Tools::CodonTable;

#locally-written modules

BEGIN {
    select(STDERR);
    $| = 1;
    select(STDOUT);
    $| = 1;
}

# edit here to log all external commands 
my $global_log_commands = 0;

# ext command failure levels
use constant {
    IGNORE_FAILURE => 0,
    WARN_ON_FAILURE => 1,
    DIE_ON_FAILURE => 2
};

# get input params and print copyright
my $global_options = checkParams();
if(!exists $global_options->{'silent'}) {
    printAtStart();
}


######################################################################
# CODE HERE
######################################################################
my %attributes_included;

foreach (split (",",$global_options->{"attributes"})){
    $attributes_included{$_}=1;
}

my %attributes_included_all=(
    'protein_id'	=>  1,
    'Name'		=>  1,
    'product'	=>  1,
);

#map {exists $attributes_included_all{$_}} keys %attributes_included;

my $upstream=$global_options->{'upstream'};
my $downstream=$global_options->{'downstream'};
my %global_gff_orfs = ();
my %global_gff_non_orfs = ();
my %global_lengths = ();
my %attributes_frame;
my $global_reject_length = overrideDefault(50,'length');
my $global_keep_non_orfs = overrideDefault(0, 'non_orfs');
my $global_protein_code =  overrideDefault(0,'protein');
if((0 != $global_keep_non_orfs) and (0 != $global_protein_code))
{
    print "**WARNING: $0 : non_orfs flag is invalid when translating into protein space --> ignoring\n";
    $global_keep_non_orfs = 0;
}
my $global_line_wrap = overrideDefault(80,'wrap');

# Read in the gff
my $current_seq_id = "";
my $gff_fh = openRead($global_options->{'gff'});
while(<$gff_fh>){
    chomp;
    next if ($_ =~ m/^#/);
    my (%attributes, %posi);
    my ($seqid, undef, $feature, $start, $end,
        undef, $strand, $frame, $attributes) = split /\t/;

    # SSW
    if (not exists $attributes_frame{$attributes}){
        $attributes_frame{$attributes} = "";
        if($strand == "+"){
            $start = $start + $frame;
        }
        else{
            $end = $end - $frame;
        }
    }
    # SSW

    if ($start =~ /[<>]/ or $end =~ /[<>]/){
        if ((not exists $global_options->{'with_ambiguous_start'}) and (not exists $global_options->{'with_ambiguous_end'})){
            ;
        }
        else{
            if (exists $global_options->{'with_ambiguous_start'}){
                $start =~ s/[<]//;
            }
            if (exists $global_options->{'with_ambiguous_end'}){
                $end =~ s/[>]//;
            }
        }
    }
    
    if (not $global_options->{'no_middle'}){
        $posi{'middle'}{'start'}=$start; # coding regions for example
        $posi{'middle'}{'end'}=$end; # coding regions for example
    }
    if ($strand eq "+"){
        $posi{'upstream'}{'end'}=$start-1 and $posi{'upstream'}{'start'}=$start-$upstream if $global_options->{'upstream'};
	$posi{'downstream'}{'start'}=$end+1 and $posi{'downstream'}{'end'}=$end+$downstream if $global_options->{'downstream'};
    } else {
        $posi{'upstream'}{'start'}=$end+1 and $posi{'upstream'}{'end'}=$end+$upstream if $global_options->{'upstream'};
	$posi{'downstream'}{'end'}=$start-1 and $posi{'downstream'}{'start'}=$start-$downstream if $global_options->{'downstream'};
    }
    $seqid =~ s/\.\d+$// if not $global_options->{'no_process_seqid'};
    while($attributes =~ /([^;=]+)=([^;]+)/g){
        next if not exists $attributes_included{$1};
        $attributes{$1}=$2;
    }
    push @{$global_gff_orfs{$seqid}}, [\%posi, $strand, \%attributes, $feature];
}
close $gff_fh;

# open the output file
my $out_fh = openWrite($global_options->{'out'});

foreach my $genome (@{$global_options->{'fasta'}}){
# Read in the fasta
my $seqio = Bio::SeqIO->new( -file => $genome, -format => 'fasta' ) or croak "**ERROR: Could not open FASTA file: $global_options->{'fasta'} $!\n";
while(my $sobj = $seqio->next_seq)
{
    my $seqid = $sobj->id;
    my $seq = $sobj->seq;
    my $seq_length = $sobj->length;

    # make sure this guy has an annotation
    if(defined($global_gff_orfs{$seqid}))
    {
	    my %combined_info;
        for(@{$global_gff_orfs{$seqid}})
        { 
	        my ($combined_seq, $strand_direction_abbr, $strand_direction);
            my ($posi_href, $strand, $attributes_href, $feature) = @$_;
	        next if (defined $global_options->{'feature'} && not exists $global_options->{feature}{$feature});
     
            # work out the length of the sub string

            foreach my $posi_type (keys %{$posi_href}){
                next if $posi_href->{$posi_type}{'start'} =~ /[<>]/ or $posi_href->{$posi_type}{'end'} =~ /[<>]/;
                my $start = max($posi_href->{$posi_type}{'start'}, 0);
                my $end = min($posi_href->{$posi_type}{'end'}, length($seq));
                my $length = $end - $start + 1;
                # check if he's long enough
                if($length < $global_reject_length){
                    if(!exists $global_options->{'silent'}) {
                        print "Rejecting: $seqid -> ($start, $end) on $strand. $length is shorter than cutoff!\n";
                    }
                    next;
                }
                my $this_seq = substr($seq, $start-1, $length); #|| print $seqid."\t".$start,$start+$length, "\n";
                $strand_direction = ($strand eq "+") ? "FORWARD" : "REVERSE";
                $strand_direction_abbr = ($strand eq "+") ? "F" : "R";
                my $seq=$strand eq "+" ? fasta_cut($this_seq, $global_protein_code) : fasta_cut(revcompl($this_seq), $global_protein_code);
                $combined_seq=$seq;

                if ($global_options->{'combine_exons'}){
                    my $ID;
                    for my $attr (keys %attributes_included){
                        $ID=$attributes_href->{$attr} if exists $attributes_href->{$attr}
                    }
                    #my $ID= exists $attributes_href->{"ID"} ? $attributes_href->{"ID"} : $attributes_href->{"Parent"};
                    #map {print $_."\t".$attributes_href->{$_}."\n"} keys %$attributes_href; exit;
                    $combined_info{$ID}{'strand'} = $strand;
                    if (defined $combined_seq){
                        push @{$combined_info{$ID}{'seqs'}}, $combined_seq;
                        push @{$combined_info{$ID}{'seq_starts'}}, $start;
                    }
                    $combined_info{$ID}{'attributes'}=$attributes_href;
                }
            
                else{
                    next if not defined $combined_seq;
                    if ($global_options->{'seq-desc'}) {
                        #print values %{$attributes_href}; print "\n";
                        print $out_fh ">".join('|',values %{$attributes_href}). "| $seqid: $strand_direction\n";
                    }
                    else{
                        print $out_fh ">$seqid"."_$strand_direction_abbr\n";
                    }
                    print $out_fh $combined_seq."\n";
                }
            }
	    }

    # ------------------------------------------------------ #
	if ($global_options->{'combine_exons'}){
	    while(my($ID, $href) = each %combined_info){
            my $title = join ("|", map {$href->{'attributes'}{$_} if exists $attributes_included{$_}} keys %{$href->{'attributes'}});
            if (exists $href->{'seqs'}){
                my @starts = @{$href->{'seq_starts'}};
                #for my $index (sort {$starts[$a]<=>$starts[$b]} 0..$#starts){
                my @ori_starts = @starts;
                @starts = sort {$a<=>$b} @starts;
                my %start_posi;
                for my $i (0..$#ori_starts){
                    my $start_2 = $ori_starts[$i];
                    my ( $j )= grep { $starts[$_] eq $start_2 } 0..$#starts;
                    $start_posi{$j} = $i;
                }
                $start_posi{-1} = $start_posi{$#starts};

                for my $index (0..$#starts){
                    #print join("\t", $index, $starts[$index])."\n";
                    given($href->{'strand'}){
                        when ('+'){
                            $href->{'seq'} .= $href->{'seqs'}->[$index];
                            my $index;
                            $index = 0 if defined($upstream);
                            $index = -1 if defined($downstream);
                            my $index2 = $start_posi{$index};
                            $href->{'first_seq'} = $href->{'seqs'}->[$index2];
                        }
                        when ('-'){
                            $href->{'seq'} = $href->{'seqs'}->[$index] . $href->{'seq'};
                            my $index;
                            $index = -1 if defined($upstream);
                            $index = 0 if defined($downstream);
                            my $index2 = $start_posi{$index};
                            $href->{'first_seq'} = $href->{'seqs'}->[$index2];
                        }
                    }
                }
            }

            my $seq = $href->{'seq'};
            if ($global_options->{"first_feature"}){
                # only first
                next if not defined $href->{'first_seq'};
                print $out_fh ">$title\n$href->{'first_seq'}\n";
            }
            else{
                next if not defined $seq;
                print $out_fh ">$title\n$seq\n";
            }
	    }
	}
    }

    elsif(exists $global_options->{'include-nulls'}){
        # include anyway
        print $out_fh ">$seqid"."_1_$seq_length"."_X\n".fasta_cut($seq, $global_protein_code)."\n";
    }
}
}

close $out_fh;

######################################################################
# CUSTOM SUBS
######################################################################
sub revcompl {
    my ($seq) = @_;
    $seq =~ tr/ACGTacgt/TGCAtgca/;
    return scalar reverse $seq;
}


sub fasta_cut {
    #-----
    # Cut up a fasta sequence 
    #
    my ($string, $prot) =  @_;
    
    # translate if need be
    if(0 != $prot) 
    { 
        my $codon_table  = Bio::Tools::CodonTable -> new ( -id => $prot );
        $string  = $codon_table->translate($string); 
    }
    
    # wrap the line if need be
    if(0 != $global_line_wrap)
    {
        my $return_str = "";
        my $len = length $string;
        my $start = 0;
        while($start < $len)
        {
            $return_str .= substr $string, $start, $global_line_wrap;
            #$return_str .="\n";
            $start += $global_line_wrap;
        }
        return $return_str;
    }
    return "$string\n";
}

######################################################################
# TEMPLATE SUBS

######################################################################
# PARAMETERS

sub checkParams {
    my %options;
    my (@features, @genomes);
    %options=(
	"protein"=>0,
	"non_orfs"=>0,
	"wrap"=>80,
	"length"=>50,
    );
    GetOptions ( "help|h+"	        =>	\$options{'help'},
         "gff|g:s"  	        =>	\$options{'gff'},
		 "fasta|f=s"	        =>  \@genomes,
         "out|o:s"	            =>	\$options{'out'},
         "protein|p:i"	        =>	\$options{'protein'},
         "length|l:i"	        =>	\$options{'length'},
         "wrap|w:i"	            =>	\$options{'wrap'},
         "non_orfs+"	        =>	\$options{'non_orfs'},
         "include_nulls+"       =>	\$options{'include_nulls'},
         "seq-desc|d"	        =>	\$options{'seq-desc'},
		 "feature=s"	        =>	\@features,
		 "attributes=s"	        =>	\$options{'attributes'},
		 "upstream=s"	        =>	\$options{'upstream'},
		 "downstream=s"	        =>	\$options{'downstream'},
		 "no_middle!"	        =>	\$options{'no_middle'},
		 "combine_exons!"       =>	\$options{'combine_exons'},
		 "first_feature!"       =>	\$options{'first_feature'},
		 "no_process_seqid!"    =>	\$options{'no_process_seqid'},
         "with_ambiguous_posi!" =>  \$options{'with_ambiguous_posi'},
         "with_ambiguous_start!"=>  \$options{'with_ambiguous_start'},
         "with_ambiguous_end!"  =>  \$options{'with_ambiguous_end'},
	         ) || die "illegal params!";
    map{delete $options{$_} if ! $options{$_}} keys %options;
    @{$options{'feature'}}{@features} = (1) x scalar(@features);
    $options{'fasta'}=\@genomes;
    # Add any other command line options, and the code to handle them
    # 
    #GetOptions( \%options, @standard_options );

    # if no arguments supplied print the usage and exit
    #
    exec("pod2usage $0") if (0 == (keys (%options) ));

    # If the -help option is set, print the usage and exit
    #
    exec("pod2usage $0") if $options{'help'};

    # Compulsory items
    if(!exists $options{'gff'} ) { printParamError ("You MUST supply a GFF file to parse"); }
    if(!exists $options{'fasta'} ) { printParamError ("You MUST supply a FASTA file to parse"); }
    if(!exists $options{'out'} ) { printParamError ("Please specify the name of the fasta file you'd like to create"); }
    #if(!exists $options{''} ) { printParamError (""); }
    if(exists $options{'combine_exons'}){
        $options{"attributes"} = $options{"attributes"} ? $options{"attributes"} : "ID,Parent";
    }
    if(exists $options{'with_ambiguous_posi'}){
        delete $options{'with_ambiguous_posi'};
        @options{qw(with_ambiguous_start with_ambiguous_end)}=(1,1);
    }
    print $options{"attributes"}."\n";
    return \%options;
}

sub printParamError
{
    #-----
    # What to do if there's something wrong with a parameter
    #  
    my ($error) = @_;  
    print "**ERROR: $0 : $error\n"; exec("pod2usage $0");
}

sub overrideDefault
{
    #-----
    # Set and override default values for parameters
    #
    my ($default_value, $option_name) = @_;
    if(exists $global_options->{$option_name}) 
    {
        return $global_options->{$option_name};
    }
    return $default_value;
}

######################################################################
# FILE IO

sub openWrite
{
    #-----
    # Open a file for writing
    #
    my ($fn) = @_;
    open my $fh, ">", $fn or croak "**ERROR: could not open file: $fn for writing $!\n";
    return $fh;
}

sub openRead
{   
    #-----
    # Open a file for reading
    #
    my ($fn) = @_;
    open my $fh, "<", $fn or croak "**ERROR: could not open file: $fn for reading $!\n";
    return $fh;
}

######################################################################
# EXTERNAL COMMANDS
#
# checkAndRunCommand("ls", {
#                          -a => ""
#                          }, 
#                          WARN_ON_FAILURE);

sub checkFileExists {
    #-----
    # Does a file exists?
    #
    my ($file) = @_;
    unless(-e $file) {
        croak "**ERROR: $0 : Cannot find:\n$file\n";
    }
}

sub logExternalCommand
{
    #-----
    # Log a command line command to the command line!
    #
    if(1 == $global_log_commands) {
        print $_[0], "\n";
    }
}

sub isCommandInPath
{
    #-----
    # Is this command in the path?
    #
    my ($cmd, $failure_type) = @_;
    if (system("which $cmd |> /dev/null")) {
        handleCommandFailure($cmd, $failure_type);
    }
}

sub runExternalCommand
{
    #-----
    # Run a command line command on the command line!
    #
    my ($cmd) = @_;
    logExternalCommand($cmd);
    system($cmd);
}

sub checkAndRunCommand
{
    #-----
    # Run external commands more sanelier
    #
    my ($cmd, $params, $failure_type) = @_;
    
    isCommandInPath($cmd, $failure_type);
    
    # join the parameters to the command
    my $param_str = join " ", map {formatParams($_)} @{$params};
    
    my $cmd_str = $cmd . " " . $param_str;
    
    logExternalCommand($cmd_str);

    # make sure that all went well
    if (system($cmd_str)) {
         handleCommandFailure($cmd_str, $failure_type)
    }
}

sub formatParams {

    #---------
    # Handles and formats the different ways of passing parameters to 
    # checkAndRunCommand
    #
    my $ref = shift;
    
    if (ref($ref) eq "ARRAY") {
        return join(" ", @{$ref});
    } elsif (ref($ref) eq "HASH") {
        return join(" ", map { $_ . " " . $ref->{$_}} keys %{$ref});
    }
    croak 'The elements of the $params argument in checkAndRunCommand can ' .
        'only contain references to arrays or hashes\n';
}


sub handleCommandFailure {
    #-----
    # What to do when all goes bad!
    #
    my ($cmd, $failure_type) = @_;
    if (defined($failure_type)) {
        if ($failure_type == DIE_ON_FAILURE) {
            croak "**ERROR: $0 : " . $! . "\n";
        } elsif ($failure_type == WARN_ON_FAILURE) {
            carp "**WARNING: $0 : " . $! . "\n";
        }
    }
}

######################################################################
# MISC

sub printAtStart {
print<<"EOF";
---------------------------------------------------------------- 
 $0
 Copyright (C) Michael Imelfort
 Edited: Kranti Konganti
    
 This program comes with ABSOLUTELY NO WARRANTY;
 This is free software, and you are welcome to redistribute it
 under certain conditions: See the source for more details.
---------------------------------------------------------------- 
EOF
}

__DATA__

=head1 NAME

    gff2fasta.pl

=head1 COPYRIGHT

   copyright (C) Michael Imelfort

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.

=head1 CHANGES

=item * 02/04/2013

   Fasta sequence headers now include either ID or Alias tags from gff files.
    
=head1 DESCRIPTION

   Convert a GFF file to fasta

=head1 SYNOPSIS

    gff2fasta.pl -gff|g GFF -fasta|f FASTA -out|o FASTA

      -gff -g GFF                  GFF3 file to parse
      -fasta -f FASTA              Original fasta file
      -out -o FASTA                Multi fasta file to create
      [-protein -p CODON_CODE]     Output protein sequences --> see below
      [-non_orfs]                  Process non-ORF regions [default: false]
      [-include_nulls]             Transparently write through contigs with no genes [default: false] 
      [-seq-desc|d]                Use either ID or Alias attributes of GFF file as sequence identifiers
      [-feature]           Print fasta for only mentioned features. Example: gene or CDS or exon etc...
      [-wrap -w LEN]               Line wrap at LEN chars [default: 80] Set to 0 for no wrap
      [-length -l LENGTH]          Reject any orfs shorter than this length [default: 50]
      [-help -h]                   Displays basic usage information
      [--first_feature]            Only displays the first feature of a gene
      [--upstream No.]             Display No. of base pairs upstream
      [--downtream No.]            Display No. of base pairs downstream
      [--no_middle]                Not display the region corresponding to the feature. It is useful when you just want to retrieve sequences corresponding to the flanking region
      [--with_ambiguous_posi]      start and end including '<' or '>' will be considered
                                   default: off
      [--with_ambiguous_start]     
      [--with_ambiguous_end]
      [--combine_exons]            combine features belonging to the same gene
      
      CODON_CODE
      
      Specify a number from the following list (Uses: Bio::Tools::CodonTable)
      
      1 Standard
      2 Vertebrate Mitochondrial
      3 Yeast Mitochondrial
      4 Mold, Protozoan,_and_CoelenterateMitochondrial_and_Mycoplasma/Spiroplasma
      5 Invertebrate Mitochondrial
      6 Ciliate, Dasycladacean_and_Hexamita_Nuclear
      9 Echinoderm Mitochondrial
      10 Euplotid Nuclear
      11 Bacterial
      12 Alternative Yeast_Nuclear
      13 Ascidian Mitochondrial
      14 Flatworm Mitochondrial
      15 Blepharisma Nuclear
      16 Chlorophycean Mitochondrial
      21 Trematode Mitochondrial
      22 Scenedesmus obliquus_Mitochondrial
      23 Thraustochytrium Mitochondrial
         
=cut
>>

