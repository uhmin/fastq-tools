#!/usr/bin/perl
use strict;
&main;

sub main{
    my %options;
    my $i;
    my @filelist;
    my $file;
    my %data;
    my $name;
    %options=&interface(@ARGV);
    for($i=0; $i<2; $i++){
	$name=sprintf("-%d", $i+1);
	@filelist=@{$options{$name}};
	foreach $file (@filelist){
	    &readOneData($i, $file, \%data);
	}
    }
    &compare(\%data);

    &output(\%data, \%options);
}

sub interface{
    my $i;
    my %options;

    if(scalar(@_) % 2 == 1){
        &help(1);
    }

    for($i=0; $i<@_; $i+=2){
        if($_[$i] eq "-1"
           || $_[$i] eq "-2"){
            push(@{$options{$_[$i]}}, $_[$i+1]);
        }elsif($_[$i] eq "-o1" 
	       || $_[$i] eq "-o2"){
            $options{$_[$i]}=$_[$i+1];
	}else{
            print stderr "Unknown option: $_[$i]\n";
            &help(1);
        }
    }
    if(!exists($options{-1})){
	print stderr "Option -1 must be specified.\n";
	&help(1);
    }
    if(!exists($options{-2})){
	print stderr "Option -2 must be specified.\n";
	&help(1);
    }

    return %options;
}

sub help{
    my $file=__FILE__;
    print << "EOF"
$file
     Compare two fastq sequence files, delete broken pairs and output to file(s).

usage:
  -1 [input file name]  Input file name for pair 1.
  -2 [input file name]  Input file name for pair 2.
Optional
  -o1 [output file name]  Outout file name for pair 1.
  -o2 [output file name]  Outout file name for pair 2. 

-1 and -2 can be used repeadedly. 
Such as:
  $file -i file1_for_pair1 -2 fileA_for_pair2 -1 file2_for_pair1 \
       -1 fileB_for_pair2 -o1 outfile.fastq

If -o2 were not specified, Interleaved FastQ 
  will be output in the file specified in -o1 option or STDOUT.

EOF
;
    exit($_[0]);
}

sub readOneData{#($i, $file, \%data);
    my $i=$_[0];
    my $file=$_[1];
    my $data=$_[2];
    my $line;
    my $count;
    my @name;
    my $sequence;
    open FIN, "$file" or die("Could not open infile: $file\n");
    for($count=0; $line=<FIN>; $count++){
	if($count % 4 == 0){
	    $data->{$name[0]}[$name[1]]=$sequence;
	  
	    # initialize sequence #
	    $sequence=$line;
	    chomp($line);
	    @name=split(/[\#\s]/, $line);
	    $name[1]=~s/.*\///;
	    $name[1]=~s/:.*//;
#	    printf("   -- %s\n   => %s %d\n", $line, $name[0], $name[1]);
	}else{
	    $sequence.=$line;
	}
    }
    $data->{$name[0]}[$name[1]]=$sequence;
    delete($data->{""});
    close FIN;
}

sub compare{#(\%data);
    my $key;
    my $data=$_[0];
    foreach $key (keys %{$data}){
	if(!exists($data->{$key}[1]) ||  !exists($data->{$key}[2])){
	    print stderr "Delete key $key\n";
	    delete($data->{$key});
	}
    }
}

sub output{#(\%data, \%options);
    my @file;
    my $key;
    my $sequence;
    if(exists($_[1]->{-o1})){
	$file[0]=$_[1]->{-o1};
	open FOUT, ">$file[0]" or die("Could not open outfile: $file[0]\n");
	if(exists($_[1]->{-o2})){
	    $file[1]=$_[1]->{-o2};
	    open FOUT2, ">$file[1]" or die("Could not open outfile: $file[1]\n");
	}else{
	    *FOUT2=*FOUT;
	}
    }else{
	*FOUT=*STDOUT;
	*FOUT2=*STDOUT;
    }
    foreach $key (keys %{$_[0]}){
	print FOUT $_[0]->{$key}[1];
	print FOUT2 $_[0]->{$key}[2];
    }
    
    close FOUT;
}
