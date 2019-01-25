#!/usr/bin/perl

use strict;
use warnings;

use autodie;
use Try::Tiny;

my $HOME_DIR = "/scratch/asl47/PDB";
my $WATER_EXEC = "/rscratch/asl47/water";
my $FREESASA_EXEC = "/rscratch/asl47/freesasa";

system("mkdir -p $HOME_DIR/results");

my ($pdb_raw, $HOMOMER) = @ARGV;
my ($pdb_id, $BA_id) = split(/_/, $pdb_raw, 2);


print "Running DIMER analysis for $pdb_id with assembly $BA_id \n";
if (-f "$HOME_DIR/${pdb_id}.pdb${BA_id}") {
    #print "Existing pdb file\n"; 
}
else {
    #print "Downloading pdb file\n";
    system("wget -q \"https://files.rcsb.org/download/${pdb_id}.pdb${BA_id}.gz\" -P $HOME_DIR");
    #print "Downloaded fasta file\n";
    system("wget -q \"https://www.rcsb.org/pdb/download/downloadFastaFiles.do\?structureIdList=${pdb_id}&compressionType=uncompressed\" -O ${HOME_DIR}/${pdb_id}.fasta ");
    #print "Extracting pdb\n";	    
    system("gunzip -f $HOME_DIR/${pdb_id}.pdb${BA_id}.gz");
}

open (my $fasta, "<", "${HOME_DIR}/${pdb_id}.fasta");
try {
    my @chains;
    {
        local $/ = '>';
        while (<$fasta>) {
            s/^>//g; # strip out '>' from beginning
            s/>$//g; # and end of line
            next if !length($_);           # ignore empty lines

            my ($header_info) = /^(.*)\n/; # capture the header
            s/^(.*)\n//;                   # and strip it out
            my $chain_id=substr($header_info,5,1);
            push @chains, $chain_id;
           
            s/\n//mg;

            open (my $fasta_chain, ">", "${HOME_DIR}/${pdb_id}.fasta${chain_id}");
            print $fasta_chain "$_";
            close $fasta_chain;
        }
    }

    my $results = "${HOME_DIR}/results/${pdb_id}.results";
    system("echo ${HOMOMER} > ${results}"); 
    if($HOMOMER) {
        if($#chains) { 
            #print "Multi chain homomer\n";
            system("${WATER_EXEC} ${HOME_DIR}/${pdb_id}.fasta${chains[0]} ${HOME_DIR}/${pdb_id}.fasta${chains[1]} -gapopen 10 -gapextend .5 -nobrief -stdout -auto | grep Longest >> ${results}");
            system("${FREESASA_EXEC} ${HOME_DIR}/${pdb_id}.pdb${BA_id} --chain-group=${chains[0]}${chains[1]}+${chains[0]}+${chains[1]} | grep Total | tail +2 >> ${results}");   

        }
        else {
            #print "One chain homomer\n";
            system("${WATER_EXEC} ${HOME_DIR}/${pdb_id}.fasta${chains[0]} ${HOME_DIR}/${pdb_id}.fasta${chains[0]} -gapopen 10 -gapextend .5 -nobrief -stdout -auto | grep Longest >> ${results}");
            system("${FREESASA_EXEC} ${HOME_DIR}/${pdb_id}.pdb${BA_id} --join-models | grep Total >> ${results}");  
            system("${FREESASA_EXEC} ${HOME_DIR}/${pdb_id}.pdb${BA_id} --separate-models | grep Total >> ${results}"); 
        }
    }
    else {
        #print "Heteromer\n";
        system("${WATER_EXEC} ${HOME_DIR}/${pdb_id}.fasta${chains[0]} ${HOME_DIR}/${pdb_id}.fasta${chains[1]} -gapopen 10 -gapextend .5 -nobrief -stdout -auto | grep Longest > ${results}");
        system("${FREESASA_EXEC} ${HOME_DIR}/${pdb_id}.pdb${BA_id} --chain-group=${chains[0]}${chains[1]}+${chains[0]}+${chains[1]} | grep Total | tail +2 >> ${results}");   
    }
} catch {
    warn "ID ${pdb_id} caught error: $_";
};

for my $file_d (grep {/fasta\S/} <${HOME_DIR}/*>) {
    #print "unlinking $file_d\n";
    unlink $file_d;
}

print "finished analysis of : $pdb_id\n";
