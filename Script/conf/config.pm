package conf::config;
require Exporter;
use warnings;

@ISA = qw(Exporter);
@EXPORT = qw($genome_reference $CCDS_file $CCDS_NM_file $samtools $java $clin_var $vcftools);


#Reference ############

$genome_reference= "/FILE_PATH/*.fasta";

$CCDS_file="/FILE_PATH/CCDS.current.txt";

$CCDS_NM_file="/FILE_PATH/CCDS2Sequence.current.txt";

$clin_var="/FILE_PATH/clinvar_*.vcf.gz";

#Tools ############

$samtools="/DIR_PATH/samtools";

$vcftools="/DIR_PATH/vcftools";

$java="java";

1;
