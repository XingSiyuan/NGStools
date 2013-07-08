use strict;
use warnings;
open(STDERR, ">myprogram.error") or die "cannot open error file: myprogram.error:$!\n";
use DBI;
use Digest::MD5::File qw( file_md5_hex );

my $database = 'ABGSAschema';
my $server = 'localhost';
my $user = 'root';
my $passwd = 'XXXXXX';

my $db = DBI->connect("dbi:mysql:$database:$server", $user, $passwd);
my $basedir='/srv/mds01/shared/Sus/ABGSA';
my $query = "select archive_name,lane_names_orig,md5sum_gzip from ABGSAschema_main;";

my $sql = $db->prepare($query);
$sql ->execute();
while (my $row = $sql->fetchrow_arrayref) {
	
	
	my $row2 = join("\t", @$row);
	my @row2array = split ("\t",$row2);
	my $ABGSA = $row2array[0];
	my $file=$row2array[1];
	my $md5sum_db=$row2array[2];
	my $pathfile = $basedir.'/'.$ABGSA.'/'.$file;
	my $filesize = -s $pathfile;
	my $md5 = file_md5_hex( $pathfile );
	my $md5ok='notok';
#       my $filesizecheck='notok';
#       my $checkfilesize_again=check($db,$file);
	if ($md5 eq $md5sum_db){
		$md5ok='ok';
	        update($db,$filesize,$file,$md5ok,$md5);
	}
#        if ($checkfilesize_again eq $filesize){
#                $filesizecheck='ok';
#	}
#       print "$ABGSA\t$file\t$checkfilesize_again\t$filesize\t$filesizecheck\t$md5\t$md5sum_db\t$md5ok\n";
        print "$ABGSA\t$file\t$filesize\t$md5\t$md5sum_db\t$md5ok\n";
	# $gelukt = update($kvl,$gene_adaptor,$hgnc_name, $ref);
	#if ($gelukt) {++$tweede;}
	
}

$db->disconnect;

print "\n\nRun ends at " . scalar (localtime) ."\n\n";

exit;

sub check {
	my ($db, $file) = @_;
	my $query = "select filesize.filesize from filesize where lane_names_orig = ('$file');";

	my $sql = $db->prepare($query);
	$sql ->execute();
	my $size;
	while (my $row = $sql->fetchrow_arrayref) {
		my $row2 = join("\t", @$row);
		my @row2array = split ("\t",$row2);
		$size = $row2array[0];
		unless ($size>10000){

			print "file too small!  ".$file."\t".$size."\n";
		}
	}
        return $size;
}

sub update {
	my ($db, $size, $file,$md5ok,$md5)= @_;

	#my $query = "update filesize set filesize.filesize = ? where lane_names_orig = ('$file');";
	#my $sql = $db ->prepare($query);
	#$sql -> execute($size);
        my $query = "update filesize set filesize.md5check = ? where lane_names_orig = ('$file');";
        my $sql = $db ->prepare($query);
        $sql -> execute($md5);
        my $query = "update filesize set filesize.filesize = ? where lane_names_orig = ('$file');";
        my $sql = $db ->prepare($query);
        $sql -> execute($size);

	$size = check($db, $file);
}

sub insert {
	my $gelukt = "";
	my ($db,$gene_adaptor,$hgnc_id)= @_;

	if (($hgnc_id)){
	
		my $query = "insert into hugo_ensembl values (?,?,NULL)";
				
		my $sql = $db->prepare($query);
	
		my @genes = @{$gene_adaptor->fetch_all_by_external_name($hgnc_id)};
	
		if ((@genes)){
			foreach my $gene (@genes){
				my $gene_id = $gene -> stable_id();
				print $gene_id."\n";
				my $sql = $db->prepare($query);
						
				$sql -> execute($hgnc_id,$gene_id);
				check($db,$hgnc_id);
				my $gelukt = 'yes';

			}
		}
		else {
			my $gene_id;
			my $sql = $db->prepare($query);
						
			$sql -> execute($hgnc_id,$gene_id);
			check($db,$hgnc_id);
			my $gelukt = 'yes';
		}
	}

	return $gelukt;
}
