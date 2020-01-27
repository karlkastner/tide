#!/usr/bin/perl
# 2013-05-19 16:31:59
# Karl Kästner, Berlin

# generate tpxo input

use warnings;
use strict;

use Text::CSV;
my $csv = Text::CSV->new({ sep_char => ';' });
my $file = 'lat_lon.csv';
my $folder;
my $lat;
my $lon;
my $y0;
my $ye;

# TODO no magic path
#my $root="/home/pia/phd/src/lib/tide";
my $root="/home/pia/large/phd/gis/tide-tpxo/OTPS1_tpxo7/";

if ($#ARGV > 0)
{
#	$file   = $ARGV[0];
	$folder = $ARGV[1];
	$lat    = $ARGV[2];
	$lon    = $ARGV[3];
	$y0     = $ARGV[4];
	$ye     = $ARGV[5];
	for (my $i=0; $i<=$#ARGV; $i++)
	{
		print $i." ".$ARGV[$i]."\n";
	}
} else {
	die "Invalid input arguments";
}
#open(my $data, '<', $file) or die "Could not open '$file' $!\n";
my $llt_name = $folder."/lat_lon_time";

#my $lat;
#my $lon;
#$lat = 0;
#$lon = 109;
#my $line = <$data>;
#chomp $line;
#if ($csv->parse($line))
#{
#	my @fields = $csv->fields();
#	$lat = $fields[0];
#	$lon = $fields[1];
#}
#close($data);
if (!defined($lat))
{
	die "Lat undefined.";
} else {
	print $lat."\n";
}

if (!defined($lon))
{
	die "Lon undefined.";
}

# export target times
print "Generating TPXO lat/lon/time file\n";
open(my $llt, '>', $llt_name) or die "Could not open ".$llt_name."\n";

my @dom = (31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31);
#@mih = ( 0, 12, 24, 36, 48 );
my @mih = ( 0, 30 );
#@mih = ( 0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55 );

# year
for (my $y = $y0; $y <= $ye; $y++)
{

	# month
	for (my $m = 1; $m <=12; $m++)
	{
		# day
		my $nd = $dom[$m-1];
		# leap year rules
                if ( 2 == $m && 0 == $y % 4)                                    
                {                                                               
                        if (0 != $y % 100 || 0 != $y % 400)
			{
				$nd = 29;
			}
		}

		for (my $d = 1; $d <= $nd; $d++)
		{
			# hour
			for (my $h=0; $h<24; $h++)
			{
				# minutes
				for (my $mi = 0; $mi < scalar(@mih); $mi++)
				{
					print $llt "$lat $lon $y $m $d $h $mih[$mi] 0\n";
					#n = n+1;
				}
			}
		}
	}
}
close($llt);


my $hc_name;
$hc_name = $folder."/hc.inp";
open(my $hc_fid, '>', $hc_name) or die "Could not open '$hc_name' $!\n";
print $hc_fid $root."/DATA/Model_tpxo7           ! 1. tidal model control file\n"
.$folder."/lat_lon_time       ! 2. latitude/longitude/<time> file\n"
."z                          ! 3. z/U/V/u/v\n"
."                           ! 4. tidal constituents to include\n"
."AP                         ! 5. AP : Amplitude,Phase; RI : Real,Imaginary\n"
."oce                        ! 6. oce/geo		???\n"
."1                          ! 7. 1/0 correct for minor constituents\n"
.$folder."/hc.out          ! 8. output file (ASCII)\n";
close ($hc_fid);

my $level_name;
$level_name = $folder."/level.inp";
open(my $level_fid, '>', $level_name) or die "Could not open '$level_name' $!\n";
print $level_fid $root."/DATA/Model_tpxo7           ! 1. tidal model control file\n"
.$folder."/lat_lon_time       ! 2. latitude/longitude/<time> file\n"
."z                          ! 3. z/U/V/u/v\n"
."                           ! 4. tidal constituents to include\n"
."AP                         ! 5. AP : Amplitude,Phase; RI : Real,Imaginary\n"
."oce                        ! 6. oce/geo		???\n"
."1                          ! 7. 1/0 correct for minor constituents\n"
.$folder."/level.out          ! 8. output file (ASCII)\n";
close ($level_fid);

my $current_name;
$current_name = $folder."/current.inp";
open(my $current_fid, '>', $current_name) or die "Could not open '$current_name' $!\n";
print $current_fid $root."/DATA/Model_tpxo7           ! 1. tidal model control file\n"
.$folder."/lat_lon_time       ! 2. latitude/longitude/<time> file\n"
."u                          ! 3. z/U/V/u/v\n"
."                           ! 4. tidal constituents to include\n"
."AP                         ! 5. AP : Amplitude,Phase; RI : Real,Imaginary\n"
."oce                        ! 6. oce/geo		???\n"
."1                          ! 7. 1/0 correct for minor constituents\n"
.$folder."/current.out          ! 8. output file (ASCII)\n";
close ($current_fid);

