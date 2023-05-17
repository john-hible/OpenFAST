# RunTurbSim.pl
# See Marshall Buhl's RunNTM Perl script for help with this script.
# Run a series of TurbSim tests, changing the seed for each case.

# NOTICE
# ------
# Look for the string "USER" for lines that you may want to change.


#-------------------------------------------------------------------
#  User-specified paths to programs and files.  USER: YOU WILL DEFINITELY NEED TO CHANGE THESE.
#-------------------------------------------------------------------

$turbsim    = "TurbSim";

# The TurbSim input/output file name:
$input_file = "InputFiles/WindPACT/WPact_W08.inp";

# The Directory to place results in:
$rslt_dir   = "J:/TurbSim_results/TestCases/WindPACT/FullField";

# The name of the resulting file in the results directory (will be appended with the seed number)
$rslt_file  = "WPact_TurbS_CTS_W08";

# The extensions of TurbSim output to copy to the results directory:
#@out_ext      = ("inp","sum","wnd","cts");
@out_ext       = ("sum","wnd","cts");

#-------------------------------------------------------------------
#  Set the number of seeds to run. Initialize the first seed.  USER: YOU MAY WANT TO CHANGE THESE
#-------------------------------------------------------------------

$seed_1  = 816416316;
$n_seeds = 1;

#-------------------------------------------------------------------
#  Get and save the start time. Initialize other variables
#-------------------------------------------------------------------

$StartTime  = date_time();
$n_runs     = 0;
$seed_line  = 4;

srand($seed_1);
#===================================================================
# Loop through all seeds.
#===================================================================

for $i_seed ( 1..$n_seeds )
{

    $ts_root = "${rslt_file}";

    #-------------------------------------------------------------------
    #     Create TurbSim Input File
    #-------------------------------------------------------------------

      $in_file = "${ts_root}.inp";

      open( IN_FILE , ">$in_file" )   or die( "Can't open '$in_file'.\n"  );
      open( IN_FILE1, "$input_file" ) or die( "Can't open '$input_file'.\n" );


      $cnt = 1;

      while ( <IN_FILE1> )
      {
          if ( $cnt == $seed_line )
          {
              printf ( IN_FILE "%s\t%s\n", "$seed_1", "The first seed" );
          }
          else
          {
                  print IN_FILE;
          }

          $cnt++;

      } # while <IN_FILE1>

      close( IN_FILE  );
      close( IN_FILE1 );


  #-------------------------------------------------------------------
  #        Run TurbSim
  #-------------------------------------------------------------------

      print "\n ==============================================\n";
      print   " TurbSim Seed $i_seed  = $seed_1 \n";
      print   " ==============================================\n";
      print   " Starting TurbSim on ", date_time(), ".\n";
      system( "$turbsim $ts_root.inp > NUL" );

  #-------------------------------------------------------------------
  # Rename output
  #-------------------------------------------------------------------

    $root = sprintf( "%s/%s_%s", $rslt_dir, $ts_root, $i_seed);

    for $i_r ( 0..$#out_ext )
    {
      rename( "$ts_root.$out_ext[$i_r]" , "$root.$out_ext[$i_r]" );
    }

  #-------------------------------------------------------------------
  # Get next seed
  #-------------------------------------------------------------------

    $seed_1 = int( rand(2147483646) ) + 1;

    $n_runs ++;


#  End of $i_seed loop
}


#-------------------------------------------------------------------
#     Done.
#-------------------------------------------------------------------

print "\n ============================================\n";
if ( $n_runs > 1 )
{
   print " After starting on $StartTime, the $n_runs batch runs\n";
   print " completed on ", date_time(), ".\n\n";
}
else
{
   print " The single batch run completed on ", date_time(), ".\n\n";
}

#-------------------------------------------------------------------
# End of main function.
#-------------------------------------------------------------------

#*******************************************************************
#*******************************************************************
# This routine returns the date and time in the form "dd-Mon-ccyy at hh:mm:ss".

sub date_time
{
   my( $sec, $min, $hour, $mday, $mon, $year, $mon_str );

   ($sec, $min, $hour, $mday, $mon, $year) = localtime;

   $mon_str = (Jan,Feb,Mar,Apr,May,Jun,Jul,Aug,Sep,Oct,Nov,Dec)[$mon];

   return sprintf( "%2.2d-%s-%4.4d at %2.2d:%2.2d:%2.2d", $mday, $mon_str, $year+1900, $hour, $min, $sec );
}

__END__

