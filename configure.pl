#!/usr/bin/perl
use warnings;

# TODO:
# - implement automatic check for dependencies?

# create the Jamfile of a top level project
sub createTopLevelJamfile {
  if (@_) {
    my $TOP = uc($_[0]);
    my $subproject;
    my $dependency;
    local *OUT;
    
    open( OUT, '+>./'.$_[0].'/Jamfile' ) or die $!;

    foreach my $project ( @{$_[1]} ) {
      $subproject = $subproject . "build-project $project ;\n";
      $dependency = $dependency . "        \$($TOP)/$project//$project\n";
    }

    print OUT "########### $_[0] Jamfile #############\n" .
        "path-constant $TOP : . ;\n" .
        "\n" .
        "# lib$_[0].so will consist of these packages:\n" .
        $subproject .
        "\n" .
        "# Build lib$_[0].so\n" .
        "lib $_[0] : \n" .
        "        [ glob *.cpp ]\n" .
        $dependency .
        "        ;\n" .
        "\n" .
        "install dist : $_[0] : <location>\$(TOP)/lib ;\n";
    close( OUT );
  }
}

# open Jamfile of a project and read its description
# enclosed in @brief...@end
sub readProjectDescription {
  if (@_) {
    local *CF;
    open( CF,'<'.$_[2].'/Jamfile' ) || die "Open $_[2]: $!";
    read( CF, my $data, -s $_[2] );
    close( CF );

    if( $data =~ /\@brief(.*)\@end/s ) {
      my $result = $1;
      $result =~ s/^\s*#\s*//mg;
      $result =~ s/[\n,\r]/ /mg;
      print "\nPackage $_[0]->$_[1]:\n" . $result . "\n".
          "Include this package in build [Y|n|?]? ";        
    } else {
      print "No description for package $_[0]->$_[1] available!\n".
          "Include this package in build [Y|n|?]? ";        
    }
  }
}

sub readYN {
  while( 1 ) {
    chomp( my $input = <STDIN> );
    if( $input =~ m/^[Y]?$/i ) {
      return 1;
    } elsif ( $input =~ m/^[N]$/i ) {
      return 0;
    } elsif ( $input =~ m/\?/ ) {
      readProjectDescription( @_ );
    } else {
      print "Please enter [Y|n|?]! ";
    }
  }
}

# get list of all projects
my @projects = `find ./ -name Jamfile -exec dirname "{}" ";"`;
# Sort list according to ASCII Numeric standards.
@projects = sort(@projects);

my @build_subprojects = ();
my @build_projects = ();

foreach (@projects) {
  # split dirnames into top level projects and sub level projects
  # omit leading './'
  chomp($_);
  ($toplevel, $sublevel) = split('/', substr( $_, 2 ), 2 );

  next if ( $toplevel eq 'Core' ); # skip 'Core' package
  next if ( $toplevel eq 'test' ); # skip 'test' package

  if ( !defined $sublevel ) {
    print "------------------------------------------\n" .
        "Configuring " . $toplevel . " package:\n";
    push @build_projects, $toplevel;
  } else {
    print "Build sub project '" . $sublevel . "' [Y|n|?]? ";
    if ( readYN( $toplevel, $sublevel, $_ ) ) {
      push @{ $buid_subprojects{$toplevel} }, $sublevel
    }
  }
}

# create/overwrite Jamfiles
foreach (@build_projects) {
  createTopLevelJamfile( $_, \@{ $buid_subprojects{$_} } );
}
