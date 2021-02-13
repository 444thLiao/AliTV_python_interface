with import <nixpkgs> {};
stdenv.mkDerivation rec {
  name = "reinforced-nethack";
  env = buildEnv { name = name; paths = buildInputs; };
  gffutils = python3.pkgs.buildPythonPackage rec {
    pname = "gffutils";
    version = "0.10.1";
    src = python3.pkgs.fetchPypi {
      inherit pname version;
      sha256 = "1jsiwd0r9z73df4rkjs9g27ic2r1w900c5iqf8a578vsdl03kz58";
    };
    propagatedBuildInputs = [
      python3Packages.argh
      python3Packages.six
      python3Packages.simplejson
      python3Packages.pyfaidx
      python3Packages.argcomplete
    ];
    doCheck = false;
  };
  ete3 = python3.pkgs.buildPythonPackage rec {
    pname = "ete3";
    version = "3.1.1";
    src = python3.pkgs.fetchPypi {
      inherit pname version;
      sha256 = "11f3p5zkgjvsxi8kbscsaf286qx7vywwdiqk9gdgndka955ks2l7";
    };
    doCheck = false;
  };
  buildInputs = [
    python3
    python3Packages.biopython
    python3Packages.pandas
    python3Packages.click
    python3Packages.tqdm
    ete3
    gffutils
  ];
}
