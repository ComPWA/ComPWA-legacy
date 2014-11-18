#!/bin/bash

function check_install_dir {
  if [ ! -d "${install_path}" ]; then
    echo "Installation directory does not exist! Please set installation path correctly and rerun."
      exit 1
  else 
    echo "Installation directory exists!"
  fi
}

function check_compiler {
if [[ "`gcc --version | head -n 1 | awk '{print $3}'`" > "4.5" ]]; then
  echo "Compatible gcc version is active! Everything is fine."
else
  echo "Please make a version of gcc 4.5 or higher available."
  exit 1
fi

compiler_path=`which c++`
}

function check_dependencies {
  hash wget 2>/dev/null || { echo >&2 "wget required, but it's not installed.  Aborting."; exit 1; }
  hash tar 2>/dev/null || { echo >&2 "tar required, but it's not installed.  Aborting."; exit 1; }
  hash cmake 2>/dev/null || { echo >&2 "cmake required, but it's not installed.  Aborting."; exit 1; }
}
