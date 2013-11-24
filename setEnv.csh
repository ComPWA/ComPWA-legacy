setenv GENEVA           /usr/local/geneva
setenv ROOTSYS          /opt/root/5.32-00
setenv ROOTINCLUDE      $ROOTSYS/include

if($?LD_LIBRARY_PATH) then
  setenv LD_LIBRARY_PATH  ${LD_LIBRARY_PATH}:$GENEVA/lib:`pwd`/lib
else
  setenv LD_LIBRARY_PATH  $GENEVA/lib:`pwd`/lib
endif
