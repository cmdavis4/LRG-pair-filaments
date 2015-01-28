#!/bin/bash

gitdir=$(git rev-parse --show-toplevel);
cd ${gitdir};
cd ..;
topdir=$(pwd);
cd ${gitdir};

#sed_str="s|{gitdir}|${gitdir}|"
#sed $sed_str <$gitdir/src/mysql/load_CAMEO.raw >$gitdir/src/mysql/load_CAMEO.mysql

pathfile=${gitdir}'/user_settings.py';

echo "Building pathfile: ${pathfile}";

############ PATHS ##########################
echo "##### PATHS  #####">$pathfile;
echo "project_path='$gitdir'">>$pathfile;

exit 0;

