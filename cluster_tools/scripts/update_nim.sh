#!/bin/bash
set -e

txtund=$(tput sgr 0 1)    # Underline
txtbld=$(tput bold)       # Bold
txtred=$(tput setaf 1)    # Red
txtgrn=$(tput setaf 2)    # Green
txtylw=$(tput setaf 3)    # Yellow
txtblu=$(tput setaf 4)    # Blue
txtpur=$(tput setaf 5)    # Purple
txtcyn=$(tput setaf 6)    # Cyan
txtwht=$(tput setaf 7)    # White
txtrst=$(tput sgr0)       # Text reset

NIM=$HOME/progs/nim
OPTION=$1
BUILT="  "

if [[ $OPTION == "-h" || $OPTION == "--help" ]]; then
  echo "Nim updating tool"
  echo "usage: update_nim [csources|withdocs]"
  echo "       'csources' will rebuild everything from csources"
  exit 0
fi;

echo -e "${txtylw}\nChecking for changes...${txtrst}"
sleep 2
cd $NIM
git fetch
GIT=$(git pull)

if [[ $GIT == "Already up-to-date." ]]; then
  echo -e "\nNothing to do."
  if [[ $OPTION != "csources"  ]]; then
    exit 0
  fi
fi

echo -e "${txtylw}\nChecking the csources for changes...${txtrst}"
sleep 2
cd csources
GIT=$(git pull)
cd ..

if [[ $GIT == "Already up-to-date." ]]; then
  echo "No changes in csources."
else
  echo "Csources were changed."
  OPTION="csources"
fi

if [[ $OPTION == "csources" ]]; then
  echo -e "${txtylw}\nBuilding from Csources...${txtrst}"
  sleep 2
  cd csources
  sh build.sh
  cd ..
  BUILT="$BUILT Csources"
fi

echo -e "${txtylw}\nBuilding Nim...${txtrst}"
sleep 2
./koch boot -d:release
BUILT="$BUILT Nim"

echo -e "${txtylw}\nBuilding Koch...${txtrst}"
sleep 2
bin/nim c koch
BUILT="$BUILT Koch"

echo -e "${txtylw}\nBuilding Tools...${txtrst}"
sleep 2
./koch tools
BUILT="$BUILT Tools (NimSuggest Nimble)"

if [[ $OPTION == "csources" || $OPTION == "withdocs" ]]; then
  echo -e "${txtylw}\nBuilding Docs...${txtrst}"
  sleep 2
  ./koch doc
  BUILT="$BUILT Docs"
fi

echo -e "${txtgrn}\nSuccessfully built:${txtrst}"
echo $BUILT
