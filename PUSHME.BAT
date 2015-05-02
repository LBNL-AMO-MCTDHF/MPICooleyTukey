
if [[ $1 == "tags" ]]
then

  git push --tags https://github.com/LBNL-AMO-MCTDHF/MPICooleyTukey master

elif [[ $1 == "force" ]]
then

  git push -f https://github.com/LBNL-AMO-MCTDHF/MPICooleyTukey master

else

  git push https://github.com/LBNL-AMO-MCTDHF/MPICooleyTukey master

fi

