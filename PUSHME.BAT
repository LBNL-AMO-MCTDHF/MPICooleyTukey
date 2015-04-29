
if [[ $1 == "tags" ]]
then

  git push --tags https://github.com/LBNL-AMO-MCTDHF/CooleyTukeyMPI master

elif [[ $1 == "force" ]]
then

  git push -f https://github.com/LBNL-AMO-MCTDHF/CooleyTukeyMPI master

else

  git push https://github.com/LBNL-AMO-MCTDHF/CooleyTukeyMPI master

fi

