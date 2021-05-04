allThreads=(1 2 4 8 16 32 64 128)
allRuntimes=()
for t in ${allThreads[@]}; do
  echo $t
done

read -p 'Continue? (Y/n): ' cont

if [[ $cont == 'Y' ]]
then
  echo 'running script'
else
  echo 'process stopping'
fi