for i in `ls | grep c1`
do
 #echo $i
 cd $i
 for j in `ls`
 do
   #echo $j
   if [ -d $j ]
   then
   cd $j
   for k in `ls`
   do
     #echo $k
     if [ -d $k ]
     then
     cd $k
     ls *.fa
     cd ..
     else
     echo $k
     fi
   done
   cd ..
   else
   echo $j
   fi
 done
 cd ..
done
