#!/bin/bash


TITLE="Adela's cfMeDIP-seq Technology"

OUTPUT="_index.Rmd"    ####  will build index.Rmd

echo '---' > $OUTPUT
echo 'title: '$TITLE >> $OUTPUT
echo '---' >> $OUTPUT
echo >> $OUTPUT

files="`ls *.html`"
for f in $files
do
  if [ $f != "index.html" ]
  then
    echo '<li><a href="'$f'">'${f}'</a></li>' >> $OUTPUT
  fi
done

###     echo '<li><a href="'$f'">'${f#_}'</a></li>' >> $OUTPUT

