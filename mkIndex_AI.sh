#!/bin/bash


TITLE="Adela's cfMeDIP-seq Technology - AI Reports"

OUTPUT="_index_AI.Rmd"    ####  will build index.Rmd

echo '---' > $OUTPUT
echo 'title: '$TITLE >> $OUTPUT
echo '---' >> $OUTPUT
echo >> $OUTPUT

files="`ls AI/*.html`"
for f in $files
do
    echo '<li><a href="'$f'">'${f}'</a></li>' >> $OUTPUT
done

###     echo '<li><a href="'$f'">'${f#_}'</a></li>' >> $OUTPUT

