#!/bin/bash


# Makes file _site.yml 


OUTPUT="_site.yml"
TITLE="Adela's cfMeDIP-seq Technology - AI Reports"

#DESCRIPTION="Put subtitle here"

echo 'name: "reAdela"' > $OUTPUT
echo 'output_dir: "_site"' >> $OUTPUT
echo 'navbar:' >> $OUTPUT
echo '  title: '$TITLE >> $OUTPUT
echo '  left:' >> $OUTPUT
echo '    - text: "Home"' >> $OUTPUT
echo '      href: ../index.html' >> $OUTPUT


