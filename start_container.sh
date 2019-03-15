#!/bin/bash

#----------------------------------
# Start YLoc Daemon Container
#----------------------------------

contact_email="abi-services@informatik.uni-tuebingen.de"
imprint_url="https://www-abi.informatik.uni-tuebingen.de/imprint"
gdpr_url="https://www-abi.informatik.uni-tuebingen.de/gdpr"

docker run --rm -it -d -p 28010:80 \
           -e YL_CONTACT_EMAIL="$contact_email" \
           -e YL_IMPRINT_URL="$imprint_url" \
           -e YL_GDPR_URL="$gdpr_url" \
           --name abi_webservice_yloc yloc
