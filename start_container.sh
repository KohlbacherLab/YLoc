#!/bin/bash

#----------------------------------
# Important Settings
#----------------------------------

if [ -z "$ABI_SERVICES_CONTACT_MAIL" ]
then
  contact_email="abi-services@informatik.uni-tuebingen.de"
else
  contact_email="$ABI_SERVICES_CONTACT_MAIL"
fi

if [ -z "$ABI_SERVICES_IMPRINT_URL" ]
then
  imprint_url="https://www-abi.informatik.uni-tuebingen.de/imprint"
else
  imprint_url="$ABI_SERVICES_IMPRINT_URL"
fi

if [ -z "$ABI_SERVICES_GDPR_URL" ]
then
  gdpr_url="https://www-abi.informatik.uni-tuebingen.de/gdpr"
else
  gdpr_url="$ABI_SERVICES_GDPR_URL"
fi

if [ -z "$ABI_SERVICES_YLOC_PORT" ]
then
  yloc_port="28010"
else
  yloc_port="$ABI_SERVICES_YLOC_PORT"
fi


#----------------------------------
# Start YLoc Daemon Container
#----------------------------------

docker run --rm -it -d -p $yloc_port:80 \
           -e YL_CONTACT_EMAIL="$contact_email" \
           -e YL_IMPRINT_URL="$imprint_url" \
           -e YL_GDPR_URL="$gdpr_url" \
           --name abi_webservice_yloc yloc
