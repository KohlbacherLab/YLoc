#!/bin/bash

#----------------------------------
# Important Settings
#----------------------------------

# Enter a valid eMail address that allows to contact the responsible colleague for this webservice
if [ -z "$ABI_SERVICES_CONTACT_MAIL" ]
then
  contact_email="abi-services@informatik.uni-tuebingen.de"
else
  contact_email="$ABI_SERVICES_CONTACT_MAIL"
fi

# Enter a valid URL that leads to the appropriate imprint page
if [ -z "$ABI_SERVICES_IMPRINT_URL" ]
then
  imprint_url="https://www-abi.informatik.uni-tuebingen.de/imprint"
else
  imprint_url="$ABI_SERVICES_IMPRINT_URL"
fi

# Enter a valid URL that leads to the appropriate GDPR declaration page
if [ -z "$ABI_SERVICES_GDPR_URL" ]
then
  gdpr_url="https://www-abi.informatik.uni-tuebingen.de/gdpr"
else
  gdpr_url="$ABI_SERVICES_GDPR_URL"
fi

# Here you can set an upper limit for the number of sequences allowed to be submitted
if [ -z "$ABI_SERVICES_YLOC_MAX_SEQ" ]
then
  yloc_max_seq="25"
else
  yloc_max_seq="$ABI_SERVICES_YLOC_MAX_SEQ"
fi

# Here you can specify th host port that is bound to port 80 from the container
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
           -e YL_MAX_SEQ="$yloc_max_seq" \
           --name abi_webservice_yloc yloc
