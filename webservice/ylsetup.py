import os

f = file("/var/www/html/cgi-bin/ylconfig.py", "w")

if os.environ.get("YL_CONTACT_EMAIL") != None:
    f.write("contact_email = '" + os.environ.get("YL_CONTACT_EMAIL") + "'\n")
else:
    f.write("contact_email = 'abi-services@informatik.uni-tuebingen.de'\n")

if os.environ.get("YL_IMPRINT_URL") != None:
    f.write("imprint_url = '" + os.environ.get("YL_IMPRINT_URL") + "'\n")
else:
    f.write("imprint_url = 'https://abi.inf.uni-tuebingen.de/impressum'\n")

if os.environ.get("YL_GDPR_URL") != None:
    f.write("gdpr_url = '" + os.environ.get("YL_GDPR_URL") + "'\n")
else:
    f.write("gdpr_url = 'https://abi.inf.uni-tuebingen.de/datenschutzerklaerung'\n")

f.close()
