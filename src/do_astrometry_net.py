

import requests
import json
import init
import secrets

from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from email.mime.application  import MIMEApplication
from email.encoders import encode_noop

R = requests.post('http://nova.astrometry.net/api/login', data={'request-json': json.dumps({"apikey": secrets.astrometrynet_api})})
sessionkey=R.json().get("session")
print(sessionkey)
filename=init.fitsdir + init.reference_frame
files={'file': open(filename, 'rb')}
jsondata={"session": sessionkey, "publicly_visible": "n", "allow_commercial_use": "n"}
#print(files, jsondataå)
m1 = MIMEBase('text', 'plain')
m1.add_header('Content-disposition',
                'form-data; name="request-json"')
m1.set_payload(jsondata)
m2 = MIMEApplication(filename,'octet-stream',encode_noop)
m2.add_header('Content-disposition',
                'form-data; name="file"; filename="%s"'%filename)
mp = MIMEMultipart('form-data', None, [m1, m2])
print(mp)
R2 = requests.post('http://nova.astrometry.net/api/url_upload', data=mp)

print(R2.text)